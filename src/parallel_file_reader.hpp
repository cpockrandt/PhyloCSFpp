#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include <string>
#include <vector>

#include "translation.hpp"

#include <cassert>

struct alignment_t
{
    uint64_t start_pos = 0;
    std::string chrom;
    std::vector<std::string> ids;
    std::vector<std::string> seqs;
    std::vector<std::vector<uint8_t> > peptides;

    alignment_t(const std::vector<newick_elem> & phylo_array)
    {
        const uint16_t leaves = (phylo_array.size() + 1)/2;
        ids.resize(leaves);
        seqs.resize(leaves, "");
        peptides.resize(leaves, {});
        for (uint16_t i = 0; i < leaves; ++i)
            ids[i] = phylo_array[i].label;
    }
};

class parallel_maf_reader
{
    unsigned threads;
    int fd = -1;
    off_t file_size;
    size_t page_size;

    size_t pages_per_thread;
    size_t last_thread_with_extra_page;

    std::unordered_map<std::string, uint16_t> *fastaid_to_alnid;

    void *file_mem = nullptr;
    size_t *file_range_pos = nullptr;
    size_t *file_range_end = nullptr;

public:

    unsigned get_jobs() const noexcept { return this->threads; }

    parallel_maf_reader(char *file_path, const unsigned threads, std::unordered_map<std::string, uint16_t> *fastaid_to_alnid)
    {
        fd = open(file_path, O_RDONLY);
        if (fd < 0)
        {
            printf("Error: Cannot open %s\n", file_path);
            exit(-1);
        }

        file_size = lseek(fd, 0, SEEK_END);
        if (file_size == (off_t) -1)
        {
            printf("Error: Cannot seek %s\n", file_path);
            close(fd);
            exit(-1);
        }

        file_mem = mmap(nullptr, file_size, PROT_READ, MAP_SHARED, fd, 0);
        if (!file_mem)
        {
            printf("Error: Cannot map %s\n", file_path);
            close(fd);
            exit(-1);
        }

        page_size = sysconf(_SC_PAGE_SIZE);

        const size_t pages = (file_size + page_size - 1) / page_size; // equals ceil(file_size / page_size)
        if (threads > pages)
        {
            // there have to be at least as many pages as threads
            this->threads = pages;
            printf("Info: Using %ld jobs instead of %d.\n", pages, threads);
        }
        else
        {
            this->threads = threads;
        }
//        printf("File size: %ld\nPage size: %ld\nPages: %ld\n\n", file_size, page_size, pages);

        // every thread will process at least floor(pages / this->threads) pages.
        // the first (pages % this->threads) threads will also process an additional page.
        pages_per_thread = pages / this->threads;
        last_thread_with_extra_page = pages % this->threads; // excl
//        printf("Pages per thread: %ld\nLast page (excl.) with excess page: %ld\n\n", pages_per_thread, last_thread_with_extra_page);

        file_range_pos = (size_t*) calloc(sizeof(size_t), this->threads);
        file_range_end = (size_t*) calloc(sizeof(size_t), this->threads);

        size_t file_range_from = 0;
        for (unsigned i = 0; i < this->threads; ++i)
        {
            size_t range_size = page_size * (pages_per_thread + (last_thread_with_extra_page > i));
            size_t file_range_to = file_range_from + range_size;

            if (i == this->threads - 1)
            {
                file_range_to = file_size;
                range_size = file_range_to - file_range_from;
            }

            file_range_pos[i] = file_range_from;
            file_range_end[i] = file_range_to;

//            printf("Thread %2d starting from %10ld (incl.) to %10ld (excl.). Size: %10ld\n", i, file_range_from, file_range_to, range_size);

            // go to beginning of alignment (i.e., go to first page with beginning of a new alignment block and set offset_in_page accordingly)
            // - alignment could start at the very beginning of the current page (or next page). note that except the very first thread, the first page
            //   has to be skipped and cannot be processed (only used for finding the beginning of the next alignment starting from the next page)
            // - there could be no new alignment at all
            // - alignment could start somewhere in the middle of the page
            // - multiple pages could be necessary to be skipped

            {
                // in case the alignment starts at the very first position of the starting page, we need to check whether the character before is a newline
                // but don't do that for the first thread, since the start position is 0 and would trigger an underflow.

                if (!(i == 0 && *((char*) file_mem + file_range_pos[i] + 0) == 'a' && *((char*) file_mem + file_range_pos[i] + 1) == ' '))
                {
                    char *aln_begin = strstr((char*) file_mem + file_range_pos[i], "\na ");

                    if (aln_begin == nullptr)
                    {
                        file_range_pos[i] = file_size;
                    }
                    else
                    {
                        file_range_pos[i] = aln_begin - (char*) file_mem + 1; // skip \n
                    }
                }
                // else
                //    file_range_pos[i] += 0; // no offset
            }
//            printf("%3ld, %3ld\n", file_range_pos[i], file_range_end[i]);

            file_range_from = file_range_to;
        }

        this->fastaid_to_alnid = fastaid_to_alnid;
    }

    // skip a line
    void skip(const unsigned thread_id)
    {
        char *newline_begin = strstr((char*) file_mem + file_range_pos[thread_id], "\n");

        if (newline_begin == nullptr)
            file_range_pos[thread_id] = file_size;
        else
            file_range_pos[thread_id] = newline_begin - (char*) file_mem + 1; // skip \n
    }

    // this function is only called on lines starting with "s "
    void get_line(const unsigned thread_id, char ** id, char ** seq, uint64_t & start_pos, uint64_t & len_wo_ref_gaps, char & strand, const bool immutable = false)
    {
        assert(file_range_pos[thread_id] < (size_t) file_size);
        assert(*((char*) file_mem + file_range_pos[thread_id]) == 's');

        char *newline_begin = strstr((char*) file_mem + file_range_pos[thread_id], "\n"); // TODO: replace all occurrences with strchr()

        size_t new_file_range_pos;
        if (newline_begin == nullptr)
            new_file_range_pos = file_size;
        else
            new_file_range_pos = newline_begin - (char*) file_mem;

        // TODO: this malloc and memcpy can be avoided by reimplementing strtok that does not modify the input string
        const size_t len = new_file_range_pos - file_range_pos[thread_id];
        char *res = (char*) malloc(len + 1);
        memcpy(res, (char*) file_mem + file_range_pos[thread_id], len + 1);
        res[len] = 0;

        char *strtok_state;
        char *token = strtok_r(res, " ", &strtok_state);
        uint16_t token_id = 0;
        while (token != nullptr)
        {
            if (token_id == 1)
            {
                const size_t token_len = strlen(token);
                *id = (char*) malloc(token_len + 1);
                memcpy(*id, token, token_len);
                (*id)[token_len] = 0;
            }
            else if (token_id == 2)
            {
                start_pos = atoi(token);
            }
            else if (token_id == 3)
            {
                len_wo_ref_gaps = atoi(token);
            }
            else if (token_id == 4)
            {
                strand = *token;
            }
            else if (token_id == 6)
            {
                const size_t token_len = strlen(token);
                *seq = (char*) malloc(token_len + 1);
                memcpy(*seq, token, token_len);
                (*seq)[token_len] = 0;
            }
            token = strtok_r(nullptr, " ", &strtok_state);
//            printf("%d. Token: %s\n", token_id, token);
            ++token_id;
        }

        assert(id != nullptr && seq != nullptr);

        free(res);

        if (!immutable)
            file_range_pos[thread_id] = new_file_range_pos + 1; // skip \n
    }

    char get_char(const unsigned thread_id) const
    {
        assert(file_range_pos[thread_id] < (size_t)file_size);

        return *((char*) file_mem + file_range_pos[thread_id]);
    }

    // in practice aln.ids will already be set and only
    // if only \n is at the end of the range for a thread, don't process it
    // if "\na" is at the end of the range for a thread, process it
    bool get_next_alignment(alignment_t & aln, const unsigned thread_id, const unsigned frame, const char strand)
    {
        if (file_range_pos[thread_id] >= file_range_end[thread_id])
            return false;

        // points to the beginning of "a score=....\n"
        skip(thread_id);

        int64_t ref_seq_id = -1;

        char line_type;
        uint64_t prev_cumulative_len_wo_ref_gaps = 0;

        // get the next alignment
        while (file_range_pos[thread_id] < (size_t)file_size && (line_type = get_char(thread_id)) != 'a')
        {
            if (line_type == 's')
            {
                // "s ID int int char int SEQ\n"
                char *id = nullptr, *seq = nullptr;
                uint64_t start_pos;
                uint64_t tmp_len_wo_ref_gaps;
                char ref_strand;
                get_line(thread_id, &id, &seq, start_pos, tmp_len_wo_ref_gaps, ref_strand);

                // cut id starting after "."
                char *ptr = strchr(id, '.');
                assert(ptr != NULL); // expect format "species_name.chrom_name"
                *ptr = 0;
                char *chrom = ptr + 1;

                // remove leading digits (e.g., galGal6 -> galGal)
                for (size_t id_i = 0; id_i < strlen(id); ++id_i)
                {
                    if (isdigit(id[id_i]))
                    {
                        id[id_i] = 0;
                        break;
                    }
                }

                auto alnid = (*fastaid_to_alnid).find(id);
                if (alnid == (*fastaid_to_alnid).end())
                {
//                    printf("ERROR: Species %s in alignment file does not exist in model (or is not translated correctly)!\n", id);
//                    exit(-1);

                    free(id);
                    free(seq);
                    continue;
                }

                aln.seqs[alnid->second] = seq; // TODO: std::move?

                // first sequence we encounter is the reference sequence. store its length!
                if (ref_seq_id == -1)
                {
                    aln.start_pos = start_pos + 1; // maf file is to be 0-indexed
                    aln.chrom = chrom;
                    ref_seq_id = alnid->second;
                    prev_cumulative_len_wo_ref_gaps = tmp_len_wo_ref_gaps;
                    assert(ref_strand == '+'); // We assume that the alignment is always on the forward strand of the reference sequence. TODO implement fix
//                    printf("Ref_seq_id: %ld\n", ref_seq_id);
                }
                else
                {
//                    printf("%ld\t%ld\n", aln.seqs[ref_seq_id].size(), aln.seqs[alnid->second].size());
                    assert(aln.seqs[ref_seq_id].size() == aln.seqs[alnid->second].size()); // all seqs same length
                }

                free(id);
                free(seq);
            }
            else
            {
                skip(thread_id);
            }
        }

        // for species not included in the alignment assign them the sequence NNNN...NNNN.
        const size_t new_ref_seq_len2 = aln.seqs[ref_seq_id].size();
        for (uint64_t i = 0; i < aln.seqs.size(); ++i)
        {
            if (aln.seqs[i].size() == 0)
            {
                aln.seqs[i] = std::string(new_ref_seq_len2, 'N');
            }
        }

        bool abort_next_alignment = false;

        // check whether we can extend the alignment (i.e., next alignment starts on the next base where the last one ended)
        while (!abort_next_alignment && file_range_pos[thread_id] < (size_t)file_size)
        {
            const size_t old_file_range_pos = file_range_pos[thread_id];
            skip(thread_id); // skip "a score=..."

            char line_type;
            int64_t ref_seq_id2 = -1;
            // get the next alignment
            while (!abort_next_alignment && file_range_pos[thread_id] < (size_t)file_size && (line_type = get_char(thread_id)) != 'a')
            {
                if (line_type == 's')
                {
                    // "s ID int int char int SEQ\n"
                    char *id = nullptr, *seq = nullptr;
                    uint64_t start_pos;
                    uint64_t tmp_len_wo_ref_gaps;
                    char ref_strand;

                    get_line(thread_id, &id, &seq, start_pos, tmp_len_wo_ref_gaps, ref_strand);

                    // cut id starting after "."
                    char *ptr = strchr(id, '.');
                    assert(ptr != NULL); // expect format "species_name.chrom_name"
                    *ptr = 0;
                    char *chrom = ptr + 1;

                    // when the first seq is retrieved (for now it's the reference sequence), check whether it can extend the previous alignment
                    if (ref_seq_id2 == -1 && !((aln.start_pos - 1) + prev_cumulative_len_wo_ref_gaps == start_pos && (strcmp(chrom, aln.chrom.c_str()) == 0)))
                    {
                        abort_next_alignment = true;
                        break;
                    }

                    // remove leading digits (e.g., galGal6 -> galGal)
                    for (size_t id_i = 0; id_i < strlen(id); ++id_i)
                    {
                        if (isdigit(id[id_i]))
                        {
                            id[id_i] = 0;
                            break;
                        }
                    }

                    auto alnid = (*fastaid_to_alnid).find(id);
                    if (alnid == (*fastaid_to_alnid).end())
                    {
//                    printf("ERROR: Species %s in alignment file does not exist in model (or is not translated correctly)!\n", id);
//                    exit(-1);

                        free(id);
                        free(seq);
                        continue;
                    }

                    aln.seqs[alnid->second] += seq; // TODO: std::move?

                    // first sequence we encounter is the reference sequence. store its length!
                    if (ref_seq_id2 == -1)
                    {
                        ref_seq_id2 = alnid->second;
                        prev_cumulative_len_wo_ref_gaps += tmp_len_wo_ref_gaps;
                        assert(ref_strand == '+'); // We assume that the alignment is always on the forward strand of the reference sequence. TODO implement fix
//                    printf("Ref_seq_id: %ld\n", ref_seq_id);
                    }
                    else
                    {
//                    printf("%ld\t%ld\n", aln.seqs[ref_seq_id].size(), aln.seqs[alnid->second].size());
                        assert(aln.seqs[ref_seq_id2].size() == aln.seqs[alnid->second].size()); // all seqs same length
                    }

                    free(id);
                    free(seq);
                }
                else
                {
                    skip(thread_id);
                }
            }

            if (abort_next_alignment)
            {
                file_range_pos[thread_id] = old_file_range_pos;
            }
            else
            {
                // for species not included in the alignment add to them the sequence NNNN...NNNN.
                const size_t new_ref_seq_len2 = aln.seqs[ref_seq_id].size();
                for (uint64_t i = 0; i < aln.seqs.size(); ++i)
                {
                    if (aln.seqs[i].size() != new_ref_seq_len2)
                    {
                        aln.seqs[i] += std::string(new_ref_seq_len2 - aln.seqs[i].size(), 'N');
                    }
                }
            }
        }

//        for (uint64_t i = 0; i < aln.seqs.size(); ++i)
//        {
//            aln.seqs[i] = aln.seqs[i].substr(0, 200);
//            printf("%27s: %s\n", aln.ids[i].c_str(), aln.seqs[i].c_str());
//        }
//        printf("\n\n\n");

        size_t new_ref_seq_len = aln.seqs[ref_seq_id].size();
        for (uint64_t pos = 0; pos < aln.seqs[ref_seq_id].size(); ++pos)
        {
            if (aln.seqs[ref_seq_id][pos] == '-')
            {
                aln.seqs[ref_seq_id][pos] = 'X';
                --new_ref_seq_len;
            }
        }

        for (uint64_t i = 0; i < aln.seqs.size(); ++i)
        {
            if (i != (size_t) ref_seq_id && aln.seqs[i].size() != 0)
            {
                for (uint64_t pos = 0; pos < aln.seqs[i].size(); ++pos)
                {
                    if (aln.seqs[ref_seq_id][pos] == 'X')
                        aln.seqs[i][pos] = 'X';
                }
            }
        }

//        for (uint64_t i = 0; i < aln.seqs.size(); ++i)
//            printf("%27s: %s\n", aln.ids[i].c_str(), aln.seqs[i].c_str());
//        printf("\n\n\n");
//        printf("new_ref_seq_len: %ld\n", new_ref_seq_len);

        for (uint64_t i = 0; i < aln.seqs.size(); ++i)
        {
            std::string & s = aln.seqs[i];
            // for species not included in the alignment assign them the sequence NNNN...NNNN. Is this necessary? (TODO)
            size_t pos_w = 0;
            for (size_t pos_r = 0; pos_r < s.size(); ++pos_r)
            {
                if (s[pos_r] != 'X')
                {
                    s[pos_w] = s[pos_r];
                    ++pos_w;
                }
            }
            s.resize(pos_w);
        }

        // get into right frame / strand
        if (strand == '+')
        {
            // != 1 for +1
            // != 2 for +2
            // != 0 for +3
            while (aln.start_pos % 3 != (frame % 3) && aln.seqs[0].size() > 0)
            {
//                size_t size_before = aln.seqs[0].size();
                for (uint16_t i = 0; i < aln.seqs.size(); ++i)
                {
                    aln.seqs[i].erase(0, 1);
                }
                ++aln.start_pos;
//                size_t size_after = aln.seqs[0].size();
//                printf("R(%ld, %ld, %ld, %ld) ", aln.start_pos - 1, aln.start_pos, size_before, size_after);
            }
//            printf("\n");
        }
        else
        {
            // !=  for -1
            // !=  for -2
            // !=  for -3

            size_t required_mod = 0;
            if (frame == 3) // -3
                required_mod = 1;
            else if (frame == 2) // -2
                required_mod = 2;
            else if (frame == 1) // -1
                required_mod = 0;
            else
                exit(77);

//            std::cout << aln.start_pos + aln.seqs[0].size() << '\n';
            while ((aln.start_pos + aln.seqs[0].size()) % 3 != required_mod && aln.seqs[0].size() > 0)
            {
//                size_t size_before = aln.seqs[0].size();
                for (uint16_t i = 0; i < aln.seqs.size(); ++i)
                {
                    aln.seqs[i].pop_back();
                }
//                ++aln.start_pos;
//                size_t size_after = aln.seqs[0].size();
//                printf("R(%ld, %ld, %ld, %ld) ", aln.start_pos - 1, aln.start_pos, size_before, size_after);
            }
//            printf("\n");

            // compute reverse complement
            for (uint16_t i = 0; i < aln.seqs.size(); ++i)
            {
                std::reverse(aln.seqs[i].begin(), aln.seqs[i].end());
                for (uint64_t j = 0; j < aln.seqs[i].size(); ++j)
                {
                    aln.seqs[i][j] = complement(aln.seqs[i][j]);
                }
            }
        }

//        for (uint64_t i = 0; i < aln.seqs.size(); ++i)
//            printf("%27s: %s\n", aln.ids[i].c_str(), aln.seqs[i].c_str());
//        printf("\n\n\n");

        // translate nucleotides
        for (uint16_t i = 0; i < aln.seqs.size(); ++i)
        {
            // alignment.peptides.push_back(translate(alignment.seqs[i]));
            aln.peptides[i].resize(aln.seqs[i].size() / 3);
            for (uint64_t aa_pos = 0; aa_pos < aln.peptides[i].size(); ++aa_pos)
            {
                aln.peptides[i][aa_pos] = get_amino_acid_id(
                        aln.seqs[i][3 * aa_pos],
                        aln.seqs[i][3 * aa_pos + 1],
                        aln.seqs[i][3 * aa_pos + 2]);
            }
        }

        return true;
    }

    ~parallel_maf_reader()
    {
        if (fd >= 0)
            close(fd);

        if (munmap(file_mem, file_size))
            printf("Error munmap %d\n", errno);

        free(file_range_pos);
        free(file_range_end);
    }
};

//int main()
//{
//    parallel_maf_reader maf_rd;
//    maf_rd.init("/home/chris/Downloads/chr22.head.maf", 4);
//
//    #pragma omp parallel for num_threads(4)
//    for (unsigned thread_id = 0; thread_id < 4; ++thread_id)
//    {
//        alignment_t aln;
//        // TODO: verify whether file_range_pos and _end are thread-safe (cache lines)? and merge both arrays for cache locality
//        while (maf_rd.get_next_alignment(aln, thread_id))
//        {
//            #pragma omp critical
//            {
//                for (size_t i = 0; i < aln.ids.size(); ++i)
//                {
//                    printf("%20s\t%s\n", aln.ids[i].c_str(), aln.seqs[i].c_str());
//                }
//                printf("\n");
//            }
//            aln.ids.clear();
//            aln.seqs.clear();
//        }
//    }
//
//    return 0;
//}
