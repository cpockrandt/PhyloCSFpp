#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <cassert>

#include <string>
#include <vector>

#include "translation.hpp"

struct alignment_t
{
    uint64_t start_pos = 0;
    uint64_t chrom_len = 0;
    unsigned skip_bases = 0;
    char strand;

    std::string chrom;
    std::vector<std::string> seqs;
    std::vector<std::vector<uint8_t> > peptides;

    alignment_t(const uint16_t leaves)
    {
        seqs.resize(leaves, "");
        peptides.resize(leaves, {});
    }

    size_t length() const noexcept
    {
        assert(seqs.size() > 0);
        for (uint16_t i = 0; i < seqs.size(); ++i)
        {
            assert(seqs[0].size() == seqs[i].size());
        }
        return seqs[0].size() - skip_bases;
    }

    void translate() noexcept
    {
        // translate nucleotides
        for (uint16_t i = 0; i < seqs.size(); ++i)
        {
            peptides[i].resize(seqs[i].size() / 3);
            for (uint64_t aa_pos = 0; aa_pos < peptides[i].size(); ++aa_pos)
            {
                peptides[i][aa_pos] = get_amino_acid_id(seqs[i][3 * aa_pos],
                                                        seqs[i][3 * aa_pos + 1],
                                                        seqs[i][3 * aa_pos + 2]);
            }
        }
    }

    void update_seqs(const uint64_t orig_start_pos, const char strand, const unsigned frame) noexcept
    {
        skip_bases = 0;
        start_pos = orig_start_pos;

        // NOTE: start_pos is 1-indexed (1 is added after extracting it from the maf-file which is 0-indexed)

        // get into right frame / strand
        if (strand == '+')
        {
            // != 1 for +1,   != 2 for +2,   != 0 for +3
            // while (start_pos % 3 != (frame % 3) && length() > 0)
            // {
            //     ++skip_bases;
            //     ++start_pos;
            // }
            int64_t tmp_skip_bases = (static_cast<int64_t>(frame - start_pos)) % 3;
            if (tmp_skip_bases < 0)
                tmp_skip_bases += 3;
            if (static_cast<uint64_t>(tmp_skip_bases) > length())
                tmp_skip_bases = length();

            skip_bases = tmp_skip_bases;
            start_pos += tmp_skip_bases;
        }
        else
        {
            // while (length() > 0 && ((chrom_len - (start_pos + length() - 1) + 1) % 3 != (frame % 3)))
            // {
            //     ++skip_bases;
            // }
            int64_t tmp_skip_bases = (static_cast<int64_t>(frame - (chrom_len - (start_pos + length()) + 2))) % 3;
            if (tmp_skip_bases < 0)
                tmp_skip_bases += 3;
            if (static_cast<uint64_t>(tmp_skip_bases) > length())
                tmp_skip_bases = length();

            skip_bases = tmp_skip_bases;
        }

        // translate nucleotides
        for (uint16_t i = 0; i < seqs.size(); ++i)
        {
            // alignment.peptides.push_back(translate(alignment.seqs[i]));
            peptides[i].resize((seqs[i].size() - skip_bases) / 3);
            for (uint64_t aa_pos = 0; aa_pos < peptides[i].size(); ++aa_pos)
            {
                peptides[i][aa_pos] = get_amino_acid_id(seqs[i][skip_bases + 3 * aa_pos],
                                                        seqs[i][skip_bases + 3 * aa_pos + 1],
                                                        seqs[i][skip_bases + 3 * aa_pos + 2]);
            }
        }
    }
};

char* rstrstr(char *haystack_end, char *haystack_begin, const char *needle)
{
    size_t needle_length = strlen(needle);

    for (char *p = haystack_end; p >= haystack_begin; --p)
    {
        for (size_t i = 0; i < needle_length; ++i)
        {
            if (p[i] != needle[i])
                goto next;
        }
        return p;

        next:;
    }
    return NULL;
}

class parallel_maf_reader
{
    unsigned jobs;
    int fd = -1;
    off_t file_size;
    size_t page_size;

    size_t pages_per_thread;

    const std::unordered_map<std::string, uint16_t> *fastaid_to_alnid;
    std::unordered_set<std::string> unresolved_identifiers;

    const bool concatenate_alignments;

    void *file_mem = nullptr;
    size_t *file_range_pos = nullptr;
    size_t *file_range_pos_orig = nullptr;
    size_t *file_range_end = nullptr;

    uint64_t *bytes_processing = nullptr;
    uint64_t total_bytes_processed = 0;

    float progress_bar_dimensions_divisor = 1.0;
    char progress_bar_format_str[50];

public:

    unsigned get_jobs() const noexcept { return this->jobs; }

    // skip a line
    void skip(const unsigned job_id)
    {
        char *newline_begin = strchr((char*) file_mem + file_range_pos[job_id], '\n');

        size_t new_file_range_pos;
        if (newline_begin == nullptr)
            new_file_range_pos = file_size;
        else
            new_file_range_pos = newline_begin - (char*) file_mem + 1; // skip \n

        bytes_processing[job_id] += new_file_range_pos - file_range_pos[job_id];

        file_range_pos[job_id] = new_file_range_pos;
    }

    // this function is only called on lines starting with "s "
    void get_line(const unsigned job_id, char ** id, char ** seq, uint64_t & start_pos, uint64_t & len_wo_ref_gaps, char & strand, uint64_t & chrom_len)
    {
        assert(file_range_pos[job_id] < (size_t) file_size);
        assert(*((char*) file_mem + file_range_pos[job_id]) == 's');

        char *newline_begin = strchr((char*) file_mem + file_range_pos[job_id], '\n');

        size_t new_file_range_pos;
        if (newline_begin == nullptr)
            new_file_range_pos = file_size;
        else
            new_file_range_pos = newline_begin - (char*) file_mem;

        // TODO: this malloc and memcpy can be avoided by reimplementing strtok that does not modify the input string
        const size_t len = new_file_range_pos - file_range_pos[job_id];
        char *res = (char*) malloc(len + 1);
        memcpy(res, (char*) file_mem + file_range_pos[job_id], len + 1);
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
            else if (token_id == 5)
            {
                chrom_len = atoi(token);
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

        file_range_pos[job_id] = new_file_range_pos + 1; // skip \n

        bytes_processing[job_id] += len + 1; // counting the skipped \n

        assert(id != nullptr && seq != nullptr);
        free(res);
    }

    char get_char(const unsigned job_id) const
    {
        assert(file_range_pos[job_id] < (size_t)file_size);

        return *((char*) file_mem + file_range_pos[job_id]);
    }

    parallel_maf_reader(const char * file_path, const unsigned jobs, const std::unordered_map<std::string, uint16_t> *fastaid_to_alnid,
                        const bool concatenate_alignments)
        : concatenate_alignments(concatenate_alignments)
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
        if (jobs > pages)
        {
            // there have to be at least as many pages as jobs
            this->jobs = pages;
//            printf("Info: Using %ld jobs instead of %d.\n", pages, jobs);
        }
        else
        {
            this->jobs = jobs;
        }
//        printf("File size: %ld\nPage size: %ld\nPages: %ld\n\n", file_size, page_size, pages);

        // every thread will process at least floor(pages / this->jobs) pages.
        // the first (pages % this->jobs) jobs will also process an additional page.
        pages_per_thread = pages / this->jobs;
        size_t last_thread_with_extra_page = pages % this->jobs; // excl
//        printf("Pages per thread: %ld\nLast page (excl.) with excess page: %ld\n\n", pages_per_thread, last_thread_with_extra_page);

        file_range_pos = (size_t*) calloc(sizeof(size_t), this->jobs);
        file_range_end = (size_t*) calloc(sizeof(size_t), this->jobs);
        file_range_pos_orig = (size_t*) calloc(sizeof(size_t), this->jobs);

        bytes_processing = (size_t*) calloc(sizeof(uint64_t), this->jobs);

        size_t file_range_from = 0;
        for (unsigned i = 0; i < this->jobs; ++i)
        {
            size_t range_size = page_size * (pages_per_thread + (last_thread_with_extra_page > i));
            size_t file_range_to = file_range_from + range_size;

            if (i == this->jobs - 1)
            {
                file_range_to = file_size;
                range_size = file_range_to - file_range_from;
            }

            file_range_pos[i] = file_range_from;
            file_range_end[i] = file_range_to;

//            printf("STRICT SPLITTING: Job %2d starting from %10ld (incl.) to %10ld (excl.). Size: %10ld\n", i, file_range_pos[i], file_range_end[i], file_range_end[i] - file_range_pos[i]);

            // go to beginning of alignment (i.e., go to first page with beginning of a new alignment block and set offset_in_page accordingly)
            // - alignment could start at the very beginning of the current page (or next page). note that except the very first thread, the first page
            //   has to be skipped and cannot be processed (only used for finding the beginning of the next alignment starting from the next page)
            // - there could be no new alignment at all
            // - alignment could start somewhere in the middle of the page
            // - multiple pages could be necessary to be skipped

            // in case the alignment starts at the very first position of the starting page, we need to check whether the character before is a newline
            // but don't do that for the first thread, since the start position is 0 and would trigger an underflow.

            if (!(i == 0 && *((char*) file_mem + file_range_pos[i] + 0) == 'a' && *((char*) file_mem + file_range_pos[i] + 1) == ' '))
            {
                char *aln_begin = strstr((char*) file_mem + file_range_pos[i], "\na ");

                if (aln_begin == nullptr)
                    file_range_pos[i] = file_size;
                else
                    file_range_pos[i] = aln_begin - (char*) file_mem + 1; // skip \n
            }

            file_range_pos_orig[i] = file_range_pos[i];

//            printf("SEARCH ALN BEGIN: Job %2d starting from %10ld (incl.) to %10ld (excl.). Size: %10ld\n", i, file_range_pos[i], file_range_end[i], file_range_end[i] - file_range_pos[i]);

            file_range_from = file_range_to;
        }

        total_bytes_processed += file_range_pos[0]; // skip the bytes that are skipped by the first job
                                                    // subsequent jobs are not considered because the bytes that they skipped here,
                                                    // will be processed (and counted) by previous jobs.

        this->fastaid_to_alnid = fastaid_to_alnid;
    }

    void setup_progressbar(const uint32_t file_id, const uint32_t files)
    {
        uint8_t progress_bar_dimensions_label_index = 0;
        char const *size_dimensions_labels[5] = {"B", "KB", "MB", "GB", "TB"};

        double formatted_filesize = file_size;
        while (formatted_filesize > 1024)
        {
            progress_bar_dimensions_divisor *= 1024;
            formatted_filesize /= 1024;
            ++progress_bar_dimensions_label_index;
        }

        if (files == 1)
        {
            sprintf(progress_bar_format_str, OUT_DEL "%%.2f / %.2f %s (%%3.2f %%%%)\r",
                        formatted_filesize, size_dimensions_labels[progress_bar_dimensions_label_index]);
        }
        else
        {
            sprintf(progress_bar_format_str, OUT_DEL "File %d of %d: %%.2f / %.2f %s (%%3.2f %%%%)\r",
                        file_id,
                        files,
                        formatted_filesize,
                        size_dimensions_labels[progress_bar_dimensions_label_index]);
        }
    }

    void print_progress()
    {
        printf(progress_bar_format_str, total_bytes_processed / progress_bar_dimensions_divisor, 100.0f * total_bytes_processed / file_size);
        fflush(stdout);
    }

    // call this at the beginning of a job, before first alignment is accessed (not for very first job)
    // this will skip any alignments for this job, that will be taken by the previous job (because it might concatenate immediate alignments
    // that fall into this jobs range)
    void skip_partial_alignment(alignment_t & aln, const unsigned job_id)
    {
        if (job_id == 0)
            return;

        // does this alignment get concatenated to the previous one?
        char *prev_aln_begin = rstrstr((char*) file_mem + file_range_pos[job_id] - 2, // haystack_end
                                       (char*) file_mem + file_range_pos_orig[job_id - 1], // haysteck_begin
                                       "\na ") + 1; // pattern; +1 to skip \n

        if (prev_aln_begin == NULL) // no alignment begin in this section, so there are no alignments to process for this job
        {
            file_range_pos[job_id] = file_size;
        }
        else
        {
            // backup
            const size_t orig_file_range_pos = file_range_pos[job_id];

            // get the previous alignment (last one from previous job)
            file_range_pos[job_id] = prev_aln_begin - (char*)file_mem;
            get_next_alignment(aln, job_id);

            // make sure the "progress" (in bytes) is not counted on the next call of get_next_alignment() from reading the previous alignment
            bytes_processing[job_id] = 0;

            for (auto & seq : aln.seqs)
                seq = "";

            // it either cannot be concatenated with the fist alignment of this job
            if (orig_file_range_pos == file_range_pos[job_id])
            {

            }
            // or it "steals" some alignments from this job (which we have to skip), i.e., do nothing
            // and the first call of get_next_alignment() will retrieve the very first complete alignment
            else
            {

            }
        }

//        printf("AFTER SKIPPING  : Job %2d starting from %10ld (incl.) to %10ld (excl.). Size: %10ld\n\n", job_id, file_range_pos[job_id], file_range_end[job_id], file_range_end[job_id] - file_range_pos[job_id]);
    }

    // in practice aln.ids will already be set and only
    // if only \n is at the end of the range for a thread, don't process it
    // if "\na" is at the end of the range for a thread, process it
    bool get_next_alignment(alignment_t & aln, const unsigned job_id)
    {
        #pragma omp atomic
        total_bytes_processed += bytes_processing[job_id];

        bytes_processing[job_id] = 0;

        if (file_range_pos[job_id] >= file_range_end[job_id])
            return false;

        // points to the beginning of "a score=....\n"
        skip(job_id);

        int64_t ref_seq_id = -1;

        char line_type;
        uint64_t prev_cumulative_len_wo_ref_gaps = 0;

        // get the next alignment
        while (file_range_pos[job_id] < (size_t)file_size && (line_type = get_char(job_id)) != 'a')
        {
            if (line_type == 's')
            {
                // "s ID int int char int SEQ\n"
                char *id = nullptr, *seq = nullptr;
                uint64_t start_pos = 0;
                uint64_t chrom_len = 0;
                uint64_t tmp_len_wo_ref_gaps = 0;
                char ref_strand = '.';
                get_line(job_id, &id, &seq, start_pos, tmp_len_wo_ref_gaps, ref_strand, chrom_len);

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
                    #pragma omp critical(unresolved_identifier)
                    {
                        if (unresolved_identifiers.find(id) == unresolved_identifiers.end())
                        {
                            unresolved_identifiers.emplace(id);
                            printf(OUT_INFO "WARNING: Not able to match species %s in alignment file to model (Use `--mapping` to fix it)!\n" OUT_RESET, id);
                        }
                    }

                    free(id);
                    free(seq);
                    continue;
                }

                assert(seq != NULL); // check because of a former bug
                aln.seqs[alnid->second] = seq; // TODO: std::move?

                // first sequence we encounter is the reference sequence. store its length!
                if (ref_seq_id == -1)
                {
                    aln.start_pos = start_pos + 1; // maf file is to be 0-indexed
                    aln.chrom = chrom;
                    aln.chrom_len = chrom_len;
                    aln.strand = ref_strand;
                    ref_seq_id = alnid->second;
                    prev_cumulative_len_wo_ref_gaps = tmp_len_wo_ref_gaps;
                    //assert(ref_strand == '+'); // We assume that the alignment is always on the forward strand of the reference sequence. TODO implement fix
                }
                else
                {
                    assert(aln.seqs[ref_seq_id].size() == aln.seqs[alnid->second].size()); // all seqs same length
                }

                free(id);
                free(seq);
            }
            else
            {
                skip(job_id);
            }
        }

        // TODO: load next alignment if the current alignment was empty (because no species was successfully mapped to the model)

        // for species not included in the alignment assign them the sequence NNNN...NNNN.
        const size_t new_ref_seq_len2 = aln.seqs[ref_seq_id].size();
        for (uint64_t i = 0; i < aln.seqs.size(); ++i)
        {
            if (aln.seqs[i].size() == 0)
            {
                aln.seqs[i] = std::string(new_ref_seq_len2, 'N');
            }
        }

        bool abort_next_alignment = !this->concatenate_alignments;

        // check whether we can extend the alignment (i.e., next alignment starts on the next base where the last one ended)
        while (!abort_next_alignment && file_range_pos[job_id] < (size_t)file_size)
        {
            const size_t old_file_range_pos = file_range_pos[job_id];
            skip(job_id); // skip "a score=..."

            char line_type;
            int64_t ref_seq_id2 = -1;
            // get the next alignment
            while (!abort_next_alignment && file_range_pos[job_id] < (size_t)file_size && (line_type = get_char(job_id)) != 'a')
            {
                if (line_type == 's')
                {
                    // "s ID int int char int SEQ\n"
                    char *id = nullptr, *seq = nullptr;
                    uint64_t start_pos = 0;
                    uint64_t chrom_len = 0;
                    uint64_t tmp_len_wo_ref_gaps = 0;
                    char ref_strand;

                    get_line(job_id, &id, &seq, start_pos, tmp_len_wo_ref_gaps, ref_strand, chrom_len);

                    // cut id starting after "."
                    char *ptr = strchr(id, '.');
                    assert(ptr != NULL); // expect format "species_name.chrom_name"
                    *ptr = 0;
                    char *chrom = ptr + 1;

                    // when the first seq is retrieved (for now it's the reference sequence), check whether it can extend the previous alignment
                    if (ref_seq_id2 == -1 && !((aln.start_pos - 1) + prev_cumulative_len_wo_ref_gaps == start_pos && (strcmp(chrom, aln.chrom.c_str()) == 0)))
                    {
                        abort_next_alignment = true;
                        free(id);
                        free(seq);
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

                    assert(seq != NULL); // check because of a former bug
                    aln.seqs[alnid->second] += seq; // TODO: std::move?

                    // first sequence we encounter is the reference sequence. store its length!
                    if (ref_seq_id2 == -1)
                    {
                        ref_seq_id2 = alnid->second;
                        prev_cumulative_len_wo_ref_gaps += tmp_len_wo_ref_gaps;
                        //assert(ref_strand == '+'); // We assume that the alignment is always on the forward strand of the reference sequence. TODO implement fix
//                    printf("Ref_seq_id: %ld\n", ref_seq_id);
                    }
                    else
                    {
                        assert(aln.seqs[ref_seq_id2].size() == aln.seqs[alnid->second].size()); // all seqs same length
                    }

                    free(id);
                    free(seq);
                }
                else
                {
                    skip(job_id);
                }
            }

            if (abort_next_alignment)
            {
                file_range_pos[job_id] = old_file_range_pos;
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

        // replace gaps in reference sequence with X
        size_t new_ref_seq_len = aln.seqs[ref_seq_id].size();
        for (uint64_t pos = 0; pos < aln.seqs[ref_seq_id].size(); ++pos)
        {
            if (aln.seqs[ref_seq_id][pos] == '-')
            {
                aln.seqs[ref_seq_id][pos] = 'X';
                --new_ref_seq_len;
            }
        }

        // replace chars in other sequences with X that are gaps/X in reference
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

        // delete all X in all sequences
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

        aln.skip_bases = 0;

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
        free(file_range_pos_orig);
        free(bytes_processing);
    }
};
