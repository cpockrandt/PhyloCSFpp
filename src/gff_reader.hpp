#pragma once

#include <sys/mman.h>
#include <cassert>
#include <fcntl.h>


#include <string>
#include <vector>

struct cds_entry
{
    uint64_t begin;
    uint64_t end;
    uint8_t phase;

    cds_entry(const uint64_t begin, const uint64_t end, const uint8_t phase): begin(begin), end(end), phase(phase) { }
};

struct gff_transcript
{
    uint64_t begin = 0;
    uint64_t end = 0;
    std::string chr = "";
    char strand = '.';

    std::vector<cds_entry> CDS;
};

class gff_reader
{
    unsigned jobs;
    int fd = -1;
    off_t file_size;
    size_t page_size;

    size_t pages_per_thread;

    void *file_mem = nullptr;
    size_t *file_range_pos = nullptr;
    size_t *file_range_end = nullptr;

    uint64_t *bytes_processing = nullptr;
    uint64_t total_bytes_processed = 0;

    float progress_bar_dimensions_divisor = 1.0;
    char progress_bar_format_str[50];

public:

    unsigned get_jobs() const noexcept { return this->jobs; }

//    // skip a line
//    void skip(const unsigned job_id)
//    {
//        char *newline_begin = strchr((char*) file_mem + file_range_pos[job_id], '\n');
//
//        size_t new_file_range_pos;
//        if (newline_begin == nullptr)
//            new_file_range_pos = file_size;
//        else
//            new_file_range_pos = newline_begin - (char*) file_mem + 1; // skip \n
//
//        bytes_processing[job_id] += new_file_range_pos - file_range_pos[job_id];
//
//        file_range_pos[job_id] = new_file_range_pos;
//    }

    // this function is only called on lines starting with "s "
    void get_line(const unsigned job_id, char ** id, char ** seq, uint64_t & start_pos, uint64_t & len_wo_ref_gaps, char & strand)
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

    gff_reader(const char * file_path, const unsigned jobs)
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
        }
        else
        {
            this->jobs = jobs;
        }

        // every thread will process at least floor(pages / this->jobs) pages.
        // the first (pages % this->jobs) jobs will also process an additional page.
        pages_per_thread = pages / this->jobs;
        size_t last_thread_with_extra_page = pages % this->jobs; // excl

        file_range_pos = (size_t*) calloc(sizeof(size_t), this->jobs);
        file_range_end = (size_t*) calloc(sizeof(size_t), this->jobs);
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

//            printf("%ld, %ld\n", file_range_pos[i], file_range_end[i]);

//            if (!(i == 0 && *((char*) file_mem + file_range_pos[i] + 0) == 'a' && *((char*) file_mem + file_range_pos[i] + 1) == ' '))
//            {
//                char *aln_begin = strstr((char*) file_mem + file_range_pos[i], "\na ");
//
//                if (aln_begin == nullptr)
//                    file_range_pos[i] = file_size;
//                else
//                    file_range_pos[i] = aln_begin - (char*) file_mem + 1; // skip \n
//            }

            file_range_from = file_range_to;
        }

//        total_bytes_processed += file_range_pos[0]; // skip the bytes that are skipped by the first job
//                                                    // subsequent jobs are not considered because the bytes that they skipped here,
//                                                    // will be processed (and counted) by previous jobs.
    }

//    void setup_progressbar(const uint32_t file_id, const uint32_t files)
//    {
//        uint8_t progress_bar_dimensions_label_index = 0;
//        char const *size_dimensions_labels[5] = {"B", "KB", "MB", "GB", "TB"};
//
//        double formatted_filesize = file_size;
//        while (formatted_filesize > 1024)
//        {
//            progress_bar_dimensions_divisor *= 1024;
//            formatted_filesize /= 1024;
//            ++progress_bar_dimensions_label_index;
//        }
//
//        if (files == 1)
//        {
//            sprintf(progress_bar_format_str, "\x1b[K%%.2f / %.2f %s (%%3.2f %%%%)\r",
//                        formatted_filesize, size_dimensions_labels[progress_bar_dimensions_label_index]);
//        }
//        else
//        {
//            sprintf(progress_bar_format_str, "\x1b[KFile %d of %d: %%.2f / %.2f %s (%%3.2f %%%%)\r",
//                        file_id,
//                        files,
//                        formatted_filesize,
//                        size_dimensions_labels[progress_bar_dimensions_label_index]);
//        }
//    }
//
//    void print_progress()
//    {
//        printf(progress_bar_format_str, total_bytes_processed / progress_bar_dimensions_divisor, 100.0f * total_bytes_processed / file_size);
//        fflush(stdout);
//    }

    bool get_next_transcript(gff_transcript & transcript, const unsigned job_id)
    {
        #pragma omp atomic
        total_bytes_processed += bytes_processing[job_id];

        bytes_processing[job_id] = 0;

        if (file_range_pos[job_id] >= file_range_end[job_id])
            return false;

        // points to the beginning of "a score=....\n"
//        skip(job_id);

//        char line_type;

        uint8_t transcript_occs = 0;
        transcript.CDS.clear();

        // get the next alignment
        while (file_range_pos[job_id] < (size_t)file_size)
        {
            std::string chr = "";
            std::string feature = "";
            uint64_t begin = 0;
            uint64_t end = 0;
            char strand = '.';
            char phase = '.';

            uint8_t col_id = 1;
            size_t col_size;
            char buf[BUFSIZ];
            char *cur_pos = (char *) file_mem + file_range_pos[job_id];
            char *newline_pos = strchr(cur_pos, '\n');
            while (cur_pos < newline_pos)
            {
                char *next_col = strchr(cur_pos, '\t');
                if (next_col == NULL || next_col > newline_pos)
                    next_col = newline_pos;
                col_size = next_col - cur_pos;
                memcpy(buf, cur_pos, col_size);
                buf[col_size] = 0;
//                printf("xxx: %s\n", buf);
                switch(col_id)
                {
                    case 1:
                        chr = buf; break;
                    case 3:
                        feature = buf; break;
                    case 4:
                        begin = std::stoull(buf); break;
                    case 5:
                        end = std::stoull(buf); break;
                    case 7:
                        strand = buf[0]; break;
                    case 8:
                        phase = buf[0]; break;
                    default:
                        break;
                }
                cur_pos = next_col + 1;
                ++col_id;
            }

//            printf("\n");

            if (feature == "transcript")
                ++transcript_occs;

            if (transcript_occs > 1)
                break;
            else
                file_range_pos[job_id] = cur_pos - (char*) file_mem;

            if (feature == "transcript")
            {
                transcript.chr = chr;
                transcript.begin = begin;
                transcript.end = end;
                transcript.strand = strand;
            }
            else if (feature == "CDS")
            {
                transcript.CDS.emplace_back(begin, end, phase - '0');
            }
        }
//            if (line_type == 's')
//            {
//                // "s ID int int char int SEQ\n"
//                char *id = nullptr, *seq = nullptr;
//                uint64_t start_pos = 0;
//                uint64_t tmp_len_wo_ref_gaps = 0;
//                char ref_strand = '.';
//                get_line(job_id, &id, &seq, start_pos, tmp_len_wo_ref_gaps, ref_strand);
//
//                // cut id starting after "."
//                char *ptr = strchr(id, '.');
//                assert(ptr != NULL); // expect format "species_name.chrom_name"
//                *ptr = 0;
//                char *chrom = ptr + 1;
//
//                // remove leading digits (e.g., galGal6 -> galGal)
//                for (size_t id_i = 0; id_i < strlen(id); ++id_i)
//                {
//                    if (isdigit(id[id_i]))
//                    {
//                        id[id_i] = 0;
//                        break;
//                    }
//                }
//
//                aln.seqs[alnid->second] = seq; // TODO: std::move?
//
//                // first sequence we encounter is the reference sequence. store its length!
//                if (ref_seq_id == -1)
//                {
//                    aln.start_pos = start_pos + 1; // maf file is to be 0-indexed
//                    aln.chrom = chrom;
//                    aln.strand = ref_strand;
//                    ref_seq_id = alnid->second;
//                    prev_cumulative_len_wo_ref_gaps = tmp_len_wo_ref_gaps;
//                    assert(ref_strand == '+'); // We assume that the alignment is always on the forward strand of the reference sequence. TODO implement fix
//                }
//
//                free(id);
//                free(seq);
//            }
//            else
//            {
//                skip(job_id);
//            }

        return true;
    }

    ~gff_reader()
    {
        if (fd >= 0)
            close(fd);

        if (munmap(file_mem, file_size))
            printf("Error munmap %d\n", errno);

        free(file_range_pos);
        free(file_range_end);
        free(bytes_processing);
    }
};
