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
    float phylo_score = NAN; // if score cannot be computed/looked up, it will automatically output NAN
    float phylo_power = NAN;

    cds_entry(const uint64_t begin, const uint64_t end, const uint8_t phase): begin(begin), end(end), phase(phase) { }
};

enum feature_t
{
    TRANSCRIPT,
    CDS,
    OTHER
};

struct gff_transcript
{
    std::string chr = "";
    uint64_t begin = 0;
    uint64_t end = 0;
    char strand = '.';
    float phylo_score = NAN; // if score cannot be computed/looked up, it will automatically output NAN
    float phylo_power = NAN;

    std::vector<cds_entry> CDS;
    std::vector<std::tuple<feature_t, std::string> > lines;
};

class gff_reader
{
    int fd = -1;
    off_t file_size;

    void *file_mem = nullptr;
    size_t file_range_pos = 0;

    float progress_bar_dimensions_divisor = 1.0;
    char progress_bar_format_str[50];

public:

    gff_reader(const char * file_path)
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
        printf(progress_bar_format_str, file_range_pos / progress_bar_dimensions_divisor, 100.0f * file_range_pos / file_size);
        fflush(stdout);
    }

    template <bool copy_lines>
    bool get_next_transcript(gff_transcript & transcript)
    {
        if (file_range_pos >= (size_t)file_size)
            return false;

        uint8_t transcript_occs = 0;
        transcript.CDS.clear();
        transcript.lines.clear();

        while (file_range_pos < (size_t)file_size)
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
            char *cur_pos = (char *) file_mem + file_range_pos;
            char *newline_pos = strchr(cur_pos, '\n');

            const size_t line_length = newline_pos - cur_pos;
            char *line = (char*)malloc(line_length + 1);
            memcpy(line, cur_pos, line_length);
            line[line_length] = 0;

            while (cur_pos < newline_pos)
            {
                char *next_col = strchr(cur_pos, '\t');
                if (next_col == NULL || next_col > newline_pos)
                    next_col = newline_pos;
                col_size = next_col - cur_pos;
                memcpy(buf, cur_pos, col_size);
                buf[col_size] = 0;
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

            file_range_pos = cur_pos - (char*) file_mem;

            feature_t f = OTHER;

            if (feature == "transcript")
            {
                f = TRANSCRIPT;
                transcript.chr = chr;
                transcript.begin = begin;
                transcript.end = end;
                transcript.strand = strand;
            }
            else if (feature == "CDS")
            {
                f = CDS;
                transcript.CDS.emplace_back(begin, end, phase - '0');
            }

            if (copy_lines)
            {
                transcript.lines.emplace_back(f, line);
//                printf("xxx: %s\n", line);
            }
        }

        return true;
    }

    ~gff_reader()
    {
        if (fd >= 0)
            close(fd);

        if (munmap(file_mem, file_size))
            printf("Error munmap %d\n", errno);
    }
};
