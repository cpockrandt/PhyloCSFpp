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

class wig_reader
{
    int fd = -1;
    off_t file_size;

    void *file_mem = nullptr;
    size_t file_range_pos = 0;

    uint64_t bytes_processing = 0;
    uint64_t total_bytes_processed = 0;

    float progress_bar_dimensions_divisor = 1.0;
    char progress_bar_format_str[50];

public:

    wig_reader(const char * file_path)
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
        printf(progress_bar_format_str, total_bytes_processed / progress_bar_dimensions_divisor, 100.0f * total_bytes_processed / file_size);
        fflush(stdout);
    }

    bool get_next_scores(std::vector<double> & scores, std::string & chr, uint64_t & start_pos)
    {
        char buf[512];

        total_bytes_processed += bytes_processing;
        bytes_processing = 0;

        while (file_range_pos < (size_t)file_size)
        {
            char *newline_begin = strchr((char*) file_mem + file_range_pos, '\n');
            if (newline_begin == nullptr)
                // did we read in scores or were we already at  the end of the file when calling this function?
                return !scores.empty();

            const uint64_t new_file_range_pos = newline_begin - (char*) file_mem;

            const size_t len = new_file_range_pos - file_range_pos;
            memcpy(buf, (char*) file_mem + file_range_pos, len); // does not copy \n
            buf[len] = 0;

            if (buf[0] == 'f')
            {
                char c_chr[100];
                uint64_t tmp_start_pos;
                // fixedStep chrom=%s start=%" PRIu64 " step=3 span=3\n
                sscanf(buf + 16, "%s %*6s%" PRIu64, c_chr, &tmp_start_pos);
                if (scores.size() == 0)
                {
                    chr = c_chr;
                    start_pos = tmp_start_pos;
                }
                else if (chr != std::string(c_chr) || start_pos + 3 * scores.size() != tmp_start_pos)
                {
                    // do not update file_range_pos, so we parse it again when calling this function the next time
                    return true;
                }
            }
            else
            {
                double d;
                sscanf(buf, "%lf", &d);
                scores.push_back(d);
            }

            file_range_pos = new_file_range_pos + 1; // skip \n
        }
        return !scores.empty();
    }

    ~wig_reader()
    {
        if (fd >= 0)
            close(fd);

        if (munmap(file_mem, file_size))
            printf("Error munmap %d\n", errno);
    }
};
