#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include <string>
#include <vector>

#include <cassert>

#ifdef OPENMP
#include <omp.h>
#endif

struct alignment_t
{
    std::vector<std::string> ids;
    std::vector<std::string> seqs;
};

class parallel_maf_reader
{
    unsigned threads;
    int fd = -1;
    off_t file_size;
    size_t page_size;

    size_t pages_per_thread;
    size_t last_thread_with_extra_page;

    void *file_mem = NULL;
    size_t *file_range_pos = NULL;
    size_t *file_range_end = NULL;

public:

    void init(char *file_path, const unsigned threads)
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

        file_mem = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
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
            printf("Info: Using %ld threads instead of %d.\n", pages, threads);
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
                    char *aln_begin = strstr((char*) (file_mem + file_range_pos[i]), "\na ");

                    if (aln_begin == NULL)
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

            file_range_from = file_range_to;
        }
    }

    // skip a line
    void skip(const unsigned thread_id)
    {
        char *newline_begin = strstr((char*) (file_mem + file_range_pos[thread_id]), "\n");

        if (newline_begin == NULL)
            file_range_pos[thread_id] = file_size;
        else
            file_range_pos[thread_id] = newline_begin - (char*) file_mem + 1; // skip \n
    }

    // this function is only called on lines starting with "s "
    void get_line(const unsigned thread_id, char ** id, char ** seq)
    {
        assert(file_range_pos[thread_id] < (size_t) file_size);
        assert(*((char*) file_mem + file_range_pos[thread_id]) == 's');

        char *newline_begin = strstr((char*) (file_mem + file_range_pos[thread_id]), "\n");

        size_t new_file_range_pos;
        if (newline_begin == NULL)
            new_file_range_pos = file_size;
        else
            new_file_range_pos = newline_begin - (char*) file_mem;

        // TODO: this malloc and memcpy can be avoided by reimplementing strtok that does not modify the input string
        const size_t len = new_file_range_pos - file_range_pos[thread_id];
        char *res = (char*) malloc(len + 1);
        memcpy(res, file_mem + file_range_pos[thread_id], len + 1);
        res[len] = 0;

//        printf("%s\n", res);

        char *strtok_state;
        char *token = strtok_r(res, " ", &strtok_state);
        uint16_t token_id = 0;
        while (token != NULL)
        {
            if (token_id == 1)
            {
                const size_t token_len = strlen(token);
                *id = (char*) malloc(token_len + 1);
                memcpy(*id, token, token_len);
                (*id)[token_len] = 0;
//                printf("%d. Token: %s (len: %ld)\n", token_id, token, token_len);
            }
            else if (token_id == 6)
            {
                const size_t token_len = strlen(token);
                *seq = (char*) malloc(token_len + 1);
                memcpy(*seq, token, token_len);
                (*seq)[token_len] = 0;
//                printf("%d. Token: %s (len: %ld)\n", token_id, token, token_len);
            }
            token = strtok_r(NULL, " ", &strtok_state);
//            printf("%d. Token: %s\n", token_id, token);
            ++token_id;
        }

        assert(id != NULL && seq != NULL);

        free(res);

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
    bool get_next_alignment(alignment_t & aln, const unsigned thread_id)
    {
        if (file_range_pos[thread_id] >= file_range_end[thread_id])
            return false;

        // points to the beginning of "a score=....\n"
        skip(thread_id);

        char line_type;
        while (file_range_pos[thread_id] < (size_t)file_size && (line_type = get_char(thread_id)) != 'a')
        {
            if (line_type == 's')
            {
                // "s ID int int char int SEQ\n"
                char *id = NULL, *seq = NULL;
                get_line(thread_id, &id, &seq);
                aln.ids.emplace_back(id);
                aln.seqs.emplace_back(seq);
                free(id);
                free(seq);
            }
            else
            {
                skip(thread_id);
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

int main()
{
    parallel_maf_reader maf_rd;
    maf_rd.init("/home/chris/Downloads/chr22.head.maf", 4);

    #pragma omp parallel for num_threads(4)
    for (unsigned thread_id = 0; thread_id < 4; ++thread_id)
    {
        alignment_t aln;
        // TODO: verify whether file_range_pos and _end are thread-safe (cache lines)? and merge both arrays for cache locality
        while (maf_rd.get_next_alignment(aln, thread_id))
        {
            #pragma omp critical
            {
                for (size_t i = 0; i < aln.ids.size(); ++i)
                {
                    printf("%20s\t%s\n", aln.ids[i].c_str(), aln.seqs[i].c_str());
                }
                printf("\n");
            }
            aln.ids.clear();
            aln.seqs.clear();
        }
    }

    return 0;
}
