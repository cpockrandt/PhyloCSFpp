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

    void *file_mem = NULL; // "fileArea"
    size_t *file_range_pos = NULL;
    size_t *file_range_end = NULL;
    size_t *buffer_pos = NULL;

    size_t buf_size;
    char **buf = NULL; // one buffer / "localArea" for each thread

public:

    size_t memcpy_cnt = 0;

    void init(char *file_path, const unsigned threads)
    {
        fd = open(file_path, O_RDONLY);
        if (fd < 0)
        {
            printf("Error: Cannot open %s\n", file_path);
            exit(-1);
        }

        file_size = lseek(fd, 0, SEEK_END);
        if (file_size == (off_t) - 1)
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
        buf_size = page_size; // TODO: is this clever? implemented this way, to reduce edge cases for now

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
        buffer_pos = (size_t*) calloc(sizeof(size_t), this->threads);
        buf = (char**) calloc(sizeof(char*), this->threads);

        size_t file_range_from = 0;
        for (unsigned i = 0; i < this->threads; ++i)
        {
            buf[i] = (char*) calloc(sizeof(char), buf_size + 1);
            buf[buf_size] = 0; // always make buf a zero-terminated strings

            size_t range_size = page_size * (pages_per_thread + (last_thread_with_extra_page > i));
            size_t file_range_to = file_range_from + range_size;

            if (i == this->threads - 1)
            {
                file_range_to = file_size;
                range_size = file_range_to - file_range_from;
            }

            file_range_pos[i] = file_range_from;
            file_range_end[i] = file_range_to;

            char *local_buf = buf[i];

            printf("Thread %2d starting from %10ld (incl.) to %10ld (excl.). Size: %10ld\n", i, file_range_from, file_range_to, range_size);

            // go to beginning of alignment (i.e., go to first page with beginning of a new alignment block and set offset_in_page accordingly)
            // - alignment could start at the very beginning of the current page (or next page). note that except the very first thread, the first page
            //   has to be skipped and cannot be processed (only used for finding the beginning of the next alignment starting from the next page)
            // - there could be no new alignment at all
            // - alignment could start somewhere in the middle of the page
            // - multiple pages could be necessary to be skipped

            {
                // in case the alignment starts at the very first position of the starting page, we need to check whether the character before is a newline
                // but don't do that for the first thread, since the start position is 0 and would trigger an underflow.
                if (i > 0)
                    --file_range_pos[i];
                memcpy(local_buf, file_mem + file_range_pos[i], buf_size); ++memcpy_cnt; // TODO: last page could be smaller than buf_size

//                printf("yyy --- %ld\n%s\n", file_range_pos[i], local_buf);

                if (i == 0 && local_buf[0] == 'a' && local_buf[1] == ' ')
                {
                    file_range_pos[i] += 0; // no offset
                }
                else
                {
                    char *alignment_block_beginning = NULL;
                    while ( (alignment_block_beginning = strstr(local_buf, "\na ")) == NULL &&
                            (file_range_pos[i] += buf_size - 2) < file_range_to // we are searching for a pattern of length 3. it could lie
                    ) // alignment block starts with "s " preceded by a newline
                    {
                        // TODO: this page could already be the last page and have less than buf_size chars!
                        memcpy(local_buf, file_mem + file_range_pos[i], buf_size); ++memcpy_cnt;
//                        printf("%s", local_buf);
//                        offset_in_page[i] += buf_size;
                    }

                    // TODO adding -1 and (buf_size - 2) when accessing memory mapped file could lead to slow-down, align read access!

                    if (alignment_block_beginning != NULL)
                    {
                        buffer_pos[i] = alignment_block_beginning - local_buf + 1; // + 1 for newline
//                        file_range_pos[i] += alignment_block_beginning - local_buf + 1; // + 1 for newline
                    }
                }

//                memcpy(local_buf, file_mem + file_range_pos[i], buf_size); ++memcpy_cnt;
//                printf("%.15s\n", local_buf);
            }

            file_range_from = file_range_to;
        }
//        exit(13);
    }

    void skip(const unsigned thread_id)
    {
        size_t & in_buffer_pos = buffer_pos[thread_id];
        char *local_buf = buf[thread_id];

        char *alignment_block_beginning = NULL;
        while ( (alignment_block_beginning = strstr(local_buf + in_buffer_pos, "\n")) == NULL &&
                (file_range_pos[thread_id] + buf_size) < (size_t)file_size // we are searching for a pattern of length 3. it could lie
                ) // alignment block starts with "s " preceded by a newline
        {
            in_buffer_pos = 0;
            file_range_pos[thread_id] += buf_size;
            size_t buf_size2 = buf_size;
            if (file_range_pos[thread_id] + buf_size > (size_t)file_size)
                buf_size2 = file_size - file_range_pos[thread_id];

            memcpy(local_buf, file_mem + file_range_pos[thread_id], buf_size2); ++memcpy_cnt;
            local_buf[buf_size2] = 0;
        }

        if (alignment_block_beginning == NULL && file_range_pos[thread_id] + buf_size >= (size_t)file_size)
        {
            assert(thread_id == this->threads - 1);
            file_range_pos[thread_id] = file_size;
            return;
        }

        in_buffer_pos += alignment_block_beginning - (local_buf + in_buffer_pos) + 1; // +1 to skip \n

        if (in_buffer_pos == buf_size)
        {
            in_buffer_pos = 0;
            file_range_pos[thread_id] += buf_size; // TODO: could be the last (partial) page
            memcpy(local_buf, file_mem + file_range_pos[thread_id], buf_size); ++memcpy_cnt;
        }
    }

    // this function is only called on lines starting with "s "
    std::string get_line(const unsigned thread_id)
    {
        size_t & in_buffer_pos = buffer_pos[thread_id];

        std::string res = "";

        char *local_buf = buf[thread_id];

        size_t buf_size2 = buf_size; // TODO: is this correct? can we already be on the last partial page when calling this function (again)?

//        printf("\n\n\n%s\n\n\n", local_buf);
        char *alignment_block_beginning = NULL;
        while ( (alignment_block_beginning = strstr(local_buf + in_buffer_pos, "\n")) == NULL &&
                (file_range_pos[thread_id] + buf_size) < (size_t)file_size // we are searching for a pattern of length 3. it could lie
                ) // alignment block starts with "s " preceded by a newline
        {
            res += (local_buf + in_buffer_pos);
            in_buffer_pos = 0;
            file_range_pos[thread_id] += buf_size;

            buf_size2 = buf_size;
            if (file_range_pos[thread_id] + buf_size > (size_t)file_size)
            {
                buf_size2 = file_size - file_range_pos[thread_id];
                local_buf[buf_size2] = 0;
            }

            memcpy(local_buf, file_mem + file_range_pos[thread_id], buf_size); ++memcpy_cnt;
        }

        size_t infix_size;
        if (alignment_block_beginning)
        {
            infix_size = alignment_block_beginning - (local_buf + in_buffer_pos);
        }
        else
        {
            // this case should only happen if we are on the very last page
            assert(thread_id == this->threads - 1);
            assert(file_range_pos[thread_id] + page_size >= file_size);
            infix_size = buf_size2 - in_buffer_pos;
            file_range_pos[thread_id] = file_size;
        }

        char *buf_prefix = (char*) malloc(sizeof(char) * (infix_size + 1));
        memcpy(buf_prefix, local_buf + in_buffer_pos, infix_size); // ++memcpy_cnt;
        buf_prefix[infix_size] = 0;
        res += buf_prefix;
        free(buf_prefix);

        in_buffer_pos = in_buffer_pos + infix_size + 1;

        return res;
    }

    char get_char(const unsigned thread_id)
    {
        size_t & in_buffer_pos = buffer_pos[thread_id];

        assert(in_buffer_pos < 2*buf_size);

        if (in_buffer_pos >= buf_size)
        {
            file_range_pos[thread_id] += buf_size;
            assert(file_range_pos[thread_id] < (size_t)file_size);
            in_buffer_pos -= buf_size;
            memcpy(buf[thread_id], file_mem + file_range_pos[thread_id], buf_size); ++memcpy_cnt;
        }
        assert(in_buffer_pos < buf_size);

        return buf[thread_id][in_buffer_pos];
    }

    // in practice aln.ids will already be set and only
    // if only \n is at the end of the range for a thread, don't process it
    // if "\na" is at the end of the range for a thread, process it
    // TODO: check whether this works!
    bool get_next_alignment(alignment_t & aln, const unsigned thread_id)
    {
        if (file_range_pos[thread_id] >= file_range_end[thread_id])
            return false;




//        char *local_buf = buf[thread_id];
//        memcpy(local_buf, file_mem + 456536063, buf_size);
//        printf("XXX: %s\n", local_buf);
//        memcpy(local_buf, file_mem + 456540159, buf_size);
//        printf("XXX: %s\n", local_buf);
//        exit(1);





        // points to the beginning of "a score=....\n"

        size_t & in_buffer_pos = buffer_pos[thread_id];
//        in_buffer_pos = 0;

        skip(thread_id);

        char line_type;
        while (file_range_pos[thread_id] < (size_t)file_size && (line_type = get_char(thread_id)) != 'a')
        {
//            if (file_range_pos[thread_id] == 608705766 && in_buffer_pos >= 4074)
//                printf("%ld\t%ld\n", file_range_pos[thread_id], in_buffer_pos);
            if (line_type == 's')
            {
                // "s ID int int char int SEQ\n"
                const size_t old_file_range_pos = file_range_pos[thread_id];
                const size_t old_in_buffer_pos = in_buffer_pos;
                std::string line = get_line(thread_id);
//                printf("%s\n", line.c_str());
//                if (thread_id == 3 && file_range_pos[thread_id] > 608701000)
                    printf("%ld\t%4ld\t%s\n", old_file_range_pos, old_in_buffer_pos, line.c_str());
            }
            else
            {
//                std::string line = get_line(thread_id);
//                printf("%s\n", line.c_str());

//                if (thread_id == 3 && file_range_pos[thread_id] > 608701000 /*&& in_buffer_pos == 4074*/)
//                    printf("Skip: %ld\t%4ld\n", file_range_pos[thread_id], in_buffer_pos);
                skip(thread_id);
            }
        }

//        if (file_range_pos[thread_id] < (size_t)file_size)
//        {
//            file_range_pos[thread_id] += in_buffer_pos;
//            size_t buf_size2 = buf_size;
//            if (file_range_pos[thread_id] + buf_size > (size_t)file_size)
//                buf_size2 = file_size - file_range_pos[thread_id];
//            memcpy(buf[thread_id], file_mem + file_range_pos[thread_id], buf_size2); ++memcpy_cnt; // TODO: avoid this by making pos_in_buffer a member of the class
//        }

        return true;
    }

    ~parallel_maf_reader()
    {
        if (fd >= 0)
            close(fd);

        for (unsigned i = 0; i < this->threads; ++i)
            free(buf[i]);

        if (munmap(file_mem, file_size))
            printf("Error munmap %d\n", errno);

        free(buf);
        free(file_range_pos);
        free(file_range_end);
    }
};

int main()
{
    parallel_maf_reader maf_rd;
    maf_rd.init("/home/chris/Downloads/chr22.head.maf", 4);

    alignment_t aln;
    for (unsigned i = 2; i < 3; ++i)
    {
//        aln.seqs.clear();
// TODO: verify whether file_range_pos and _end are thread-safe (cache lines)? and merge both arrays for cache locality
        while (maf_rd.get_next_alignment(aln, i))
        {
//            printf("\n\n---------------------\n\n\n");
//            break;
//            printf("a\n");
//            for (auto & s : aln.seqs)
//            {
//                printf("%s\n", s.c_str());
//            }
        }
    }

    printf("memcpy count: %ld\n", maf_rd.memcpy_cnt);

    return 0;
}

//static int clean(
//        int retVal,     // value to return.
//        char *err,      // error/NULL, allows optional %s for strerror(errno).
//        int fd,         // fd/-1 to close.
//        void *filMem,   // memory/NULL to unmap.
//        off_t sz,       // size if unmapping.
//        void *locMem    // memory/NULL to free.
//) {
//    if (err)     printf (err, strerror(errno));
//    if (locMem)  free(locMem);
//    if (filMem)  munmap(filMem, sz);
//    if (fd >= 0) close(fd);
//    return retVal;
//}
//
//int main(void)
//{
//    const int fd = open("/home/chris/Downloads/chr22.head.maf", O_RDONLY);
//    if (fd < 0)
//        return clean(-1, "Can't open: %s\n", -1, NULL, 0, NULL);
////    printf("File opened okay, fd = %d.\n", fd);
//
//    const off_t sz = lseek(fd, 0, SEEK_END);
//    if (sz == (off_t) - 1)
//        return clean(-1, "Can't seek: %s\n", fd, NULL, 0, NULL);
//    printf("File size: %ld\n", sz);
//
//    size_t page_size = sysconf(_SC_PAGE_SIZE);
//    printf("Page size: %ld\n", page_size);
//
//    // sz bytes starting at 0
//    size_t from = 4 * page_size; // 0
//    size_t mapsize = 100; // sz
//    size_t cpysize = 20;
//    void *fileArea = mmap(NULL, mapsize, PROT_READ, MAP_SHARED, fd, from);
//    if (!fileArea)
//        return clean(-1, "Can't map: %s\n", fd, NULL, mapsize, NULL);
//    printf("File mapped to address %p.\n", fileArea);
//
//    char *localArea = (char*)calloc(1, sz + 1);
//    if (!localArea)
//        return clean(-1, "Can't allocate\n", fd, fileArea, mapsize, NULL);
//
////    memcpy(localArea, fileArea, mapsize);
////    printf("Data copied to %p, result is [\n%s]\n", localArea, localArea);
//    for (size_t i = 0; i < mapsize; i += cpysize)
//    {
//        memcpy(localArea, fileArea + i, cpysize);
//        printf("%s", localArea, localArea);
//    }
//
//    munmap(fileArea, mapsize);
//    free(localArea);
//    close(fd);
//
//    return 0;
////    return clean(0, NULL, fd, fileArea, sz, localArea);
//}