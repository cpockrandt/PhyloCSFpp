//#include <sys/stat.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include <string>
#include <vector>

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
    size_t *offset_in_page = NULL;
//    size_t *offset_pages = NULL;

    size_t buf_size;
    char **buf = NULL; // one buffer / "localArea" for each thread

public:
    void init(char *file_path, unsigned threads)
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

//        offset_pages = (size_t *) calloc(sizeof(size_t), this->threads);
//        file_start_pos = (void **) calloc(sizeof(size_t*), this->threads);
        offset_in_page = (size_t*) calloc(sizeof(size_t), this->threads);
        buf = (char**) calloc(sizeof(char*), this->threads);

        // NOTE: (except for the very first thread), every threads starts one page before the first one it can be processed.
        // this is necessary in case an alignment block starts at the very beginning of a (by chance). then we need to check whether
        // the last character on the previous page is a newline!
        size_t file_range_from = 0;
        for (unsigned i = 0; i < this->threads; ++i)
        {
            buf[i] = (char*) calloc(sizeof(char), buf_size);

            size_t range_size = page_size * (pages_per_thread + (last_thread_with_extra_page > i));
            size_t file_range_to = file_range_from + range_size;

            if (i == this->threads - 1)
            {
                file_range_to = file_size;
                range_size = file_range_to - file_range_from;
            }

            if (i > 0)
            {
                file_range_from -= page_size;
                range_size += page_size;
            }

            offset_in_page[i] = file_range_from;

            char *local_buf = buf[i];
//            void *local_file_pos = file_start_pos[i];

            printf("Thread %2d starting from %10ld (incl.) to %10ld (excl.). Size: %10ld\n", i, file_range_from, file_range_to, range_size);
//            local_file_pos = mmap(NULL, range_size, PROT_READ, MAP_SHARED, fd, file_range_from);
//            if (!local_file_pos)
//            {
//                printf("Error: Cannot map %s\n", file_path);
//                exit(-1); // TODO: free memory, etc
//            }

            // go to beginning of alignment (i.e., go to first page with beginning of a new alignment block and set offset_in_page accordingly)
            // - alignment could start at the very beginning of the current page (or next page). note that except the very first thread, the first page
            //   has to be skipped and cannot be processed (only used for finding the beginning of the next alignment starting from the next page)
            // - there could be no new alignment at all
            // - alignment could start somewhere in the middle of the page
            // - multiple pages could be necessary to be skipped
//            exit(0);
//            if (i == 1)
            {
                size_t j = 0;
                char *alignment_block_beginning = NULL;
                do {
//                    memcpy(local_buf, local_file_pos + j, buf_size);
                    memcpy(local_buf, file_mem + offset_in_page[i] + j, buf_size);
                    printf("%s", local_buf);
                    j += buf_size;
//                    printf("%ld\n", j);
//                    exit(13);
                } while (j < range_size && (alignment_block_beginning = strstr(local_buf, "\na ")) == NULL); // alignment block starts with "s " preceded by a newline
                // TODO: do not go past last page!

                offset_in_page[i] = alignment_block_beginning - local_buf;

                printf("\n\n%ld\n%s\n", offset_in_page[i], local_buf + offset_in_page[i]);
                printf("\n%c\n", local_buf[0]);
            }
            exit(0);

            file_range_from = file_range_to;
        }
    }

    // in practice aln.ids will already be set and only
    void get_next_alignment(alignment_t & aln, unsigned thread_id)
    {
//        char *local_buf = buf[thread_id];
//        void *local_file_pos = file_start_pos[thread_id];

//        memcpy(localArea, fileArea + i, cpysize);
//        printf("%s", localArea, localArea);
    }

    ~parallel_maf_reader()
    {
        if (fd >= 0)
            close(fd);

        for (unsigned i = 0; i < this->threads; ++i)
        {
            free(buf[i]);
//            if (!munmap((file_start_pos + i), pages_per_thread + (last_thread_with_extra_page > i)))
//                printf("Error munmap %d\n", errno);
        }

        if (!munmap(file_mem, file_size))
            printf("Error munmap %d\n", errno);

        free(buf);
        free(offset_in_page);
//        free(file_start_pos);
    }
};

int main()
{
    parallel_maf_reader maf_rd;
    maf_rd.init("/home/chris/Downloads/chr22.head.maf", 4);
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