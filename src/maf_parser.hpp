#pragma once

#include <stdio.h>

#include <string>
#include <vector>

int parse_maf(const char * const file_path, alignment_t & alignment)
{
    char * linebuf;
    size_t linesiz = 0;
    ssize_t linelen = 0;

    FILE *fptr;
    if ((fptr = fopen(file_path, "r")) == NULL)
    {
        printf("Error! opening file");
        return -1;
    }

    while ((linelen = getline(&linebuf, &linesiz, fptr)) > 0)
    {
        assert(linelen > 0); // linebuf cannot be empty (at least has a newline character)
        if (linebuf[0] == '#' || linebuf[0] == 'i' || linebuf[0] == 'q' || isspace(linebuf[0]))
        {
            continue;
        }
        else if (linebuf[0] == 'a')
        {
            printf("\n");
        }
        else if (linebuf[0] == 's')
        {
            printf("%s", linebuf);
        }
    }
    free(linebuf);

    fclose(fptr);

    return 0;
}
