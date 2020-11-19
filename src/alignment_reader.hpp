#include <stdio.h>

#include <string>
#include <vector>

int read_alignment(const char * const file_path, std::vector<std::string>& ids, std::vector<std::string>& seqs)
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

    int16_t i = -1;
    while ((linelen = getline(&linebuf, &linesiz, fptr)) > 0)
    {
        if (linebuf[0] == '>') // TODO: remove space after ">" and  suffix after first space after species name
        {
            ++i;
            ids.emplace_back(std::string(linebuf + 1, linelen - 2)); // remove first char, i.e. '>', and remove last char, i.e., '\n'
//            std::cout << "id: " << ids[i] << '\n';
            seqs.emplace_back();
        } else {
            seqs[i] += std::string(linebuf, linelen - 1); // remove last char, i.e., '\n'
//            std::cout << "seqs: " << seqs[i] << '\n';
        }

        if (i == -1) // first line does not start with ">"
        {
            return -1;
        }
    }
    free(linebuf);
    linebuf = NULL;

    fclose(fptr);
    return 0;
}
