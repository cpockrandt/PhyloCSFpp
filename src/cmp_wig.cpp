#include <stdio.h>
#include <string>

#include <regex>

#include "common.hpp"

// TODO: get rid of regex and use sscanf
const std::regex pieces_regex("fixedStep chrom=([A-Za-z0-9_]+) start=([0-9]+) step=3 span=3");
std::smatch pieces_match;

bool get_next_value(FILE *f, std::string & chr, uint64_t & pos, float & val)
{
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    bool new_header = false;

    while (true)
    {
        do {
            read = getline(&line, &len, f);
        } while (read == 0); // skip empty lines

        if (read == -1)
        {
            free(line);
            return false;
        }

        // remove newline
        if (line[read - 1] == '\n')
            line[read - 1] = 0;

        if (line[0] == 'f') // fixedStep
        {
            const std::string str_line = line;
            if (std::regex_match(str_line, pieces_match, pieces_regex))
            {
                chr = pieces_match[1].str();
                pos = std::stoi(pieces_match[2].str());
//                printf("New info: %s, %ld\n", chr.c_str(), pos);
                new_header = true;
            }
            else
            {
                printf(OUT_ERROR "Regex did not match!\n" OUT_RESET);
                exit(-1);
            }
        }
        else
        {
            val = std::atof(line);
//            printf("New float: %f\n", val);
            if (!new_header)
                pos += 3;

            free(line);
            return true;
        }
    }
}

int main(int argc, char ** argv)
{
    if (argc != 4)
    {
        printf(OUT_ERROR "ERROR: Please pass two pathes as arguments and a threshold!\n" OUT_RESET);
        return -1;
    }

    const char *p1 = argv[1];
    const char *p2 = argv[2];
    const float threshold = atof(argv[3]);

    FILE *f1 = NULL, *f2 = NULL;

    if ((f1 = fopen(p1, "r")) == NULL)
    {
        printf(OUT_ERROR "Cannot open file for reading: %s\n" OUT_RESET, p1);
        return -1;
    }

    if ((f2 = fopen(p2, "r")) == NULL)
    {
        printf(OUT_ERROR "Cannot open file for reading: %s\n" OUT_RESET, p2);
        return -1;
    }

    std::string chr1, chr2;
    uint64_t pos1, pos2;
    float val1, val2;

    bool not_eof1 = get_next_value(f1, chr1, pos1, val1);
    bool not_eof2 = get_next_value(f2, chr2, pos2, val2);

    while (not_eof1 || not_eof2)
    {
        if (chr1 == chr2)
        {
            if (pos1 == pos2)
            {
                if (std::abs(val1 - val2) > threshold)
                    printf("DIFF: %s %ld: %f vs. %f\n", chr1.c_str(), pos1, val1, val2);
                not_eof1 = get_next_value(f1, chr1, pos1, val1);
                not_eof2 = get_next_value(f2, chr2, pos2, val2);
            }
            else
            {
                if (pos1 < pos2)
                {
                    if (val1 != 0.0)
                        printf("MISSING in file %d: %s %ld (exp. val: %f)\n", 2, chr1.c_str(), pos1, val1);
                    not_eof1 = get_next_value(f1, chr1, pos1, val1);
                }
                else
                {
                    if (val2 != 0.0)
                        printf("MISSING in file %d: %s %ld (exp. val: %f)\n", 1, chr1.c_str(), pos2, val2);
                    not_eof2 = get_next_value(f2, chr2, pos2, val2);
                }
            }
        }
        else
        {
            exit(1);
        }
    }

    fclose(f1);
    fclose(f2);

    return 0;
}