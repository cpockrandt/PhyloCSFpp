#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct empirical_codon_model
{
    float matrix[64][64]; // TODO: only store triangle and flatten array for cache locality, also don't store diag since it's 0
    float codon_freq[64];

    // makes matrix symmetric and set diagonal to 0
    bool open(const char* path)
    {
        FILE* stream = fopen(path, "r");

        char line[1024];

        // read file
        uint16_t line_id = 1; // the first line would be a zeros, that's why it is skipped
        matrix[0][0] = .0f;

        while (fgets(line, 1024, stream))
        {
            uint16_t field_id = 0;

            const char *line2 = line;
            char *end;
            if (line_id <= 63) // codon matrix
            {
                matrix[line_id][line_id] = .0f; // diagonal is 0
                for (float val = strtof(line2, &end); line2 != end; val = strtof(line2, &end))
                {
                    line2 = end;
                    matrix[field_id][line_id] = val; // symmetric matrix
                    matrix[line_id][field_id] = val;
                    ++field_id;
                }
                assert(field_id == line_id);
            }
            else if (line_id == 65) // codon frequencies
            {
                for (float val = strtof(line2, &end); line2 != end; val = strtof(line2, &end))
                {
                    line2 = end;
                    codon_freq[field_id] = val;
                    ++field_id;
                }
                assert(field_id == 64);
            }

            ++line_id;
        }

        return true; // TODO: error handling if file does not exist
    }

    void print_model() const noexcept
    {
        for (uint16_t i = 0; i < 64; ++i)
        {
            for (uint16_t j = 0; j < 64; ++j)
            {
                printf("%.3f\t", matrix[i][j]);
            }
            printf("\n");
        }

        printf("\n");

        for (uint16_t i = 0; i < 64; ++i)
        {
            printf("%.3f\t", codon_freq[i]);
        }
        printf("\n");
    }
};
