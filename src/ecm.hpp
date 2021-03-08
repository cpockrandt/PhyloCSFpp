#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct my_model {
    std::string * tree;
    double * coding_matrix; // 63*64/2
    double * noncoding_matrix; // 63*64/2
    double * coding_codon_freq; // 64
    double * noncoding_codon_freq; // 64
};

struct empirical_codon_model
{
    double matrix[64][64]; // TODO: only store triangle and flatten array for cache locality, also don't store diag since it's 0
    double codon_freq[64];

    // makes matrix symmetric and set diagonal to 0
    bool open(const char* path)
    {
        FILE* stream = fopen(path, "r");

        char line[4096];

        // read file
        uint16_t line_id = 1; // the first line would be a zeros, that's why it is skipped
        matrix[0][0] = .0f;

        while (fgets(line, 4096, stream))
        {
            uint16_t field_id = 0;

            const char *line2 = line;
            char *end;
            if (line_id <= 63) // codon matrix
            {
                matrix[line_id][line_id] = .0f; // diagonal is 0
                for (double val = strtod(line2, &end); line2 != end; val = strtod(line2, &end))
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
                for (double val = strtod(line2, &end); line2 != end; val = strtod(line2, &end))
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

    // makes matrix symmetric and set diagonal to 0
    void load(const my_model & model, const bool coding_mode)
    {
        double *m_matrix;
        double *m_codon_freq;
        if (coding_mode)
        {
            m_matrix = model.coding_matrix;
            m_codon_freq = model.coding_codon_freq;
        } else {
            m_matrix = model.noncoding_matrix;
            m_codon_freq = model.noncoding_codon_freq;
        }

        uint32_t i = 0, j = 1;
        matrix[0][0] = 0.0;
        for (uint32_t a = 0; a < 63*64/2; ++a)
        {
            if (i == j)
            {
                matrix[i][j] = 0.0;
                matrix[j][i] = 0.0;
                i = 0;
                ++j;
            }
            matrix[i][j] = m_matrix[a];
            matrix[j][i] = m_matrix[a];
            ++i;
        }
        matrix[63][63] = 0.0;

        memcpy(codon_freq, m_codon_freq, 64 * sizeof(m_codon_freq));
    }
};
