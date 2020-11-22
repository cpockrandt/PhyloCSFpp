#include <stdio.h>

#include <string>
#include <vector>

#include "translation.h"

struct alignment_t
{
    std::vector<std::string> ids;
    std::vector<std::string> seqs;
    std::vector<std::vector<uint8_t> > peptides;
};

//inline std::string translate(const std::string & s)
//{
//    const seqan::Dna5String s2(s.c_str());
//    seqan::Peptide p2;
//    seqan::translate(p2, s2, seqan::SINGLE_FRAME, seqan::Serial());
//    seqan::CharString c2 = p2;
//    return seqan::toCString(c2);
//}

int read_alignment(const char * const file_path, alignment_t & alignment)
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
            alignment.ids.emplace_back(std::string(linebuf + 1, linelen - 2)); // remove first char, i.e. '>', and remove last char, i.e., '\n'
//            std::cout << "id: " << ids[i] << '\n';
            alignment.seqs.emplace_back();
        } else {
            alignment.seqs[i] += std::string(linebuf, linelen - 1); // remove last char, i.e., '\n'
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

    // translate nucleotides
    alignment.peptides.resize(alignment.seqs.size());
    for (uint16_t i = 0; i < alignment.seqs.size(); ++i)
    {
        assert(alignment.seqs[0].size() == alignment.seqs[i].size());

        // alignment.peptides.push_back(translate(alignment.seqs[i]));
        alignment.peptides[i].resize(alignment.seqs[i].size() / 3);
        for (uint64_t aa_pos = 0; aa_pos < alignment.peptides[i].size(); ++aa_pos)
        {
            alignment.peptides[i][aa_pos] = get_amino_acid_id(
                    alignment.seqs[i][3 * aa_pos],
                    alignment.seqs[i][3 * aa_pos + 1],
                    alignment.seqs[i][3 * aa_pos + 2]);
        }
    }

    return 0;
}
