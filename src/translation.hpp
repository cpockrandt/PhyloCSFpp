#include <cassert>

char translation_table[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',

    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',

    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',

    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    '*', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',

    '-' // Marginalize "code"
};

uint8_t get_dna_id(const char n)
{
    switch (n)
    {
        case 'A': case 'a':
            return 0;
        case 'C': case 'c':
            return 1;
        case 'G': case 'g':
            return 2;
        case 'T': case 't':
            return 3;
        case '.':
        case '-':
        case 'N': // TODO: some amino acids could be called even with a single N!
            return 4;
        default:
            // TODO: error handling?
            assert(false);
            printf("AHHHH: #%c#\n", n);
            exit(37);
            return 99;
    }
}

uint8_t get_amino_acid_id(const char n1, const char n2, const char n3)
{
    const uint8_t id1 = get_dna_id(n1);
    const uint8_t id2 = get_dna_id(n2);
    const uint8_t id3 = get_dna_id(n3);
    if (id1 == 4 || id2 == 4 || id3 == 4)
        return 64; // Marginalize "code"
    return 16 * id1 + 4 * id2 + id3;
}

void from_amino_acid_id(const uint8_t codon_id, uint8_t & i1, uint8_t & i2, uint8_t & i3)
{
    assert(0 <= codon_id && codon_id < 64);
    i1 = codon_id / 16;
    i2 = (codon_id - 16 * i1) / 4;
    i3 = (codon_id - 16 * i1 - 4 * i2);
    assert(codon_id == get_amino_acid_id(i1, i2, i3))
}

char get_amino_acid(const uint8_t amino_acid_id)
{
    assert(0 <= amino_acid_id && amino_acid_id <= 64);
    return translation_table[amino_acid_id];
}

std::string print_peptide(const std::vector<uint8_t> & p)
{
    std::string res;
    res.resize(p.size());
    for (uint32_t aa_pos = 0; aa_pos < p.size(); ++aa_pos)
    {
        res[aa_pos] = translation_table[p[aa_pos]];
    }
    return res;
}