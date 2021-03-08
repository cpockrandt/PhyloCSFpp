#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>

#include <stdio.h>
#include <stdlib.h>

#include "src/newick.hpp"
#include "src/fit.hpp"
#include "src/models.hpp"
#include "phylocsf++.cpp"

void newick_get_species(newick_node * n, std::vector<std::string> & v)
{
    if (newick_is_leaf(n))
        v.push_back(n->label);
    else
    {
        newick_get_species(n->left, v);
        newick_get_species(n->right, v);
    }
}

std::string random_seq(const size_t len, const bool with_gaps = false)
{
    std::string s;
    s.reserve(len);
    for (size_t i = 0; i < len; ++i)
    {
        const uint8_t r = rand() % (5 + with_gaps);
        switch (r)
        {
            case 0: s.push_back('A'); break;
            case 1: s.push_back('C'); break;
            case 2: s.push_back('G'); break;
            case 3: s.push_back('T'); break;
            case 4: s.push_back('N'); break;
            case 5: s.push_back('-'); break;
            default: s.push_back('-'); break;
        }
    }
    return s;
}

int main(int argc, char ** argv)
{
    auto seed = 1607119489; // time(NULL); // 1607113429; //
    std::cout << "SEED: " << seed << '\n';
    srand(seed);

    std::string dataset = "100vertebrates";
    auto model = models.find(dataset)->second;
    // get species:
    std::vector<std::string> species;
    {
        newick_node * tree_node = new newick_node;
        newick_parse(*model.tree, tree_node);
        newick_get_species(tree_node, species);
        newick_free(tree_node);
    }
//    species.resize(5);

    while (true)
    {
        std::unordered_set<std::string> species_aln;
        const uint16_t nbr_species_aln = rand()%((species.size() + 1) - 2) + 2; // min 2, max species.size()
        while (species_aln.size() < nbr_species_aln)
        {
            const uint16_t rand_id = rand() % species.size();
            if (species_aln.count(species[rand_id]) == 0)
                species_aln.insert(species[rand_id]);
        }

        std::unordered_set<std::string> species_selection = species_aln;
        const uint16_t nbr_species_selection = rand() % ((species.size() + 1) - nbr_species_aln) + nbr_species_aln;
        while (species_selection.size() < nbr_species_selection)
        {
            const uint16_t rand_id = rand() % species.size();
            if (species_selection.count(species[rand_id]) == 0)
                species_selection.insert(species[rand_id]);
        }

        std::string species_selection_str = "";
        for (auto & s : species_selection)
            species_selection_str += s + ",";
        species_selection_str.pop_back();

//        std::cout << nbr_species_aln << " ... " << nbr_species_selection << '\n';
//    for (auto & s : species_aln)
//        std::cout << s << '\n';

        // create random alignment:
        uint32_t aln_len = rand() % ((20 + 1) - 1) + 1;
        FILE *fptr = fopen("/tmp/test","w");
        if(fptr == NULL)
        {
            printf("Assignment writing to file error!\n");
            exit(1);
        }

        for (auto & s : species_aln)
        {
            std::string seq = random_seq(aln_len, false);
//        std::cout << ">" << s << '\n';
//        std::cout << seq << '\n';
            fprintf(fptr,">%s\n%s\n", s.c_str(), seq.c_str());
        }
        fclose(fptr);


        // OVERWRITE
//        {
//            species_selection_str = "Zebrafish,Cape_golden_mole,Horse,Zebra_finch";
//            FILE *fptr = fopen("/tmp/test","w");
//            fprintf(fptr,
//                ">Zebrafish\n"
//                "AGATGNGGNNAATCTC\n"
//                ">Zebra_finch\n"
//                "TTGNATCTGACNAGGG\n"
//                ">Horse\n"
//                "CGGANCGTANANATAA\n"
//                ">Cape_golden_mole\n"
//                "NANNTTGGAANGGCGC\n");
//            fclose(fptr);
//        }
        // OVERWRITE END



//        std::cout << "phylocsf\tbls score\tanc score\n";

        char path[4096];
        std::string phylocsf_string = "/home/chris/dev-uni/PhyloCSF/PhyloCSF.Linux.x86_64 --bls --ancComp"
                                      " --strategy=omega"
                                      " --species=" + species_selection_str +
                                      " ~/dev-uni/PhyloCSF_vm/PhyloCSF_Parameters/" + dataset +
                                      " /tmp/test";
        std::cout << phylocsf_string << '\n';
        FILE *fp = popen(phylocsf_string.c_str(), "r");
        if (fp == NULL)
        {
            printf("Failed to run command\n");
            exit(1);
        }

        double orig_phylo, orig_bls, orig_anc;
        while (fgets(path, sizeof(path), fp) != NULL)
        {
            printf("orig: %s", path);
            char delim[] = "\t";

            char *ptr = strtok(path, delim);
            orig_phylo = strtod(ptr, &ptr);
            ptr = strtok(NULL, delim);
            orig_bls = strtod(ptr, &ptr);
            ptr = strtok(NULL, delim);
            orig_anc = strtod(ptr, &ptr);
        }
        pclose(fp);

        char* species_c = (char*)species_selection_str.c_str();
        auto new_results = run("/tmp/test", "100vertebrates", species_c, algorithm_t::OMEGA);
        double new_phylo = std::get<0>(new_results);
        double new_bls   = std::get<1>(new_results);
        double new_anc   = std::get<2>(new_results);
        printf("new : %f\t%f\t%f\n", new_phylo, new_bls, new_anc);
        printf("\n");

        double diff_phylo = fabs(new_phylo - orig_phylo);
        double diff_bls = fabs(new_bls - orig_bls);
        double diff_anc = fabs(new_anc - orig_anc);
        if (diff_phylo > 0.00001 || diff_bls > 0.00001 || diff_anc > 0.00001)
        {
            printf("%f\n", diff_phylo);
            printf("%f\n", diff_bls);
            printf("%f\n", diff_anc);
            exit(1);
        }
    }
    return 0;
}

