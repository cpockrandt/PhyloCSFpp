#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>

#include "src/newick.hpp"
#include "src/models.hpp"

std::ostream& operator<<(std::ostream & os, const newick_elem & e)
{
    os  << e.id << '\t'
        << e.label << '\t'
        << e.parent_id << '\t'
        << '(' << e.child1_id << ", " << e.child2_id << ")\t"
        << e.sibling_id << '\t'
        << e.branch_length;
    return os;
}

std::ostream& operator<<(std::ostream & os, const std::vector<newick_elem> & v)
{
    for (uint32_t i = 0; i < v.size(); ++i)
    {
        os << v[i] << '\n';
    }
    return os;
}

std::ostream& operator<<(std::ostream & os, const std::vector<double> & v)
{
    for (uint32_t i = 0; i < v.size(); ++i)
    {
        char buf[10];
        sprintf (buf, "%f\t", v[i]);
        os << buf;
    }
    std::cout << '\n';
    return os;
}

std::ostream& operator<<(std::ostream & os, const std::vector<std::vector<double> > & v)
{
    for (uint32_t i = 0; i < v.size(); ++i)
    {
        for (uint32_t j = 0; j < v[i].size(); ++j)
        {
            char buf[10];
            sprintf (buf, "%.3f ", v[i][j]);
            os << buf << ' ';
        }
        std::cout << '\n';
    }
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream & os, const std::vector<T> & v)
{
    std::cout << '[';
    for (uint32_t i = 0; i < v.size(); ++i)
    {
        os << v[i] << ' ';
    }
    std::cout << ']' << '\n';
    return os;
}

#include "src/alignment_reader.hpp"
#include "src/ecm.hpp"
#include "src/instance.hpp"
#include "src/fixed_lik.hpp"
#include "src/additional_scores.hpp"

enum algorithm_t
{
    MLE,
    FIXED,
    OMEGA
};

struct Options
{
    algorithm_t algorithm = algorithm_t::FIXED;

    bool remove_gaps = false;
    bool compute_ancestor_score = false;
    bool compute_branch_length_score = false;

    size_t threads = 1; // TODO: default value
};

struct Data
{
    empirical_codon_model c_model;
    empirical_codon_model nc_model;
    instance_t c_instance;
    instance_t nc_instance;
    newick_node * tree = NULL;

    std::vector<newick_elem> newick_flattened;

    std::unordered_set<std::string> selected_species;

    bool load_model(const char * dataset_or_path)
    {
        auto dataset = models.find(dataset_or_path);
        // load model from included dataset
        if (dataset != models.end())
        {
            c_model.load(dataset->second, 1);
            nc_model.load(dataset->second, 0);

            tree = new newick_node;
            tree->parent = NULL;

            newick_parse(*dataset->second.tree, tree);
        }
        // load model from files
        else // TODO: check and validate files
        {
            c_model.open((std::string(dataset_or_path) + "_coding.ECM").c_str());
            nc_model.open((std::string(dataset_or_path) + "_noncoding.ECM").c_str());
            tree = newick_open((std::string(dataset_or_path) + ".nh").c_str());
        }
        // print error message
//        else
//        {
//            printf("ERROR: The model \"23flies\" does not exist. Please choose one from the following set or give a path to your model:\n");
//            for (const auto & m : models)
//                printf("\t- %s\n", m.first.c_str());
//            exit(-1);
//        }
        assert(tree->branch_length == 0.0);

        // reduce tree when --species is passed
        // std::cout << newick_print(root) << '\n';
        if (selected_species.size() > 0)
        {
            std::unordered_set<std::string> missing_species = selected_species; // merge into one set with bool
            newick_check_missing_species(tree, missing_species);
            if (missing_species.size() > 0)
            {
                printf("ERROR: The following selected species are missing in the phylogenetic tree (TODO: output filename):\n");
                for (const std::string & sp : selected_species)
                    printf("\t- %s\n", sp.c_str());
                return -1;
            }

            newick_reduce(tree, selected_species);
            // std::cout << newick_print(tree) << '\n';
            assert(tree->branch_length == 0.0);
        }

        // get array representation of tree
        newick_flatten(tree, newick_flattened);
        //    for (const auto & elem : newick_flattened)
        //        std::cout << elem << '\n';
    }

    ~Data()
    {
        if (tree != NULL)
            newick_free(tree);
    }
};

//void run(const Options & opt, Data & data)
//{
//
//}

int main(int argc, char ** argv)
{
    Data data;
    Options opt;

//    std::string aln_path = argv[1];
//    std::string model_str = argv[2];
//    char* selected_species_str = argv[3];
    char delim[] = ",";

//    std::string aln_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/ALDH2.exon5.fa";
//    std::string model_str = "100vertebrates";
//    char selected_species_str[] = "Human,Mouse,Rat,Cow,Horse,Dog";

    std::string aln_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny3.fa";
    std::string model_str = "23flies";
    char selected_species_str[] = "dmel,dsim,dsec";

    char *ptr = strtok(selected_species_str, delim);
    while(ptr != NULL)
    {
        data.selected_species.emplace(ptr);
        ptr = strtok(NULL, delim);
    }

    data.load_model(model_str.c_str());

    // read alignment
    alignment_t alignment;
    const uint16_t nbr_leaves = (data.newick_flattened.size() + 1)/2;
    alignment.ids.resize(nbr_leaves);
    alignment.seqs.resize(nbr_leaves, "");
    alignment.peptides.resize(nbr_leaves, {});
    for (uint16_t i = 0; i < nbr_leaves; ++i)
    {
        alignment.ids[i] = data.newick_flattened[i].label;
    }

    read_alignment(aln_path.c_str(), alignment);
//    for (uint16_t i = 0; i < alignment.seqs.size(); ++i)
//        std::cout << alignment.ids[i] << '\t' << alignment.seqs[i] << '\t' << print_peptide(alignment.peptides[i]) << '\n';

    // initialize model
    const uint16_t nbr_leaves_in_tree = newick_count_leaves(data.tree);
    PhyloCSFModel_make(data.c_instance, data.c_model, data.newick_flattened);
    PhyloCSFModel_make(data.nc_instance, data.nc_model, data.newick_flattened);

//    std::cout << "Coding instance:\n" << instance_coding;
//    std::cout << "NonCoding instance:\n" << instance_noncoding;

    double lpr_c, lpr_nc, elpr_anc_c, elpr_anc_nc;
    lpr_leaves(data.c_instance, alignment, 1.0, nbr_leaves_in_tree, lpr_c, elpr_anc_c);
    lpr_leaves(data.nc_instance, alignment, 1.0, nbr_leaves_in_tree, lpr_nc, elpr_anc_nc);

    const double decibans_score = 10.0 * (lpr_c - lpr_nc) / log(10.0);
    const double anchestral_score = 10.0 * (elpr_anc_c - elpr_anc_nc) / log(10.0);
    const double bls_score = compute_bls_score(data.tree, alignment);
    printf("%f\t%f\t%f", decibans_score, bls_score, anchestral_score);

    return 0;
}
