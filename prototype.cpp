#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>

#include "src/newick.hpp"
#include "src/fit.hpp"
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

std::ostream& operator<<(std::ostream & os, const std::tuple<double, double, double> & t)
{
    std::cout << '(' << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ')';
    return os;
}

#include "src/alignment_reader.hpp"
#include "src/ecm.hpp"
#include "src/instance.hpp"
#include "src/fixed_lik.hpp"
#include "src/additional_scores.hpp"
#include "src/omega.hpp"
#include "src/maf_parser.hpp"

std::ostream& operator<<(std::ostream & os, const gsl_matrix& m)
{
    for (uint8_t i = 0; i < m.size1; ++i)
    {
        for (uint8_t j = 0; j < m.size2; ++j)
        {
            os << gsl_matrix_get(&m, i, j) << '\t';
        }
        os << '\n';
    }
    return os;
}

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
    newick_node * phylo_tree = NULL;

    std::vector<newick_elem> phylo_array;

    std::unordered_set<std::string> selected_species;

    void load_model(const char * dataset_or_path, const bool only_phylo_tree)
    {
        auto dataset = models.find(dataset_or_path);
        // load model from included dataset
        if (dataset != models.end())
        {
            if (!only_phylo_tree)
            {
                c_model.load(dataset->second, 1);
                nc_model.load(dataset->second, 0);
            }

            phylo_tree = new newick_node;
            phylo_tree->parent = NULL;

            newick_parse(*dataset->second.tree, phylo_tree);
        }
        // load model from files
        else // TODO: check and validate files
        {
            if (!only_phylo_tree)
            {
                c_model.open((std::string(dataset_or_path) + "_coding.ECM").c_str());
                nc_model.open((std::string(dataset_or_path) + "_noncoding.ECM").c_str());
            }
            phylo_tree = newick_open((std::string(dataset_or_path) + ".nh").c_str());
        }
        // print error message
//        else
//        {
//            printf("ERROR: The model \"23flies\" does not exist. Please choose one from the following set or give a path to your model:\n");
//            for (const auto & m : models)
//                printf("\t- %s\n", m.first.c_str());
//            exit(-1);
//        }
        assert(phylo_tree->branch_length == 0.0);

        // reduce tree when --species is passed
        // std::cout << newick_print(root) << '\n';
        if (selected_species.size() > 0)
        {
            std::unordered_set<std::string> missing_species = selected_species; // merge into one set with bool
            newick_check_missing_species(phylo_tree, missing_species);
            if (missing_species.size() > 0)
            {
                printf("ERROR: The following selected species are missing in the phylogenetic tree (TODO: output filename):\n");
                for (const std::string & sp : missing_species)
                    printf("\t- %s\n", sp.c_str());
                exit(-1);
            }

            newick_reduce(phylo_tree, selected_species);
            // std::cout << newick_print(tree) << '\n';
            assert(phylo_tree->branch_length == 0.0);
        }

        // get array representation of tree
        newick_flatten(phylo_tree, phylo_array);
        //    for (const auto & elem : newick_flattened)
        //        std::cout << elem << '\n';
    }

//    Data() = default;
//    Data(const Data&) = default;
//    Data(Data&&) = default;
//    Data& operator=(const Data&) = default;
//    Data& operator=(Data&&) = default;
    ~Data()
    {
        if (phylo_tree != NULL)
            newick_free(phylo_tree);
    }
};

void run(char aln_path[], char model_str[], char selected_species_str[], algorithm_t algo)
{
    Data data;
    char delim[] = ",";

    char *ptr = strtok(selected_species_str, delim);
    while(ptr != NULL)
    {
        data.selected_species.emplace(ptr);
        ptr = strtok(NULL, delim);
    }

    data.load_model(model_str, algo == algorithm_t::OMEGA); // TODO: check whether selected species exist!

    // prepare alignment
    alignment_t alignment(data.phylo_array);
    // read alignment
    read_alignment(aln_path, alignment);
//     for (uint16_t i = 0; i < alignment.seqs.size(); ++i)
//        std::cout << alignment.ids[i] << '\t' << alignment.seqs[i] << '\t' << print_peptide(alignment.peptides[i]) << '\n';

    if (algo == algorithm_t::OMEGA)
    {
        auto & inst = data.c_instance;
        // let new_instance ?(kappa:2.5):2.5 ?(omega=1.0) ?(sigma=1.0) ?(tree_scale=1.0) tree_shape =
        {
            // let q_settings = Array.concat [ [| kappa; omega; sigma |]; (Array.make 9 1.0) ]
            inst.q_settings = {2.5 /*kappa*/, 1.0 /*omega*/, 1.0/*sigma*/, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
            // let tree_settings = [| tree_scale |]
            inst.tree_settings = 1.0;

            // let p14n = make_p14n tree_shape
            {
                inst.p14n.tree_shape = data.phylo_array;
                inst.p14n.q_domains = std::vector<domain>(12, domain::NonNeg);
                inst.p14n.tree_domains = std::vector<domain>(1, domain::Pos);
                // skipped these, because they only build expressions. they are evaluated in instantiate_tree and instantiate_qs (in Ocaml). That's where we will start with them in PhyloCSF++
                // q_p14ns = [| q_p14n |];
                // q_scale_p14ns = [| q_scale |];
                // tree_p14n = make_tree_p14n tree_shape;
            }

            // instantiate_tree
            inst.instantiate_tree(inst.tree_settings); // TODO: tree_p14n = Array.init (T.size tree_shape - 1) (fun br -> Mul (Var 0,Val (T.branch tree_shape br)))'

            // instantiate_qs
            compute_q_p14ns_and_q_scale_p14ns_omega(inst);
            inst.instantiate_qs();
            // PM.P14n.instantiate p14n ~q_settings:q_settings ~tree_settings:tree_settings

            //let update_f3x4 inst leaves = // REVIEW: update_f3x4 seems to be correct!
            gsl_matrix * counts = gsl_matrix_alloc(3, 4);
            for (uint8_t i = 0; i < 12; ++i)
                counts->data[i] = 1.0;

            for (const auto peptide : alignment.peptides)
            {
                for (uint32_t i = 0; i < peptide.size(); ++i)
                {
                    if (peptide[i] != 64)
                    {
                        uint8_t i1, i2, i3;
                        from_amino_acid_id(peptide[i], i1, i2, i3);
//                        std::cout << "X: " << (unsigned)peptide[i]
//                                << ' ' << (unsigned)i1
//                                << ' ' << (unsigned)i2
//                                << ' ' << (unsigned)i3
//                        << '\n';

                        gsl_matrix_set(counts, 0, i1, gsl_matrix_get(counts, 0, i1) + 1); // counts[0][i1]++
                        gsl_matrix_set(counts, 1, i2, gsl_matrix_get(counts, 1, i2) + 1); // counts[1][i2]++
                        gsl_matrix_set(counts, 2, i3, gsl_matrix_get(counts, 2, i3) + 1); // counts[2][i3]++
                    }
                }
            }

            std::cout << *counts << '\n';

            for (uint8_t i = 0; i < 3; ++i)
            {
                for (uint8_t j = 0; j < 3; ++j)
                {
                    inst.q_settings[3 + (3 * i + j)] = gsl_matrix_get(counts, i, j) / gsl_matrix_get(counts, i, 3);
                }
            }

            //    PM.P14n.update ~q_settings:qs inst
            gsl_matrix_free(counts);

            exit(1);
        }

        const double phylocsf_score = 0.0;
        const double anchestral_score = NAN;
        const double bls_score = 0.0;
        printf("%f\t%f\t%f\n", phylocsf_score, bls_score, anchestral_score);
    }
    else
    {
        // initialize model
        PhyloCSFModel_make(data.c_instance, data.c_model, data.phylo_array);
        PhyloCSFModel_make(data.nc_instance, data.nc_model, data.phylo_array);

        // std::cout << "Coding instance:\n" << instance_coding;
        // std::cout << "NonCoding instance:\n" << instance_noncoding;

        double lpr_c, lpr_nc, elpr_anc_c, elpr_anc_nc;

        if (algo == algorithm_t::MLE)
        {
            max_lik_lpr_leaves(data.c_instance, alignment, 1.0, lpr_c, elpr_anc_c);
            max_lik_lpr_leaves(data.nc_instance, alignment, 1.0, lpr_nc, elpr_anc_nc);
        }
        else // if (algo == algorithm_t::FIXED)
        {
            lpr_leaves(data.c_instance, alignment, 1.0, lpr_c, elpr_anc_c);
            lpr_leaves(data.nc_instance, alignment, 1.0, lpr_nc, elpr_anc_nc);
        }

        const double phylocsf_score = 10.0 * (lpr_c - lpr_nc) / log(10.0);
        const double anchestral_score = 10.0 * (elpr_anc_c - elpr_anc_nc) / log(10.0);
        const double bls_score = compute_bls_score(data.phylo_tree, alignment);
        printf("%f\t%f\t%f\n", phylocsf_score, bls_score, anchestral_score);
    }
}

int main(int argc, char ** argv)
{
//    parse_maf("/home/chris/dev-uni/PhyloCSF++/test/test.maf");
//    exit(1);

//    Options opt;

//    char* aln_path = argv[1];
//    char* model_str = argv[2];
//    char* selected_species_str = argv[3];

//    char aln_path[] = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny3.fa";
//    char model_str[] = "23flies";
//    char selected_species_str[] = "";

//    char aln_path[] = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/ALDH2.exon5.fa";
//    char model_str[] = "100vertebrates";
//    char selected_species_str[] = "Dog,Cow,Horse,Human,Mouse,Rat";

//    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny3.fa", "23flies", "", algorithm_t::OMEGA);
    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA.fa", "12flies", "", algorithm_t::FIXED);
    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA.fa", "12flies", "", algorithm_t::MLE);

    // run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny3.fa", "23flies", "", algorithm_t::OMEGA);
    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/ALDH2.exon5.fa", "100vertebrates", "", algorithm_t::FIXED);
    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/ALDH2.exon5.fa", "100vertebrates", "", algorithm_t::MLE);

    printf("\nSOLL:\n"
           "361.687566\t0.731391\t48.024101\n"
           "297.623468\t0.731391\t48.257442\n"
           "-176.807976\t0.058448\t-33.030802\n"
           "-179.110842\t0.058448\t-32.878297\n");

    return 0;
}

