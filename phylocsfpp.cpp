#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>

#include "src/newick.hpp"
#include "src/fit.hpp"
#include "src/models.hpp"
#include "src/parallel_file_reader.hpp"

//#ifdef OPENMP
//#include <omp.h>
//#endif

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
    os << '\n';
    return os;
}

//std::ostream& operator<<(std::ostream & os, const std::vector<std::vector<double> > & v)
//{
//    for (uint32_t i = 0; i < v.size(); ++i)
//    {
//        for (uint32_t j = 0; j < v[i].size(); ++j)
//        {
//            char buf[10];
//            sprintf (buf, "%.3f ", v[i][j]);
//            os << buf << ' ';
//        }
//        std::cout << '\n';
//    }
//    return os;
//}

template <typename T>
std::ostream& operator<<(std::ostream & os, const std::vector<T> & v)
{
    os << '[';
    for (uint32_t i = 0; i < v.size(); ++i)
    {
        os << v[i] << ' ';
    }
    os << ']';
    return os;
}

std::ostream& operator<<(std::ostream & os, const std::tuple<double, double, double> & t)
{
    os << '(' << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << ')';
    return os;
}

#include "src/ecm.hpp"
#include "src/instance.hpp"
#include "src/fixed_lik.hpp"
#include "src/additional_scores.hpp"
#include "src/omega.hpp"

//std::ostream& operator<<(std::ostream & os, const gsl_matrix& m)
//{
//    for (uint8_t i = 0; i < m.size1; ++i)
//    {
//        for (uint8_t j = 0; j < m.size2; ++j)
//        {
//            os << gsl_matrix_get(&m, i, j) << '\t';
//        }
//        os << '\n';
//    }
//    return os;
//}

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

    void clear()
    {
        c_instance.clear();
        nc_instance.clear();
    }

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

auto run(Data & data, alignment_t & alignment, algorithm_t algo)
{
    if (algo == algorithm_t::OMEGA)
    {
        auto & inst = data.c_instance;
        // let new_instance ?(kappa:2.5):2.5 ?(omega=1.0) ?(sigma=1.0) ?(tree_scale=1.0) tree_shape =
        {
            // let q_settings = Array.concat [ [| kappa; omega; sigma |]; (Array.make 9 1.0) ]
            inst.q_settings = gsl_vector_alloc(12);
            gsl_vector_set(inst.q_settings, 0, 2.5); // kappa
            gsl_vector_set(inst.q_settings, 1, 1.0); // omega
            gsl_vector_set(inst.q_settings, 2, 1.0); // sigma
            for (uint8_t i = 3; i < 12; ++i)
                gsl_vector_set(inst.q_settings, i, 1.0);

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
            inst.model.qms.reserve(inst.p14n.tree_p14n.size() - 1);
            compute_q_p14ns_and_q_scale_p14ns_omega(inst);
//            std::cout << "instance.p14n.q_p14ns:\n" << *inst.p14n.q_p14ns << '\n';
            inst.instantiate_qs();
            // PM.P14n.instantiate p14n ~q_settings:q_settings ~tree_settings:tree_settings
            PhyloModel_make(inst, inst.q_settings, true);
//            print(inst);
//            exit(13);
        }

//        std::cout << *inst.model.qms[0].q << '\n';
//        std::cout << *inst.model.qms[0].eig.r_l << '\n';
//        exit(17);
        {
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

//            std::cout << *counts << '\n';

            for (uint8_t i = 0; i < 3; ++i)
            {
                for (uint8_t j = 0; j < 3; ++j)
                {
                    gsl_vector_set(inst.q_settings, 3 + (3 * i + j), gsl_matrix_get(counts, i, j) / gsl_matrix_get(counts, i, 3));
                }
            }
            gsl_matrix_free(counts);

            //    PM.P14n.update ~q_settings:qs inst
            compute_q_p14ns_and_q_scale_p14ns_omega(inst);

//            std::cout << *inst.q_settings << '\n';

//            std::cout << "----------------------------------\n";
            inst.instantiate_qs();
            PhyloModel_make(inst, NULL, true);

//            print(inst);
//            exit(19);
        }

        // STATUS 2020-12-03 02:40 pms and qms are identical here with Ocaml version!!!! :)
//        std::cout.precision(6);
        // kr_map leaves inst
        double lpr_H0 = 0.0, lpr_H1 = 0.0;
        {
            double elpr_anc = 0.0;
//            f_roh = (lpr_rho rho +. lpr_leaves inst_rho leaves);

            {
                for (uint8_t i = 0; i < 3; ++i) // TODO: 3 iterations
                {
                    const double init_rho = inst.tree_settings;
//                    std::cout << "init rho: " << init_rho << '\n';
                    max_lik_lpr_leaves(inst, alignment, lpr_H0, elpr_anc, init_rho, 0.001, 10.0, &minimizer_lpr_leaves_rho);
//                    std::cout << "lpr after min_rho: " << lpr_H0 << " ---\n";
//            print(inst);
//            exit(24);

                    const double init_kappa = gsl_vector_get(inst.q_settings, 0);
                    max_lik_lpr_leaves(inst, alignment, lpr_H0, elpr_anc, init_kappa, 1.0, 10.0, &minimizer_lpr_leaves_kappa);
//                    std::cout << "lpr after min_kappa: " << lpr_H0 << " ---\n";
//                print(inst);
//                exit(106);
                }
//                printf("lpr_H0: %f\n", lpr_H0);
            }
//            exit(19);


            {
                gsl_vector_set(inst.q_settings, 1, 0.2); // omega_H1
                gsl_vector_set(inst.q_settings, 2, 0.01); // sigma_H1
                compute_q_p14ns_and_q_scale_p14ns_omega(inst);
                // std::cout << *inst.q_settings << '\n';
                inst.instantiate_qs();
                PhyloModel_make(inst, NULL, true);

                for (uint8_t i = 0; i < 3; ++i)
                {
                    const double init_rho = inst.tree_settings;
//                    std::cout << "init rho: " << init_rho << '\n';
                    max_lik_lpr_leaves(inst, alignment, lpr_H1, elpr_anc, init_rho, 0.001, 10.0, &minimizer_lpr_leaves_rho);
//                    std::cout << "lpr after min_rho: " << lpr_H1 << " ---\n";
//            print(inst);
//            exit(24);

                    const double init_kappa = gsl_vector_get(inst.q_settings, 0);
                    max_lik_lpr_leaves(inst, alignment, lpr_H1, elpr_anc, init_kappa, 1.0, 10.0, &minimizer_lpr_leaves_kappa);
//                    std::cout << "lpr after min_kappa: " << lpr_H1 << " ---\n";
//                print(inst);
//                exit(106);
                }
//                std::cout << "lpr_H1: " << lpr_H1 << '\n';
            }
        }

//         exit(1);

        const double phylocsf_score = (10.0 * (lpr_H1 - lpr_H0) / log(10.0));
        const double anchestral_score = NAN;
        std::vector<double> TODO_deleteme;
        const double bls_score = compute_bls_score(data.phylo_tree, alignment, TODO_deleteme);
//        printf("%f\t%f\t%f\n", phylocsf_score, bls_score, anchestral_score);
        return std::make_tuple(phylocsf_score, bls_score, anchestral_score);
    }
    else
    {
        // initialize model
        PhyloCSFModel_make(data.c_instance, data.c_model, data.phylo_array);
        PhyloCSFModel_make(data.nc_instance, data.nc_model, data.phylo_array);

        // std::cout << "Coding instance:\n" << instance_coding;
        // std::cout << "NonCoding instance:\n" << instance_noncoding;

        double lpr_c, lpr_nc, elpr_anc_c, elpr_anc_nc;

        std::vector<double> c_lpr_per_codon, nc_lpr_per_codon, bls_per_codon;

        if (algo == algorithm_t::MLE)
        {
            max_lik_lpr_leaves(data.c_instance, alignment, lpr_c, elpr_anc_c, 1.0, 1e-2, 10.0, &minimizer_lpr_leaves);
            max_lik_lpr_leaves(data.nc_instance, alignment, lpr_nc, elpr_anc_nc, 1.0, 1e-2, 10.0, &minimizer_lpr_leaves);
        }
        else // if (algo == algorithm_t::FIXED)
        {
            lpr_leaves(data.c_instance, alignment, 1.0, lpr_c, elpr_anc_c, c_lpr_per_codon);
            lpr_leaves(data.nc_instance, alignment, 1.0, lpr_nc, elpr_anc_nc, nc_lpr_per_codon);
        }

        const double phylocsf_score = 10.0 * (lpr_c - lpr_nc) / log(10.0);
        const double anchestral_score = 10.0 * (elpr_anc_c - elpr_anc_nc) / log(10.0);
        const double bls_score = compute_bls_score(data.phylo_tree, alignment, bls_per_codon);

//        for (uint32_t xx = 0; xx < c_lpr_per_codon.size(); ++xx)
//        {
//            const double _score = 10.0 * (c_lpr_per_codon[xx] - nc_lpr_per_codon[xx]) / log(10.0);
//            printf("Codon %d\t%f\t%f\n", xx, _score, bls_per_codon[xx]);
//        }

//        printf("%f\t%f\t%f\n", phylocsf_score, bls_score, anchestral_score);
        return std::make_tuple(phylocsf_score, bls_score, anchestral_score);
    }
}

struct comp_result
{
    double fixed_score, fixed_anc, mle_score, mle_anc, omega_score, bls;

    comp_result(double fixed_score, double fixed_anc, double mle_score, double mle_anc, double omega_score, double bls)
    {
        this->fixed_score = fixed_score;
        this->fixed_anc = fixed_anc;
        this->mle_score = mle_score;
        this->mle_anc = mle_anc;
        this->omega_score = omega_score;
        this->bls = bls;
    }

    void print() const noexcept
    {
        printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", fixed_score, fixed_anc, mle_score, mle_anc, omega_score, bls);
    }
};

int main(int /*argc*/, char ** /*argv*/)
{
    unsigned threads = 1;
    unsigned jobs = 1;
    char model_str[] = "100vertebrates";
    char selected_species_str[] = ""; // "Dog,Cow,Horse,Human,Mouse,Rat";
    char aln_path[] = "/home/chris/Downloads/chr22.head.maf";
//    char aln_path[] = "/home/chris/dev-uni/PhyloCSF++/ALDH2.exon5.maf";
    char scores_path[] = "/home/chris/dev-uni/PhyloCSF++/chr22.head.orig.results.formatted";

    Data data_omega, data_fixed_mle;
    char delim[] = ",";

    char *ptr = strtok(selected_species_str, delim);
    while (ptr != NULL)
    {
        data_omega.selected_species.emplace(ptr);
        data_fixed_mle.selected_species.emplace(ptr);
        ptr = strtok(NULL, delim);
    }

    data_omega.load_model(model_str, true); // TODO: check whether selected species exist!
    data_fixed_mle.load_model(model_str, false); // TODO: check whether selected species exist!

    std::unordered_map<std::string, uint16_t> fastaid_to_alnid;
    // we use the order of data.phylo_array
    {
        FILE *file = fopen("/home/chris/dev-uni/PhyloCSF++/commonNames_assemblies.txt", "r");
        if (file != NULL)
        {
            char line[BUFSIZ];
            char c1[BUFSIZ];
            char c2[BUFSIZ];
            while (fgets(line, sizeof line, file) != NULL)
            {
                sscanf(line, "%s\t%s", c1, c2); // PhyloCSF identifier => maf identifier

                bool found = false;
                for (uint16_t i = 0; i < data_omega.phylo_array.size(); ++i)
                {
                    if (strcmp(data_omega.phylo_array[i].label.c_str(), c1) == 0)
                    {
                        fastaid_to_alnid.emplace(c2, i); // maf identifier => id in vector for ids and seqs
                        fastaid_to_alnid.emplace(c1, i); // maf identifier => id in vector for ids and seqs
                        found = true;
                        break;
                    }
                }
                if (!found)
                    printf("ERROR: %20s\t%10s\tmapping missing!\n", c1, c2);
            }
        }
        fclose(file);

        // TODO: check whether there are mappings missing (not sure) or superfluous mappings (also not sure)
    }

    // prepare alignment
    std::vector<alignment_t> alignments;
    for (unsigned i = 0; i < threads; ++i)
    {
        alignments.emplace_back(data_omega.phylo_array); // TODO: remove id's outside of aln struct because they are read-only
    }

    // open correct values from orig PhyloCSF++
    // awk 'BEGIN{ OFS="\t"; print "ALN-ID", "FIXED", "FIXED-ANC", "MLE", "MLE-ANC", "OMEGA", "BLS" } ($0 != "" && !($0 ~ /^\/tmp/)) { if ($3 ~ /^\/tmp/) { $3 = "nan" } if ($2 == "Fixed:") { bls[$1] = $4 } anc[$1][$2] = $5; f[$1][$2] = $3 } END { for (id in f) print id, f[id]["Fixed:"], anc[id]["Fixed:"], f[id]["MLE:"], anc[id]["MLE:"], f[id]["Omega:"], bls[id] }' chr22.head.orig.results > chr22.head.orig.results.formatted
    std::vector<comp_result> results;
    {
        FILE *file = fopen(scores_path, "r");
        if (file != NULL)
        {
            char line[BUFSIZ];
            size_t id;
            double fixed_score, fixed_anc, mle_score, mle_anc, omega_score, bls;
            fgets(line, sizeof line, file); // skip header line
            while (fgets(line, sizeof line, file) != NULL)
            {
                sscanf(line, "%lu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &id, &fixed_score, &fixed_anc, &mle_score, &mle_anc, &omega_score, &bls);
                assert(id == results.size());
                results.emplace_back(fixed_score, fixed_anc, mle_score, mle_anc, omega_score, bls);
            }
        }
        fclose(file);
    }
//    printf("%ld\n%lf\n%lf\n%lf\n%lf\n%lf\n%lf\n\n", 0,  results[0].fixed_score, results[0].fixed_anc, results[0].mle_score, results[0].mle_anc, results[0].omega_score, results[0].bls);

    parallel_maf_reader maf_rd(aln_path, jobs, &fastaid_to_alnid);

//    #pragma omp parallel for num_threads(threads) default(none) shared(jobs, alignments, maf_rd, data, mode)
    for (unsigned job_id = 0; job_id < jobs; ++job_id) // TODO: split it more parts than threads
    {
        auto & aln = alignments[0/*omp_get_thread_num()*/];
        size_t a_id = 0;
        // TODO: verify whether file_range_pos and _end are thread-safe (cache lines)? and merge both arrays for cache locality

        printf("ALN-ID\tFIXED\tFIXED-ANC\tMLE\tMLE-ANC\tOMEGA\tBLS\n");

        maf_rd.get_next_alignment(aln, job_id); // skip 1st (which is big) alignment
//        maf_rd.get_next_alignment(aln, job_id);

        while (a_id < 2 /*results.size()*/)
        {
            auto & r = results[0/*a_id*/];

            std::tuple<double, double, double> results_fixed, results_mle, results_omega;

//            {
//                results_fixed = run(data_fixed_mle, aln, algorithm_t::FIXED);
//                data_fixed_mle.clear();
//                const double score = fabs(std::get<0>(results_fixed) - r.fixed_score);
//                const double anc = fabs(std::get<2>(results_fixed) - r.fixed_anc);
//                const double bls = fabs(std::get<1>(results_fixed) - r.bls);
//                if (score >= 0.000001 || anc >= 0.000001 || bls >= 0.000001)
//                    printf("%ld\t%s:\t%f\t%f\t%f\t*****\n", a_id, "Fixed", score, anc, bls);
//            }
//
//            {
//                results_mle = run(data_fixed_mle, aln, algorithm_t::MLE);
//                data_fixed_mle.clear();
//                const double score = fabs(std::get<0>(results_mle) - r.mle_score);
//                const double anc = fabs(std::get<2>(results_mle) - r.mle_anc);
//                if (score >= 0.000001 || anc >= 0.000001)
//                    printf("%ld\t%s:\t%f\t%f\t--------\t*****\n", a_id, "MLE", score, anc);
//            }

            {
                results_omega = run(data_fixed_mle, aln, algorithm_t::OMEGA);
                data_omega.clear();
                const double score = fabs(std::get<0>(results_omega) - r.omega_score);
                if (score >= 0.000001)
                    printf("%ld\t%s:\t%f\t--------\t--------\t*****\n", a_id, "Omega", score);
            }

            printf("%ld\t%f\t%f\t%f\t%f\t%f\t%f\n",
                       a_id,
                       std::get<0>(results_fixed), std::get<2>(results_fixed),
                       std::get<0>(results_mle), std::get<2>(results_mle),
                       std::get<0>(results_omega), std::get<1>(results_fixed)
            );
            fflush(stdout);

//            for (auto & seq : aln.seqs)
//                seq = "";
            ++a_id;
        }
    }

    return 0;
}

//    char selected_species_str[] = "Dog,Cow,Horse,Human,Mouse,Rat";
//    auto new_results = run("/home/chris/dev-uni/chr4_100vert_alignment.fa", "100vertebrates", "", algorithm_t::FIXED);
//    auto new_results = run("/tmp/test", "100vertebrates", "Dog,Cow,Horse,Human,Mouse,Rat", algorithm_t::MLE);
//    char aln_path[] = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny3.fa";
//    char model_str[] = "23flies";
//    "Dog,Cow,Horse,Human,Mouse,Rat";
//    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny3.fa", "23flies", "", algorithm_t::OMEGA);
//    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA.fa", "12flies", "", algorithm_t::MLE);
//    run("/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/ALDH2.exon5.fa", "100vertebrates", "", algorithm_t::MLE);
//    printf("\nSOLL:\n"
//           "361.687566\t0.731391\t48.024101\n"
//           "297.623468\t0.731391\t48.257442\n"
//           "-176.807976\t0.058448\t-33.030802\n"   // FIXED
//           "-179.110842\t0.058448\t-32.878297\n"); // MLE
