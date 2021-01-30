#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>

#include "newick.hpp"
#include "models.hpp"
#include "parallel_file_reader.hpp"

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

#include "ecm.hpp"
#include "instance.hpp"
#include "fixed_lik.hpp"
#include "additional_scores.hpp"
#include "omega.hpp"

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
    instance_t c_instance;
    instance_t nc_instance;

    void clear()
    {
        c_instance.clear();
        nc_instance.clear();
    }
};

void run_tracks(Data & data, const Model & model, const alignment_t & alignment, std::vector<double> & lpr_per_codon/*, std::vector<double> & bls_per_codon*/)
{
    // initialize model
    PhyloCSFModel_make(data.c_instance, model.c_model, model.phylo_array);
    PhyloCSFModel_make(data.nc_instance, model.nc_model, model.phylo_array);

    double lpr_c, lpr_nc, elpr_anc_c, elpr_anc_nc;

    std::vector<double> nc_lpr_per_codon;
    nc_lpr_per_codon.reserve(alignment.length() / 3);

    // TODO: skip triplets with PhyloPower < 0.1 to speed up computation!
    lpr_leaves(data.c_instance, alignment, 1.0, lpr_c, elpr_anc_c, lpr_per_codon);
    lpr_leaves(data.nc_instance, alignment, 1.0, lpr_nc, elpr_anc_nc, nc_lpr_per_codon);

    for (uint32_t i = 0; i < lpr_per_codon.size(); ++i)
    {
        lpr_per_codon[i] = 10.0 * (lpr_per_codon[i] - nc_lpr_per_codon[i]) / log(10.0);
    }
}

std::tuple<float, float> run(Data & data, const Model & model, const alignment_t & alignment, const algorithm_t algo)
{
    if (algo == algorithm_t::OMEGA)
    {
        auto & inst = data.c_instance;
        // let new_instance ?(kappa:2.5):2.5 ?(omega=1.0) ?(sigma=1.0) ?(tree_scale=1.0) tree_shape =
        {
            // let q_settings = Array.concat [ [| kappa; omega; sigma |]; (Array.make 9 1.0) ]x
            if (inst.q_settings != NULL && inst.q_settings->size != 12)
            {
                gsl_vector_free(inst.q_settings);
                inst.q_settings = gsl_vector_alloc(12);
            }
            else if (inst.q_settings == NULL)
            {
                inst.q_settings = gsl_vector_alloc(12);
            }

            gsl_vector_set(inst.q_settings, 0, 2.5); // kappa
            gsl_vector_set(inst.q_settings, 1, 1.0); // omega
            gsl_vector_set(inst.q_settings, 2, 1.0); // sigma
            for (uint8_t i = 3; i < 12; ++i)
                gsl_vector_set(inst.q_settings, i, 1.0);

            // let tree_settings = [| tree_scale |]
            inst.tree_settings = 1.0;

            // let p14n = make_p14n tree_shape
            {
                inst.p14n.tree_shape = model.phylo_array;
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
        }

//        std::cout << *inst.model.qms[0].q << '\n';
//        std::cout << *inst.model.qms[0].eig.r_l << '\n';
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

                    const double init_kappa = gsl_vector_get(inst.q_settings, 0);
                    max_lik_lpr_leaves(inst, alignment, lpr_H0, elpr_anc, init_kappa, 1.0, 10.0, &minimizer_lpr_leaves_kappa);
//                    std::cout << "lpr after min_kappa: " << lpr_H0 << " ---\n";
                }
//                printf("lpr_H0: %f\n", lpr_H0);
            }


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

                    const double init_kappa = gsl_vector_get(inst.q_settings, 0);
                    max_lik_lpr_leaves(inst, alignment, lpr_H1, elpr_anc, init_kappa, 1.0, 10.0, &minimizer_lpr_leaves_kappa);
//                    std::cout << "lpr after min_kappa: " << lpr_H1 << " ---\n";
//                print(inst);
                }
//                std::cout << "lpr_H1: " << lpr_H1 << '\n';
            }
        }

        const double phylocsf_score = (10.0 * (lpr_H1 - lpr_H0) / log(10.0));
        const double anchestral_score = NAN;
//        printf("%f\t%f\t%f\n", phylocsf_score, bls_score, anchestral_score);
        return std::make_tuple(phylocsf_score, anchestral_score);
    }
    else
    {
        // initialize model
        PhyloCSFModel_make(data.c_instance, model.c_model, model.phylo_array);
        PhyloCSFModel_make(data.nc_instance, model.nc_model, model.phylo_array);

        // std::cout << "Coding instance:\n" << instance_coding;
        // std::cout << "NonCoding instance:\n" << instance_noncoding;

        double lpr_c, lpr_nc, elpr_anc_c, elpr_anc_nc;

        std::vector<double> c_lpr_per_codon, nc_lpr_per_codon;

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

//        printf("%f\t%f\t%f\n", phylocsf_score, bls_score, anchestral_score);
        return std::make_tuple(phylocsf_score, anchestral_score);
    }
}

int print_model_info(const std::string & model_name)
{
    Model model;
    load_model(model, model_name, "", false, 0, "");

    printf("The model %s contains the following species.\n\n", model_name.c_str());
    printf("%35s\t%s\n", "Species name", "Alternative name(s)");
    for (const auto & elem : model.phylo_array)
    {
        if (elem.label != "")
        {
            std::string scientific_names;
            for (const auto & sn : sequence_name_mapping[elem.label])
                scientific_names += sn + " ";
            printf("%35s\t%s\n", elem.label.c_str(), scientific_names.c_str());
        }
    }

    return 0;
}