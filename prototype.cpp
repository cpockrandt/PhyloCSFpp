#include <iostream>
#include <vector>
#include <string>
#include <gsl/gsl_matrix.h>

#include "src/alignment_reader.hpp"
#include "src/newick.hpp"
#include "src/ecm.hpp"

// #include <seqan/sequence.h>
//
// using namespace seqan;

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

enum domain : uint8_t { // this is a C++11 feature
    Real,
    Pos,
    Neg,
    NonPos,
    NonNeg,
    Probability
};

struct instance_t
{
    struct model_t {

        struct q_diag_t {
            gsl_matrix q;

            gsl_matrix r_s;
            gsl_matrix r_s2;
            gsl_matrix r_l;

            gsl_matrix_complex nr_s; // don't duplicate here! maybe templatize?
            gsl_matrix_complex nr_s2;
            gsl_matrix_complex nr_l;
//        eig : [`r of eig_r | `nr of eig_nr] = {
//            r_s : Gsl.Matrix.matrix;			(* S = right eigenvectors (in the columns) *)
//            r_s' : Gsl.Matrix.matrix;			(* S' = left eigenvectors (in the rows) *)
//            r_l : Gsl.Vector.vector 			(* diag(L) = eigenvalues *)
//                          --- OR ---
//            nr_s : Gsl.Matrix_complex.matrix;	(* S = right eigenvectors (in the columns) *)
//            nr_s' : Gsl.Matrix_complex.matrix;	(* S' = left eigenvectors (in the rows) *)
//            nr_l : Gsl.Vector_complex.vector;	(* diag(L) = eigenvalues *)
//        }
            gsl_vector  pi;

            bool have_pi; // mutable
//        mutable memoized_to_Pt : (float -> Gsl.Matrix.matrix) option;
            float tol; // mutable
        };

        std::vector<newick_elem> tree;
        std::vector<q_diag_t> qms;
        std::vector<gsl_matrix*> pms;
        std::vector<float> prior; // this used to be of type "float array option", i.e., float array or nothing
    } model;
    struct p14n_t {
//        q_p14ns : q_p14n array = Expr.t array array array
//        q_scale_p14ns : Expr.t array;
        std::vector<domain> q_domains;

        std::vector<newick_elem> tree_shape;
//        tree_p14n : Expr.t array;
        std::vector<domain> tree_domains;
    } p14n;

    std::vector<float>  q_settings;
    std::vector<float>  tree_settings;
};

int main(int argc, char ** argv)
{
//    gsl_matrix *test = gsl_matrix_alloc(64, 64);
//    std::cout << gsl_matrix_get(test, 0, 0) << '\n';

    std::string aln_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Examples/tal-AA-tiny.fa";
    std::string model_path = "/home/chris/dev-uni/PhyloCSF_vm/PhyloCSF_Parameters/12flies";

    empirical_codon_model ecm_coding, ecm_noncoding;
    ecm_coding.open(std::string(model_path + "_coding.ECM").c_str());
    ecm_noncoding.open(std::string(model_path + "_noncoding.ECM").c_str());
//    ecm_coding.print_model();

    // read alignment
    std::vector<std::string> ids;
    std::vector<std::string> seqs;
    read_alignment(aln_path.c_str(), ids, seqs);

//    for (auto const & id : ids)
//        std::cout << id << '\n';
//    for (auto const & s : seqs)
//        std::cout << s << '\n';

    // open newick (and ignore spaces)
    newick_node* root = newick_open(std::string(model_path + ".nh").c_str());
    std::cout << newick_print(root) << '\n';
    // transform tree into array notation
    std::vector<newick_elem> newick_flattened;
    newick_flatten(root, newick_flattened);
    for (const auto & elem : newick_flattened)
        std::cout << elem << '\n';

    // initialize model
    instance_t instance;

    return 0;
}
