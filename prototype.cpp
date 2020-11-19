#include <iostream>
#include <vector>
#include <string>
#include <gsl/gsl_matrix.h>

#include "src/alignment_reader.hpp"
#include "src/newick.hpp"
#include "src/ecm.hpp"

#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/translation.h>

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

inline std::string translate(const std::string & s)
{
    const seqan::Dna5String s2(s.c_str());
    seqan::Peptide p2;
    seqan::translate(p2, s2, seqan::SINGLE_FRAME, seqan::Serial());
    seqan::CharString c2 = p2;
    return seqan::toCString(c2);
}

struct instance_t
{
    struct model_t {

        struct q_diag_t {
            gsl_matrix q;

            struct eig_t {
                gsl_matrix r_s;  // S = right eigenvectors (in the columns)
                gsl_matrix r_s2; // S' = left eigenvectors (in the rows)
                gsl_matrix r_l;  // diag(L) = eigenvalues

                gsl_matrix_complex nr_s; // TODO: don't duplicate here! maybe templatize?
                gsl_matrix_complex nr_s2;
                gsl_matrix_complex nr_l;
            } eig;

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
    std::vector<std::string> peptides;
    read_alignment(aln_path.c_str(), ids, seqs);

    for (uint16_t i = 0; i < seqs.size(); ++i)
    {
        peptides.push_back(translate(seqs[i]));
        std::cout << ids[i] << '\t' << seqs[i] << '\t' << peptides[i] << '\n';
    }

    // open newick (and ignore spaces)
    newick_node* root = newick_open(std::string(model_path + ".nh").c_str());
//    std::cout << newick_print(root) << '\n';
    // transform tree into array notation
    std::vector<newick_elem> newick_flattened;
    newick_flatten(root, newick_flattened);
//    for (const auto & elem : newick_flattened)
//        std::cout << elem << '\n';

    // TODO: initialize model
    instance_t instance;

    return 0;
}
