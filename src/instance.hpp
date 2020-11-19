#include <vector>

//#include "expr.hpp"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>

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
            gsl_matrix * q;

            struct eig_t {
                gsl_matrix * r_s;  // S = right eigenvectors (in the columns)
                gsl_matrix * r_s2; // S' = left eigenvectors (in the rows)
                gsl_matrix * r_l;  // diag(L) = eigenvalues

                gsl_matrix_complex * nr_s; // TODO: don't duplicate here! maybe templatize?
                gsl_matrix_complex * nr_s2;
                gsl_matrix_complex * nr_l;
            } eig;

            gsl_vector * pi;

            bool have_pi; // mutable
//        mutable memoized_to_Pt : (float -> Gsl.Matrix.matrix) option; // TODO
            double tol; // mutable
        };

        std::vector<newick_elem> tree;
        std::vector<q_diag_t> qms; // size = numer of tree node without root (tree.size() - 1)
        std::vector<gsl_matrix*> pms; // size = numer of tree node without root (tree.size() - 1)
        std::vector<double> prior; // this used to be of type "float array option", i.e., float array or nothing
    } model;

    struct p14n_t {
        std::vector<std::vector<double> > q; // already evaluated
        double q_scale; // already evaluated

//        std::vector<std::vector<std::vector<expr_node *> > > q_p14ns; // one dim of vectors can be reduced. there is only one outer-most element in PhylocCSF
//        std::vector<expr_node *> q_scale_p14ns; // one dim of vectors can be reduced. there is only one outer-most element in PhylocCSF
        std::vector<domain> q_domains;

        std::vector<newick_elem> tree_shape;
        std::vector<newick_elem> tree; // already computed
        std::vector<domain> tree_domains;
    } p14n;

    std::vector<double>  q_settings;
    std::vector<double>  tree_settings; // this only stores the tree scale (a single float/double), no need for a vector!

    void compute_q_p14ns_and_q_scale_p14ns( // see instantiate_q and instantiate_qs,
                                            // memorization "from previous branches" is happening in Ocaml. Optimize here?
            const std::vector<std::vector<double> > & _q, // original from ECM model!!!
            const std::vector<double> & variables) noexcept
    {
        // - multiply every element q[i][j] with a variable j
        // - then set diagonal to: q[i][i] = sum_(i <> j) (-q[i][j])
        this->p14n.q = _q;
        this->p14n.q_scale = 0.0;

        for (uint8_t i = 0; i < 64; ++i)
        {
            assert(_q[i][i] == 0.0); // otherwise we need an IF before multiplying and summing up below
            double sum = .0;
            for (uint8_t j = 0; j < 64; ++j)
            {
                this->p14n.q[i][j] *= variables[j];
                sum -= this->p14n.q[i][j];
            }
            this->p14n.q[i][i] = sum; // set diagonal such that all rows and cols sum up to 0
            this->p14n.q_scale -= this->p14n.q[i][i] * variables[i];
        }

        for (uint8_t i = 0; i < 64; ++i)
        {
            for (uint8_t j = 0; j < 64; ++j)
            {
                this->p14n.q[i][j] /= this->p14n.q_scale;
            }
        }
    }

    void compute_tree_p14n(const double factor) noexcept // see instantiate_tree
    {
        // multiply all branch lengths with "factor"
        this->p14n.tree = this->p14n.tree_shape;
        for (newick_elem & elem : this->p14n.tree)
        {
            elem.branch_length *= factor;
        }
    }
};

inline void make(instance_t & instance, const empirical_codon_model & ecm, std::vector<newick_elem> & tree_array)
{
//    instance.p14n.q_p14ns;
//    std::vector<std::vector<expr_node *> > q_p14n(64);
//    for (uint16_t i = 0; i < 64; ++i)
//    {
//        q_p14n[i].resize(64);
//        for (uint16_t j = 0; j < 64; ++j)
//        {
//
//        }
//    }
//    instance.p14n.q_p14ns.push_back(std::move(q_p14n));

//    instance.p14n.q_scale_p14ns;

    instance.p14n.q_domains = {};

    instance.p14n.tree_shape = tree_array;
    instance.p14n.tree = tree_array;

//    instance.p14n.tree_p14n;
//    for (const newick_elem & elem : tree_array)
//    {
//        if (elem.parent_id != -1) // not the root! TODO: maybe remove the root from tree_array at the first place!
//        {
//            instance.p14n.tree_p14n.push_back(
//                    new expr_node(expr_op::MUL, // MUL
//                                  new expr_node(expr_op::VAR, NULL, NULL, 0, .0f), // VAR 0
//                                  new expr_node(expr_op::VAL, NULL, NULL, 0, elem.branch_length), // VAL branch_length
//                                  0,
//                                  .0f
//                    )
//            );
//            std::cout << expr_to_string(instance.p14n.tree_p14n.back()) << '\n';
//        }
//    }

    instance.p14n.tree_domains = {domain::Pos};



    instance.q_settings.assign(ecm.codon_freq, ecm.codon_freq + 64);
    instance.tree_settings = {1.0};

    std::vector<std::vector<double> > matrix;
    std::vector<double> codon_freq;
    codon_freq.resize(64);
    matrix.resize(64);
    for (int i = 0; i < 64; ++i)
    {
        codon_freq[i] = ecm.codon_freq[i];
        matrix[i].resize(64);
        for (int j = 0; j < 64; ++j)
        {
            matrix[i][j] = ecm.matrix[i][j];
        }
    }
    assert(instance.q_settings == codon_freq);
    instance.compute_q_p14ns_and_q_scale_p14ns(matrix, codon_freq); // "evaluate" formulas in q_expr and q_scale_expr

//    std::cout << "Q (after evaluation of expr):\n";
//    for (int i = 0; i < 64; ++i)
//    {
//        for (int j = 0; j < 64; ++j)
//            printf("%.3f\t", instance.p14n.q[i][j]);
//        std::cout << '\n';
//    }

    // TODO: somewhere also Q.of_Q is computed on result. Eigen computations, etc.
    instance.model.qms.resize(instance.p14n.tree.size() - 1);
    instance.model.qms[0].q = gsl_matrix_alloc(64, 64);
    for (uint8_t i = 0; i < 64; ++i)
    {
        for (uint8_t j = 0; j < 64; ++j)
        {
            gsl_matrix_set(instance.model.qms[0].q, i, j, instance.p14n.q[i][j]);
//            gsl_eigen_nonsymmv()
        }
    }
//    std::cout << gsl_matrix_get(test, 0, 0) << '\n';

    // TODO: down in make also include PhyloModel.make that does some (yet unknown) stuff! Make is also called after each update (with new tree and old q)!


}