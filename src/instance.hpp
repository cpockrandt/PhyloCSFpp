#include <vector>

//#include "expr.hpp"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_blas.h>

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
                gsl_vector * r_l;  // diag(L) = eigenvalues

                gsl_matrix_complex * nr_s; // TODO: don't duplicate here! maybe templatize?
                gsl_matrix_complex * nr_s2;
                gsl_vector_complex * nr_l;
            } eig;

            gsl_vector * pi;

            bool have_pi; // mutable
//        mutable memoized_to_Pt : (float -> Gsl.Matrix.matrix) option; // TODO
            double tol; // mutable

            gsl_matrix * _deep_copy_matrix(const gsl_matrix * src)
            {
                if (src == NULL)
                    return NULL;
                gsl_matrix * dest = gsl_matrix_alloc(src->size1, src->size2);
                gsl_matrix_memcpy(dest, src);
                return dest;
            }

            gsl_vector * _deep_copy_vector(const gsl_vector * src)
            {
                if (src == NULL)
                    return NULL;
                gsl_vector * dest = gsl_vector_alloc(src->size);
                gsl_vector_memcpy(dest, src);
                return dest;
            }

            gsl_matrix_complex * _deep_copy_matrix_complex(const gsl_matrix_complex * src)
            {
                if (src == NULL)
                    return NULL;
                gsl_matrix_complex * dest = gsl_matrix_complex_alloc(src->size1, src->size2);
                gsl_matrix_complex_memcpy(dest, src);
                return dest;
            }

            gsl_vector_complex * _deep_copy_vector_complex(const gsl_vector_complex * src)
            {
                if (src == NULL)
                    return NULL;
                gsl_vector_complex * dest = gsl_vector_complex_alloc(src->size);
                gsl_vector_complex_memcpy(dest, src);
                return dest;
            }

            q_diag_t() = default;
            q_diag_t(const q_diag_t&) = default; // copy constructor
            q_diag_t(q_diag_t&&) = default; // move constructor
            q_diag_t& operator=(const q_diag_t& other) // copy assignment
            {
                this->q = _deep_copy_matrix(other.q);

                this->eig.r_l = _deep_copy_vector(other.eig.r_l);
                this->eig.r_s = _deep_copy_matrix(other.eig.r_s);
                this->eig.r_s2 = _deep_copy_matrix(other.eig.r_s2);

                this->eig.nr_l = _deep_copy_vector_complex(other.eig.nr_l);
                this->eig.nr_s = _deep_copy_matrix_complex(other.eig.nr_s);
                this->eig.nr_s2 = _deep_copy_matrix_complex(other.eig.nr_s2);

                this->pi = _deep_copy_vector(other.pi);

                this->have_pi = other.have_pi;
                // mutable memoized_to_Pt
                this->tol = other.tol;
                return *this;
            };
            q_diag_t& operator=(q_diag_t&&) = default; // move assignment
            virtual ~q_diag_t() = default; // destructor
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

        std::vector<newick_elem> tree_shape; // original tree
        std::vector<newick_elem> tree; // tree with evaluated branch lengths (already computed)
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

    instance.compute_tree_p14n(instance.tree_settings.back());

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

    instance.model.qms.resize(instance.p14n.tree.size() - 1);
    instance.model.qms[0].q = gsl_matrix_alloc(64, 64);
    for (uint8_t i = 0; i < 64; ++i)
    {
        for (uint8_t j = 0; j < 64; ++j)
        {
            gsl_matrix_set(instance.model.qms[0].q, i, j, instance.p14n.q[i][j]);
//            printf("%.3f\t", instance.p14n.q[i][j]);
        }
//        printf("\n");
    }

    gsl_vector_complex *l = gsl_vector_complex_alloc(64);
    gsl_matrix_complex *s = gsl_matrix_complex_alloc(64, 64);
    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(64);

    gsl_eigen_nonsymmv(instance.model.qms[0].q, l, s, w);
    gsl_eigen_nonsymmv_free(w);
    // NOTE: it seems that order does not seem to matter! can sort in Ocaml and does not change the result
    // NOTE: even though we sort, some values are neg. instead of pos and vice versa
    gsl_eigen_nonsymmv_sort (l, s, GSL_EIGEN_SORT_ABS_ASC); // TODO: not necessary, just for comparing results with ocaml helpful

//    let p = Gsl.Permut.make n (***)
//    let lu = Gsl.Vectmat.cmat_convert (`CM (Gsl.Matrix_complex.copy m)) (** REVIEW: it seems cmat_convert doesn't do anything according to source: https://github.com/mmottl/gsl-ocaml/blob/master/src/vectmat.ml#L73 *)
//    ignore (Gsl.Linalg.complex_LU_decomp lu p) (** REVIEW: ignore might have a side effect ??? *)
//    let m' = Gsl.Vectmat.cmat_convert (`CM (Gsl.Matrix_complex.create n n)) (** REVIEW: same *)
//    Gsl.Linalg.complex_LU_invert lu p m'

    int signum;
    gsl_permutation * permut = gsl_permutation_alloc(64);
    gsl_matrix_complex *lu = gsl_matrix_complex_alloc(64, 64);
    gsl_matrix_complex_memcpy(lu, s);
    gsl_matrix_complex *s2 = gsl_matrix_complex_alloc(64, 64);
    gsl_linalg_complex_LU_decomp(lu, permut, &signum);
    gsl_linalg_complex_LU_invert(lu, permut, s2 /* inverse */);

    for (uint8_t i = 0; i < 64; ++i)
    {
//        const gsl_complex &x = gsl_vector_complex_get(l, i);
//        printf("(%.3f, %.3f)\t", x.dat[0], x.dat[1]);
        for (uint8_t j = 0; j < 64; ++j)
        {
//            const gsl_complex &x = gsl_matrix_complex_get(s2, i, j);
//            printf("(%.3f, %.3f)\t", x.dat[0], x.dat[1]);
        }
//        printf("\n");
    }
//    printf("\n");

    constexpr double tol = 1e-6;

    // im = 0. || (abs_float im) *. 1000. < (Complex.norm z) || (abs_float re < tol && abs_float im < tol)
    // check whether l can be a real vector
    bool is_l_complex = false;
    for (uint8_t i = 0; i < 64 && !is_l_complex; ++i)
    {
        const gsl_complex &x = gsl_vector_complex_get(l, i);
        if (!(
                x.dat[1] == 0.0 ||
                (fabs(x.dat[1]) * 1000. < gsl_complex_abs(x)) ||
                (fabs(x.dat[0]) < tol && fabs(x.dat[1]) < tol)
            )) // NOTE: fabs() is for doubles, fabsf() is for floats
        {
            is_l_complex = true;
        }
    }

//    std::cout << "is complex? " << is_l_complex << '\n';
    if (is_l_complex)
    {
        instance.model.qms[0].eig.nr_s = s;
        instance.model.qms[0].eig.nr_l = l;
        instance.model.qms[0].eig.nr_s2 = s2;

        instance.model.qms[0].eig.r_s = NULL;
        instance.model.qms[0].eig.r_l = NULL;
        instance.model.qms[0].eig.r_s2 = NULL;
    }
    else
    {
        // transform complex matrices and vector into real ones
        instance.model.qms[0].eig.r_s = gsl_matrix_alloc(64, 64);
        instance.model.qms[0].eig.r_l = gsl_vector_alloc(64);
        instance.model.qms[0].eig.r_s2 = gsl_matrix_alloc(64, 64);

        instance.model.qms[0].eig.nr_s = NULL;
        instance.model.qms[0].eig.nr_l = NULL;
        instance.model.qms[0].eig.nr_s2 = NULL;

        for (uint8_t i = 0; i < 64; ++i)
        {
            const gsl_complex &l_elem = gsl_vector_complex_get(l, i);
            gsl_vector_set(instance.model.qms[0].eig.r_l, i, l_elem.dat[0]);
            for (uint8_t j = 0; j < 64; ++j)
            {
                const gsl_complex &s_elem = gsl_matrix_complex_get(s, i, j);
                const gsl_complex &s2_elem = gsl_matrix_complex_get(s2, i, j);
                gsl_matrix_set(instance.model.qms[0].eig.r_s, i, j, s_elem.dat[0]);
                gsl_matrix_set(instance.model.qms[0].eig.r_s2, i, j, s2_elem.dat[0]);
            }
        }
    }

    // TODO: PhyloModel.make (remember: this is also called in each update (with new tree and old q))
    instance.model.qms[0].tol = tol;
    instance.model.qms[0].have_pi = false;
    instance.model.qms[0].pi = gsl_vector_alloc(64); // unused empty vector? // pi = Gsl.Vector.create (fst (Gsl.Matrix.dims qm))

    // do this only in the initialization step (copy qms[0] to all other (uninitialized) members)
    instance.model.pms.resize(instance.p14n.tree.size() - 1);
    instance.model.qms.resize(instance.p14n.tree.size() - 1);
    for (uint16_t i = 1; i < instance.model.qms.size(); ++i)
    {
        instance.model.qms[i] = instance.model.qms[0]; // deep-copy
//        gsl_vector_scale(instance.model.qms[i].eig.r_l, (double)i);
    }

//    for (uint16_t i = 0; i < instance.model.qms.size(); ++i)
//        std::cout << "TEST: " << gsl_vector_get(instance.model.qms[i].eig.r_l, 0) << '\n';

//    struct q_diag_t {
//        gsl_matrix * q;
//
//        struct eig_t {
//            gsl_matrix * r_s;  // S = right eigenvectors (in the columns)
//            gsl_matrix * r_s2; // S' = left eigenvectors (in the rows)
//            gsl_vector * r_l;  // diag(L) = eigenvalues
//
//            gsl_matrix_complex * nr_s; // TODO: don't duplicate here! maybe templatize?
//            gsl_matrix_complex * nr_s2;
//            gsl_vector_complex * nr_l;
//        } eig;
//
//        gsl_vector * pi;
//
//        bool have_pi; // mutable
////        mutable memoized_to_Pt : (float -> Gsl.Matrix.matrix) option; // TODO
//        double tol; // mutable
//    };

    instance.model.tree = instance.p14n.tree; // tree with evaluated expressions (i.e., multiplied)
    instance.model.prior = codon_freq;
    // pms = Array.init (T.size t - 1) (fun br -> Q.Diag.to_Pt qms.(br) (T.branch t br))
    for (uint16_t i = 0; i < instance.p14n.tree_shape.size() - 1; ++i)
    {
        // to_Pt:
        // TODO: here is memoization happening with q.memoized_to_Pt
        // especially inefficient since for initialization all qms[i] are identical!
        instance_t::model_t::q_diag_t &q = instance.model.qms[i]; // TODO: make sure whether tree[i] matches qms[i] (but should)
        gsl_matrix *&p = instance.model.pms[i];
        const double t = instance.p14n.tree[i].branch_length;

        p = gsl_matrix_alloc(64, 64); // p is the i-th matrix in pms, and will have first gemm_c (variable name from Ocaml), and afterwards sm (substitution matrix of gemm_c)
        // real part
        if (q.eig.nr_l == NULL)
        {
            assert(q.eig.r_l != NULL && q.eig.r_s != NULL && q.eig.r_s2 != NULL);
            assert(q.eig.nr_l == NULL && q.eig.nr_s == NULL && q.eig.nr_s2 == NULL);

            gsl_vector * expLt = gsl_vector_alloc(64);
            gsl_vector_memcpy(expLt, q.eig.r_l); // remove memcopy and insert below!
            for (uint8_t expLt_i = 0; expLt_i < 64; ++expLt_i)
                gsl_vector_set(expLt, expLt_i, gsl_sf_exp(gsl_vector_get(expLt, expLt_i) * t)); // gsl_sf_exp might to error checking! maybe avoid in the end

            // gemm r_s (diagm expLt r_s')
            gsl_matrix * diagm_c = gsl_matrix_alloc(64, 64);
            for (uint8_t diagm_i = 0; diagm_i < 64; ++diagm_i)
            {
                const double expLt_i = gsl_vector_get(expLt, diagm_i);
                for (uint8_t diagm_j = 0; diagm_j < 64; ++diagm_j)
                {
                    gsl_matrix_set(diagm_c, diagm_i, diagm_j, gsl_matrix_get(q.eig.r_s2, diagm_i, diagm_j) * expLt_i);
                }
            }

            gsl_blas_dgemm(CBLAS_TRANSPOSE::CblasNoTrans,
                           CBLAS_TRANSPOSE::CblasNoTrans,
                           1.0 /* alpha */,
                           q.eig.r_s /* A */,
                           diagm_c /* B */,
                           0.0 /* beta */,
                           p /* C result */);

            // results are equal with ocaml until here!!!! :)
//            for (uint8_t tmp_i = 0; tmp_i < 64; ++tmp_i)
//            {
//                for (uint8_t tmp_j = 0; tmp_j < 64; ++tmp_j)
//                {
//                    printf("%.3f\t", gsl_matrix_get(gemm_c, tmp_i, tmp_j));
//                }
//                printf("\n");
//            }
        }
        // non-real-part
        else
        {
            // NOTE: we just create all q's (exept the very first one). the arrays should all have 0s!
            //assert(q.eig.r_l == NULL && q.eig.r_s == NULL && q.eig.r_s2 == NULL);
            //assert(q.eig.nr_l != NULL && q.eig.nr_s != NULL && q.eig.nr_s2 != NULL);

            gsl_vector_complex * expLt = gsl_vector_complex_alloc(64);
            for (uint8_t expLt_i = 0; expLt_i < 64; ++expLt_i)
            {
                const gsl_complex complex_t = gsl_complex_rect(t, 0.0); // TODO: test whether this works!
                // val = exp((q.eig.nr_l[expLt_i] * complex_t))
                const gsl_complex expLt_val = gsl_complex_exp(gsl_complex_mul(gsl_vector_complex_get(q.eig.nr_l, expLt_i), complex_t));
                gsl_vector_complex_set(expLt, expLt_i, expLt_val);
            }

            // OLD: gemm r_s (diagm expLt r_s')
            // NEW: m_of_cm (zgemm nr_s (zdiagm expLt nr_s'))
            gsl_matrix_complex * zdiagm_c = gsl_matrix_complex_alloc(64, 64);
            for (uint8_t diagm_i = 0; diagm_i < 64; ++diagm_i)
            {
                const gsl_complex expLt_i = gsl_vector_complex_get(expLt, diagm_i);
                for (uint8_t diagm_j = 0; diagm_j < 64; ++diagm_j)
                {
                    gsl_matrix_complex_set(zdiagm_c, diagm_i, diagm_j, gsl_complex_mul(gsl_matrix_complex_get(q.eig.nr_s2, diagm_i, diagm_j), expLt_i));
                }
            }

            gsl_matrix_complex * zgemm_c = gsl_matrix_complex_alloc(64, 64);
            gsl_blas_zgemm(CBLAS_TRANSPOSE::CblasNoTrans,
                           CBLAS_TRANSPOSE::CblasNoTrans,
                           GSL_COMPLEX_ONE /* alpha */,
                           q.eig.nr_s /* A */,
                           zdiagm_c /* B */,
                           GSL_COMPLEX_ZERO /* beta */,
                           zgemm_c /* C result */);

            // extract real part of matrix into gemm_c
            for (uint8_t diagm_i = 0; diagm_i < 64; ++diagm_i)
            {
                for (uint8_t diagm_j = 0; diagm_j < 64; ++diagm_j)
                {
                    gsl_matrix_set(p, diagm_i, diagm_j, GSL_REAL(gsl_matrix_complex_get(zgemm_c, diagm_i, diagm_j)));
                }
            }
        }

        // now doing something on gemm_c ... gemm_c is referred to as matrix "sm" (substitution matrix) in the ocaml code
        for (uint8_t index_i = 0; index_i < 64; ++index_i)
        {
            double total = 0.0;
            double smii = 1.0;

            for (uint8_t index_j = 0; index_j < 64; ++index_j)
            {
                double cell_value = gsl_matrix_get(p, index_i, index_j); // sm[i][j]

                total += cell_value;
                if (cell_value < 0.0)
                {
                    if (fabs(cell_value) > q.tol) // fabs is for doubles, fabsf is for floats
                    {
                        printf("CamlPaml.Q.substition_matrix: expm(%.2e*Q)[%d,%d] = %e < 0", t, index_i, index_j, cell_value);
                        exit(1);
                    }
                    else
                    {
                        gsl_matrix_set(p, index_i, index_j, 0.0);
                    }
                }
                if (index_i != index_j)
                {
                    smii -= gsl_matrix_get(p, index_i, index_j);
                }
            }

            if (fabs(total - 1.0) > q.tol)
            {
                printf("CamlPaml.Q.substitution matrix: sum(expm(%.2e*Q)[%d,] = %e > 1", t, index_i, total);
                exit(2);
            }
            assert(smii <= 1.0 && smii > 0.0);

            gsl_matrix_set(p, index_i, index_i, smii);
        }
        // p now is the i-th gemm_c matrix / substitution matrix sm
    }

    // results are equal with ocaml until here!!!! :) (except complex case - didn't occur yet!)
    for (auto & p : instance.model.pms)
    {
        for (uint8_t tmp_i = 0; tmp_i < 64; ++tmp_i)
        {
            for (uint8_t tmp_j = 0; tmp_j < 64; ++tmp_j)
            {
                printf("%.3f\t", gsl_matrix_get(p, tmp_i, tmp_j));
            }
            printf("\n");
        }
        printf("---\n");
    }
}