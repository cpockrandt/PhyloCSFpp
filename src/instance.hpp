#pragma once

#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_blas.h>

enum domain {
    Pos,
    NonNeg
};

std::ostream& operator<<(std::ostream & os, const gsl_vector v)
{
    os << '[';
    for (uint32_t i = 0; i < v.size; ++i)
    {
        os << gsl_vector_get(&v, i) << ' ';
    }
    os << ']';
    return os;
}

std::ostream& operator<<(std::ostream & os, const gsl_vector_complex v)
{
    os << '[';
    for (uint32_t i = 0; i < v.size; ++i)
    {
        auto c = gsl_vector_complex_get(&v, i);
        os << c.dat[0] << ',' << c.dat[1] << ' ';
    }
    os << ']';
    return os;
}

std::ostream& operator<<(std::ostream & os, const gsl_matrix m)
{
    for (uint32_t i = 0; i < m.size1; ++i)
    {
        for (uint32_t j = 0; j < m.size2; ++j)
        {
//            os << gsl_matrix_get(&m, i, j) << ' ';
            char buf[10];
            sprintf(buf, "%f ", gsl_matrix_get(&m, i, j));
            os << buf;
        }
        os << '\n';
    }
    return os;
}

std::ostream& operator<<(std::ostream & os, const gsl_matrix_complex m)
{
    for (uint32_t i = 0; i < m.size1; ++i)
    {
        for (uint32_t j = 0; j < m.size2; ++j)
        {
            auto c = gsl_matrix_complex_get(&m, i, j);
//            os << gsl_matrix_get(&m, i, j) << ' ';
            char buf[20];
            sprintf(buf, "%f,%f ", c.dat[0], c.dat[1]);
            os << buf;
        }
        os << '\n';
    }
    return os;
}

std::ostream& operator<<(std::ostream & os, const domain d)
{
    switch (d)
    {
        case domain::Pos:
            os << "Pos";
            break;
        case domain::NonNeg:
            os << "NonNeg";
            break;
        default:
            os << "UNKNOWN!";
    }
    return os;
}

inline bool check_real(const gsl_complex& c, const double tol)
{
    return
        c.dat[1] == 0.0 ||
        (fabs(c.dat[1]) * 1000. < gsl_complex_abs(c)) ||
        (fabs(c.dat[0]) < tol && fabs(c.dat[1]) < tol);
}

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

struct instance_t
{
    struct model_t {

        struct q_diag_t {
            gsl_matrix * q = NULL;

            struct eig_t {
                gsl_matrix * r_s = NULL;  // S = right eigenvectors (in the columns)
                gsl_matrix * r_s2 = NULL; // S' = left eigenvectors (in the rows)
                gsl_vector * r_l = NULL;  // diag(L) = eigenvalues

                gsl_matrix_complex * nr_s = NULL; // TODO: don't duplicate here! maybe templatize?
                gsl_matrix_complex * nr_s2 = NULL;
                gsl_vector_complex * nr_l = NULL;

                ~eig_t()
                {
                    gsl_matrix_free(r_s);
                    gsl_matrix_free(r_s2);
                    gsl_vector_free(r_l);
                    gsl_matrix_complex_free(nr_s);
                    gsl_matrix_complex_free(nr_s2);
                    gsl_vector_complex_free(nr_l);
                }
            } eig;

            gsl_vector * pi = NULL;

            bool have_pi; // mutable
//        mutable memoized_to_Pt : (float -> Gsl.Matrix.matrix) option; // TODO
            double tol; // mutable

            q_diag_t() = default;
            q_diag_t(const q_diag_t& other)// copy constructor
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
            }
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
            virtual ~q_diag_t()
            {
                gsl_matrix_free(q);
                gsl_vector_free(pi);
            }; // destructor
        };

        std::vector<newick_elem> tree;
        std::vector<q_diag_t> qms; // size = numer of tree node without root (tree.size() - 1)
        std::vector<gsl_matrix*> pms; // size = numer of tree node without root (tree.size() - 1)
        gsl_vector * prior = NULL; // this used to be of type "float array option", i.e., float array or nothing

        model_t() = default;
        model_t(const model_t& other) = delete; // copy constructor: TODO: deep-copy of prior
//        {
//            this->tree = other.tree;
//            this->qms = other.qms;
//            this->pms.resize(other.pms.size());
//            for (uint16_t i = 0; i < other.pms.size(); ++i)
//                this->pms[i] = _deep_copy_matrix(other.pms[i]);
//            this->prior = other.prior;
//        };
        model_t(model_t&&) = default; // move constructor
        model_t& operator=(const model_t& other) = delete; // copy assignment TODO: deep-copy of prior
//        {
//            this->tree = other.tree;
//            this->qms = other.qms;
//            this->pms.resize(other.pms.size());
//            for (uint16_t i = 0; i < other.pms.size(); ++i)
//                this->pms[i] = _deep_copy_matrix(other.pms[i]);
//            this->prior = other.prior;
//            return *this;
//        };
        model_t& operator=(model_t&&) = default; // move assignment
        ~model_t()
        {
            if (prior != NULL)
                gsl_vector_free(prior);
            for (auto p : pms)
                gsl_matrix_free(p);
        };
    } model;

    struct p14n_t {
        gsl_matrix * q_p14ns = NULL; // already evaluated
        double q_scale_p14ns; // already evaluated

        std::vector<domain> q_domains;

        std::vector<newick_elem> tree_shape; // original tree
        std::vector<newick_elem> tree_p14n; // tree with evaluated branch lengths (already computed)
        std::vector<domain> tree_domains;

        ~p14n_t()
        {
            if (q_p14ns != NULL)
                gsl_matrix_free(q_p14ns);
        };
    } p14n;

    gsl_vector * q_settings = NULL;
    double tree_settings; // this only stores the tree scale (a single float/double), no need for a vector!

    ~instance_t()
    {
        if (q_settings != NULL)
            gsl_vector_free(q_settings);
    }

    void instantiate_tree(const double factor) noexcept // see instantiate_tree
    {
        // multiply all branch lengths with "factor"
        this->p14n.tree_p14n = this->p14n.tree_shape;
        for (newick_elem & elem : this->p14n.tree_p14n)
        {
            elem.branch_length *= factor;
        }
    }

    void instantiate_qs() noexcept
    {
        if (this->model.qms.size() == 0)
            this->model.qms.resize(1);

        this->model.qms[0].q = gsl_matrix_alloc(64, 64);
        gsl_matrix_memcpy(this->model.qms[0].q, this->p14n.q_p14ns);
//        std::cout << *this->model.qms[0].q << '\n';

        // -- of_Q
        // ---- Gsl.Eigen.nonsymmv
        gsl_vector_complex *l = gsl_vector_complex_alloc(64);
        gsl_matrix_complex *s = gsl_matrix_complex_alloc(64, 64);
        gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(64);

        gsl_eigen_nonsymmv(this->model.qms[0].q, l, s, w);
        gsl_eigen_nonsymmv_free(w);
        // NOTE: it seems that order does not seem to matter! can sort in Ocaml and does not change the result
        // NOTE: even though we sort, some values are neg. instead of pos and vice versa
        // TODO: not necessary, just for comparing results with ocaml helpful
        gsl_eigen_nonsymmv_sort(l, s, GSL_EIGEN_SORT_ABS_ASC);

        // ---- zinvm
        int signum;
        gsl_permutation * permut = gsl_permutation_alloc(64);
        gsl_matrix_complex *lu = gsl_matrix_complex_alloc(64, 64);
        gsl_matrix_complex_memcpy(lu, s);
        gsl_matrix_complex *s2 = gsl_matrix_complex_alloc(64, 64);
        gsl_linalg_complex_LU_decomp(lu, permut, &signum);
        gsl_linalg_complex_LU_invert(lu, permut, s2 /* inverse */);
        gsl_permutation_free(permut);
        gsl_matrix_complex_free(lu);

        // im = 0. || (abs_float im) *. 1000. < (Complex.norm z) || (abs_float re < tol && abs_float im < tol)
        // check whether l can be a real vector
        bool is_l_complex = false;
        constexpr double tol = 1e-6;
        for (uint8_t i = 0; i < 64 && !is_l_complex; ++i)
        {
            const gsl_complex &x = gsl_vector_complex_get(l, i);
            if (!check_real(x, tol)) // NOTE: fabs() is for doubles, fabsf() is for floats
            {
                is_l_complex = true;
            }
        }

        // -- still in of_Q
        if (is_l_complex)
        {
            this->model.qms[0].eig.nr_s = s;
            this->model.qms[0].eig.nr_l = l;
            this->model.qms[0].eig.nr_s2 = s2;

            this->model.qms[0].eig.r_s = NULL;
            this->model.qms[0].eig.r_l = NULL;
            this->model.qms[0].eig.r_s2 = NULL;
        }
        else
        {
            // transform complex matrices and vector into real ones
            this->model.qms[0].eig.r_s = gsl_matrix_alloc(64, 64);
            this->model.qms[0].eig.r_l = gsl_vector_alloc(64);
            this->model.qms[0].eig.r_s2 = gsl_matrix_alloc(64, 64);

            this->model.qms[0].eig.nr_s = NULL;
            this->model.qms[0].eig.nr_l = NULL;
            this->model.qms[0].eig.nr_s2 = NULL;

            for (uint8_t i = 0; i < 64; ++i)
            {
                const gsl_complex &l_elem = gsl_vector_complex_get(l, i);
                gsl_vector_set(this->model.qms[0].eig.r_l, i, l_elem.dat[0]);
                for (uint8_t j = 0; j < 64; ++j)
                {
                    const gsl_complex &s_elem = gsl_matrix_complex_get(s, i, j);
                    const gsl_complex &s2_elem = gsl_matrix_complex_get(s2, i, j);
                    gsl_matrix_set(this->model.qms[0].eig.r_s, i, j, s_elem.dat[0]);
                    gsl_matrix_set(this->model.qms[0].eig.r_s2, i, j, s2_elem.dat[0]);
                }
            }

            gsl_matrix_complex_free(s);
            gsl_matrix_complex_free(s2);
            gsl_vector_complex_free(l);
        }

        this->model.qms[0].tol = tol;
        this->model.qms[0].have_pi = false;
        this->model.qms[0].pi = gsl_vector_alloc(64); // TODO: unused empty vector? // pi = Gsl.Vector.create (fst (Gsl.Matrix.dims qm))
    }
};

std::ostream& operator<<(std::ostream & os, const instance_t::p14n_t & x)
{
//    os << "q_p14ns\n" << x.q_p14ns << '\n';
//    char buf[10];
//    sprintf(buf, "%.3f", x.q_scale_p14ns);
//    os << "q_scale_p14ns: " << buf << '\n';
//    os << "q_domains: " << x.q_domains << '\n';
//    os << "tree_shape:\n" << x.tree_shape << '\n';
//    os << "tree_p14n:\n" << x.tree_p14n << '\n';
//    os << "tree_domains: " << x.tree_domains << '\n';
    return os;
}

std::ostream& operator<<(std::ostream & os, const instance_t::model_t & x)
{
//    os << "tree\n" << x.tree << '\n';
//    os << "qms\n" << x.qms << '\n';
//    os << "pms\n" << x.pms << '\n';
//    os << "prior: " << *x.prior << '\n';
//    char buf[10];
//    sprintf(buf, "%.3f", x.q_scale_p14ns);

    return os;
}

std::ostream& operator<<(std::ostream & os, const instance_t & x)
{
//    os << x.p14n << '\n';
//    os << x.model;
//    os << "tree_settings: " << x.tree_settings << '\n';
//    os << "q_settings: " << x.q_settings << '\n';
    return os;
}

// TODO: let make ?prior t qms
// the new tree and the new qms (created by instantiate_tree or instantiate_qs) are passed be reference in Ocaml
// the results are stored (or computed) in copies.
// on initialization it is okay, when we change them in place. See "let instantiate " in PhyloModel.ml:
// the model, which is part of instance, is constructed the first time on return, i.e., inst or model is not passed by argument
// (because it does not exist yet)
// one subsequent call (for fixedLik, through lpr_leaves):
// calls update() which then eventually calls instantiate_tree, and afterwards make on that tree.
// Then (as described above) qms makes a copy of the reference, but when it returns the copy, it's incorporate into the
// new returned instance "inst" that overwrites the old one in lpr_leaves
// TODO: since in make() in PhyloModel.ml the argument variable (qms) is shadowed by the copy of qms,
//  also named qms, we cannot access old values while computing in qms anyway.
// So I think there is (for now) no need to copy "instance". not sure though how it works for MLE mode where lpr_leaves performs multiple iterations.
// I would expect them in each iteration to work on the latest "instance" but I don't see it in the code yet.
void PhyloModel_make(instance_t & instance, gsl_vector * prior)
{

    // NOTE: taken from the beginning of PhyloModel_make()
    instance.model.pms.resize(instance.p14n.tree_p14n.size() - 1);


    // do this only in the initialization step (copy qms[0] to all other (uninitialized) members)
    if (instance.model.qms.size() == 1)
    {
        instance.model.qms.resize(instance.p14n.tree_p14n.size() - 1);
        for (uint16_t i = 1; i < instance.model.qms.size(); ++i)
            instance.model.qms[i] = instance.model.qms[0]; // deep-copy of object diag_t
        // NOTE: we cannot resize here, because resizing qms will call destructors of qms.model.eig
        // which will delete the matrices (actually it should perform deep-copies of qms[0] members.
        // TODO: figure out why it's not working. maybe a pointer existing to something in qms[0] that is stored outside qms[0]?
        // for simplicity we moved this to the end of PhyloCSFModel_make()
//        instance.model.pms.resize(instance.p14n.tree_p14n.size() - 1);
//        instance.model.qms.resize(instance.p14n.tree_p14n.size() - 1/*, instance.model.qms[0]*/);
//        for (uint16_t i = 1; i < instance.model.qms.size(); ++i)
//            instance.model.qms[i] = instance.model.qms[0]; // deep-copy of object diag_t
    }

    // tree_p14n is the result of instantiate_tree that is passed to make() in PhyloModel.ml
    // it is then returned as a copy in the "model" struct: { tree = T.copy t; qms; pms; prior }
    instance.model.tree = instance.p14n.tree_p14n; // tree with evaluated expressions (i.e., multiplied)
    if (instance.model.prior != NULL)
    {
        gsl_vector_free(instance.model.prior);
        instance.model.prior = NULL;
    }
    if (prior != NULL)
    {
        // this is codon_freq on very first initialization
        instance.model.prior = _deep_copy_vector(prior); // NOTE: this makes a deep-copy
    }

    // pms = Array.init (T.size t - 1) (fun br -> Q.Diag.to_Pt qms.(br) (T.branch t br))
    for (uint16_t i = 0; i < instance.p14n.tree_shape.size() - 1; ++i)
    {
        // to_Pt:
        // TODO: here is memoization happening with q.memoized_to_Pt
        // especially inefficient since for initialization all qms[i] are identical!
        instance_t::model_t::q_diag_t &q = instance.model.qms[i]; // TODO: make sure whether tree[i] matches qms[i] (but should)
        gsl_matrix *&p = instance.model.pms[i];
        const double t = instance.p14n.tree_p14n[i].branch_length;

        if (p != NULL)
            gsl_matrix_free(p);

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

            gsl_vector_free(expLt);

            gsl_blas_dgemm(CBLAS_TRANSPOSE::CblasNoTrans,
                           CBLAS_TRANSPOSE::CblasNoTrans,
                           1.0 /* alpha */,
                           q.eig.r_s /* A */,
                           diagm_c /* B */,
                           0.0 /* beta */,
                           p /* C result */);

            gsl_matrix_free(diagm_c);

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
            std::cout << "NON_REAL!!!\n";
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
            gsl_matrix_complex_free(zdiagm_c);
            gsl_matrix_complex_free(zgemm_c);
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
}

void compute_q_p14ns_and_q_scale_p14ns_fixed_mle( // see instantiate_q and instantiate_qs,
        // memorization "from previous branches" is happening in Ocaml. Optimize here?
        instance_t & instance,
        gsl_matrix * _q, // original from ECM model!!!
        const gsl_vector * const & variables) noexcept
{
    // - multiply every element q[i][j] with a variable j
    // - then set diagonal to: q[i][i] = sum_(i <> j) (-q[i][j])
    instance.p14n.q_p14ns = _q;
    instance.p14n.q_scale_p14ns = 0.0;

    for (uint8_t i = 0; i < 64; ++i)
    {
        assert(gsl_matrix_get(_q, i, i) == 0.0); // otherwise we need an IF before multiplying and summing up below
        double sum = .0;
        for (uint8_t j = 0; j < 64; ++j)
        {
            const double val_ij = gsl_matrix_get(instance.p14n.q_p14ns, i, j) * gsl_vector_get(variables, j);
            gsl_matrix_set(instance.p14n.q_p14ns, i, j, val_ij);
            sum -= val_ij;
        }
        gsl_matrix_set(instance.p14n.q_p14ns, i, i, sum); // set diagonal such that all rows and cols sum up to 0
        instance.p14n.q_scale_p14ns -= sum * gsl_vector_get(variables, i);
    }

//    gsl_matrix_scale(instance.p14n.q_p14ns, 1/instance.p14n.q_scale_p14ns);
    for (uint8_t i = 0; i < 64; ++i)
    {
        for (uint8_t j = 0; j < 64; ++j)
        {
            gsl_matrix_set(instance.p14n.q_p14ns, i, j, gsl_matrix_get(instance.p14n.q_p14ns, i, j) / instance.p14n.q_scale_p14ns);
//            instance.p14n.q_p14ns[i][j] /= instance.p14n.q_scale_p14ns;
        }
    }
}

void PhyloCSFModel_make(instance_t & instance, const empirical_codon_model & ecm, std::vector<newick_elem> & tree_array)
{
    instance.p14n.q_domains = {};
    instance.p14n.tree_shape = tree_array;
//    instance.p14n.tree_p14n = tree_array; // redundant work
    instance.p14n.tree_domains = {domain::Pos};
    instance.q_settings = gsl_vector_alloc(64);
    for (uint8_t i = 0; i < 64; ++i)
        gsl_vector_set(instance.q_settings, i, ecm.codon_freq[i]);
    instance.tree_settings = 1.0;

    // instantiate_tree
    instance.instantiate_tree(instance.tree_settings);

    // instantiate_qs
    gsl_matrix * substitution_matrix = gsl_matrix_alloc(64, 64);
    for (int i = 0; i < 64; ++i)
        for (int j = 0; j < 64; ++j)
            gsl_matrix_set(substitution_matrix, i, j, ecm.matrix[i][j]);

    instance.model.qms.reserve(instance.p14n.tree_p14n.size() - 1);
    compute_q_p14ns_and_q_scale_p14ns_fixed_mle(instance, substitution_matrix, instance.q_settings /*codon_freq*/);

    instance.instantiate_qs();

    // -- make (remember: this is also called in each update (with new tree and old q))
    PhyloModel_make(instance, instance.q_settings/*&codon_freq*/);
}
