#pragma once

#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

static constexpr int64_t MY_MIN_INT = -(1ULL << 62);
static constexpr int64_t MY_MAX_INT = (1ULL << 62) - 1;

struct workspace_t
{
    int64_t workspace_generation;
    gsl_matrix * workspace_data;

    gsl_matrix_view alpha;
    gsl_matrix_view beta;

    bool have_alpha;
    bool have_beta;

    double z;

    ~workspace_t()
    {
        gsl_matrix_free(this->workspace_data);
    }
};

// for now we create a new vector and don't use views, because we don't know what else we will use this function for and how we will handle everything, and scope, etc.
gsl_vector * alpha_get(const instance_t & instance, workspace_t & workspace, const int16_t br,
                       const alignment_t & alignment, const uint32_t codon_pos) // don't think & is necessary for alpha
{
    // Inside algo (aka Felsenstein algo)
    const uint16_t k = workspace.alpha.matrix.size2; // snd (Gsl.Matrix.dims x.alpha)
    const uint16_t nl = (instance.model.tree.size() + 1)/2; // let nl = T.leaves tree
    if (br < nl)
    {
        // get_leaf k alignment.peptides[br][codon_pos]
        gsl_vector * result = gsl_vector_alloc(k);
        if (alignment.peptides[br][codon_pos] == 64) // Marginalize code
        {
            // Gsl.Vector.of_array (Array.make k 1.)
            gsl_vector_set_all(result, 1.0);
        }
        else // Certain i
        {
            // raw_leaf (k, i = alignment.peptides[br][codon_pos] i.e. AA_id)
            // Gsl.Vector.of_array (Array.init k (fun j -> if i = j then 1. else 0.))
            for (uint16_t tmp_j = 0; tmp_j < k; ++tmp_j)
            {
                if (tmp_j == alignment.peptides[br][codon_pos])
                    gsl_vector_set(result, tmp_j, 1.0);
                else
                    gsl_vector_set(result, tmp_j, 0.0);
            }
        }
        return result;
    }
    else
    {
        auto row = gsl_matrix_row(&workspace.alpha.matrix, br - nl);
        gsl_vector * result = NULL;
        _deeep_copy_vector(&result, &row.vector);
        return result;
//        return gsl_matrix_row(alpha, br - nl); // Gsl.Matrix.row x.alpha (br-nl)
    }
}

struct alpha_result
{
    gsl_vector * alpha = gsl_vector_alloc(64);
    uint8_t aa = 0;
    bool is_distribution = false; // TODO: is the naming correct? bl < nl

    ~alpha_result()
    {
        gsl_vector_free(alpha);
        alpha = NULL;
    }
};

// for now we create a new vector and don't use views, because we don't know what else we will use this function for and how we will handle everything, and scope, etc.
void alpha_get(const instance_t & instance, workspace_t & workspace, const int16_t br,
               const alignment_t & alignment, const uint32_t codon_pos, alpha_result & result) // don't think & is necessary for alpha
{
    // Inside algo (aka Felsenstein algo)
    assert(workspace.alpha.matrix.size2 == 64); // = k
    const uint16_t nl = (instance.model.tree.size() + 1)/2; // let nl = T.leaves tree
    if (br < nl)
    {
        result.is_distribution = false;
        result.aa = alignment.peptides[br][codon_pos];
    }
    else
    {
        result.is_distribution = true;
        auto row = gsl_matrix_row(&workspace.alpha.matrix, br - nl);
        _deeep_copy_vector(&result.alpha, &row.vector);
    }
}

void dot_with_alpha(gsl_vector * v, alpha_result & result, double & d)
{
    if (result.is_distribution)
    {
        gsl_blas_ddot(v, result.alpha, &d);
    }
    else if (result.aa == 64)
    {
        // gsl_blas_ddot(v, [1 1 1 1 ... 1], &d);
        d = 0.0;
        for (size_t i = 0; i < v->size; ++i)
            d += gsl_vector_get(v, i);
    }
    else // if (result.aa != 64)
    {
        // gsl_blas_ddot(v, [0 0 ... 0 1 0 ... 0 0], &d); // only at position result.alpha is a 1
        d = gsl_vector_get(v, result.aa);
    }
}

void ensure_alpha(const instance_t & instance, workspace_t & workspace, const alignment_t & alignment, const uint32_t codon_pos) // don't think & is necessary for alpha
{
    if (workspace.have_alpha)
        return;

    const uint16_t n = instance.model.tree.size();
    const uint16_t nl = (instance.model.tree.size() + 1)/2; // let nl = T.leaves tree
    const uint16_t k = workspace.alpha.matrix.size2; // let k = snd (Gsl.Matrix.dims x.alpha)

    alpha_result alc, arc;
    for (uint16_t i = nl; i < n; ++i)
    {
        const int16_t lc = instance.model.tree[i].child1_id;
        const int16_t rc = instance.model.tree[i].child2_id;
        assert (lc >= 0 && rc >= 0);
        assert (lc < i && rc < i);
        gsl_matrix * ls = instance.model.pms[lc]; // let ls = x.pms.(lc)
        gsl_matrix * rs = instance.model.pms[rc]; // let rs = x.pms.(rc)

        alpha_get(instance, workspace, lc, alignment, codon_pos, alc);
        alpha_get(instance, workspace, rc, alignment, codon_pos, arc);

        for (uint16_t a = 0; a < k; ++a)
        {
            auto lsa = gsl_matrix_row(ls, a);
            auto rsa = gsl_matrix_row(rs, a);
            double res1, res2;
            dot_with_alpha(&lsa.vector, alc, res1);
            dot_with_alpha(&rsa.vector, arc, res2);

            gsl_matrix_set(&workspace.alpha.matrix, i - nl, a, res1 * res2); // x.alpha.{i-nl,a} <- (Gsl.Blas.dot lsa alc) *. (Gsl.Blas.dot rsa arc)
        }
    }

    alpha_get(instance, workspace, n - 1, alignment, codon_pos, alc);
    auto x2 = gsl_matrix_row(&workspace.beta.matrix, n - 1);
    dot_with_alpha(&x2.vector, alc, workspace.z);

    workspace.have_alpha = true;
}

void ensure_beta(const instance_t & instance, workspace_t & workspace, const alignment_t & alignment, const uint32_t codon_pos)
{
    ensure_alpha(instance, workspace, alignment, codon_pos);
    if (!workspace.have_beta)
    {
        const uint16_t k = workspace.beta.matrix.size2;
        gsl_vector * inter = gsl_vector_alloc(k);
        gsl_vector * ps_colb = gsl_vector_alloc(k);

        alpha_result a_res;

        for (int16_t i = ((instance.model.tree.size() + 1) / 2) - 1; i >= 0; --i) // for i = (T.root x.tree)-1 downto 0 do
        {
            const int16_t p = instance.model.tree[i].parent_id;
            assert(p > i);
            // TODO: update siblings on reduce!
            const int16_t s = instance.model.tree[i].sibling_id;
            gsl_matrix * ps = instance.model.pms[i];
            gsl_matrix * ss = instance.model.pms[s];
            auto bp = gsl_matrix_row(&workspace.beta.matrix, p); // TODO: compare output
            alpha_get(instance, workspace, s, alignment, codon_pos, a_res);

            for (uint16_t a = 0; a < k; ++a)
            {
                // inter.{a} <- bp.{a} *. (Gsl.Blas.dot (Gsl.Matrix.row ss a) xas)
                double dot_result;
                auto ss_row = gsl_matrix_row(ss, a);
                dot_with_alpha(&ss_row.vector, a_res, dot_result);
                dot_result *= gsl_vector_get(&bp.vector, a);
                gsl_vector_set(inter, a, dot_result);
            }

            for (uint16_t b = 0; b < k; ++b)
            {
                for (uint16_t a = 0; a < k; ++a) // instead of this loop it'd probably be a little faster to transpose ps
                {
                    gsl_vector_set(ps_colb, a, gsl_matrix_get(ps, a, b));
//                    Bigarray.Array1.unsafe_set ps_colb a (Bigarray.Array2.unsafe_get ps a b)
                }
                double dot_result;
                gsl_blas_ddot(inter, ps_colb, &dot_result);
                gsl_matrix_set(&workspace.beta.matrix, i, b, dot_result);
            }
        }

        gsl_vector_free(ps_colb);
        workspace.have_beta = true;
    }
}

gsl_vector * node_posterior(instance_t & instance, workspace_t & workspace, const alignment_t & alignment, const uint32_t codon_pos,
                            const uint16_t node_id)
{
    // assert(workspace_generation == inferred.my_generation) // invalidated workspace
    ensure_alpha(instance, workspace, alignment, codon_pos);
    const uint16_t k = workspace.alpha.matrix.size2;

    if (workspace.z == 0.0) // the data were impossible by the model
    {
        gsl_vector * pr_root = gsl_vector_alloc(k);
        gsl_vector_set_zero(pr_root);
        return pr_root;
    }
    else if (node_id < (instance.model.tree.size() + 1) / 2) // (root_node_id <  T.leaves tree) leaf: usually certain
    {
        return alpha_get(instance, workspace, node_id, alignment, codon_pos); // alpha_get(inferred node_id)
    }
    else // calculate the betas (except for the root, whose betas get prefilled in prepare)
    {
        gsl_vector * pr_root = gsl_vector_alloc(k);
        if (node_id != instance.model.tree.size() - 1)
            ensure_beta(instance, workspace, alignment, codon_pos);
        gsl_vector * alpha = alpha_get(instance, workspace, node_id, alignment, codon_pos); // alpha_get(inferred node_id)
        for (uint16_t x = 0; x < k; ++x)
        {
            const double val = gsl_vector_get(alpha, x) * gsl_matrix_get(&workspace.beta.matrix, node_id, x) / workspace.z;
            gsl_vector_set(pr_root, x, val);
        }
        gsl_vector_free(alpha);
        return pr_root;
    }
}

gsl_vector * complex_to_real(const gsl_vector_complex * const cv, const double tol = 1e-6)
{
    gsl_vector * v = gsl_vector_alloc(cv->size);
    for (uint16_t i = 0; i < cv->size; ++i)
    {
        const gsl_complex & c = gsl_vector_complex_get(cv, i);
        if (!check_real(c, tol))
        {
            throw std::runtime_error("real_of_complex " + std::to_string(GSL_REAL(c)) + "+" + std::to_string(GSL_IMAG(c)) + "i");
        }
        gsl_vector_set(v, i, GSL_REAL(c));
    }
    return v;
}

gsl_matrix * complex_to_real(const gsl_matrix_complex * const cm, const double tol = 1e-6)
{
    gsl_matrix * m = gsl_matrix_alloc(cm->size1, cm->size2);
    for (uint16_t i = 0; i < cm->size1; ++i)
    {
        for (uint16_t j = 0; j < cm->size2; ++j)
        {
            const gsl_complex &c = gsl_matrix_complex_get(cm, i, j);
            if (!check_real(c, tol))
            {
                throw std::runtime_error("real_of_complex " + std::to_string(GSL_REAL(c)) + "+" + std::to_string(GSL_IMAG(c)) + "i");
            }
            gsl_matrix_set(m, i, j, GSL_REAL(c));
        }
    }
    return m;
}

gsl_vector * get_prior(instance_t & instance)
{
    if (instance.model.prior != NULL)
    {
        return instance.model.prior;
    }
    else
    {
        // tmp_prior = Q.Diag.equilibrium qms.(snd (T.children instance.model.tree (T.root instance.model.tree)))
        // qms is instance.model.qms
        const uint16_t root = instance.model.tree.size() - 1;
        const uint16_t child2 = instance.model.tree[root].child2_id;
        auto & q = instance.model.qms[child2];

        // equilibrium
        if (!q.have_pi)
        {
            gsl_vector * real_eig_l;
            gsl_matrix * real_eig_s2;
            assert((
                           q.eig.r_l == NULL && q.eig.r_s == NULL && q.eig.r_s2 == NULL &&
                           q.eig.nr_l != NULL && q.eig.nr_s != NULL && q.eig.nr_s2 != NULL
                   ) ||
                   (
                           q.eig.r_l != NULL && q.eig.r_s != NULL && q.eig.r_s2 != NULL &&
                           q.eig.nr_l == NULL && q.eig.nr_s == NULL && q.eig.nr_s2 == NULL
                   )
            );
            if (q.eig.r_l == NULL)
            {
                // todo: implement "reversible q" (not important)
                // if no (reversible q) then failwith "CamlPaml.Q.equilibrium: non-reversible model"

                real_eig_l = complex_to_real(q.eig.nr_l);
                real_eig_s2 = complex_to_real(q.eig.nr_s2);
            }
            else
            {
                real_eig_l = q.eig.r_l;
                real_eig_s2 = q.eig.r_s2;
            }

            const uint16_t n = q.eig.r_l->size;
            double min_L = fabs(gsl_vector_get(q.eig.r_l, 0)); // infinity
            int64_t min_Lp = 0;
            for (uint16_t i = 1; i < n; ++i)
            {
                const double mag_i = fabs(gsl_vector_get(q.eig.r_l, i));
                if (mag_i < min_L)
                {
                    min_L = mag_i;
                    min_Lp = i;
                }
            }
//                assert (min_Lp >= 0);
            if (fabs(min_L) > q.tol)
                printf("equilibrium: smallest-magnitude eigenvalue %e is unacceptably large; check rate matrix validity or increase tol", min_L);
//                auto lev = gsl_matrix_row(q.eig.r_s2, min_Lp);
            double mass = 0.0;
            for (uint16_t j = 0; j < n; ++j)
            {
                mass += gsl_matrix_get(q.eig.r_s2, min_Lp, j); // equiv. to get val j from vector lev
            }
            for (uint16_t j = 0; j < n; ++j)
            {
                gsl_vector_set(q.pi, j, gsl_matrix_get(q.eig.r_s2, min_Lp, j) / mass); // equiv. to get val j from vector lev
            }

            q.have_pi = true;

            if (q.eig.r_l == NULL)
            {
                // todo: implement assert with "reversible q" (not important)
                gsl_vector_free(real_eig_l);
                gsl_matrix_free(real_eig_s2);
            }
        }

        return q.pi;
    }
}

void lpr_leaves(instance_t & instance, const alignment_t & alignment, const double t, double & lpr, double & elpr_anc, std::vector<double> & lpr_per_codon)
{
    // don't think a deep-copy is necessary. there's a
    // 1.  copy instance (not sure whether it's necessary), but it doesn't seem that
//    instance_t instance_new = instance; // verify whether it's a deep-copy

    // 2.  let inst = PM.P14n.update ~tree_settings:[t] inst
    // 2a. instantiate_tree inst.p14n.tree_shape inst.p14n.tree_p14n t:double // updates instance_new.tree_settings
    instance.instantiate_tree(t);
    // 2b. make ?prior:NULL newtree newq // updates instance_new.model
    PhyloModel_make(instance, NULL, false);

    // let workspace = PhyloLik.new_workspace (PM.tree (PM.P14n.model inst)) Codon.dim
    workspace_t workspace;
    const uint16_t nbr_leaves_in_tree = (instance.model.tree.size() + 1) / 2;
    const uint16_t rows = 2 * instance.p14n.tree_shape.size() - nbr_leaves_in_tree;
    workspace.workspace_generation = MY_MIN_INT; // TODO: min_int from Ocaml, but should be 1ULL << 63??? std::numeric_limits<int64_t>::min()
    workspace.workspace_data = gsl_matrix_alloc(rows, 64);

// BEGIN - THIS IS FOR ANC SCORE
    gsl_vector * anc_lprior = gsl_vector_alloc(64);
    gsl_vector_memcpy(anc_lprior, get_prior(instance));
    for (uint8_t i = 0; i < 64; ++i)
        gsl_vector_set(anc_lprior, i, log(gsl_vector_get(anc_lprior, i)));
    const uint16_t root_node_id = instance.model.tree.size() - 1;
// END - THIS IS FOR ANC SCORE

    gsl_vector * tmp_prior = get_prior(instance); // TODO: move this out of the for loop below!

    lpr = 0.0;
    elpr_anc = 0.0;
    for (uint32_t aa_pos = 0; aa_pos < alignment.peptides[0].size(); ++aa_pos) // lvs = peptides.(...)(pos)
    {
        // let info = PM.prepare_lik workspace instance.model lvs
        // let info = PhyloLik.prepare workspace instance.model.tree instance.model.pms (prior instance.model) lvs
        // prior is a function that either returns instance.model.prior or computes the equilibrium
//        auto & m = instance.model;

        // now since (prior instance.model) is computed, we can evaluate the actual value:
        // let info = PhyloLik.prepare workspace instance.model.tree instance.model.pms (prior instance.model) lvs
        const uint16_t k = tmp_prior->size;
        const uint16_t n = instance.model.tree.size();
        const uint16_t nl = (instance.model.tree.size() + 1)/2; // let nl = T.leaves tree

//        std::cout << nl << ' ' << alignment.peptides.size() << '\n';
        assert(nl == alignment.peptides.size()); // if nl <> Array.length leaves then invalid_arg "CamlPaml.Infer.prepare: length(leaves) != leaves(t)"
        assert(instance.model.pms.size() >= (size_t)(n - 1)); // if Array.length pms < n-1 then invalid_arg "CamlPaml.Infer.prepare: not enough P matrices"

        // NOTE: ignored, since until now, we always passed a workspace in Ocaml
        // let workspace = match workspace with Some x -> x | None -> new_workspace tree k

        workspace.workspace_generation = (workspace.workspace_generation == MY_MAX_INT) ? MY_MIN_INT : (workspace.workspace_generation + 1);

        assert(workspace.workspace_data->size1 >= (size_t)(2*n-nl)); // if rows < (2*n-nl) || cols <> k then invalid_arg "CamlPaml.Infer.prepare: inappropriate workspace dimensions"
        assert(workspace.workspace_data->size2 == k);

        // let alpha = Bigarray.Array2.sub_left workspace.data 0 (n-nl)
        // let beta = Bigarray.Array2.sub_left workspace.data (n-nl) n
        // I don't think submatrices are really necessary here. we can do the calculations ourselves
        workspace.have_alpha = false;
        workspace.have_beta = false;
        workspace.alpha = gsl_matrix_submatrix(workspace.workspace_data, 0, 0, n - nl, 64); // upper "half"
        workspace.beta  = gsl_matrix_submatrix(workspace.workspace_data, n - nl, 0, n, 64); // lower "half"

        for (uint8_t a = 0; a < k; ++a)
        {
            gsl_matrix_set(&workspace.beta.matrix, n - 1, a, gsl_vector_get(tmp_prior, a)); // TODO: check whether submatrix works and setting values here!
        }

        ensure_alpha(instance, workspace, alignment, aa_pos); // eventually return more than just the info.z score, because it might be used for the ancestor computation
        const double log_z = log(workspace.z);
        lpr += log_z;
        lpr_per_codon.push_back(log_z);

// BEGIN - THIS IS FOR ANC SCORE
        gsl_vector * pr_root = node_posterior(instance, workspace, alignment, aa_pos, root_node_id);
        assert(pr_root->size == anc_lprior->size);
        for (uint16_t xx = 0; xx < anc_lprior->size; ++xx)
        {
            elpr_anc += gsl_vector_get(anc_lprior, xx) * gsl_vector_get(pr_root, xx);
        }
        gsl_vector_free(pr_root);
// END - THIS IS FOR ANC SCORE
    }
// BEGIN - THIS IS FOR ANC SCORE
    gsl_vector_free(anc_lprior);
// END - THIS IS FOR ANC SCORE
}

struct minimizer_params_t {
    instance_t& instance;
    const alignment_t& alignment;
    double x;
    double lpr;
    double elpr_anc;
};

double minimizer_lpr_leaves(const double x, void * params)
{
    minimizer_params_t * min_params = (minimizer_params_t*) params;
    min_params->x = x;
    std::vector<double> TODO_remove;
    lpr_leaves(min_params->instance, min_params->alignment, x, min_params->lpr, min_params->elpr_anc, TODO_remove);
    return (-1) * min_params->lpr;
}

template <typename function_t>
void fit_find_init(const uint32_t max_tries, const double init, const double lo, const double hi, function_t * f, minimizer_params_t * params, std::mt19937 & gen)
{
    assert(lo < hi && lo > 0.0);

    const double width = log(hi) - log(lo);
    // because the function f is also used for the MLE step (using the GSL minimizer), it negates the results.
    // So we have to "un-negate" it in the init-step
    const double flo = (-1) * (*f)(lo, params);
    const double fhi = (-1) * (*f)(hi, params);
    double x = init;
    double fx = (-1) * (*f)(init, params);

    std::uniform_real_distribution<> dis(0.0, width);

    uint32_t i = 0;
    while (i < max_tries && (fx <= flo || fx <= fhi))
    {
        const double r = dis(gen); // width * ((float) rand()/RAND_MAX);
        x = exp(log(lo) + r);
        fx = (-1) * (*f)(x, params);
        ++i;
    }

    if (i == max_tries)
    {
        if (flo > fhi)
        {
            params->x = lo;
            params->lpr = flo;
        }
        else
        {
            params->x = hi;
            params->lpr = fhi;
        }
    }
    // NOTE: otherwise x and lpr in params already set by call to (*f)(x, params)

    (*f)(params->x, params); // NOTE: this is suuuuper improtant! (maybe!) we need to set stuff like tree_settings and pms to the values for the x that we eventually chose!
}

template <typename function_t>
void max_lik_lpr_leaves(instance_t &instance, const alignment_t &alignment, double &lpr, double &elpr_anc, const double init, double lo, double hi, function_t * min_func, std::mt19937 & gen)
{
    const double accuracy = 0.01;

    minimizer_params_t params {instance, alignment, 0.0, 0.0, 0.0};
    fit_find_init(250/*max_tries*/, init, lo, hi, &min_func, &params, gen);
    if (lo < params.x && params.x < hi)
    {
        const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
        instance.minimizer = gsl_min_fminimizer_alloc(T);
        gsl_function F{min_func, &params};

        gsl_min_fminimizer_set(instance.minimizer, &F, params.x, lo, hi);

        int64_t max_iter = 250;
        do {
            gsl_min_fminimizer_iterate(instance.minimizer);

            const double x = gsl_min_fminimizer_x_minimum(instance.minimizer);
            const double lb = gsl_min_fminimizer_x_lower(instance.minimizer);
            const double ub = gsl_min_fminimizer_x_upper(instance.minimizer);

            if (((ub - lb) / x) <= accuracy)
                break;
            --max_iter;
        } while (max_iter > 0);
        gsl_min_fminimizer_free(instance.minimizer);
        instance.minimizer = NULL;
    }

    lpr = params.lpr;
    elpr_anc = params.elpr_anc;
}
