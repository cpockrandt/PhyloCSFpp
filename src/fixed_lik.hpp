#pragma once

#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

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
        gsl_vector * result = _deep_copy_vector(&row.vector);
        return result;
//        return gsl_matrix_row(alpha, br - nl); // Gsl.Matrix.row x.alpha (br-nl)
    }
}

void ensure_alpha(const instance_t & instance, workspace_t & workspace, const alignment_t & alignment, const uint32_t codon_pos) // don't think & is necessary for alpha
{
    if (workspace.have_alpha)
        return;

    const uint16_t n = instance.model.tree.size();
    const uint16_t nl = (instance.model.tree.size() + 1)/2; // let nl = T.leaves tree
    const uint16_t k = workspace.alpha.matrix.size2; // let k = snd (Gsl.Matrix.dims x.alpha)

    for (uint16_t i = nl; i < n; ++i)
    {
        const int16_t lc = instance.model.tree[i].child1_id;
        const int16_t rc = instance.model.tree[i].child2_id;
        assert (lc >= 0 && rc >= 0);
        assert (lc < i && rc < i);
        gsl_matrix * ls = instance.model.pms[lc]; // let ls = x.pms.(lc)
        gsl_matrix * rs = instance.model.pms[rc]; // let rs = x.pms.(rc)

        gsl_vector * alc = alpha_get(instance, workspace, lc, alignment, codon_pos);
        gsl_vector * arc = alpha_get(instance, workspace, rc, alignment, codon_pos);

        for (uint16_t a = 0; a < k; ++a)
        {
            auto lsa = gsl_matrix_row(ls, a);
            auto rsa = gsl_matrix_row(rs, a);
            double res1, res2;
            gsl_blas_ddot(&lsa.vector, alc, &res1);
            gsl_blas_ddot(&rsa.vector, arc, &res2);
            gsl_matrix_set(&workspace.alpha.matrix, i - nl, a, res1 * res2); // x.alpha.{i-nl,a} <- (Gsl.Blas.dot lsa alc) *. (Gsl.Blas.dot rsa arc)
        }

        gsl_vector_free(alc);
        gsl_vector_free(arc);
    }

    auto x1 = alpha_get(instance, workspace, n - 1, alignment, codon_pos);
    auto x2 = gsl_matrix_row(&workspace.beta.matrix, n - 1);
    gsl_blas_ddot(x1, &x2.vector, &workspace.z);
    gsl_vector_free(x1);

    workspace.have_alpha = true;
}

void ensure_beta(const instance_t & instance, workspace_t & workspace, const alignment_t & alignment, const uint32_t codon_pos)
{
    ensure_alpha(instance, workspace, alignment, codon_pos);
//    std::cout << "zzzzz\n";
    if (!workspace.have_beta)
    {
        const uint16_t k = workspace.beta.matrix.size2;
        gsl_vector * inter = gsl_vector_alloc(k);
        gsl_vector * ps_colb = gsl_vector_alloc(k);
        for (int16_t i = ((instance.model.tree.size() + 1) / 2) - 1; i >= 0; --i) // for i = (T.root x.tree)-1 downto 0 do
        {
            std::cout << "XXXXXXXXX\n";
            const int16_t p = instance.model.tree[i].parent_id;
            assert(p > i);
            // TODO: update siblings on reduce!
            const int16_t s = instance.model.tree[i].sibling_id;
            gsl_matrix * ps = instance.model.pms[i];
            gsl_matrix * ss = instance.model.pms[s];
            auto bp = gsl_matrix_row(&workspace.beta.matrix, p); // TODO: compare output
            gsl_vector * xas = alpha_get(instance, workspace, s, alignment, codon_pos);

            for (uint16_t a = 0; a < k; ++a)
            {
                // inter.{a} <- bp.{a} *. (Gsl.Blas.dot (Gsl.Matrix.row ss a) xas)
                double dot_result;
                auto ss_row = gsl_matrix_row(ss, a);
                gsl_blas_ddot(&ss_row.vector, xas, &dot_result);
                dot_result *= gsl_vector_get(&bp.vector, a);
                gsl_vector_set(inter, a, dot_result);
            }

            gsl_vector_free(inter);

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
            printf("real_of_complex %g+%gi", GSL_REAL(c), GSL_IMAG(c));
            exit(19);
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
            if (!check_real(c, tol)) {
                printf("real_of_complex %g+%gi", GSL_REAL(c), GSL_IMAG(c));
                exit(20);
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
//        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!! equilibrium " << child2 << "\n";
        auto & q = instance.model.qms[child2];

        // equilibrium
        if (!q.have_pi)
        {
//            std::cout << "Has no PI!\n";
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
//        else
//            std::cout << "Has PI already!\n";

        return q.pi;
    }
}

void lpr_leaves(instance_t & instance, const alignment_t & alignment, const double t, double & lpr, double & elpr_anc)
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

    gsl_vector * anc_lprior = gsl_vector_alloc(64);
    gsl_vector_memcpy(anc_lprior, get_prior(instance));
    for (uint8_t i = 0; i < 64; ++i)
        gsl_vector_set(anc_lprior, i, log(gsl_vector_get(anc_lprior, i)));
//    for (const auto a : anc_lprior)
//        printf("%f ", a);
//    std::cout << '\n';

    const uint16_t root_node_id = instance.model.tree.size() - 1;
//    std::cout << root_node_id << '\n';

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

//        for (uint16_t x = 0; x < 7; ++x)
//        {
//            double test = 100 * (x + 1);
//            for (uint16_t y = 0; y < 64; ++y)
//            {
//                gsl_matrix_set(workspace_data, x, y, test);
//                ++test;
//            }
//        }

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

//        for (uint16_t x = 0; x < 7; ++x)
//        {
//            for (uint16_t y = 0; y < 64; ++y)
//            {
//                std::cout << gsl_matrix_get(workspace_data, x, y) << ' ';
//            }
//            std::cout << '\n';
//        }
//        std::cout << "--------------------------\n";
//        for (uint16_t x = 0; x < alpha.matrix.size1; ++x)
//        {
//            for (uint16_t y = 0; y < alpha.matrix.size2; ++y)
//            {
//                std::cout << gsl_matrix_get(&alpha.matrix, x, y) << ' ';
//            }
//            std::cout << '\n';
//        }
//        std::cout << "--------------------------\n";
//        for (uint16_t x = 0; x < beta.matrix.size1; ++x)
//        {
//            for (uint16_t y = 0; y < beta.matrix.size2; ++y)
//            {
//                std::cout << gsl_matrix_get(&beta.matrix, x, y) << ' ';
//            }
//            std::cout << '\n';
//        }
//        exit(10);

//        { tree = tree; pms = pms; leaves = leaves; workspace = workspace; my_generation = workspace.generation;
//            alpha = alpha; z = nan; have_alpha = false; beta = beta; have_beta = false }

        ensure_alpha(instance, workspace, alignment, aa_pos); // eventually return more than just the info.z score, because it might be used for the ancestor computation
        lpr += log(workspace.z);

        gsl_vector * pr_root = node_posterior(instance, workspace, alignment, aa_pos, root_node_id);
        assert(pr_root->size == anc_lprior->size);
        for (uint16_t xx = 0; xx < anc_lprior->size; ++xx)
        {
            elpr_anc += gsl_vector_get(anc_lprior, xx) * gsl_vector_get(pr_root, xx);
        }
        gsl_vector_free(pr_root);
    }
    gsl_vector_free(anc_lprior);
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
    lpr_leaves(min_params->instance, min_params->alignment, x, min_params->lpr, min_params->elpr_anc);
    return (-1) * min_params->lpr;
}

static size_t rand_id = 0;

double my_ocaml_rand()
{
    static std::vector<double> r = {7.822293,5.798030,3.989764,1.063207,7.482271,3.649649,7.579218,8.293643,8.523592,1.206956,7.282246,6.466290,7.484271,7.682103,1.596537,0.096310,7.516328,5.537377,6.502269,3.281553,7.537951,2.965443,7.916212,5.738046,6.264473,6.752350,0.721391,1.842304,8.722071,3.473594,3.147176,7.485667,7.230315,5.263190,1.426113,5.758685,3.995294,3.889733,5.984792,1.424801,1.125528,1.457110,1.403364,0.995787,2.035756,5.209547,6.550507,7.332020,3.486653,7.001480,3.795608,5.089984,1.844483,8.669368,3.934994,6.944483,5.152360,0.013000,4.518018,6.008967,2.224518,2.423278,5.651355,1.834266,6.497965,2.748319,9.175865,4.997588,3.257924,3.930131,4.350889,3.074516,6.199744,5.137806,8.659628,1.674099,2.633948,5.021263,2.960110,7.187632,6.283709,3.795092,8.775025,1.028682,2.599440,7.846708,3.859253,7.798536,3.020941,3.120367,2.393945,8.650188,1.898054,6.251291,0.036876,5.909511,6.828188,7.889084,2.448866,8.290341,5.576253,2.819084,9.111765,7.796012,2.002562,7.485182,2.942386,2.505075,5.119834,0.060475,1.227442,3.477496,9.045787,2.327140,8.632030,4.078893,0.282804,3.812583,3.577886,0.917564,5.036360,8.611727,3.335440,2.435461,2.733632,6.425171,7.890829,4.100519,6.142246,7.774043,2.159690,3.282000,4.082643,5.339580,0.845743,0.992296,7.467267,7.825164,8.303643,2.511173,7.487288,4.901596,6.340900,1.717617,5.518605,8.779051,8.273968,6.638098,0.281870,7.001361,3.806238,8.168235,9.096416,3.666193,2.846685,6.097120,2.442748,3.787557,6.661030,4.339549,2.024124,4.163124,6.909457,9.041447,8.482296,8.848310,1.103360,3.997171,4.047571,5.790041,5.321769,0.587898,7.167071,8.497480,8.703287,5.331802,8.926357,6.158044,8.032398,0.883934,7.308699,5.625973,1.977815,8.144576,5.990840,8.687504,9.195913,2.201845,7.994647,1.106824,7.358382,7.284752,3.272319,2.523356,3.706573,5.722759,0.064958,2.004703,6.396779,6.086434,1.215683,6.022801,6.589996,3.382713,7.354740,8.171676,1.668836,8.606837,2.842373,1.521058,8.238068,1.437507,6.188858,1.217426,3.629684,0.885047,2.643184,0.905827,6.448448,4.642211,6.273793,3.663685,5.659357,6.019781,6.798340,5.451232,6.357655,7.090366,3.069555,6.267293,2.444180,1.090158,8.972063,7.278150,0.409440,2.000894,4.711551,8.070772,1.602930,3.876541,8.216012,2.271986,3.324750,8.811580,5.222365,7.553890,6.197910,8.879333,7.640081,8.587510};
    return r[rand_id++];
}

template <typename function_t>
void fit_find_init(const uint32_t max_tries, const double init, const double lo, const double hi, function_t * f, minimizer_params_t * params)
{
    assert(lo < hi && lo > 0.0);

    srand(0);

    const double width = log(hi) - log(lo);
    // because the function f is also used for the MLE step (using the GSL minimizer), it negates the results.
    // So we have to "un-negate" it in the init-step
    const double flo = (-1) * (*f)(lo, params);
    const double fhi = (-1) * (*f)(hi, params);
    double x = init;
    double fx = (-1) * (*f)(init, params);
//    std::cout << "xxx: " << x << " ||| " << fx << '\n';

    uint32_t i = 0;
    while (i < max_tries && (fx <= flo || fx <= fhi))
    {
        const double r = width * ((float) rand()/RAND_MAX);
//        const double r = my_ocaml_rand();
//        std::cout << "RANDOM: " << r << '\n';
        x = exp(log(lo) + r);
        fx = (-1) * (*f)(x, params);
//        std::cout << "xxx: " << x << " ||| " << fx << '\n';
        ++i;
    }
//    std::cout << "------------------------------------------------\n";
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
//        return (flo > fhi) ? flo : fhi;
    }
    else
    {
//        return fx; // NOTE: x and lpr in params already set by call to (*f)(x, params)
    }
    (*f)(params->x, params); // NOTE: this is suuuuper improtant! (maybe!) we need to set stuff like tree_settings and pms to the values for the x that we eventually chose!
}

template <typename function_t>
void max_lik_lpr_leaves(instance_t &instance, const alignment_t &alignment, double &lpr, double &elpr_anc, const double init, double lo, double hi, function_t * min_func)
{
    const double accuracy = 0.01;

    minimizer_params_t params {instance, alignment};
    fit_find_init(250/*max_tries*/, init, lo, hi, &min_func, &params);
//    std::cout << "good init: " << params.x << " => " << params.lpr << '\n';
//    exit(44);
    if (lo < params.x && params.x < hi)
    {
        const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
        gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);
        gsl_function F{min_func, &params};

        gsl_min_fminimizer_set(s, &F, params.x, lo, hi);

        int64_t max_iter = 250;
        do {
            gsl_min_fminimizer_iterate(s); // int status =

            const double x = gsl_min_fminimizer_x_minimum(s);
            const double lb = gsl_min_fminimizer_x_lower(s);
            const double ub = gsl_min_fminimizer_x_upper(s);
//            printf("[%.7f, %.7f] %.7f %.7f\n", lb, ub, x, params.lpr);

            if (((ub - lb) / x) <= accuracy)
                break;
            --max_iter;
        } while (max_iter > 0);

        gsl_min_fminimizer_free(s);
    }

    lpr = params.lpr;
    elpr_anc = params.elpr_anc;
}















