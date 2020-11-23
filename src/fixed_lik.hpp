#pragma once

#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

static constexpr int64_t MY_MIN_INT = -(1ULL << 62);
static constexpr int64_t MY_MAX_INT = (1ULL << 62) - 1;

// for now we create a new vector and don't use views, because we don't know what else we will use this function for and how we will handle everything, and scope, etc.
template <typename alpha_t>
gsl_vector * alpha_get(const instance_t & instance, alpha_t & alpha, const int16_t br, const alignment_t & alignment, const uint32_t codon_pos) // don't think & is necessary for alpha
{
    // Inside algo (aka Felsenstein algo)
    const uint16_t k = alpha.matrix.size2; // snd (Gsl.Matrix.dims x.alpha)
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
        auto row = gsl_matrix_row(&alpha.matrix, br - nl);
        gsl_vector * result = _deep_copy_vector(&row.vector);
        return result;
//        return gsl_matrix_row(alpha, br - nl); // Gsl.Matrix.row x.alpha (br-nl)
    }
}

template <typename alpha_t, typename beta_t>
double ensure_alpha(const instance_t & instance, alpha_t & alpha, beta_t & beta, const alignment_t & alignment, const uint32_t codon_pos) // don't think & is necessary for alpha
{
    // if not x.have_alpha then

    const uint16_t n = instance.model.tree.size();
    const uint16_t nl = (instance.model.tree.size() + 1)/2; // let nl = T.leaves tree
    const uint16_t k = alpha.matrix.size2; // let k = snd (Gsl.Matrix.dims x.alpha)

    for (uint16_t i = nl; i < n; ++i)
    {
        const int16_t lc = instance.model.tree[i].child1_id;
        const int16_t rc = instance.model.tree[i].child2_id;
        assert (lc >= 0 && rc >= 0);
        assert (lc < i && rc < i);
        gsl_matrix * ls = instance.model.pms[lc]; // let ls = x.pms.(lc)
        gsl_matrix * rs = instance.model.pms[rc]; // let rs = x.pms.(rc)

        gsl_vector * alc = alpha_get(instance, alpha, lc, alignment, codon_pos);
        gsl_vector * arc = alpha_get(instance, alpha, rc, alignment, codon_pos);

        for (uint16_t a = 0; a < k; ++a)
        {
            auto lsa = gsl_matrix_row(ls, a);
            auto rsa = gsl_matrix_row(rs, a);
            double res1, res2;
            gsl_blas_ddot(&lsa.vector, alc, &res1);
            gsl_blas_ddot(&rsa.vector, arc, &res2);
            gsl_matrix_set(&alpha.matrix, i - nl, a, res1 * res2); // x.alpha.{i-nl,a} <- (Gsl.Blas.dot lsa alc) *. (Gsl.Blas.dot rsa arc)
        }

        gsl_vector_free(alc);
        gsl_vector_free(arc);
    }

    auto x1 = alpha_get(instance, alpha, n - 1, alignment, codon_pos);
    auto x2 = gsl_matrix_row(&beta.matrix, n - 1);
    double final_result;
    gsl_blas_ddot(x1, &x2.vector, &final_result);
    gsl_vector_free(x1);
    return final_result;

    // x.have_alpha <- true
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

double lpr_leaves(instance_t & instance, const alignment_t & alignment, const double t, const uint16_t nbr_leaves_in_tree)
{
    // don't think a deep-copy is necessary. there's a
    // 1.  copy instance (not sure whether it's necessary), but it doesn't seem that
//    instance_t instance_new = instance; // verify whether it's a deep-copy

    // 2.  let inst = PM.P14n.update ~tree_settings:[t] inst
    // 2a. instantiate_tree inst.p14n.tree_shape inst.p14n.tree_p14n t:double // updates instance_new.tree_settings
    instance.compute_tree_p14n(t);
    // 2b. make ?prior:NULL newtree newq // updates instance_new.model
    PhyloModel_make(instance, NULL);

    // let workspace = PhyloLik.new_workspace (PM.tree (PM.P14n.model inst)) Codon.dim
    const uint16_t rows = 2 * instance.p14n.tree_shape.size() - nbr_leaves_in_tree; // 7
    int64_t workspace_generation = MY_MIN_INT; // TODO: min_int from Ocaml, but should be 1ULL << 63??? std::numeric_limits<int64_t>::min()
    gsl_matrix * workspace_data = gsl_matrix_alloc(rows, 64);

    double lpr = 0.0;
    for (uint32_t aa_pos = 0; aa_pos < alignment.peptides[0].size(); ++aa_pos) // lvs = peptides.(...)(pos)
    {
        // let info = PM.prepare_lik workspace instance.model lvs
        // let info = PhyloLik.prepare workspace instance.model.tree instance.model.pms (prior instance.model) lvs
        // prior is a function that either returns instance.model.prior or computes the equilibrium
        auto & m = instance.model;

        std::vector<double> tmp_prior; // don't want to overwrite instance.model.prior (Ocaml code doesn't seem to do that either)
        if (instance.model.prior != NULL)
        {
            tmp_prior = *instance.model.prior; // deep copy
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
                // TODO: implement
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
                    // todo: implement "reversible q" (not important)
                    gsl_vector_free(real_eig_l);
                    gsl_matrix_free(real_eig_s2);
                }
            }

            tmp_prior.resize(q.pi->size);
            for (uint16_t i = 0; i < q.pi->size; ++i)
                tmp_prior[i] = gsl_vector_get(q.pi, i);
        }

//        std::cout << "tmp_prior: " << tmp_prior << '\n';

        // now since (prior instance.model) is computed, we can evaluate the actual value:
        // let info = PhyloLik.prepare workspace instance.model.tree instance.model.pms (prior instance.model) lvs
        const uint16_t k = tmp_prior.size();
        const uint16_t n = instance.model.tree.size();
        const uint16_t nl = (instance.model.tree.size() + 1)/2; // let nl = T.leaves tree

//        std::cout << nl << ' ' << alignment.peptides.size() << '\n';
        assert(nl == alignment.peptides.size()); // if nl <> Array.length leaves then invalid_arg "CamlPaml.Infer.prepare: length(leaves) != leaves(t)"
        assert(instance.model.pms.size() >= n - 1); // if Array.length pms < n-1 then invalid_arg "CamlPaml.Infer.prepare: not enough P matrices"

        // NOTE: ignored, since until now, we always passed a workspace in Ocaml
        // let workspace = match workspace with Some x -> x | None -> new_workspace tree k

        workspace_generation = (workspace_generation == MY_MAX_INT) ? MY_MIN_INT : (workspace_generation + 1);

        assert(workspace_data->size1 >= (2*n-nl)); // if rows < (2*n-nl) || cols <> k then invalid_arg "CamlPaml.Infer.prepare: inappropriate workspace dimensions"
        assert(workspace_data->size2 == k);

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
        auto alpha = gsl_matrix_submatrix(workspace_data, 0, 0, n - nl, 64); // upper "half"
        auto beta  = gsl_matrix_submatrix(workspace_data, n - nl, 0, n, 64); // lower "half"

        for (uint8_t a = 0; a < k; ++a)
        {
            gsl_matrix_set(&beta.matrix, n - 1, a, tmp_prior[a]); // TODO: check whether submatrix works and setting values here!
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

        const double info_z = ensure_alpha(instance, alpha, beta, alignment, aa_pos); // eventually return more than just the info.z score, because it might be used for the ancestor computation
        lpr += log(info_z);
    }

    gsl_matrix_free(workspace_data);

    return lpr;
}