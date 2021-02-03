#pragma once

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>

// this is used in q_p14n and q_scale, but the expressions seem to be evaluated with the same input (I think).
double pi_expr_sc(const gsl_vector * const variables, const uint8_t codon_id) noexcept
{
    uint8_t i1, i2, i3;
    // TODO: what happens with N's or gaps???
    from_amino_acid_id(codon_id, i1, i2, i3);

    const double f1 = ((i1 == 3) ? 1.0 : gsl_vector_get(variables, 3 + i1)) / (1.0 + gsl_vector_get(variables, 3) + gsl_vector_get(variables,  4) + gsl_vector_get(variables,  5));
    const double f2 = ((i2 == 3) ? 1.0 : gsl_vector_get(variables, 6 + i2)) / (1.0 + gsl_vector_get(variables, 6) + gsl_vector_get(variables,  7) + gsl_vector_get(variables,  8));
    const double f3 = ((i3 == 3) ? 1.0 : gsl_vector_get(variables, 9 + i3)) / (1.0 + gsl_vector_get(variables, 9) + gsl_vector_get(variables, 10) + gsl_vector_get(variables, 11));

    return f1 * f2 * f3;
}

gsl_vector * pi_expr(const gsl_vector * const variables) noexcept
{
    gsl_vector * pi = gsl_vector_alloc(64);

    const double sigma = gsl_vector_get(variables, 2);
    const double denom = 1.0 - ((1.0 - sigma) * (pi_expr_sc(variables, 3*16 + 0*4 + 0 /*TAA*/) +
                                                 pi_expr_sc(variables, 3*16 + 0*4 + 2 /*TAG*/) +
                                                 pi_expr_sc(variables, 3*16 + 2*4 + 0 /*TGA*/)));

    for (uint8_t i = 0; i < 64; ++i)
    {
        gsl_vector_set(pi, i, pi_expr_sc(variables, i) / denom);
    }

    return pi;
}

gsl_matrix * comp_q_p14n(const gsl_vector * const variables, const gsl_vector * const pi) noexcept
{
    gsl_matrix * sq = gsl_matrix_alloc(64, 64);

    const double kappa = gsl_vector_get(variables, 0);
    const double omega = gsl_vector_get(variables, 1);

    for (uint8_t i = 0; i < 64; ++i)
    {
        uint8_t i1, i2, i3;
        from_amino_acid_id(i, i1, i2, i3);
        const char i_aa = get_amino_acid(i);

        for (uint8_t j = 0; j < 64; ++j)
        {
            // TODO: continue if i==j because we overwrite diagonal below anyway

            uint8_t j1, j2, j3;
            from_amino_acid_id(j, j1, j2, j3);

            double val = 0.0;

            const uint8_t nbr_nucl_changes = (i1 != j1) + (i2 != j2) + (i3 != j3);
            if (nbr_nucl_changes == 1)
            {
                bool transition = false;
                if (i1 != j1 && (i1 + j1 == 2 || i1 + j1 == 4)) // TODO: might get problems with Ns or gaps! // 2 == A and G (0 and 2), 4 == C and T (1 and 3)
                    transition = true;
                if (i2 != j2 && (i2 + j2 == 2 || i2 + j2 == 4))
                    transition = true;
                if (i3 != j3 && (i3 + j3 == 2 || i3 + j3 == 4))
                    transition = true;
                val = (transition) ? kappa : 1.0;

                const char j_aa = get_amino_acid(j);
                val *= (i_aa != '*' && j_aa != '*' && i_aa != j_aa) ? omega : 1.0;
                val *= gsl_vector_get(pi, j);
            }

            gsl_matrix_set(sq, i, j, val);
        }
    }

    // fill_q_diagonal
    for (uint8_t i = 0; i < 64; ++i)
    {
        double val = 0.0;
        for (uint8_t j = 0; j < 64; ++j)
        {
            if (i == j)
                continue;
            val -= gsl_matrix_get(sq, i, j);
        }
        gsl_matrix_set(sq, i, i, val);
    }

    return sq;
}

//void print(instance_t & i)
//{
//
//    for (uint16_t k = 0; k < i.model.qms.size(); ++k)
//    {
//        std::cout << "q[" << k << "]: " << (gsl_matrix_row(i.model.qms[k].q, 0)).vector << '\n'; // TODO: not identical!!!
////        std::cout << "q[" << k << "]:\n" << (gsl_matrix_row(i.model.qms[k].q, 0)).vector << '\n'; // TODO: not identical!!!
//        std::cout << "eig[" << k << "] r_s: " << (gsl_matrix_row(i.model.qms[k].eig.r_s, 0)).vector << '\n';
//        std::cout << "eig[" << k << "] r_s': " << (gsl_matrix_row(i.model.qms[k].eig.r_s2, 0)).vector << '\n';
//        std::cout << "eig[" << k << "] r_l: " << *i.model.qms[k].eig.r_l << '\n';
//        std::cout << "pi[" << k << "]: " << *i.model.qms[k].pi << '\n';
//    }
//
//    for (uint16_t k = 0; k < i.model.pms.size(); ++k)
//        std::cout << "pms[" << k << "]: " << (gsl_matrix_row(i.model.pms[k], 0)).vector << '\n'; // TODO: not identical!!!
//
////    // identical
////    std::cout << i.model.tree << '\n';
//    if (i.model.prior != NULL)
//        std::cout << "PRIOR: " << *i.model.prior << '\n';
//    else
//        std::cout << "PRIOR: No inst.model.prior" << '\n';
//
////    // identical
////    std::cout << *i.p14n.q_p14ns << '\n'; // is correct (after new_instance and after initial update)
////    std::cout << i.p14n.q_scale_p14ns << '\n'; // is correct (after new_instance and after initial update)
////    std::cout << i.p14n.q_domains << '\n';
////    std::cout << i.p14n.tree_shape << '\n';
////    std::cout << i.p14n.tree_p14n << '\n';
////    std::cout << i.p14n.tree_domains << '\n';
////
////    // identical
//    std::cout << "q_settings: " << *i.q_settings << '\n';
//    char buf[10];
//    sprintf(buf, "%f ", i.tree_settings);
//    std::cout << "tree_settings: [" << buf << ']' << '\n';
//}

double comp_q_scale(const gsl_vector * const pi, const gsl_matrix * const q_p14n) noexcept
{
    double factor = 0.0;
    for (uint8_t i = 0; i < 64; ++i)
        factor -= gsl_vector_get(pi, i) * gsl_matrix_get(q_p14n, i, i);
    return factor;
}

// see instantiate_q and instantiate_qs,
// memorization "from previous branches" is happening in Ocaml. Optimize here?
void compute_q_p14ns_and_q_scale_p14ns_omega(instance_t & instance)
{
    gsl_vector * pi_evaluated = pi_expr(instance.q_settings);

    gsl_matrix * q_p14n_tmp = comp_q_p14n(instance.q_settings, pi_evaluated);

    instance.p14n.q_scale_p14ns = comp_q_scale(pi_evaluated, q_p14n_tmp);

    if (instance.p14n.q_p14ns == NULL)
        instance.p14n.q_p14ns = gsl_matrix_alloc(64, 64);

    for (uint8_t i = 0; i < 64; ++i)
    {
        for (uint8_t j = 0; j < 64; ++j)
        {
            gsl_matrix_set(instance.p14n.q_p14ns, i, j, gsl_matrix_get(q_p14n_tmp, i, j) / instance.p14n.q_scale_p14ns);
        }
    }
//    std::cout << "q_scale_p14ns: " << instance.p14n.q_scale_p14ns << '\n';
    gsl_matrix_free(q_p14n_tmp);

    gsl_vector_free(pi_evaluated);
}

double get_lpr_rho(const double rho) noexcept
{
    constexpr double mode = 1.0;
    constexpr double scale = 0.5;
    // half_cauchy_lpdf
    assert(!(rho < 0.0 || scale <= 0.0 || mode < 0.0));
    const double numer = 1.0 / (M_PI * scale * (1.0 + pow(((rho - mode) / scale), 2.0)));
    const double cauchy_cdf = atan ((0.0 - mode) / scale) / M_PI + 0.5;
    const double denom = 1.0 - cauchy_cdf;
    assert(denom > 0.0);
    return log(numer) - log(denom);
}

double get_lpr_kappa(const double kappa) noexcept
{
    const double epsilon_float = DBL_EPSILON;
//    printf("%f\n", epsilon_float);
    const double k = kappa - 1.0 + epsilon_float;
    const double gamma = gsl_ran_gamma_pdf(k, 7.0, 0.25);
    return log(gamma);
}



void lpr_leaves_omega(instance_t & instance, const alignment_t & alignment, const double /*t*/, double & lpr)
{
    // let workspace = PhyloLik.new_workspace (PM.tree (PM.P14n.model inst)) Codon.dim
    workspace_t workspace;
    const uint16_t nbr_leaves_in_tree = (instance.model.tree.size() + 1) / 2;
    const uint16_t rows = 2 * instance.p14n.tree_shape.size() - nbr_leaves_in_tree;
    workspace.workspace_generation = MY_MIN_INT; // TODO: min_int from Ocaml, but should be 1ULL << 63??? std::numeric_limits<int64_t>::min()
    workspace.workspace_data = gsl_matrix_alloc(rows, 64);

//    const uint16_t root_node_id = instance.model.tree.size() - 1; // TODO: why did we remove this???
//    std::cout << root_node_id << '\n';

    lpr = 0.0;
    for (uint32_t aa_pos = 0; aa_pos < alignment.peptides[0].size(); ++aa_pos) // lvs = peptides.(...)(pos)
    {
        // let info = PM.prepare_lik workspace instance.model lvs
        // let info = PhyloLik.prepare workspace instance.model.tree instance.model.pms (prior instance.model) lvs
        // prior is a function that either returns instance.model.prior or computes the equilibrium
//        auto & m = instance.model;

        // now since (prior instance.model) is computed, we can evaluate the actual value:
        // let info = PhyloLik.prepare workspace instance.model.tree instance.model.pms (prior instance.model) lvs

        gsl_vector * tmp_prior = get_prior(instance); // TODO: move this out of the for loop?
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
//                gsl_matrix_set(workspace.workspace_data, x, y, test);
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

        gsl_matrix_set_zero(&workspace.alpha.matrix);

//        std::cout << workspace.alpha.matrix << '\n';

        for (uint8_t a = 0; a < k; ++a)
        {
            gsl_matrix_set(&workspace.beta.matrix, n - 1, a, gsl_vector_get(tmp_prior, a)); // TODO: check whether submatrix works and setting values here!
        }

//        for (uint16_t x = 0; x < 7; ++x)
//        {
//            for (uint16_t y = 0; y < 64; ++y)
//                std::cout << gsl_matrix_get(workspace.workspace_data, x, y) << ' ';
//            std::cout << '\n';
//        }
//        std::cout << "--------------------------\n";
//        for (uint16_t x = 0; x < workspace.alpha.matrix.size1; ++x)
//        {
//            for (uint16_t y = 0; y < workspace.alpha.matrix.size2; ++y)
//                std::cout << gsl_matrix_get(&workspace.alpha.matrix, x, y) << ' ';
//            std::cout << '\n';
//        }
//        std::cout << "--------------------------\n";
//        for (uint16_t x = 0; x < workspace.beta.matrix.size1; ++x)
//        {
//            for (uint16_t y = 0; y < workspace.beta.matrix.size2; ++y)
//                std::cout << gsl_matrix_get(&workspace.beta.matrix, x, y) << ' ';
//            std::cout << '\n';
//        }

        ensure_alpha(instance, workspace, alignment, aa_pos); // eventually return more than just the info.z score, because it might be used for the ancestor computation
//        std::cout << workspace.alpha.matrix << '\n';
        lpr += log(workspace.z);
//        std::cout << "log summand: " << log(workspace.z) << '\n';
//        printf("log summand: %f\n", log(workspace.z));
    }
}


double minimizer_lpr_leaves_rho(const double x, void * params)
{
    minimizer_params_t * min_params = (minimizer_params_t*) params;

    // PM.P14n.update ~tree_settings:ts inst
    min_params->x = x;
    min_params->instance.tree_settings = x;
    min_params->instance.instantiate_tree(x); // x is a variable that is used for maximization in maximize_lpr
    PhyloModel_make(min_params->instance, NULL, false);

    lpr_leaves_omega(min_params->instance, min_params->alignment, x, min_params->lpr/*, min_params->elpr_anc*/);

//    print(min_params->instance);
//    std::cout << "get_lpr_rho(" << x << "): " << get_lpr_rho(x) << '\n';
//    printf("lpr_leaves result (rho, without get_rho): %f => %f\n", x, min_params->lpr);
//    std::cout << "lpr_leaves result (rho, without get_rho): " << x << " => " << min_params->lpr << '\n';
    min_params->lpr += get_lpr_rho(x);
//    printf("lpr_leaves result (rho, with    get_rho): %f => %f\n", x, min_params->lpr);
//    std::cout << "lpr_leaves result (rho, with    get_rho): " << x << " => " << min_params->lpr << '\n';
    return (-1) * min_params->lpr;
}

double minimizer_lpr_leaves_kappa(const double x, void * params)
{
    minimizer_params_t * min_params = (minimizer_params_t*) params;

    // qms identical, pms is not! how can that be? pms is computed based on qms???

    // PM.P14n.update ~tree_settings:ts inst
    min_params->x = x;
    gsl_vector_set(min_params->instance.q_settings, 0, x);
    compute_q_p14ns_and_q_scale_p14ns_omega(min_params->instance);
    min_params->instance.instantiate_qs(); // x is a variable that is used for maximization in maximize_lpr
    PhyloModel_make(min_params->instance, NULL, true);

    lpr_leaves_omega(min_params->instance, min_params->alignment, x, min_params->lpr/*, min_params->elpr_anc*/);
//    double lpr_old = min_params->lpr;
    min_params->lpr += get_lpr_kappa(x);
//    printf("lpr_leaves result (kap, with    get_kap): %f => %f %f\n", x, min_params->lpr, lpr_old);
//    print(min_params->instance);
//    printf("lpr_leaves result (kap, with    get_kap): %f => %f\n", x, min_params->lpr);
//    printf("get_kappa(%f): %f\n", x, get_lpr_kappa(x));
    return (-1) * min_params->lpr;
}