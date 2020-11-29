#pragma once

#include <vector>

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

double comp_q_scale(const gsl_vector * const pi, const gsl_matrix * const q_p14n) noexcept
{
    double factor = 0.0;
    for (uint8_t i = 0; i < 64; ++i)
        factor -= gsl_vector_get(pi, i) * gsl_matrix_get(q_p14n, i, i);
    return factor;
}


gsl_vector * vector2gsl(const std::vector<double> & v)
{
    gsl_vector * gv = gsl_vector_alloc(v.size());
    for (uint8_t i = 0; i < v.size(); ++i)
    {
        gsl_vector_set(gv, i, v[i]);
    }
    return gv;
}

// see instantiate_q and instantiate_qs,
// memorization "from previous branches" is happening in Ocaml. Optimize here?
void compute_q_p14ns_and_q_scale_p14ns_omega(instance_t & instance)
{
    gsl_vector * pi_evaluated = pi_expr(instance.q_settings);

    auto q_p14n_tmp = comp_q_p14n(instance.q_settings, pi_evaluated);
    instance.p14n.q_p14ns = gsl_matrix_alloc(64, 64);
    for (uint8_t i = 0; i < 64; ++i)
    {
        for (uint8_t j = 0; j < 64; ++j)
        {
            gsl_matrix_set(instance.p14n.q_p14ns, i, j, gsl_matrix_get(q_p14n_tmp, i, j));
        }
    }
    instance.p14n.q_scale_p14ns = comp_q_scale(pi_evaluated, q_p14n_tmp);
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




// f shall return std::pair<double, double>(x, f(x))
//template <typename function_t>
//auto new_fit_find_init(const uint32_t max_tries, const double init, const double lo, const double hi, function_t f, minimizer_params_t & params)
//{
//    assert(lo < hi && lo > 0.0);
//
//    srand(0);
//
//    const double width = log(hi) - log(lo);
//    const auto flo = f(lo);
//    const auto fhi = f(hi);
//    double x = init;
//    auto fx = f(init);
////    std::cout << "xxx: " << x << " ||| " << fx << '\n';
//
//    uint32_t i = 0;
//    while (i < max_tries && (std::get<1>(fx) <= std::get<1>(flo) || std::get<1>(fx) <= std::get<1>(fhi)))
//    {
//        const double r = width * ((float) rand()/RAND_MAX);
//        x = exp(log(lo) + r);
//        fx = f(x);
////        std::cout << "xxx: " << x << " ||| " << fx << '\n';
//        ++i;
//    }
////    std::cout << "------------------------------------------------\n";
//    if (i == max_tries)
//        return (flo > fhi) ? flo : fhi;
//    return fx;
//}
//
//void new_maximize_lpr(instance_t &instance, const alignment_t &alignment, const double init, double &lpr, double &elpr_anc, double lo, double hi)
//{
//    constexpr double accuracy = 0.01;
//
//    minimizer_params_t params {instance, alignment};
//    auto good_init = new_fit_find_init(250/*max_tries*/, init /*init*/, lo, hi, f, params);
////    std::cout << good_init << '\n';
//    if (lo < std::get<0>(good_init) && std::get<0>(good_init) < hi)
//    {
//        const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
//        gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);
//        gsl_function F{&minimizer_lpr_leaves, &params};
//
//        gsl_min_fminimizer_set(s, &F, std::get<0>(good_init), lo, hi);
//
//        int64_t max_iter = 250;
//        do {
//            gsl_min_fminimizer_iterate(s); // int status =
//
//            const double x = gsl_min_fminimizer_x_minimum(s);
//            const double lb = gsl_min_fminimizer_x_lower(s);
//            const double ub = gsl_min_fminimizer_x_upper(s);
//
////            printf("[%.7f, %.7f] %.7f %.7f\n", lb, ub, x, params.lpr);
//
//            if (((ub - lb) / x) <= accuracy)
//                break;
//            --max_iter;
//        } while (max_iter > 0);
//
//        gsl_min_fminimizer_free(s);
//
//        lpr = params.lpr;
//        elpr_anc = params.elpr_anc;
//    }
//    else
//    {
//        lpr = std::get<1>(good_init);
//        elpr_anc = std::get<2>(good_init);
//    }
//}