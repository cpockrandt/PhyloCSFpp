#pragma once

#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

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


// this is almost a 1:1 copy of PhyloCSFModel_make(), just beginning and end removed
void Omega_instantiate(instance_t & instance, std::vector<newick_elem> & tree_array)
{
    // computes instance.p14n.tree_p14n
    // (which is implemented with expressions in Ocaml)
    // instantiate_tree
    instance.compute_tree_p14n(instance.tree_settings);

    // instantiate_qs

    // -- computes instance.p14n.q_p14ns and instance.p14n.q_scale_p14ns
    // -- (which are implemented with expressions in Ocaml)

    // this is just because ecm stores raw arrays, and instance.compute_q_p14ns_and_q_scale_p14ns() wants a vector
//    std::vector<std::vector<double> > matrix;
//    std::vector<double> codon_freq;
//    codon_freq.resize(64);
//    matrix.resize(64);
//    for (int i = 0; i < 64; ++i)
//    {
//        codon_freq[i] = ecm.codon_freq[i];
//        matrix[i].resize(64);
//        for (int j = 0; j < 64; ++j)
//            matrix[i][j] = ecm.matrix[i][j];
//    }
//    assert(instance.q_settings == codon_freq);
//
//    instance.compute_q_p14ns_and_q_scale_p14ns(matrix, codon_freq);

// NOVEL NOVEL NOVEL

    // q_settings from std::vector to gsl_vector
    gsl_vector * q_settings_gsl = vector2gsl(instance.q_settings);
    gsl_vector * pi_evaluated = pi_expr(q_settings_gsl);

    auto q_p14n_tmp = comp_q_p14n(q_settings_gsl, pi_evaluated);
    instance.p14n.q_p14ns.resize(64);
    for (uint8_t i = 0; i < 64; ++i)
    {
        instance.p14n.q_p14ns[i].resize(64);
        for (uint8_t j = 0; j < 64; ++j)
        {
            instance.p14n.q_p14ns[i][j] = gsl_matrix_get(q_p14n_tmp, i, j);
        }
    }
    instance.p14n.q_scale_p14ns = comp_q_scale(pi_evaluated, q_p14n_tmp);
    gsl_matrix_free(q_p14n_tmp);

    gsl_vector_free(pi_evaluated);
    gsl_vector_free(q_settings_gsl);

// NOVEL NOVEL NOVEL

//    std::cout << "codon_freq: " << codon_freq << '\n';

//    instance.model.qms.resize(1);
    instance.model.qms.resize(instance.p14n.tree_p14n.size() - 1);
    instance.model.qms[0].q = gsl_matrix_alloc(64, 64);
    for (uint8_t i = 0; i < 64; ++i)
        for (uint8_t j = 0; j < 64; ++j)
            gsl_matrix_set(instance.model.qms[0].q, i, j, instance.p14n.q_p14ns[i][j]);

    // -- of_Q
    // ---- Gsl.Eigen.nonsymmv
    gsl_vector_complex *l = gsl_vector_complex_alloc(64);
    gsl_matrix_complex *s = gsl_matrix_complex_alloc(64, 64);
    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(64);

    gsl_eigen_nonsymmv(instance.model.qms[0].q, l, s, w);
    gsl_eigen_nonsymmv_free(w);
    // NOTE: it seems that order does not seem to matter! can sort in Ocaml and does not change the result
    // NOTE: even though we sort, some values are neg. instead of pos and vice versa
    // TODO: not necessary, just for comparing results with ocaml helpful
    gsl_eigen_nonsymmv_sort (l, s, GSL_EIGEN_SORT_ABS_ASC);

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

        gsl_matrix_complex_free(s);
        gsl_matrix_complex_free(s2);
        gsl_vector_complex_free(l);
    }

    instance.model.qms[0].tol = tol;
    instance.model.qms[0].have_pi = false;
    instance.model.qms[0].pi = gsl_vector_alloc(64); // TODO: unused empty vector? // pi = Gsl.Vector.create (fst (Gsl.Matrix.dims qm))
}