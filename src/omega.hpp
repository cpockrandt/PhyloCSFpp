#pragma once

#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

// this is used in q_p14n and q_scale, but the expressions seem to be evaluated with the same input (I think).
double pi_expr_sc(const gsl_vector * const x, const uint8_t codon_id) noexcept
{
    uint8_t i1, i2, i3;
    // TODO: what happens with N's or gaps???
    from_amino_acid_id(codon_id, i1, i2, i3);

    const double f1 = ((i1 == 3) ? 1.0 : gsl_vector_get(x, 3 + i1)) / (1.0 + gsl_vector_get(x, 3) + gsl_vector_get(x, 4) + gsl_vector_get(x, 5));
    const double f2 = ((i2 == 3) ? 1.0 : gsl_vector_get(x, 6 + i2)) / (1.0 + gsl_vector_get(x, 6) + gsl_vector_get(x, 7) + gsl_vector_get(x, 8));
    const double f3 = ((i3 == 3) ? 1.0 : gsl_vector_get(x, 9 + i3)) / (1.0 + gsl_vector_get(x, 9) + gsl_vector_get(x, 10) + gsl_vector_get(x, 11));

    return f1 * f2 * f3;
}

gsl_vector * pi_expr(const gsl_vector * const x) noexcept
{
    gsl_vector * pi = gsl_vector_alloc(64);

    const double sigma = gsl_vector_get(x, 2);
    const double denom = 1.0 - ((1.0 - sigma) * (pi_expr_sc(x, 3*16 + 0*4 + 0 /*TAA*/) +
                                                 pi_expr_sc(x, 3*16 + 0*4 + 2 /*TAG*/) +
                                                 pi_expr_sc(x, 3*16 + 2*4 + 0 /*TGA*/)));

    for (uint8_t i = 0; i < 64; ++i)
    {
        gsl_vector_set(pi, i, pi_expr_sc(x, i) / denom);
    }

    return pi;
}

gsl_matrix * comp_q_p14n(const gsl_vector * const x, const gsl_vector * const pi) noexcept
{
    gsl_matrix * sq = gsl_matrix_alloc(64, 64);

    const double kappa = gsl_vector_get(x, 0);
    const double omega = gsl_vector_get(x, 1);

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
