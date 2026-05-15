/**
 * @file dp.c
 * @brief Diagonal Polynomial packing function implementation.
 */

#include <math.h>
#include "dp.h"
#include "gtp.h"
#include "arith.h"

/* ---------- static helpers ---------- */
 
/*
 * Returns the unique d >= 0 such that C(d+m-1, m) <= alpha < C(d+m, m).
 *
 * Initial estimate: d0 = floor((m! * alpha)^(1/m)), which satisfies
 * d <= d0 <= d + (m-1). At most m-1 downward corrections are needed.
 */
static pf_int_t diagonal_class(pf_int_t alpha, size_t m)
{
    pf_int_t f = 1;
    for (size_t k = 1; k <= m; k++) f *= (pf_int_t)k;
 
    pf_int_t d = pf_iroot(f * alpha, (pf_int_t)m);
 
    while (pf_binomial(d + (pf_int_t)m - 1, (pf_int_t)m) > alpha) d--;
    while (pf_binomial(d + (pf_int_t)m,     (pf_int_t)m) <= alpha) d++;
 
    return d;
}

/* ---------- direct ---------- */


/**
 * @brief Encodes a 2-D vector using the closed-form Cantor pairing formula.
 */
pf_int_t pf_direct_dp_2d(pf_int_t w1, pf_int_t w2) {
    return w1 + pf_binomial(1 + w1 + w2, 2);
}


/**
 * @brief Encodes an m-dimensional vector as an integer.
 *
 * Maintains a running sum s of all components seen so far. At each position i,
 * it adds C(i-1+s, i) to alpha, which is the contribution of the i-th diagonal
 * layer.
 */
pf_int_t pf_direct_dp(const pf_int_t *W, size_t m) {
    pf_int_t s     = 0;
    pf_int_t alpha = 0;
    for (size_t i = 1; i <= m; i++) {
        s     += W[i - 1];
        alpha += pf_binomial(i - 1 + s, i);
    }
    return alpha;
}


/* ---------- inverse closed-form cases ---------- */

/**
 *
 * Finds the diagonal index d, then derives
 * the two components from d and alpha directly without searching.
 */
void pf_inverse_dp_2d(pf_int_t alpha, pf_int_t *W) {
    pf_int_t d  = (pf_isqrt(1 + 8 * alpha) - 1) / 2;
    pf_int_t w1 = alpha - d * (d + 1) / 2;
    pf_int_t w2 = d - w1;
    W[0] = w1;
    W[1] = w2;
}

/*
 * Depressed cubic t^3 - t - 6*alpha = 0 (after substituting d = t - 1 in
 * d^3 + 3d^2 + 2d - 6*alpha = 0). Discriminant < 0 for alpha >= 1, so there
 * is exactly one real root, given by Cardano:
 *
 *   t = cbrt(3a + sqrt(9a^2 - 1/27)) + cbrt(3a - sqrt(9a^2 - 1/27))
 *   d = t - 1
 *
 * alpha = 0 is handled separately because 9a^2 - 1/27 < 0 there.
 */
void pf_inverse_dp_3d(pf_int_t alpha, pf_int_t *W)
{
    pf_int_t d;
 
    if (alpha == 0) {
        d = 0;
    } else {
        double a     = (double)alpha;
        double disc  = sqrt(9.0*a*a - 1.0/27.0);
        double t     = cbrt(3.0*a + disc) + cbrt(3.0*a - disc);
        d = (pf_int_t)floor(t - 1.0);
 
        /* correct for floating-point rounding (at most 1 step each) */
        while (pf_binomial(d + 2, 3) > alpha) d--;
        while (pf_binomial(d + 3, 3) <= alpha) d++;
    }
 
    pf_int_t alpha2 = alpha - pf_binomial(d + 2, 3);
    pf_inverse_dp_2d(alpha2, W);
    W[2] = d - W[0] - W[1];
}
 
/*
 * The quartic d^4 + 6d^3 + 11d^2 + 6d - 24*alpha = 0 depresses to
 * t^4 - (5/2)t^2 + 9/16 - 24*alpha = 0 under d = t - 3/2. Its roots are
 * symmetric about zero, making it bicuadratic. Setting u = t^2:
 *
 * Simplified: d = (sqrt(5 + 4*sqrt(1 + 24*alpha)) - 3) / 2.
 */
void pf_inverse_dp_4d(pf_int_t alpha, pf_int_t *W)
{
    double a = (double)alpha;
    double t = (sqrt(5.0 + 4.0 * sqrt(1.0 + 24.0 * a)) - 3.0) / 2.0;
    pf_int_t d = (pf_int_t)floor(t);
 
    /* correct for floating-point rounding (at most 1 step each) */
    while (pf_binomial(d + 3, 4) > alpha) d--;
    while (pf_binomial(d + 4, 4) <= alpha) d++;
 
    pf_int_t alpha3 = alpha - pf_binomial(d + 3, 4);
    pf_inverse_dp_3d(alpha3, W);
    W[3] = d - W[0] - W[1] - W[2];
}

/* ---------- inverse: general methods ---------- */
 
/*
 * At each dimension m, finds the diagonal class d via integer root, subtracts
 * the class contribution C(d+m-1, m), recurses on (alpha', m-1), then sets
 * the last component by W[m-1] = d - sum(W[0..m-2]).
 */
void pf_inverse_dp_direct(pf_int_t alpha, size_t m, pf_int_t *W)
{
    if (m == 0) return;
    if (m == 1) { W[0] = alpha; return; }
    if (m == 2) { pf_inverse_dp_2d(alpha, W); return; }
 
    pf_int_t d      = diagonal_class(alpha, m);
    pf_int_t alpha2 = alpha - pf_binomial(d + (pf_int_t)m - 1, (pf_int_t)m);
 
    pf_inverse_dp_direct(alpha2, m - 1, W);
 
    pf_int_t s = 0;
    for (size_t i = 0; i < m - 1; i++) s += W[i];
    W[m - 1] = d - s;
}
 
/**
 * Recovers the GTP representation with pf_inverse_gtp and then converts it to
 * DP with pf_gtp_to_dp. Supports up to 64 dimensions.
 */
void pf_inverse_dp(pf_int_t alpha, size_t m, pf_int_t *W) {
    pf_int_t V[64];
    if (m == 0 || m > 64) return;
    pf_inverse_gtp(alpha, m, V);
    pf_gtp_to_dp(V, m, W);
}
