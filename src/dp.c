/**
 * @file dp.c
 * @brief Diagonal Polynomial packing function implementation.
 */

#include "dp.h"
#include "gtp.h"
#include "arith.h"

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

/**
 * @brief Decodes an integer into a 2-D vector using the closed-form inverse.
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

/**
 * @brief Decodes an integer into an m-dimensional vector via the GTP bijection.
 *
 * Recovers the GTP representation with pf_inverse_gtp and then converts it to
 * DP with pf_gtp_to_dp. Supports up to 64 dimensions.
 */
void pf_inverse_dp(pf_int_t alpha, size_t m, pf_int_t *W) {
    pf_int_t V[64];
    if (m == 0 || m > 64) return;
    pf_inverse_gtp(alpha, m, V);
    pf_gtp_to_dp(V, m, W);
}
