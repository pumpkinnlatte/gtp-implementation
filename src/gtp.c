/**
 * @file gtp.c
 * @brief Greater-Than Polynomial packing function implementation.
 */

#include "gtp.h"
#include "arith.h"

/**
 * @brief Encodes a strictly-increasing vector as an integer.
 *
 * Accumulates C(V[i-1], i) for each position i from 1 to m.
 */
pf_int_t pf_direct_gtp(const pf_int_t *V, size_t m) {
    pf_int_t alpha = 0;
    for (size_t i = 1; i <= m; i++) {
        alpha += pf_binomial(V[i - 1], i);
    }
    return alpha;
}

/**
 * @brief Decodes an integer into a strictly-increasing vector.
 *
 * Processes components from largest index to smallest. For each position i,
 * the initial candidate v_i is set to floor((i! * alpha)^(1/i)), which is
 * guaranteed to be <= the true value. The candidate is then stepped up by one
 * at a time using pf_next_binomial until the next binomial would exceed the
 * remaining alpha.
 */
void pf_inverse_gtp(pf_int_t alpha, size_t m, pf_int_t *V) {
    pf_int_t f = 1;
    for (size_t k = 1; k <= m; k++) f *= k;

    for (size_t i = m; i >= 1; i--) {
        pf_int_t v_i  = pf_iroot(f * alpha, (pf_int_t)i);
        pf_int_t delta = pf_binomial(v_i, i);

        while (pf_next_binomial(delta, v_i, i) <= alpha) {
            delta = pf_next_binomial(delta, v_i, i);
            v_i++;
        }

        V[i - 1] = v_i;
        alpha -= delta;
        f /= i;

        if (i == 1) break;
    }
}

/**
 * @brief Converts a GTP vector to the equivalent DP vector.
 *
 * Each DP component measures the gap between consecutive GTP components
 * minus one. The first component is copied directly.
 */
void pf_gtp_to_dp(const pf_int_t *V, size_t m, pf_int_t *W) {
    if (m == 0) return;
    W[0] = V[0];
    for (size_t i = 1; i < m; i++) {
        W[i] = V[i] - V[i - 1] - 1;
    }
}

/**
 * @brief Converts a DP vector to the equivalent GTP vector.
 *
 * Accumulates the running sum of DP components and adds the position offset
 * to recover each strictly-increasing GTP value.
 */
void pf_dp_to_gtp(const pf_int_t *W, size_t m, pf_int_t *V) {
    pf_int_t s = 0;
    for (size_t i = 1; i <= m; i++) {
        s += W[i - 1];
        V[i - 1] = (i - 1) + s;
    }
}

/**
 * @brief Iterates over all t-element subsets of {0, ..., k-1} in GTP order.
 *
 * Maintains the current subset in V and advances it lexicographically.
 * The pointer omega tracks the leftmost position that can still be incremented,
 * allowing the loop to reset trailing positions after each carry.
 */
void pf_all_subsets(pf_int_t k, size_t t, pf_subset_cb cb, void *user) {
    if (t == 0 || k < t) return;

    pf_int_t V[64];
    if (t > 64) return;

    for (size_t i = 0; i < t; i++) V[i] = (pf_int_t)i;

    size_t omega = (k == t) ? 1 : t;

    while (V[t - 1] != k || omega != 1) {
        cb(V, t, user);

        V[t - 1] += 1;
        if (V[t - 1] == k && omega != 1) {
            V[omega - 2] += 1;
            for (size_t i = omega - 1; i < t; i++) V[i] = V[i - 1] + 1;
            if (V[omega - 2] == k - t + omega - 2) omega -= 1;
            else                                   omega = t;
        }
    }
}
