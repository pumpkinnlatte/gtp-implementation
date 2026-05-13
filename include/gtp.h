/**
 * @file gtp.h
 * @brief Greater-Than Polynomial packing function.
 *
 * Encodes strictly-increasing integer vectors into unique non-negative integers
 * using the combinatorial number system. Also exposes the bijection between GTP
 * and Diagonal Polynomial vectors, and an iterator over all t-element subsets.
 *
 * Reference: Torres-Jimenez et al., Applied Mathematics and Information
 * Sciences 9(6), 2757-2766, 2015.
 */

#ifndef PACKING_GTP_H
#define PACKING_GTP_H

#include "packing_types.h"

/**
 * @brief Encodes a strictly-increasing vector as an integer.
 *
 * The caller must ensure V[0] < V[1] < ... < V[m-1]; no domain check is
 * performed internally.
 *
 * @param[in] V Array of m non-negative integers in strictly increasing order.
 * @param[in] m Length of V; must be >= 1.
 * @return Unique non-negative integer corresponding to V.
 */
pf_int_t pf_direct_gtp(const pf_int_t *V, size_t m);

/**
 * @brief Decodes an integer into a strictly-increasing vector.
 *
 * Recovers the unique V such that pf_direct_gtp(V, m) == alpha. Each component
 * is approximated with an integer root and then corrected upward one step at a
 * time using pf_next_binomial, avoiding floating-point arithmetic entirely.
 *
 * @param[in]  alpha Integer to decode.
 * @param[in]  m     Number of vector components; must be >= 1 and <= 64.
 * @param[out] V     Output array of length m filled with the recovered vector.
 */
void pf_inverse_gtp(pf_int_t alpha, size_t m, pf_int_t *V);

/**
 * @brief Converts a GTP vector to the equivalent DP vector.
 *
 * W[0] = V[0] and W[i] = V[i] - V[i-1] - 1 for i >= 1.
 * Both vectors map to the same integer under their respective polynomials.
 *
 * @param[in]  V GTP vector of length m.
 * @param[in]  m Length of V; returns immediately if 0.
 * @param[out] W Output DP vector of length m.
 */
void pf_gtp_to_dp(const pf_int_t *V, size_t m, pf_int_t *W);

/**
 * @brief Converts a DP vector to the equivalent GTP vector.
 *
 * V[i] = i + sum of W[0..i].
 * Both vectors map to the same integer under their respective polynomials.
 *
 * @param[in]  W DP vector of length m.
 * @param[in]  m Length of W.
 * @param[out] V Output GTP vector of length m.
 */
void pf_dp_to_gtp(const pf_int_t *W, size_t m, pf_int_t *V);

/**
 * @brief Callback invoked by pf_all_subsets for each enumerated subset.
 *
 * @param[in] V    Current GTP vector representing a t-element subset.
 * @param[in] t    Size of the subset.
 * @param[in] user Caller-supplied context pointer.
 */
typedef void (*pf_subset_cb)(const pf_int_t *V, size_t t, void *user);

/**
 * @brief Iterates over all t-element subsets of {0, ..., k-1} in GTP order.
 *
 * Delivers C(k, t) subsets to cb without computing their GTP rank explicitly.
 * Subsets arrive in the order matching GTP encoding 0, 1, ..., C(k,t)-1.
 *
 * @param[in] k    Upper bound on component values; must satisfy k >= t.
 * @param[in] t    Subset size; must satisfy 1 <= t <= 64.
 * @param[in] cb   Function called once per subset; must not be NULL.
 * @param[in] user Opaque pointer forwarded to every cb invocation.
 */
void pf_all_subsets(pf_int_t k, size_t t, pf_subset_cb cb, void *user);

#endif /* PACKING_GTP_H */
