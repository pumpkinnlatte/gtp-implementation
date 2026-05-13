/**
 * @file dp.h
 * @brief Diagonal Polynomial packing function.
 *
 * Encodes unrestricted non-negative integer vectors into unique non-negative
 * integers using the diagonal enumeration scheme. The 2-D case corresponds to
 * Cantor's classical pairing function. The general inverse is computed via the
 * GTP bijection.
 *
 * Reference: Torres-Jimenez et al., Applied Mathematics and Information
 * Sciences 9(6), 2757-2766, 2015.
 */

#ifndef PACKING_DP_H
#define PACKING_DP_H

#include "packing_types.h"

/**
 * @brief Encodes an unrestricted non-negative vector as an integer.
 *
 * @param[in] W Array of m non-negative integers.
 * @param[in] m Length of W; must be >= 1.
 * @return Unique non-negative integer corresponding to W.
 */
pf_int_t pf_direct_dp(const pf_int_t *W, size_t m);

/**
 * @brief Encodes a 2-D vector as an integer using the closed-form expression.
 *
 * Faster than pf_direct_dp for the 2-D case because it avoids the loop.
 *
 * @param[in] w1 First component.
 * @param[in] w2 Second component.
 * @return Unique non-negative integer corresponding to the pair.
 */
pf_int_t pf_direct_dp_2d(pf_int_t w1, pf_int_t w2);

/**
 * @brief Decodes an integer into a 2-D vector using the closed-form inverse.
 *
 * Computes the diagonal index d = floor((sqrt(1 + 8*alpha) - 1) / 2),
 * then derives both components from d and alpha directly.
 *
 * @param[in]  alpha Integer to decode.
 * @param[out] W     Output array of length 2.
 */
void pf_inverse_dp_2d(pf_int_t alpha, pf_int_t *W);

/**
 * @brief Decodes an integer into an m-dimensional vector via the GTP bijection.
 *
 * Calls pf_inverse_gtp to obtain the GTP representation, then converts it to
 * DP with pf_gtp_to_dp. Supports up to 64 dimensions.
 *
 * @param[in]  alpha Integer to decode.
 * @param[in]  m     Number of vector components; must be >= 1 and <= 64.
 * @param[out] W     Output array of length m filled with the recovered vector.
 */
void pf_inverse_dp(pf_int_t alpha, size_t m, pf_int_t *W);

#endif /* PACKING_DP_H */
