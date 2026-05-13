/**
 * @file bp.h
 * @brief Box Polynomial packing function.
 *
 * Encodes unrestricted non-negative integer vectors by classifying them
 * according to their maximum component and the index where that maximum first
 * appears. This classification yields a mixed-radix encoding that is a factor
 * of m faster than GTP and DP for both direct and inverse computations.
 *
 * Reference: Torres-Jimenez et al., Applied Mathematics and Information
 * Sciences 9(6), 2757-2766, 2015.
 */

#ifndef PACKING_BP_H
#define PACKING_BP_H

#include "packing_types.h"

/**
 * @brief Encodes a non-negative integer vector as an integer.
 *
 * Identifies the maximum component beta and its first position iota, computes
 * the class contribution from beta and iota, then encodes the remaining
 * components in mixed radix base beta / base beta+1.
 *
 * @param[in] B Array of m non-negative integers.
 * @param[in] m Length of B; must be >= 1.
 * @return Unique non-negative integer corresponding to B.
 */
pf_int_t pf_direct_bp(const pf_int_t *B, size_t m);

/**
 * @brief Decodes an integer into a non-negative integer vector.
 *
 * Estimates beta with the integer m-th root of alpha, then finds iota by
 * stepping through class boundaries. The remaining components are recovered
 * by mixed-radix division using bases beta and beta+1.
 *
 * @param[in]  alpha Integer to decode.
 * @param[in]  m     Number of vector components; must be >= 1.
 * @param[out] B     Output array of length m filled with the recovered vector.
 */
void pf_inverse_bp(pf_int_t alpha, size_t m, pf_int_t *B);

#endif /* PACKING_BP_H */
