/**
 * @file bp.c
 * @brief Box Polynomial packing function implementation.
 */

#include "bp.h"
#include "arith.h"

/**
 * @brief Encodes a non-negative integer vector as an integer.
 *
 * First pass finds the maximum component beta and the index iota of its first
 * occurrence (1-based).
 *
 * Second pass computes the class contribution phi1 + phi2 - phi3, where
 *   phi1 = beta^m,  phi2 = (beta+1)^m,  phi3 = beta^(iota-1) * (beta+1)^(m-iota+1).
 * This offset places every vector in the same beta-iota class above all vectors
 * in earlier classes.
 *
 * Third pass interprets the remaining m-1 components as a mixed-radix number
 * with base beta for positions before iota and base beta+1 for positions after.
 */
pf_int_t pf_direct_bp(const pf_int_t *B, size_t m) {
    if (m == 0) return 0;

    pf_int_t beta = B[0];
    size_t iota = 1;
    for (size_t i = 2; i <= m; i++) {
        if (B[i - 1] > beta) {
            beta = B[i - 1];
            iota = i;
        }
    }

    pf_int_t phi1 = 1, phi2 = 1, phi3 = 1;
    for (size_t i = 1; i <= m; i++) {
        phi1 *= beta;
        phi2 *= (beta + 1);
        if (i < iota) phi3 *= beta;
        else          phi3 *= (beta + 1);
    }

    pf_int_t alpha;
    size_t i;
    if (iota == 1) {
        alpha = (m >= 2) ? B[1] : 0;
        i = 3;
    } else {
        alpha = B[0];
        i = 2;
    }

    while (i <= m) {
        if (i < iota) alpha = alpha * beta      + B[i - 1];
        if (i > iota) alpha = alpha * (beta + 1) + B[i - 1];
        i++;
    }

    alpha += phi1 + phi2 - phi3;
    return alpha;
}

/**
 * @brief Decodes an integer into a non-negative integer vector.
 *
 * Estimates beta as the integer m-th root of alpha, which gives the smallest
 * beta such that beta^m <= alpha. Then steps iota upward from 1 while the
 * class contribution for iota is still less than or equal to alpha, identifying
 * the correct class. After subtracting the class offset, the remaining value is
 * decoded as a mixed-radix number with the bases determined by iota.
 */
void pf_inverse_bp(pf_int_t alpha, size_t m, pf_int_t *B) {
    if (m == 0) return;

    pf_int_t beta = pf_iroot(alpha, (pf_int_t)m);

    pf_int_t phi1 = 1, phi2 = 1;
    for (size_t i = 1; i <= m; i++) {
        phi1 *= beta;
        phi2 *= (beta + 1);
    }

    pf_int_t phi3 = phi2;
    size_t iota = 1;
    while (beta > 0 && phi1 + phi2 - phi3 * beta / (beta + 1) <= alpha) {
        iota++;
        phi3 = phi3 * beta / (beta + 1);
    }

    alpha -= (phi1 + phi2 - phi3);
    B[iota - 1] = beta;

    for (size_t i = m; i >= 1; i--) {
        if (i > iota) {
            B[i - 1] = alpha % (beta + 1);
            alpha   /= (beta + 1);
        } else if (i < iota) {
            B[i - 1] = alpha % beta;
            alpha   /= beta;
        }
        if (i == 1) break;
    }
}
