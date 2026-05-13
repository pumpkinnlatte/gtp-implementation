/**
 * @file arith.h
 * @brief Arithmetic primitives shared by all packing-function algorithms.
 *
 * Provides exact binomial-coefficient evaluation and integer root extraction
 * using only integer arithmetic, no floating-point operations.
 */

#ifndef PACKING_ARITH_H
#define PACKING_ARITH_H

#include "packing_types.h"

/**
 * @brief Computes the binomial coefficient C(k, r).
 *
 * Uses the multiplicative formula so every intermediate division is exact.
 * Returns 0 when k < r.
 *
 * @param[in] k Total number of items.
 * @param[in] r Number of items to choose.
 * @return C(k, r), or 0 if k < r.
 */
pf_int_t pf_binomial(pf_int_t k, pf_int_t r);

/**
 * @brief Advances a binomial coefficient by one: computes C(k+1, r) from C(k, r).
 *
 * Avoids recomputing from scratch in loops that increment k one at a time.
 * Uses one multiplication and one division.
 *
 * @param[in] p Current value C(k, r).
 * @param[in] k Current lower parameter.
 * @param[in] r Upper parameter, fixed across calls.
 * @return C(k+1, r), or 0 if k+1 < r.
 */
pf_int_t pf_next_binomial(pf_int_t p, pf_int_t k, pf_int_t r);

/**
 * @brief Computes the integer square root of x via binary search.
 *
 * @param[in] x Non-negative integer.
 * @return Largest integer s such that s*s <= x.
 */
pf_int_t pf_isqrt(pf_int_t x);

/**
 * @brief Computes the integer i-th root of x via binary search.
 *
 * Used by inverse packing functions to produce the initial approximation
 * for each vector component. Delegates to pf_isqrt when i equals 2.
 *
 * @param[in] x Non-negative integer.
 * @param[in] i Root degree; must be >= 1.
 * @return Largest integer r such that r^i <= x.
 */
pf_int_t pf_iroot(pf_int_t x, pf_int_t i);

#endif /* PACKING_ARITH_H */
