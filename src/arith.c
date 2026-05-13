/**
 * @file arith.c
 * @brief Arithmetic for packing algorithms.
 */

#include "arith.h"

/**
 * @brief Computes C(k, r) using the multiplicative formula.
 *
 * Divides after each multiplication to keep intermediate values small and
 * guarantee exact integer division at every step.
 */
pf_int_t pf_binomial(pf_int_t k, pf_int_t r) {
    if (k < r) return 0;

    pf_int_t b = 1;
    for (pf_int_t i = 1; i <= r; i++) {
        b = b * (k - i + 1) / i;
    }
    return b;
}

/**
 * @brief Computes C(k+1, r) from the already-known C(k, r).
 *
 * Applies the identity C(k+1, r) = C(k, r) * (k+1) / (k+1-r).
 */
pf_int_t pf_next_binomial(pf_int_t p, pf_int_t k, pf_int_t r) {
    if (k + 1 <  r) return 0;
    if (k + 1 == r) return 1;
    return p * (k + 1) / (k + 1 - r);
}

/**
 * @brief Computes floor(sqrt(x)) via binary search.
 */
pf_int_t pf_isqrt(pf_int_t x) {
    if (x < 2) return x;

    pf_int_t lo = 1;
    pf_int_t hi = (pf_int_t)1u << 32;

    while (lo < hi) {
        pf_int_t mid = lo + (hi - lo + 1) / 2;
        if (mid <= x / mid) lo = mid;
        else                hi = mid - 1;
    }
    return lo;
}

/**
 * @brief Computes base^i, returning limit+1 on overflow to keep comparisons safe.
 */
static pf_int_t pow_saturating(pf_int_t base, pf_int_t i, pf_int_t limit) {
    pf_int_t r = 1;
    for (pf_int_t k = 0; k < i; k++) {
        if (r > limit / base) return limit + 1;
        r *= base;
    }
    return r;
}

/**
 * @brief Computes floor(x^(1/i)) via binary search.
 *
 * Uses pow_saturating to avoid overflow when testing candidate roots.
 */
pf_int_t pf_iroot(pf_int_t x, pf_int_t i) {
    if (i == 1) return x;
    if (i == 2) return pf_isqrt(x);
    if (x < 2)  return x;

    pf_int_t lo = 1;
    pf_int_t hi = (pf_int_t)1u << ((64 + i - 1) / i);

    while (lo < hi) {
        pf_int_t mid = lo + (hi - lo + 1) / 2;
        if (pow_saturating(mid, i, x) <= x) lo = mid;
        else                                hi = mid - 1;
    }
    return lo;
}
