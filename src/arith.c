#include "arith.h"
// mcd / gcd
static pf_int_t gcd_u64(pf_int_t a, pf_int_t b) {
    while (b != 0) { pf_int_t t = a % b; a = b; b = t; }
    return a;
}

// Algoritmo 2. C(k, r) con r <= k/2.
pf_int_t pf_binomial(pf_int_t k, pf_int_t r) {
    if (k < r){
        return 0;
    } 
    if (r > k - r){
        r = k - r;
    }   

    pf_int_t b = 1;
    for (pf_int_t i = 1; i <= r; i++) {
        pf_int_t num = k - i + 1;
        pf_int_t den = i;

        pf_int_t g1 = gcd_u64(b, den);  b /= g1;  den /= g1;
        pf_int_t g2 = gcd_u64(num, den); num /= g2;
        // tras los dos gcd, den == 1
        b = b * num;
    }
    return b;
}

// Algoritmo 4. C(k+1, r) a partir de C(k, r).
pf_int_t pf_next_binomial(pf_int_t p, pf_int_t k, pf_int_t r) {
    pf_int_t k1 = k + 1;
    if (k1 < r)  return 0;
    if (k1 == r) return 1;
    return p * k1 / (k1 - r);
}

// Busqueda binaria entre 0 y 2^32 (cota: (2^32)^2 = 2^64 desborda).
pf_int_t pf_isqrt(pf_int_t x) {
    if (x < 2) return x;
    pf_int_t lo = 1, hi = (pf_int_t)1u << 32;
    while (lo < hi) {
        pf_int_t mid = lo + (hi - lo + 1) / 2;
        // mid*mid puede desbordar; comparamos como mid <= x/mid.
        if (mid <= x / mid) lo = mid;
        else                hi = mid - 1;
    }
    return lo;
}


// Raiz i-esima entera por busqueda binaria.
pf_int_t pf_iroot(pf_int_t x, pf_int_t i) {
    if (i == 1) return x;
    if (i == 2) return pf_isqrt(x);
    if (x < 2)  return x;

    // Cota superior: 2^ceil(64/i) basta porque (2^(64/i))^i >= 2^64 > x.
    pf_int_t hi = (pf_int_t)1u << ((64 + i - 1) / i);
    pf_int_t lo = 1;

    while (lo < hi) {
        pf_int_t mid = lo + (hi - lo + 1) / 2;
        pf_int_t p = pow_saturating(mid, i, x);
        if (p <= x) lo = mid;
        else        hi = mid - 1;
    }
    return lo;
}

static pf_int_t pow_saturating(pf_int_t base, pf_int_t i, pf_int_t limit) {
    pf_int_t r = 1;
    for (pf_int_t k = 0; k < i; k++) {
        if (r > limit / base + 1) return limit + 1;  // ya no cabe
        r *= base;
        if (r > limit) return limit + 1;
    }
    return r;
}
