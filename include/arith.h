#ifndef PACKING_ARITH_H
#define PACKING_ARITH_H

#include "packing_types.h"

// C(k, r). Devuelve 0 si k < r.
pf_int_t pf_binomial(pf_int_t k, pf_int_t r);

// Dado p = C(k, r), devuelve C(k+1, r).
pf_int_t pf_next_binomial(pf_int_t p, pf_int_t k, pf_int_t r);

// Raiz cuadrada entera: floor(sqrt(x)).
pf_int_t pf_isqrt(pf_int_t x);

// Raiz i-esima entera: floor(x^(1/i)).
pf_int_t pf_iroot(pf_int_t x, pf_int_t i);

#endif