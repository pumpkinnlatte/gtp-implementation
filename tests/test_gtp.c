// Tests del modulo GTP.

#include <stdio.h>
#include <inttypes.h>
#include "gtp.h"
#include "dp.h"

static int fails = 0;

#define CHECK(cond, msg) do {                                  \
    if (!(cond)) { printf("  FAIL %s:%d  %s\n",                \
                          __FILE__, __LINE__, msg); fails++; } \
} while (0)

#define CHECK_EQ(a, b, msg) do {                               \
    pf_int_t _a = (a), _b = (b);                               \
    if (_a != _b) {                                            \
        printf("  FAIL %s:%d  %s  got=%" PRIu64                \
               " expected=%" PRIu64 "\n",                      \
               __FILE__, __LINE__, msg, _a, _b);               \
        fails++;                                               \
    }                                                          \
} while (0)

// Tabla 5x5 para GTP. 
static const pf_int_t example_table_gtp[5][5] = {
    {999,   0,   1,   3,   6},
    {999, 999,   2,   4,   7},
    {999, 999, 999,   5,   8},
    {999, 999, 999, 999,   9},
    {999, 999, 999, 999, 999}
};

static void test_direct_gtp_example_table(void) {
    for (pf_int_t v1 = 0; v1 < 5; v1++) {
        for (pf_int_t v2 = v1 + 1; v2 < 5; v2++) {
            pf_int_t V[2] = {v1, v2};
            CHECK_EQ(pf_direct_gtp(V, 2), example_table_gtp[v1][v2],
                     "tabla GTP del paper");
        }
    }
}

static void test_direct_gtp_known(void) {
    // Casos triviales.
    pf_int_t V1[1] = {0};
    CHECK_EQ(pf_direct_gtp(V1, 1), 0, "GTP m=1: V=(0) -> 0");

    pf_int_t V1b[1] = {7};
    CHECK_EQ(pf_direct_gtp(V1b, 1), 7, "GTP m=1: V=(7) -> 7");

    // Vector creciente con valores pequenos.
    pf_int_t V3[3] = {0, 1, 2};
    CHECK_EQ(pf_direct_gtp(V3, 3), 0, "GTP m=3: V=(0,1,2) -> 0");
}

static void test_bijections_gtp_dp(void) {
    // Tomamos un DP arbitrario, lo pasamos a GTP, y de vuelta.
    // El resultado debe coincidir con el original.
    for (pf_int_t a = 0; a < 5; a++) {
        for (pf_int_t b = 0; b < 5; b++) {
            for (pf_int_t c = 0; c < 5; c++) {
                pf_int_t W[3] = {a, b, c};
                pf_int_t V[3], W2[3];
                pf_dp_to_gtp(W, 3, V);
                pf_gtp_to_dp(V, 3, W2);
                CHECK_EQ(W2[0], W[0], "DP->GTP->DP componente 0");
                CHECK_EQ(W2[1], W[1], "DP->GTP->DP componente 1");
                CHECK_EQ(W2[2], W[2], "DP->GTP->DP componente 2");
            }
        }
    }
}

static void test_bijection_preserves_alpha(void) {
    // direct_dp(W) debe ser igual a direct_gtp(dp_to_gtp(W)).
    // Es la propiedad clave de las dos representaciones.
    for (pf_int_t a = 0; a <= 6; a++) {
        for (pf_int_t b = 0; b <= 6; b++) {
            for (pf_int_t c = 0; c <= 6; c++) {
                pf_int_t W[3] = {a, b, c};
                pf_int_t V[3];
                pf_dp_to_gtp(W, 3, V);
                CHECK_EQ(pf_direct_dp(W, 3), pf_direct_gtp(V, 3),
                         "DP y GTP via biyeccion dan el mismo alpha");
            }
        }
    }
}

static void test_inverse_gtp_2d(void) {
    // Para cada celda valida de la tabla, decodificar alpha.
    for (pf_int_t v1 = 0; v1 < 5; v1++) {
        for (pf_int_t v2 = v1 + 1; v2 < 5; v2++) {
            pf_int_t alpha = example_table_gtp[v1][v2];
            pf_int_t V[2];
            pf_inverse_gtp(alpha, 2, V);
            CHECK_EQ(V[0], v1, "inverse_gtp 2D: v1");
            CHECK_EQ(V[1], v2, "inverse_gtp 2D: v2");
        }
    }
}

static void test_bijection_inverse_gtp(void) {
    // Invariante GTP. inverse_gtp(direct_gtp(V)) == V.
    // Recorremos vectores GTP de tamaño 3 con componentes pequenas.
    for (pf_int_t v1 = 0; v1 <= 8; v1++) {
        for (pf_int_t v2 = v1 + 1; v2 <= 9; v2++) {
            for (pf_int_t v3 = v2 + 1; v3 <= 10; v3++) {
                pf_int_t V[3] = {v1, v2, v3};
                pf_int_t alpha = pf_direct_gtp(V, 3);
                pf_int_t Vrec[3];
                pf_inverse_gtp(alpha, 3, Vrec);
                CHECK_EQ(Vrec[0], v1, "inverse_gtp 3D: v1");
                CHECK_EQ(Vrec[1], v2, "inverse_gtp 3D: v2");
                CHECK_EQ(Vrec[2], v3, "inverse_gtp 3D: v3");
            }
        }
    }
}

static void test_inverse_dp_general_matches_2d(void) {
    // El inverso m-dim general debe coincidir con el 2D para m=2.
    for (pf_int_t alpha = 0; alpha < 200; alpha++) {
        pf_int_t W2[2], Wgen[2];
        pf_inverse_dp_2d(alpha, W2);
        pf_inverse_dp(alpha, 2, Wgen);
        CHECK_EQ(Wgen[0], W2[0], "inverse_dp general vs 2D: w1");
        CHECK_EQ(Wgen[1], W2[1], "inverse_dp general vs 2D: w2");
    }
}

static void test_bijection_inverse_dp_higher_dim(void) {
    // Invariante DP. inverse(direct(W)) == W.
    for (pf_int_t a = 0; a <= 4; a++) {
        for (pf_int_t b = 0; b <= 4; b++) {
            for (pf_int_t c = 0; c <= 4; c++) {
                pf_int_t W[3] = {a, b, c};
                pf_int_t alpha = pf_direct_dp(W, 3);
                pf_int_t Wrec[3];
                pf_inverse_dp(alpha, 3, Wrec);
                CHECK_EQ(Wrec[0], a, "inverse_dp 3D: w1");
                CHECK_EQ(Wrec[1], b, "inverse_dp 3D: w2");
                CHECK_EQ(Wrec[2], c, "inverse_dp 3D: w3");
            }
        }
    }
}


int main(void) {
    test_direct_gtp_example_table();
    test_direct_gtp_known();
    test_bijections_gtp_dp();
    test_bijection_preserves_alpha();
    test_inverse_gtp_2d();
    test_bijection_inverse_gtp();
    test_inverse_dp_general_matches_2d();
    test_bijection_inverse_dp_higher_dim();

    if (fails == 0) printf("gtp: OK\n");
    else            printf("gtp: %d FAILS\n", fails);
    return fails == 0 ? 0 : 1;
}
