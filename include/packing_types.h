/**
 * @file packing_types.h
 * @brief Shared integer type and status codes for all packing-function modules.
 */

#ifndef PACKING_TYPES_H
#define PACKING_TYPES_H

#include <stdint.h>
#include <stddef.h>

/** @brief Unsigned 64-bit integer used for vector components and encoded values. */
typedef uint64_t pf_int_t;

/** @brief Maximum value representable by pf_int_t. */
#define PF_INT_MAX UINT64_MAX

/** @brief Status codes returned by fallible operations. */
typedef enum {
    PF_OK = 0,            /**< Operation completed successfully. */
    PF_ERR_OVERFLOW,      /**< Result does not fit in pf_int_t. */
    PF_ERR_INVALID_INPUT, /**< Arguments are outside the function domain. */
    PF_ERR_NULL_POINTER,  /**< A required output pointer is NULL. */
    PF_ERR_DOMAIN         /**< Domain constraint violated, e.g. vector not strictly increasing. */
} pf_status_t;

/**
 * @brief Fixed-length vector of non-negative integers.
 *
 * Used as both input and output for the packing functions.
 * The caller is responsible for allocating and freeing data.
 */
typedef struct {
    pf_int_t *data; /**< Array of m components. */
    size_t    m;    /**< Number of components. */
} pf_vector_t;

/**
 * @brief Returns a human-readable string for a status code.
 * @param[in] s Status code to describe.
 * @return Pointer to a static string; never NULL.
 */
const char *pf_status_str(pf_status_t s);

#endif /* PACKING_TYPES_H */
