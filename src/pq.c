/*!
 * \file pq.c
 *
 * \brief An algorithm for quantization of discrete probability distributions.
 *
 *  This is an implementatioon of a quantization algoritm described in:
 *    Y. A. Reznik, "An Algorithm for Quantization of Discrete Probability Distributions,"
 *    Proc. Data Compression Conference (DCC'11), pp. 333-343, Snowbird, UT, March 2011.
 *
 *  Copyright (c) 2011-2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 2.0
 *  \date    March 15, 2026
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "pq.h"

/* we assume that integers are at least 32 bits long */
#if INT_MAX < 0x7FFFFFFF
#error "Integers must be at least 32 bits long"
#endif

/*!
 *  \brief Limits of present implementation.
 */
#define M_MAX   20	/*!< max number of dimensions */
#define N_MAX   33	/*!< max value of common denominator/fidelity of the quantizer */

/***** 
 * Note on limits:
 *  The current choice of N_MAX is such that binomial[N_MAX][N_MAX/2] ~ 30.14 < 31, 
 *  i.e. it allows binomial coefficients to be stored in an int. 
 *  If you need to support larger values of N_MAX, you will need to change the type of 
 *  binomial coefficients to a larger integer type (e.g. long long) 
 *  and adjust the code accordingly.
 *  The same change in data type will also be needed for n_points parameter in pq_points() 
 *  and the idx parametrer in pq_quant() and pq_reconst() functions.
 ****/

/*!
 *  Local data/tables used to support construction.
 */
static int* binomial[N_MAX + 1];                       /*!< Pascal's triangle row pointers */
static int b_data[(N_MAX + 1) * (N_MAX + 2) / 2];      /*!< Pascal's triangle data storage */
static int initialized = 0;                            /*!< flag indicating whether pq_init() has been called */

/*!
 * \brief Initialize quantization module.
 *
 * \return PQ_SUCCESS on success.
 */
int pq_init(void)
{
    int n, k;
    int *b = b_data;
    for (n = 0; n <= N_MAX; n++) {
        binomial[n] = b; b += n + 1;            /* row pointer */
        binomial[n][0] = binomial[n][n] = 1;    /* 1st & last coeffs */
        for (k = 1; k < n; k++)                 /* compute coeffs in the middle */
            binomial[n][k] = binomial[n - 1][k - 1] + binomial[n - 1][k];
    }
    initialized = 1;
    return PQ_SUCCESS;
}

/*!
 * \brief Compute number of reconstruction points.
 *
 * \param[in]   m           the number of dimensions
 * \param[in]   n           precision parameter
 * \param[out]  n_points    pointer to receive the number of reconstruction points
 *
 * \return      PQ_SUCCESS if success, error code otherwise.
 */
int pq_points(int m, int n, int* n_points)
{
    /* check parameters: */
    if (m < 2 || m >= M_MAX) return PQ_INVARG;
    if (n < 1 || n >= N_MAX) return PQ_INVARG;
    if (n + m - 1 >= N_MAX)  return PQ_INVARG;
    if (n_points == NULL) return PQ_INVARG;
    if (!initialized) return PQ_NOTINIT;

    /* return number of types: */
    *n_points = binomial[n + m - 1][m - 1];
    return PQ_SUCCESS;
}

/*! structure for sorting deltas: */
typedef struct { float delta; int i; } delta_t;

/*! callback comparison function for qsort(): */
static int delta_cmp(const void* a, const void* b)
{
    const delta_t* da = (const delta_t*)a;
    const delta_t* db = (const delta_t*)b;
    return (da->delta > db->delta) ? 1 : ((da->delta < db->delta) ? -1 : 0);
}

/* tolerance on distributions: */
#define EPS 1e-6f

/*!
 * \brief Quantize probability distribution.
 *
 * \param[in]   m       the number of dimensions
 * \param[in]   n       precision parameter
 * \param[in]   p       distribution to quantize, p[0..m-1]
 * \param[out]  idx     pointer to receive index of nearest lattice point
 *
 * \return      PQ_SUCCESS if success, error code otherwise
 */
int pq_quantize(int m, int n, const float* p, int* idx)
{
    int k[M_MAX]; 
    delta_t d[M_MAX];
    float ps;
    int i, j, n_rem, n1, s, Delta;

    /* check parameters: */
    if (m < 2 || m >= M_MAX) return PQ_INVARG;
    if (n < 1 || n >= N_MAX) return PQ_INVARG;
    if (n + m - 1 >= N_MAX)  return PQ_INVARG;
    if (p == NULL || idx == NULL) return PQ_INVARG;
    if (!initialized) return PQ_NOTINIT;

    /* check if input distribution is valid: */
    for (i = 0, ps = 0; i < m; i++) {
        if (p[i] < 0.f || p[i] > 1.f) return PQ_INVARG;
        ps += p[i];
    }
    if (fabs(ps - 1.f) > EPS) return PQ_INVARG; 
    
    /* 
     * Quantize distribution p[] to the nearest type k[]/n: 
     */

    /* step 1: quantize to nearest; n1 = total: */
    for (n1 = 0, i = 0; i < m; i++) 
        n1 += (k[i] = (int)floor(p[i] * n + 0.5));     
    
    /* step 2: compute delta: */
    Delta = n1 - n;
    if (Delta != 0) 
    {
        /* step 3: compute and sort errors: */
        for (i = 0; i < m; i++) { d[i].delta = (float)k[i] - p[i] * n; d[i].i = i; } 
        qsort(d, (size_t)m, sizeof(delta_t), delta_cmp);

        /* step 4: adjust k[] for largest of smallest errors such that sum_i k[i] = n */
        if (Delta > 0) 
            { for (j = m - Delta; j < m; j++) k[d[j].i]--; }
        else 
            { for (j = 0; j < abs(Delta); j++) k[d[j].i]++; }
    }

    /* 
     * Compute lexicographic index of type k[]/n: 
     */
    n_rem = n;
    *idx = 0;
    for (i = 0; i < m - 2; i++) {
        s = 0;
        for (j = 0; j < k[i]; j++)
            s += binomial[n_rem - j + m - i - 2][m - i - 2];
        *idx += s;
        n_rem -= k[i];
    }
    *idx += k[m - 2];

    return PQ_SUCCESS;
}

/*!
 * \brief Reconstruct quantized probability distribution.
 *
 * \param[in]   m       number of dimensions
 * \param[in]   n       precision parameter (denominator)
 * \param[in]   idx     index of a reconstruction point
 * \param[out]  p       reconstructed probability distribution, p[0..m-1]
 *
 * \return      PQ_SUCCESS if success, error code otherwise
 */
int pq_reconstruct(int m, int n, int idx, float* p)
{
	int k[M_MAX];
    float n_inv;
    int i, j, s, x;

    /* check parameters: */
    if (m < 2 || m >= M_MAX) return PQ_INVARG;
    if (n < 1 || n >= N_MAX) return PQ_INVARG;
    if (n + m - 1 >= N_MAX)  return PQ_INVARG;
    if (p == NULL) return PQ_INVARG;
    if (!initialized) return PQ_NOTINIT;
    if (idx < 0 || idx >= binomial[n + m - 1][m - 1]) return PQ_INVARG;

    /* compute the inverse of n: */
    n_inv = 1.f / (float)n;

    /* decode type: idx -> k[]/n: */
    for (i = 0; i < m - 2; i++) 
    {
        s = 0;
        for (j = 0; j < n; j++) {
            x = binomial[n - j + m - i - 2][m - i - 2];
            if (idx - s < x) break;
            s += x;
        }
        k[i] = j;
        idx -= s;
        n -= j;
    }
    k[m - 2] = idx;
    k[m - 1] = n - idx;

    /* type to distribution: */
    for (j = 0; j < m; j++)
        p[j] = (float)k[j] * n_inv;

    return PQ_SUCCESS;
}
