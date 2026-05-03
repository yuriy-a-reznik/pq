/*!
 * \file pq.h
 * 
 * \brief An algorithm for quantization of discrete probability distributions.
 *
 *  This is an implementation of a quantization algorithm described in:
 *    Y. A. Reznik, "An Algorithm for Quantization of Discrete Probability Distributions," 
 *    Proc. Data Compression Conference (DCC'11), pp. 333-343, Snowbird, UT, March 2011.
 * 
 *  Main functions: 
 * 
 *    pq_init()        - initialize the quantization module (must be called before any other function)
 *    pq_points()      - compute the number of reconstruction points for given parameters m and n
 *    pq_quantize()    - quantize a given probability distribution
 *    pq_reconstruct() - reconstruct a probability distribution from its quantized index.
 * 
 *  Copyright (c) 2011-2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 2.0
 *  \date    March 15, 2026
 */

#ifndef PQ_H_
#define PQ_H_ 1

#ifdef __cplusplus
extern "C" {
#endif

/*! error codes: */
enum {
	PQ_SUCCESS = 0,     /*!< operation completed successfully */
	PQ_NOTINIT = 1,		/*!< module not initialized */
	PQ_INVARG = 2       /*!< invalid argument (e.g. NULL pointer or invalid range) */
};

/* functions: */
extern int pq_init(void);
extern int pq_points(int m, int n, int* n_points);
extern int pq_quantize(int m, int n, const float* p, int* idx);
extern int pq_reconstruct(int m, int n, int idx, float* p);

#ifdef __cplusplus
}
#endif

#endif /* PQ_H_ */