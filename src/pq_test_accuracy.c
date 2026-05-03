/*!
 *  \file pq_test_accuracy.c
 * 
 *  \brief Test program for pq (probability simplex quantization) library.
 *
 *  This program measures the worst-case L1 quantization error for various lattice parameters 
 *  and compares it against the known expression for L1 covering radius of such lattice.
 *  For details, see Proposition 2 in the paper: 
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
#include <math.h>
#include <float.h>
#include "pq.h"

#define M_MAX 20

/*!
 *  \brief L1 distance between two distributions.
 *
 *  \param[in] m     number of dimensions
 *  \param[in] p1    first distribution 
 *  \param[in] p2    second distribution
 *
 *  \return L1 distance.
 */
static double l1_dist(int m, float* p1, float* p2)
{
	double s = 0.; int i;
	for (i = 0; i < m; i++) s += fabs(p1[i] - p2[i]);
	return s;
}

/*!
 *  \brief Test program:
 */
int main(void)
{
	float p[M_MAX], q[M_MAX];   /* original & quantized distributions */
	int m, m_max = 6;			/* the number of dimensions  */
	int n, n_max = 27;			/* lattice accuracy parameters */
	int a, i, j, t, t_max, err;
	double e, e1, r;

	/* init library: */
	pq_init();

	/* print header: */
	printf("Testing quantization lattices:\n");

	/* try different dimensions: */
	for (m = 2; m < m_max; m++)
	{
		/* compute # of points in the reference (fine-resolution lattice): */
		err = pq_points(m, n_max, &t_max);
		if (err != PQ_SUCCESS) {
			fprintf(stderr, "pq_points() failed with error code %d\n", err);
			return 1;
		}

		/* try different resolutions of lattices: */
		for (n = 1; n < n_max; n++)
		{
			/* report # of points & rate: */
			err = pq_points(m, n, &t);
			if (err != PQ_SUCCESS) {
				fprintf(stderr, "pq_points() failed with error code %d\n", err);
				return 1;
			}
			printf("m=%d, n=%d, t=%d, rate=%g, ", m, n, t, log((double)t) / log(2.));

			/* scan points in the reference lattice (n = n_max): */
			e = 0.;
			for (i = 0; i < t_max; i++)
			{
				/* retrieve reference distribution: */
				err = pq_reconstruct(m, n_max, i, p);
				if (err != PQ_SUCCESS) {
					fprintf(stderr, "pq_reconstruct() failed with error code %d\n", err);
					return 1;
				}

				/* quantize it to lattice with a given parameter n: */
				j = 0;
				err = pq_quantize(m, n, p, &j);
				if (err != PQ_SUCCESS) {
					fprintf(stderr, "p_quantize() failed with error code %d\n", err);
					return 1;
				}
				/* reconstruct quantized distribution: */
				err = pq_reconstruct(m, n, j, q);
				if (err != PQ_SUCCESS) {
					fprintf(stderr, "pq_reconstruct() failed with error code %d\n", err);
					return 1;
				}

				/* compute quantization error: */
				e1 = l1_dist(m, p, q);
				if (e1 > e) e = e1;
			}

			/* compute L1 covering radius: */
			a = m / 2;
			r = (double)(2 * a * (m - a)) / (double)(m * n);

			/* report measured quantization error and compare it with L1 radius: */
			printf("L1 error = %g vs. %g\n", e, r);

			/* check if measured error exceeds L1 covering radius: */
			if (e > (1.0 + 6*FLT_EPSILON) * r) {
				printf("Test failed: measured error exceeds L1 covering radius!\n");
				return 1;
			}
		}
	}

	/* print final report and exit: */
	printf("Test completed!\n");

	return 0;
}
