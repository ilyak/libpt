/*
 * Copyright (c) 2016 Ilya Kaliman
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <err.h>
#include <stdlib.h>

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "pt.h"

#define D_OV(i, a) d_ov[i*v+a]
#define F_OV(i, a) f_ov[i*v+a]
#define I_OOOV(i, j, k, a) i_ooov[i*o*o*v+j*o*v+k*v+a]
#define I_OOVV(i, j, a, b) i_oovv[i*o*v*v+j*v*v+a*v+b]
#define I_OVVV(i, a, b, c) i_ovvv[i*v*v*v+a*v*v+b*v+c]
#define T1(i, a) t1[i*v+a]
#define T2(i, j, a, b) t2[i*o*v*v+j*v*v+a*v+b]
#define T3AIJK(a, b, c) t3a[0*v*v*v+a*v*v+b*v+c]
#define T3AJIK(a, b, c) t3a[1*v*v*v+a*v*v+b*v+c]
#define T3AKJI(a, b, c) t3a[2*v*v*v+a*v*v+b*v+c]
#define T3AIKJ(a, b, c) t3a[3*v*v*v+a*v*v+b*v+c]
#define T3AJKI(a, b, c) t3a[4*v*v*v+a*v*v+b*v+c]
#define T3AKIJ(a, b, c) t3a[5*v*v*v+a*v*v+b*v+c]
#define T3BIJK(a, b, c) t3b[0*v*v*v+a*v*v+b*v+c]
#define T3BJIK(a, b, c) t3b[1*v*v*v+a*v*v+b*v+c]
#define T3BKJI(a, b, c) t3b[2*v*v*v+a*v*v+b*v+c]
#define T3BIKJ(a, b, c) t3b[3*v*v*v+a*v*v+b*v+c]
#define T3BJKI(a, b, c) t3b[4*v*v*v+a*v*v+b*v+c]
#define T3BKIJ(a, b, c) t3b[5*v*v*v+a*v*v+b*v+c]

static void
ccsd_asymm_t3(size_t v, double *t3a)
{
	size_t a, b, c;

	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		double x;
		x = T3AIJK(a, b, c) -
		    T3AJIK(a, b, c) -
		    T3AKJI(a, b, c) -
		    T3AIKJ(a, b, c) +
		    T3AJKI(a, b, c) +
		    T3AKIJ(a, b, c);
		T3AIJK(a, b, c) = x;
	}}}

	for (a = 0; a < v; a++) {
	for (b = a+1; b < v; b++) {
	for (c = b+1; c < v; c++) {
		double x;
		x = T3AIJK(a, b, c) -
		    T3AIJK(b, a, c) -
		    T3AIJK(c, b, a) -
		    T3AIJK(a, c, b) +
		    T3AIJK(b, c, a) +
		    T3AIJK(c, a, b);
		T3AIJK(a, b, c) = x;
	}}}
}

void dgemm_(char *, char *, int *, int *, int *, double *, double *,
    int *, double *, int *, double *, double *, int *);

static void
gemm(int m, int n, int k, const double *a, const double *b, double *c)
{
	double alpha = 1.0;
	double beta = 0.0;
	int lda = m;
	int ldb = k;
	int ldc = m;
	char transa = 'N';
	char transb = 'N';

	dgemm_(&transa, &transb, &m, &n, &k, &alpha, (double *)a, &lda,
	    (double *)b, &ldb, &beta, c, &ldc);
}

static void
ccsd_t3a(size_t o, size_t v, size_t i, size_t j, size_t k, double *t3a,
    const double *t2, const double *i_ooov, const double *i_ovvv,
    double *work)
{
	double *mt, *mov, *mvv, *movv, *mvvv;
	size_t l, a, b, c, d;

	mov = work;
	mvv = mov + o*v;
	movv = mvv + v*v;
	mvvv = movv + o*v*v;

	mt = mvv;
	for (d = 0; d < v; d++) {
	for (a = 0; a < v; a++) {
		*mt++ = T2(i, j, a, d);
	}}
	mt = mvvv;
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
	for (d = 0; d < v; d++) {
		*mt++ = I_OVVV(k, d, c, b);
	}}}
	gemm(v, v*v, v, mvv, mvvv, t3a+0*v*v*v);

	mt = mvv;
	for (d = 0; d < v; d++) {
	for (a = 0; a < v; a++) {
		*mt++ = T2(k, j, a, d);
	}}
	mt = mvvv;
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
	for (d = 0; d < v; d++) {
		*mt++ = I_OVVV(i, d, c, b);
	}}}
	gemm(v, v*v, v, mvv, mvvv, t3a+1*v*v*v);

	mt = mvv;
	for (d = 0; d < v; d++) {
	for (a = 0; a < v; a++) {
		*mt++ = T2(i, k, a, d);
	}}
	mt = mvvv;
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
	for (d = 0; d < v; d++) {
		*mt++ = I_OVVV(j, d, c, b);
	}}}
	gemm(v, v*v, v, mvv, mvvv, t3a+2*v*v*v);

	mt = movv;
	for (l = 0; l < o; l++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
		*mt++ = T2(i, l, a, b);
	}}}
	mt = mov;
	for (c = 0; c < v; c++) {
	for (l = 0; l < o; l++) {
		*mt++ = I_OOOV(j, k, l, c);
	}}
	gemm(v*v, v, o, movv, mov, t3a+3*v*v*v);

	mt = movv;
	for (l = 0; l < o; l++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
		*mt++ = T2(j, l, a, b);
	}}}
	mt = mov;
	for (c = 0; c < v; c++) {
	for (l = 0; l < o; l++) {
		*mt++ = I_OOOV(i, k, l, c);
	}}
	gemm(v*v, v, o, movv, mov, t3a+4*v*v*v);

	mt = movv;
	for (l = 0; l < o; l++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
		*mt++ = T2(k, l, a, b);
	}}}
	mt = mov;
	for (c = 0; c < v; c++) {
	for (l = 0; l < o; l++) {
		*mt++ = I_OOOV(j, i, l, c);
	}}
	gemm(v*v, v, o, movv, mov, t3a+5*v*v*v);

	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		double t3a1ijk = T3AIJK(a, b, c);
		double t3a1kji = T3AKJI(a, b, c);
		double t3a1ikj = T3AIKJ(a, b, c);
		double t3a1jik = -t3a1ijk;
		double t3a1jki = -t3a1kji;
		double t3a1kij = -t3a1ikj;
		double t3a2ijk = T3AJIK(a, b, c);
		double t3a2jik = T3AJKI(a, b, c);
		double t3a2kji = T3AKIJ(a, b, c);
		double t3a2ikj = -t3a2ijk;
		double t3a2jki = -t3a2jik;
		double t3a2kij = -t3a2kji;
		T3AIJK(a, b, c) = t3a1ijk - t3a2ijk;
		T3AJIK(a, b, c) = t3a1jik - t3a2jik;
		T3AKJI(a, b, c) = t3a1kji - t3a2kji;
		T3AIKJ(a, b, c) = t3a1ikj - t3a2ikj;
		T3AJKI(a, b, c) = t3a1jki - t3a2jki;
		T3AKIJ(a, b, c) = t3a1kij - t3a2kij;
	}}}
}

static void
ccsd_t3b(size_t o, size_t v, size_t i, size_t j, size_t k,
    double *t3b, const double *t1, const double *t2,
    const double *i_oovv, const double *f_ov)
{
	size_t a, b, c;

	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		T3BIJK(a, b, c) = T1(i, a) * I_OOVV(j, k, b, c) +
		    F_OV(i, a) * T2(j, k, b, c);
		T3BJIK(a, b, c) = T1(j, a) * I_OOVV(i, k, b, c) +
		    F_OV(j, a) * T2(i, k, b, c);
		T3BKJI(a, b, c) = T1(k, a) * I_OOVV(j, i, b, c) +
		    F_OV(k, a) * T2(j, i, b, c);
		T3BIKJ(a, b, c) = T1(i, a) * I_OOVV(k, j, b, c) +
		    F_OV(i, a) * T2(k, j, b, c);
		T3BJKI(a, b, c) = T1(j, a) * I_OOVV(k, i, b, c) +
		    F_OV(j, a) * T2(k, i, b, c);
		T3BKIJ(a, b, c) = T1(k, a) * I_OOVV(i, j, b, c) +
		    F_OV(k, a) * T2(i, j, b, c);
	}}}
}

static double
ccsd_pt_energy(size_t v, size_t i, size_t j, size_t k,
    const double *t3a, const double *t3b, const double *d_ov)
{
	double dn, e_pt = 0.0;
	size_t a, b, c;

	for (a = 0; a < v; a++) {
	for (b = a+1; b < v; b++) {
	for (c = b+1; c < v; c++) {
		dn = D_OV(i, a) + D_OV(j, b) + D_OV(k, c);
		e_pt += T3AIJK(a, b, c) * T3BIJK(a, b, c) / dn;
	}}}

	return (e_pt);
}

static double
ccsd_pt_worker(int id, int nid, size_t o, size_t v, const double *d_ov,
    const double *f_ov, const double *i_ooov, const double *i_oovv,
    const double *i_ovvv, const double *t1, const double *t2)
{
	double e_pt = 0.0, *t3a, *t3b, *work;
	size_t i, j, k, n;
	int iter;

	t3a = malloc(6 * v*v*v * sizeof(double));
	if (t3a == NULL)
		err(1, "malloc");
	t3b = malloc(6 * v*v*v * sizeof(double));
	if (t3b == NULL)
		err(1, "malloc");
	work = malloc((o*v + v*v + o*v*v + v*v*v) * sizeof(double));
	if (work == NULL)
		err(1, "malloc");

	for (i = 0, iter = 0; i < o; i++) {
	for (j = i+1; j < o; j++) {
	for (k = j+1; k < o; k++, iter++) {
		if (iter % nid != id)
			continue;
		ccsd_t3a(o, v, i, j, k, t3a, t2, i_ooov, i_ovvv, work);
		ccsd_t3b(o, v, i, j, k, t3b, t1, t2, i_oovv, f_ov);
		ccsd_asymm_t3(v, t3a);
		ccsd_asymm_t3(v, t3b);

		for (n = 0; n < v*v*v; n++) t3b[n] += t3a[n];

		e_pt += ccsd_pt_energy(v, i, j, k, t3a, t3b, d_ov);
	}}}
	e_pt *= (1.0 / 16.0);

	free(t3a);
	free(t3b);
	free(work);
	return (e_pt);
}

double
ccsd_pt(size_t o, size_t v, const double *d_ov,
    const double *f_ov, const double *i_ooov, const double *i_oovv,
    const double *i_ovvv, const double *t1, const double *t2)
{
	double e_pt = 0.0;
	int pid, npid;

	if (o == 0 || v == 0)
		return (0.0);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &npid);

#ifdef _OPENMP
#pragma omp parallel reduction(+:e_pt)
#endif
	{
		int id, nid, tid = 0, ntid = 1;
#ifdef _OPENMP
		tid = omp_get_thread_num();
		ntid = omp_get_num_threads();
#endif
		id = pid * ntid + tid;
		nid = npid * ntid;

		e_pt = ccsd_pt_worker(id, nid, o, v, d_ov, f_ov, i_ooov,
		    i_oovv, i_ovvv, t1, t2);
	}

	MPI_Allreduce(MPI_IN_PLACE, &e_pt, 1, MPI_DOUBLE, MPI_SUM,
	    MPI_COMM_WORLD);
	return (e_pt);
}
