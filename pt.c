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
#include <string.h>

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
#define T3AABC(i, j, k) t3a[0*o*o*o+i*o*o+j*o+k]
#define T3ABAC(i, j, k) t3a[1*o*o*o+i*o*o+j*o+k]
#define T3ACBA(i, j, k) t3a[2*o*o*o+i*o*o+j*o+k]
#define T3AACB(i, j, k) t3a[3*o*o*o+i*o*o+j*o+k]
#define T3ABCA(i, j, k) t3a[4*o*o*o+i*o*o+j*o+k]
#define T3ACAB(i, j, k) t3a[5*o*o*o+i*o*o+j*o+k]
#define T3BABC(i, j, k) t3b[0*o*o*o+i*o*o+j*o+k]
#define T3BBAC(i, j, k) t3b[1*o*o*o+i*o*o+j*o+k]
#define T3BCBA(i, j, k) t3b[2*o*o*o+i*o*o+j*o+k]
#define T3BACB(i, j, k) t3b[3*o*o*o+i*o*o+j*o+k]
#define T3BBCA(i, j, k) t3b[4*o*o*o+i*o*o+j*o+k]
#define T3BCAB(i, j, k) t3b[5*o*o*o+i*o*o+j*o+k]
#define MVVV(a, b, c) mvvv[a*v*v+b*v+c]

static void
ccsd_asymm_t3(size_t o, double *t3a)
{
	size_t i, j, k;

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		double x;
		x = T3AABC(i, j, k) -
		    T3ABAC(i, j, k) -
		    T3ACBA(i, j, k) -
		    T3AACB(i, j, k) +
		    T3ABCA(i, j, k) +
		    T3ACAB(i, j, k);
		T3AABC(i, j, k) = x;
	}}}

	for (i = 0; i < o; i++) {
	for (j = i+1; j < o; j++) {
	for (k = j+1; k < o; k++) {
		double x;
		x = T3AABC(i, j, k) -
		    T3AABC(j, i, k) -
		    T3AABC(k, j, i) -
		    T3AABC(i, k, j) +
		    T3AABC(j, k, i) +
		    T3AABC(k, i, j);
		T3AABC(i, j, k) = x;
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
ccsd_t3a(size_t o, size_t v, size_t a, size_t b, size_t c, double *t3a,
    const double *t2, const double *i_ooov, const double *i_ovvv,
    double *work)
{
	//double *mt, *moo, *mov, *mooo, *moov;
	size_t i, j, k, l, d;

	memset(t3a, 0, 6*o*o*o*sizeof(double));

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		for (d = 0; d < v; d++) {
			T3AABC(i, j, k) -=
			    T2(i, j, a, d) * I_OVVV(k, d, b, c);
			T3ACBA(i, j, k) -=
			    T2(i, j, c, d) * I_OVVV(k, d, b, a);
//			T3AACB(i, j, k) -=
//			    T2(i, j, a, d) * I_OVVV(k, d, c, b);
			T3ABAC(i, j, k) -=
			    T2(i, j, b, d) * I_OVVV(k, d, a, c);
//			T3ABCA(i, j, k) -=
//			    T2(i, j, b, d) * I_OVVV(k, d, c, a);
//			T3ACAB(i, j, k) -=
//			    T2(i, j, c, d) * I_OVVV(k, d, a, b);
		}
		for (l = 0; l < o; l++) {
			T3ACAB(i, j, k) -=
			    T2(i, l, a, b) * I_OOOV(j, k, l, c);
			T3AACB(i, j, k) -=
			    T2(i, l, c, b) * I_OOOV(j, k, l, a);
			T3ABCA(i, j, k) -=
			    T2(i, l, a, c) * I_OOOV(j, k, l, b);
//			T3ABAC(i, j, k) -=
//			    T2(i, l, b, a) * I_OOOV(j, k, l, c);
//			T3ABCA(i, j, k) -=
//			    T2(i, l, b, c) * I_OOOV(j, k, l, a);
//			T3ACAB(i, j, k) -=
//			    T2(i, l, c, a) * I_OOOV(j, k, l, b);
		}
	}}}

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		double t3a1abc = T3AABC(i, j, k);
		double t3a1cba = T3ACBA(i, j, k);
		double t3a1bac = T3ABAC(i, j, k);
		double t3a1acb = -t3a1abc;
		double t3a1cab = -t3a1cba;
		double t3a1bca = -t3a1bac;
		double t3a2abc = T3ACAB(i, j, k);
		double t3a2cba = T3AACB(i, j, k);
		double t3a2acb = T3ABCA(i, j, k);
		double t3a2bac = -t3a2abc;
		double t3a2bca = -t3a2cba;
		double t3a2cab = -t3a2acb;
		T3AABC(i, j, k) = t3a1abc + t3a2abc;
		T3ABAC(i, j, k) = t3a1bac + t3a2bac;
		T3ACBA(i, j, k) = t3a1cba + t3a2cba;
		T3AACB(i, j, k) = t3a1acb + t3a2acb;
		T3ABCA(i, j, k) = t3a1bca + t3a2bca;
		T3ACAB(i, j, k) = t3a1cab + t3a2cab;
	}}}

//	moo = work;
//	mov = moo + o*o;
//	mooo = mov + o*v;
//	moov = mooo + o*o*o;
//
//	mt = moov;
//	for (d = 0; d < v; d++) {
//	for (i = 0; i < o; i++) {
//	for (j = 0; j < o; j++) {
//		*mt++ = T2(i, j, a, d);
//	}}}
//	mt = mov;
//	for (k = 0; k < o; k++) {
//	for (d = 0; d < v; d++) {
//		*mt++ = I_OVVV(k, d, c, b);
//	}}
////	memset(mvvv, 0, v*v*v*sizeof(double));
////	for (l = 0; l < i_ovvv->len; l++)
////		if (i_ovvv->idx[l].a == k) {
////			d = i_ovvv->idx[l].b;
////			c = i_ovvv->idx[l].c;
////			b = i_ovvv->idx[l].d;
////			MVVV(b, c, d) = i_ovvv->data[l];
////		}
//	gemm(o*o, o, v, moov, mov, t3a+0*o*o*o);
//
//	mt = moov;
//	for (d = 0; d < v; d++) {
//	for (k = 0; k < o; k++) {
//	for (j = 0; j < o; j++) {
//		*mt++ = T2(k, j, a, d);
//	}}}
//	mt = mov;
//	for (i = 0; i < o; i++) {
//	for (d = 0; d < v; d++) {
//		*mt++ = I_OVVV(i, d, c, b);
//	}}
////	memset(mvvv, 0, v*v*v*sizeof(double));
////	for (l = 0; l < i_ovvv->len; l++)
////		if (i_ovvv->idx[l].a == i) {
////			d = i_ovvv->idx[l].b;
////			c = i_ovvv->idx[l].c;
////			b = i_ovvv->idx[l].d;
////			MVVV(b, c, d) = i_ovvv->data[l];
////		}
//	gemm(o*o, o, v, moov, mov, t3a+1*o*o*o);
//
//	mt = moov;
//	for (d = 0; d < v; d++) {
//	for (i = 0; i < o; i++) {
//	for (k = 0; k < o; k++) {
//		*mt++ = T2(i, k, a, d);
//	}}}
//	mt = mov;
//	for (j = 0; j < o; j++) {
//	for (d = 0; d < v; d++) {
//		*mt++ = I_OVVV(j, d, c, b);
//	}}
////	memset(mvvv, 0, v*v*v*sizeof(double));
////	for (l = 0; l < i_ovvv->len; l++)
////		if (i_ovvv->idx[l].a == j) {
////			d = i_ovvv->idx[l].b;
////			c = i_ovvv->idx[l].c;
////			b = i_ovvv->idx[l].d;
////			MVVV(b, c, d) = i_ovvv->data[l];
////		}
//	gemm(o*o, o, v, moov, mov, t3a+2*o*o*o);
//
//	mt = moo;
//	for (l = 0; l < o; l++) {
//	for (i = 0; i < o; i++) {
//		*mt++ = T2(i, l, a, b);
//	}}
//	mt = mooo;
//	for (j = 0; j < o; j++) {
//	for (k = 0; k < o; k++) {
//	for (l = 0; l < o; l++) {
//		*mt++ = I_OOOV(j, k, l, c);
//	}}}
//	gemm(o, o*o, o, moo, mooo, t3a+3*o*o*o);
//
//	mt = moo;
//	for (l = 0; l < o; l++) {
//	for (j = 0; j < o; j++) {
//		*mt++ = T2(j, l, a, b);
//	}}
//	mt = mooo;
//	for (i = 0; i < o; i++) {
//	for (k = 0; k < o; k++) {
//	for (l = 0; l < o; l++) {
//		*mt++ = I_OOOV(i, k, l, c);
//	}}}
//	gemm(o, o*o, o, moo, mooo, t3a+4*o*o*o);
//
//	mt = moo;
//	for (l = 0; l < o; l++) {
//	for (k = 0; k < o; k++) {
//		*mt++ = T2(k, l, a, b);
//	}}
//	mt = mooo;
//	for (j = 0; j < o; j++) {
//	for (i = 0; i < o; i++) {
//	for (l = 0; l < o; l++) {
//		*mt++ = I_OOOV(j, i, l, c);
//	}}}
//	gemm(o, o*o, o, moo, mooo, t3a+5*o*o*o);
}

static void
ccsd_t3b(size_t o, size_t v, size_t a, size_t b, size_t c,
    double *t3b, const double *t1, const double *t2,
    const double *i_oovv, const double *f_ov)
{
	size_t i, j, k;

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		T3BABC(i, j, k) = T1(i, a) * I_OOVV(j, k, b, c) +
		    F_OV(i, a) * T2(j, k, b, c);
		T3BBAC(i, j, k) = T1(i, b) * I_OOVV(j, k, a, c) +
		    F_OV(i, b) * T2(j, k, a, c);
		T3BCBA(i, j, k) = T1(i, c) * I_OOVV(j, k, b, a) +
		    F_OV(i, c) * T2(j, k, b, a);
		T3BACB(i, j, k) = T1(i, a) * I_OOVV(j, k, c, b) +
		    F_OV(i, a) * T2(j, k, c, b);
		T3BBCA(i, j, k) = T1(i, b) * I_OOVV(j, k, c, a) +
		    F_OV(i, b) * T2(j, k, c, a);
		T3BCAB(i, j, k) = T1(i, c) * I_OOVV(j, k, a, b) +
		    F_OV(i, c) * T2(j, k, a, b);
	}}}
}

static double
ccsd_pt_energy(size_t o, size_t v, size_t a, size_t b, size_t c,
    const double *t3a, const double *t3b, const double *d_ov)
{
	double dn, e_pt = 0.0;
	size_t i, j, k;

	for (i = 0; i < o; i++) {
	for (j = i+1; j < o; j++) {
	for (k = j+1; k < o; k++) {
		dn = D_OV(i, a) + D_OV(j, b) + D_OV(k, c);
		e_pt += T3AABC(i, j, k) * T3BABC(i, j, k) / dn;
	}}}

	return (e_pt);
}

static double
ccsd_pt_worker(int id, int nid, size_t o, size_t v, const double *d_ov,
    const double *f_ov, const double *i_ooov, const double *i_oovv,
    const double *i_ovvv, const double *t1, const double *t2)
{
	double e_pt = 0.0, *t3a, *t3b, *work;
	size_t a, b, c, n;
	int iter;

	t3a = malloc(6*o*o*o*sizeof(double));
	if (t3a == NULL)
		err(1, "malloc");
	t3b = malloc(6*o*o*o*sizeof(double));
	if (t3b == NULL)
		err(1, "malloc");
	work = malloc((o*o + o*v + o*o*o + o*o*v) * sizeof(double));
	if (work == NULL)
		err(1, "malloc");

	for (a = 0, iter = 0; a < v; a++) {
	for (b = a+1; b < v; b++) {
	for (c = b+1; c < v; c++, iter++) {
		if (iter % nid != id)
			continue;
		ccsd_t3a(o, v, a, b, c, t3a, t2, i_ooov, i_ovvv, work);
		ccsd_t3b(o, v, a, b, c, t3b, t1, t2, i_oovv, f_ov);
		ccsd_asymm_t3(o, t3a);
		ccsd_asymm_t3(o, t3b);

		for (n = 0; n < o*o*o; n++) t3b[n] += t3a[n];

		e_pt += ccsd_pt_energy(o, v, a, b, c, t3a, t3b, d_ov);
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
