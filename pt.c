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
//#define I_OOOV(i, j, k, a) i_ooov[i*o*o*v+j*o*v+k*v+a]
//#define I_OOVV(i, j, a, b) i_oovv[i*o*v*v+j*v*v+a*v+b]
//#define I_OVVV(i, a, b, c) i_ovvv[i*v*v*v+a*v*v+b*v+c]
#define T1(i, a) t1[i*v+a]
//#define T2(i, j, a, b) t2[i*o*v*v+j*v*v+a*v+b]
#define T3AABC(i, j, k) t3a[0*o*o*o+i*o*o+j*o+k]
//#define T3ABAC(i, j, k) t3a[1*o*o*o+i*o*o+j*o+k]
#define T3ACBA(i, j, k) t3a[1*o*o*o+i*o*o+j*o+k]
#define T3AACB(i, j, k) t3a[2*o*o*o+i*o*o+j*o+k]
//#define T3ABCA(i, j, k) t3a[4*o*o*o+i*o*o+j*o+k]
//#define T3ACAB(i, j, k) t3a[5*o*o*o+i*o*o+j*o+k]
#define T3BABC(i, j, k) t3b[0*o*o*o+i*o*o+j*o+k]
#define T3BBAC(i, j, k) t3b[1*o*o*o+i*o*o+j*o+k]
#define T3BCBA(i, j, k) t3b[2*o*o*o+i*o*o+j*o+k]
//#define T3BACB(i, j, k) t3b[3*o*o*o+i*o*o+j*o+k]
//#define T3BBCA(i, j, k) t3b[4*o*o*o+i*o*o+j*o+k]
//#define T3BCAB(i, j, k) t3b[5*o*o*o+i*o*o+j*o+k]
#define MOV(i, a) mov[i*v+a]
#define MOO1(i, j) moo1[i*o+j]
#define MOO2(i, j) moo2[i*o+j]
#define MOOO(i, j, k) mooo[i*o*o+j*o+k]
#define MVOO(a, i, j) mvoo[a*o*o+i*o+j]

static void
ccsd_asymm_t3a(size_t o, double *t3a)
{
	size_t i, j, k;

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		double x;
		x = T3AABC(i, j, k) -
		    T3ACBA(i, j, k) -
		    T3AACB(i, j, k);
		T3AABC(i, j, k) = 2.0 * x;
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

static void
ccsd_asymm_t3b(size_t o, double *t3b)
{
	size_t i, j, k;

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		double x;
		x = T3BABC(i, j, k) -
		    T3BBAC(i, j, k) -
		    T3BCBA(i, j, k);
		T3BABC(i, j, k) = 2.0 * x;
	}}}

	for (i = 0; i < o; i++) {
	for (j = i+1; j < o; j++) {
	for (k = j+1; k < o; k++) {
		double x;
		x = T3BABC(i, j, k) -
		    T3BABC(j, i, k) -
		    T3BABC(k, j, i) -
		    T3BABC(i, k, j) +
		    T3BABC(j, k, i) +
		    T3BABC(k, i, j);
		T3BABC(i, j, k) = x;
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
    double *t3b, const struct st4 *t2, const struct st4 *i_ooov,
    const struct st4 *i_ovvv, double *work)
{
	double *moo1, *mov, *mooo, *mvoo;
	size_t i, j, k, l, d;

	/*memset(t3a, 0, 6*o*o*o*sizeof(double));
	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		for (d = 0; d < v; d++) {
			T3AABC(i, j, k) +=
			    T2(i, j, a, d) * I_OVVV(k, d, b, c);
			T3ACBA(i, j, k) +=
			    T2(i, j, c, d) * I_OVVV(k, d, b, a);
//			T3AACB(i, j, k) +=
//			    T2(i, j, a, d) * I_OVVV(k, d, c, b);
			T3ABAC(i, j, k) +=
			    T2(i, j, b, d) * I_OVVV(k, d, a, c);
//			T3ABCA(i, j, k) +=
//			    T2(i, j, b, d) * I_OVVV(k, d, c, a);
//			T3ACAB(i, j, k) +=
//			    T2(i, j, c, d) * I_OVVV(k, d, a, b);
		}
		for (l = 0; l < o; l++) {
			T3ACAB(i, j, k) +=
			    T2(i, l, a, b) * I_OOOV(j, k, l, c);
			T3AACB(i, j, k) +=
			    T2(i, l, c, b) * I_OOOV(j, k, l, a);
			T3ABCA(i, j, k) +=
			    T2(i, l, a, c) * I_OOOV(j, k, l, b);
//			T3ABAC(i, j, k) +=
//			    T2(i, l, b, a) * I_OOOV(j, k, l, c);
//			T3ABCA(i, j, k) +=
//			    T2(i, l, b, c) * I_OOOV(j, k, l, a);
//			T3ACAB(i, j, k) +=
//			    T2(i, l, c, a) * I_OOOV(j, k, l, b);
		}
	}}}*/

	mov = work;
	mvoo = mov + o*v;

//	mt = moov;
//	for (d = 0; d < v; d++) {
//	for (i = 0; i < o; i++) {
//	for (j = 0; j < o; j++) {
//		*mt++ = T2(i, j, a, d);
//	}}}
//	mt = mov;
//	for (k = 0; k < o; k++) {
//	for (d = 0; d < v; d++) {
//		*mt++ = I_OVVV(k, d, b, c);
//	}}
	memset(mvoo, 0, v*o*o*sizeof(double));
	for (l = 0; l < t2->len; l++)
		if (t2->idx[l].d == a) {
			i = t2->idx[l].a;
			j = t2->idx[l].b;
			d = t2->idx[l].c;
			MVOO(d, i, j) = t2->data[l];
		}
	memset(mov, 0, o*v*sizeof(double));
	for (l = 0; l < i_ovvv->len; l++)
		if (i_ovvv->idx[l].c == b && i_ovvv->idx[l].d == c) {
			k = i_ovvv->idx[l].a;
			d = i_ovvv->idx[l].b;
			MOV(k, d) = i_ovvv->data[l];
		}
	gemm(o*o, o, v, mvoo, mov, &(T3BABC(0,0,0)));

//	mt = moov;
//	for (d = 0; d < v; d++) {
//	for (k = 0; k < o; k++) {
//	for (j = 0; j < o; j++) {
//		*mt++ = T2(k, j, c, d);
//	}}}
//	mt = mov;
//	for (i = 0; i < o; i++) {
//	for (d = 0; d < v; d++) {
//		*mt++ = I_OVVV(i, d, b, a);
//	}}
	memset(mvoo, 0, v*o*o*sizeof(double));
	for (l = 0; l < t2->len; l++)
		if (t2->idx[l].d == c) {
			k = t2->idx[l].a;
			j = t2->idx[l].b;
			d = t2->idx[l].c;
			MVOO(d, k, j) = t2->data[l];
		}
	memset(mov, 0, o*v*sizeof(double));
	for (l = 0; l < i_ovvv->len; l++)
		if (i_ovvv->idx[l].c == b && i_ovvv->idx[l].d == a) {
			i = i_ovvv->idx[l].a;
			d = i_ovvv->idx[l].b;
			MOV(i, d) = i_ovvv->data[l];
		}
	gemm(o*o, o, v, mvoo, mov, &(T3BCBA(0,0,0)));

//	mt = moov;
//	for (d = 0; d < v; d++) {
//	for (i = 0; i < o; i++) {
//	for (k = 0; k < o; k++) {
//		*mt++ = T2(i, k, b, d);
//	}}}
//	mt = mov;
//	for (j = 0; j < o; j++) {
//	for (d = 0; d < v; d++) {
//		*mt++ = I_OVVV(j, d, a, c);
//	}}
	memset(mvoo, 0, v*o*o*sizeof(double));
	for (l = 0; l < t2->len; l++)
		if (t2->idx[l].d == b) {
			i = t2->idx[l].a;
			k = t2->idx[l].b;
			d = t2->idx[l].c;
			MVOO(d, i, k) = t2->data[l];
		}
	memset(mov, 0, o*v*sizeof(double));
	for (l = 0; l < i_ovvv->len; l++)
		if (i_ovvv->idx[l].c == a && i_ovvv->idx[l].d == c) {
			j = i_ovvv->idx[l].a;
			d = i_ovvv->idx[l].b;
			MOV(j, d) = i_ovvv->data[l];
		}
	gemm(o*o, o, v, mvoo, mov, &(T3BBAC(0,0,0)));

	moo1 = work;
	mooo = moo1 + o*o;

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
	memset(moo1, 0, o*o*sizeof(double));
	for (l = 0; l < t2->len; l++)
		if (t2->idx[l].c == a && t2->idx[l].d == b) {
			i = t2->idx[l].a;
			d = t2->idx[l].b;
			MOO1(d, i) = t2->data[l];
		}
	memset(mooo, 0, o*o*o*sizeof(double));
	for (l = 0; l < i_ooov->len; l++)
		if (i_ooov->idx[l].d == c) {
			j = i_ooov->idx[l].a;
			k = i_ooov->idx[l].b;
			d = i_ooov->idx[l].c;
			MOOO(j, k, d) = i_ooov->data[l];
		}
	gemm(o, o*o, o, moo1, mooo, &(T3AABC(0,0,0)));

//	mt = moo;
//	for (l = 0; l < o; l++) {
//	for (j = 0; j < o; j++) {
//		*mt++ = T2(j, l, c, b);
//	}}
//	mt = mooo;
//	for (i = 0; i < o; i++) {
//	for (k = 0; k < o; k++) {
//	for (l = 0; l < o; l++) {
//		*mt++ = I_OOOV(i, k, l, a);
//	}}}
	memset(moo1, 0, o*o*sizeof(double));
	for (l = 0; l < t2->len; l++)
		if (t2->idx[l].c == c && t2->idx[l].d == b) {
			j = t2->idx[l].a;
			d = t2->idx[l].b;
			MOO1(d, j) = t2->data[l];
		}
	memset(mooo, 0, o*o*o*sizeof(double));
	for (l = 0; l < i_ooov->len; l++)
		if (i_ooov->idx[l].d == a) {
			i = i_ooov->idx[l].a;
			k = i_ooov->idx[l].b;
			d = i_ooov->idx[l].c;
			MOOO(i, k, d) = i_ooov->data[l];
		}
	gemm(o, o*o, o, moo1, mooo, &(T3ACBA(0,0,0)));

//	mt = moo;
//	for (l = 0; l < o; l++) {
//	for (k = 0; k < o; k++) {
//		*mt++ = T2(k, l, a, c);
//	}}
//	mt = mooo;
//	for (j = 0; j < o; j++) {
//	for (i = 0; i < o; i++) {
//	for (l = 0; l < o; l++) {
//		*mt++ = I_OOOV(j, i, l, b);
//	}}}
	memset(moo1, 0, o*o*sizeof(double));
	for (l = 0; l < t2->len; l++)
		if (t2->idx[l].c == a && t2->idx[l].d == c) {
			k = t2->idx[l].a;
			d = t2->idx[l].b;
			MOO1(d, k) = t2->data[l];
		}
	memset(mooo, 0, o*o*o*sizeof(double));
	for (l = 0; l < i_ooov->len; l++)
		if (i_ooov->idx[l].d == b) {
			j = i_ooov->idx[l].a;
			i = i_ooov->idx[l].b;
			d = i_ooov->idx[l].c;
			MOOO(j, i, d) = i_ooov->data[l];
		}
	gemm(o, o*o, o, moo1, mooo, &(T3AACB(0,0,0)));

//	for (i = 0; i < o; i++) {
//	for (j = 0; j < o; j++) {
//	for (k = 0; k < o; k++) {
//		double t3a1abc = T3AABC(i, j, k);
//		double t3a1cba = T3ACBA(i, j, k);
//		double t3a1bac = T3ABAC(i, j, k);
//		double t3a1acb = -t3a1abc;
//		double t3a1cab = -t3a1cba;
//		double t3a1bca = -t3a1bac;
//		double t3a2abc = T3ACAB(i, j, k);
//		double t3a2cba = T3AACB(i, j, k);
//		double t3a2acb = T3ABCA(i, j, k);
//		double t3a2bac = -t3a2abc;
//		double t3a2bca = -t3a2cba;
//		double t3a2cab = -t3a2acb;
//		T3AABC(i, j, k) = -(t3a1abc + t3a2abc);
//		T3ABAC(i, j, k) = -(t3a1bac + t3a2bac);
//		T3ACBA(i, j, k) = -(t3a1cba + t3a2cba);
//		T3AACB(i, j, k) = -(t3a1acb + t3a2acb);
//		T3ABCA(i, j, k) = -(t3a1bca + t3a2bca);
//		T3ACAB(i, j, k) = -(t3a1cab + t3a2cab);
//	}}}
}

static void
ccsd_t3b(size_t o, size_t v, size_t a, size_t b, size_t c,
    double *t3b, const double *f_ov, const double *t1,
    const struct st4 *t2, const struct st4 *i_oovv, double *work)
{
	double *moo1, *moo2;
	size_t i, j, k, l;

	moo1 = work;
	moo2 = moo1 + o*o;

	memset(moo1, 0, o*o*sizeof(double));
	for (l = 0; l < i_oovv->len; l++)
		if (i_oovv->idx[l].c == b && i_oovv->idx[l].d == c) {
			j = i_oovv->idx[l].a;
			k = i_oovv->idx[l].b;
			MOO1(j, k) = i_oovv->data[l];
		}
	memset(moo2, 0, o*o*sizeof(double));
	for (l = 0; l < t2->len; l++)
		if (t2->idx[l].c == b && t2->idx[l].d == c) {
			j = t2->idx[l].a;
			k = t2->idx[l].b;
			MOO2(j, k) = t2->data[l];
		}
	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		T3BABC(i, j, k) = T1(i, a) * MOO1(j, k) +
		    F_OV(i, a) * MOO2(j, k);
//		T3BACB(i, j, k) = -T3BABC(i, j, k);
	}}}

	memset(moo1, 0, o*o*sizeof(double));
	for (l = 0; l < i_oovv->len; l++)
		if (i_oovv->idx[l].c == a && i_oovv->idx[l].d == c) {
			j = i_oovv->idx[l].a;
			k = i_oovv->idx[l].b;
			MOO1(j, k) = i_oovv->data[l];
		}
	memset(moo2, 0, o*o*sizeof(double));
	for (l = 0; l < t2->len; l++)
		if (t2->idx[l].c == a && t2->idx[l].d == c) {
			j = t2->idx[l].a;
			k = t2->idx[l].b;
			MOO2(j, k) = t2->data[l];
		}
	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		T3BBAC(i, j, k) = T1(i, b) * MOO1(j, k) +
		    F_OV(i, b) * MOO2(j, k);
//		T3BBCA(i, j, k) = -T3BBAC(i, j, k);
	}}}

	memset(moo1, 0, o*o*sizeof(double));
	for (l = 0; l < i_oovv->len; l++)
		if (i_oovv->idx[l].c == b && i_oovv->idx[l].d == a) {
			j = i_oovv->idx[l].a;
			k = i_oovv->idx[l].b;
			MOO1(j, k) = i_oovv->data[l];
		}
	memset(moo2, 0, o*o*sizeof(double));
	for (l = 0; l < t2->len; l++)
		if (t2->idx[l].c == b && t2->idx[l].d == a) {
			j = t2->idx[l].a;
			k = t2->idx[l].b;
			MOO2(j, k) = t2->data[l];
		}
	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		T3BCBA(i, j, k) = T1(i, c) * MOO1(j, k) +
		    F_OV(i, c) * MOO2(j, k);
//		T3BCAB(i, j, k) = -T3BCBA(i, j, k);
	}}}

//	for (i = 0; i < o; i++) {
//	for (j = 0; j < o; j++) {
//	for (k = 0; k < o; k++) {
//		T3BABC(i, j, k) = T1(i, a) * I_OOVV(j, k, b, c) +
//		    F_OV(i, a) * T2(j, k, b, c);
//		T3BBAC(i, j, k) = T1(i, b) * I_OOVV(j, k, a, c) +
//		    F_OV(i, b) * T2(j, k, a, c);
//		T3BCBA(i, j, k) = T1(i, c) * I_OOVV(j, k, b, a) +
//		    F_OV(i, c) * T2(j, k, b, a);
////		T3BACB(i, j, k) = T1(i, a) * I_OOVV(j, k, c, b) +
////		    F_OV(i, a) * T2(j, k, c, b);
////		T3BBCA(i, j, k) = T1(i, b) * I_OOVV(j, k, c, a) +
////		    F_OV(i, b) * T2(j, k, c, a);
////		T3BCAB(i, j, k) = T1(i, c) * I_OOVV(j, k, a, b) +
////		    F_OV(i, c) * T2(j, k, a, b);
//	}}}
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
    const double *f_ov, const double *t1, const struct st4 *t2,
    const struct st4 *i_ooov, const struct st4 *i_oovv,
    const struct st4 *i_ovvv)
{
	double e_pt = 0.0, *t3a, *t3b, *work;
	size_t a, b, c, n, iter;

	t3a = malloc(3*o*o*o*sizeof(double));
	if (t3a == NULL)
		err(1, "malloc");
	t3b = malloc(3*o*o*o*sizeof(double));
	if (t3b == NULL)
		err(1, "malloc");
	n = o > v ? o : v;
	work = malloc((o*n + o*o*n) * sizeof(double));
	if (work == NULL)
		err(1, "malloc");

	for (a = 0, iter = 0; a < v; a++) {
	for (b = a+1; b < v; b++) {
	for (c = b+1; c < v; c++, iter++) {
		if (iter % nid != id)
			continue;
		ccsd_t3a(o, v, a, b, c, t3a, t3b, t2, i_ooov, i_ovvv, work);
		ccsd_asymm_t3a(o, t3a);
		ccsd_asymm_t3b(o, t3b);

		for (n = 0; n < o*o*o; n++) t3a[n] = t3b[n]-t3a[n];

		ccsd_t3b(o, v, a, b, c, t3b, f_ov, t1, t2, i_oovv, work);
		ccsd_asymm_t3b(o, t3b);

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
ccsd_pt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *t1, const struct st4 *t2, const struct st4 *i_ooov,
    const struct st4 *i_oovv, const struct st4 *i_ovvv)
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

		e_pt = ccsd_pt_worker(id, nid, o, v, d_ov, f_ov, t1, t2,
		    i_ooov, i_oovv, i_ovvv);
	}

	MPI_Allreduce(MPI_IN_PLACE, &e_pt, 1, MPI_DOUBLE, MPI_SUM,
	    MPI_COMM_WORLD);
	return (e_pt);
}
