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
#define T3ACBA(i, j, k) t3a[1*o*o*o+i*o*o+j*o+k]
#define T3AACB(i, j, k) t3a[2*o*o*o+i*o*o+j*o+k]
#define T3BABC(i, j, k) t3b[0*o*o*o+i*o*o+j*o+k]
#define T3BBAC(i, j, k) t3b[1*o*o*o+i*o*o+j*o+k]
#define T3BCBA(i, j, k) t3b[2*o*o*o+i*o*o+j*o+k]
#define OVX(i, j, k) ovx[i*v*x+j*x+k]
#define VVX(i, j, k) vvx[i*v*x+j*x+k]
#define MOV(i, a) mov[i*v+a]
#define MOO1(i, j) moo1[i*o+j]
#define MOO2(i, j) moo2[i*o+j]
#define MOOO(i, j, k) mooo[i*o*o+j*o+k]
#define MVOO(a, i, j) mvoo[a*o*o+i*o+j]
#define MOX(i, j) mox[i*x+j]
#define MXV(i, j) mxv[i*v+j]

static void
ccsd_asymm_t3a(size_t o, double *t3a)
{
	size_t i, j, k;

//	for (i = 0; i < o; i++) {
//	for (j = 0; j < o; j++) {
//	for (k = 0; k < o; k++) {
//		double x;
//		x = T3AABC(i, j, k) -
//		    T3ACBA(i, j, k) -
//		    T3AACB(i, j, k);
//		T3AABC(i, j, k) = x;
//	}}}

	for (i = 0; i < o/2; i++) {
	for (j = i+1; j < o/2; j++) {
	for (k = j+1; k < o; k++) {
		double abcijk, abckji, abcikj;
		abcijk = T3AABC(i, j, k) -
			 T3ACBA(i, j, k) -
			 T3AACB(i, j, k);
		abckji = T3AABC(k, j, i) -
			 T3ACBA(k, j, i) -
			 T3AACB(k, j, i);
		abcikj = T3AABC(i, k, j) -
			 T3ACBA(i, k, j) -
			 T3AACB(i, k, j);
		T3AABC(i, j, k) = abcijk - abckji - abcikj;
//		x = T3AABC(i, j, k) -
//		    T3AABC(k, j, i) -
//		    T3AABC(i, k, j);
//		T3AABC(i, j, k) = x;
	}}}
}

static void
ccsd_asymm_t3b(size_t o, double *t3b)
{
	size_t i, j, k;

//	for (i = 0; i < o; i++) {
//	for (j = 0; j < o; j++) {
//	for (k = 0; k < o; k++) {
//		double x;
//		x = T3BABC(i, j, k) -
//		    T3BBAC(i, j, k) -
//		    T3BCBA(i, j, k);
//		T3BABC(i, j, k) = x;
//	}}}

	for (i = 0; i < o/2; i++) {
	for (j = i+1; j < o/2; j++) {
	for (k = j+1; k < o; k++) {
		double abcijk, abcjik, abckji;
		abcijk = T3BABC(i, j, k) -
			 T3BBAC(i, j, k) -
			 T3BCBA(i, j, k);
		abcjik = T3BABC(j, i, k) -
			 T3BBAC(j, i, k) -
			 T3BCBA(j, i, k);
		abckji = T3BABC(k, j, i) -
			 T3BBAC(k, j, i) -
			 T3BCBA(k, j, i);
		T3BABC(i, j, k) = abcijk - abcjik - abckji;
//		x = T3BABC(i, j, k) -
//		    T3BABC(j, i, k) -
//		    T3BABC(k, j, i);
//		T3BABC(i, j, k) = x;
	}}}
}

void dgemm_(char *, char *, int *, int *, int *, double *, double *,
    int *, double *, int *, double *, double *, int *);

static void
gemm(int m, int n, int k, double alpha, const double *a, const double *b,
    double beta, double *c)
{
	int lda = m;
	int ldb = k;
	int ldc = m;
	char transa = 'N';
	char transb = 'N';

	dgemm_(&transa, &transb, &m, &n, &k, &alpha, (double *)a, &lda,
	    (double *)b, &ldb, &beta, c, &ldc);
}

static void
ccsd_t3a(size_t o, size_t v, size_t x, size_t a, size_t b, size_t c,
    double *t3a, double *t3b, const struct st4 *t2, const struct st4 *i_ooov,
    const struct st4 *i_ovvv, const double *ovx, const double *vvx,
    double *work)
{
	double *moo1, *mov, *mooo, *mvoo, *mox, *mxv;
	size_t i, j, k, l, d;

	/* t3a(i,j,k,a,b,c) = contract(l, t2(i,l,a,b), i_ooov(j,k,l,c)) */
	/* t3b(i,j,k,a,b,c) = contract(d, t2(i,j,d,a), i_ovvv(k,d,b,c)) */

	moo1 = work;
	mooo = moo1 + o*o;

	memset(moo1, 0, o*o*sizeof(double));
	for (l = t2->offset[b]; l < t2->offset[b+1]; l++) {
		if (t2->idx[l].c == a) {
			i = t2->idx[l].a;
			d = t2->idx[l].b;
			MOO1(d, i) = t2->data[l];
		}
	}
	memset(mooo, 0, o*o*o*sizeof(double));
	for (l = i_ooov->offset[c]; l < i_ooov->offset[c+1]; l++) {
		j = i_ooov->idx[l].a;
		k = i_ooov->idx[l].b;
		d = i_ooov->idx[l].c;
		MOOO(j, k, d) = i_ooov->data[l];
	}
	gemm(o, o*o, o, 1.0, moo1, mooo, 0.0, &(T3AABC(0,0,0)));

	memset(moo1, 0, o*o*sizeof(double));
	for (l = t2->offset[b]; l < t2->offset[b+1]; l++) {
		if (t2->idx[l].c == c) {
			j = t2->idx[l].a;
			d = t2->idx[l].b;
			MOO1(d, j) = t2->data[l];
		}
	}
	memset(mooo, 0, o*o*o*sizeof(double));
	for (l = i_ooov->offset[a]; l < i_ooov->offset[a+1]; l++) {
		i = i_ooov->idx[l].a;
		k = i_ooov->idx[l].b;
		d = i_ooov->idx[l].c;
		MOOO(i, k, d) = i_ooov->data[l];
	}
	gemm(o, o*o, o, 1.0, moo1, mooo, 0.0, &(T3ACBA(0,0,0)));

	memset(moo1, 0, o*o*sizeof(double));
	for (l = t2->offset[c]; l < t2->offset[c+1]; l++) {
		if (t2->idx[l].c == a) {
			k = t2->idx[l].a;
			d = t2->idx[l].b;
			MOO1(d, k) = t2->data[l];
		}
	}
	memset(mooo, 0, o*o*o*sizeof(double));
	for (l = i_ooov->offset[b]; l < i_ooov->offset[b+1]; l++) {
		j = i_ooov->idx[l].a;
		i = i_ooov->idx[l].b;
		d = i_ooov->idx[l].c;
		MOOO(j, i, d) = i_ooov->data[l];
	}
	gemm(o, o*o, o, 1.0, moo1, mooo, 0.0, &(T3AACB(0,0,0)));

	mov = work;
	mvoo = mov + o*v;
	mox = mvoo + v*o*o;
	mxv = mox + o*x;

	memset(mvoo, 0, v*o*o*sizeof(double));
	for (l = t2->offset[a]; l < t2->offset[a+1]; l++) {
		i = t2->idx[l].a;
		j = t2->idx[l].b;
		d = t2->idx[l].c;
		MVOO(d, i, j) = t2->data[l];
	}
	if (x == 0) {
		memset(mov, 0, o*v*sizeof(double));
		for (l = i_ovvv->offset[c]; l < i_ovvv->offset[c+1]; l++) {
			if (i_ovvv->idx[l].c == b) {
				k = i_ovvv->idx[l].a;
				d = i_ovvv->idx[l].b;
				MOV(k, d) = i_ovvv->data[l];
			}
		}
	} else {
//		for (k = 0; k < o; k++) {
//		for (d = 0; d < v; d++) {
//		for (l = 0; l < x; l++) {
//			MOV(k, d) += OVX(k, b, l) * VVX(d, c, l) -
//				     OVX(k, c, l) * VVX(d, b, l);
//		}}}

		for (k = 0; k < o; k++) {
		for (l = 0; l < x; l++) {
			MOX(k, l) = OVX(k, b, l);
		}}
		for (l = 0; l < x; l++) {
		for (d = 0; d < v; d++) {
			MXV(l, d) = VVX(d, c, l);
		}}
		gemm(v, o, x, 1.0, mxv, mox, 0.0, mov);
		for (k = 0; k < o; k++) {
		for (l = 0; l < x; l++) {
			MOX(k, l) = OVX(k, c, l);
		}}
		for (l = 0; l < x; l++) {
		for (d = 0; d < v; d++) {
			MXV(l, d) = VVX(d, b, l);
		}}
		gemm(v, o, x, -1.0, mxv, mox, 1.0, mov);
	}
	gemm(o*o, o, v, 1.0, mvoo, mov, 0.0, &(T3BABC(0,0,0)));

	memset(mvoo, 0, v*o*o*sizeof(double));
	for (l = t2->offset[c]; l < t2->offset[c+1]; l++) {
		k = t2->idx[l].a;
		j = t2->idx[l].b;
		d = t2->idx[l].c;
		MVOO(d, k, j) = t2->data[l];
	}
	if (x == 0) {
		memset(mov, 0, o*v*sizeof(double));
		for (l = i_ovvv->offset[a]; l < i_ovvv->offset[a+1]; l++) {
			if (i_ovvv->idx[l].c == b) {
				i = i_ovvv->idx[l].a;
				d = i_ovvv->idx[l].b;
				MOV(i, d) = i_ovvv->data[l];
			}
		}
	} else {
//		for (k = 0; k < o; k++) {
//		for (d = 0; d < v; d++) {
//		for (l = 0; l < x; l++) {
//			MOV(k, d) += OVX(k, b, l) * VVX(d, a, l) -
//				     OVX(k, a, l) * VVX(d, b, l);
//		}}}

		for (k = 0; k < o; k++) {
		for (l = 0; l < x; l++) {
			MOX(k, l) = OVX(k, b, l);
		}}
		for (l = 0; l < x; l++) {
		for (d = 0; d < v; d++) {
			MXV(l, d) = VVX(d, a, l);
		}}
		gemm(v, o, x, 1.0, mxv, mox, 0.0, mov);
		for (k = 0; k < o; k++) {
		for (l = 0; l < x; l++) {
			MOX(k, l) = OVX(k, a, l);
		}}
		for (l = 0; l < x; l++) {
		for (d = 0; d < v; d++) {
			MXV(l, d) = VVX(d, b, l);
		}}
		gemm(v, o, x, -1.0, mxv, mox, 1.0, mov);
	}
	gemm(o*o, o, v, 1.0, mvoo, mov, 0.0, &(T3BCBA(0,0,0)));

	memset(mvoo, 0, v*o*o*sizeof(double));
	for (l = t2->offset[b]; l < t2->offset[b+1]; l++) {
		i = t2->idx[l].a;
		k = t2->idx[l].b;
		d = t2->idx[l].c;
		MVOO(d, i, k) = t2->data[l];
	}
	if (x == 0) {
		memset(mov, 0, o*v*sizeof(double));
		for (l = i_ovvv->offset[c]; l < i_ovvv->offset[c+1]; l++) {
			if (i_ovvv->idx[l].c == a) {
				j = i_ovvv->idx[l].a;
				d = i_ovvv->idx[l].b;
				MOV(j, d) = i_ovvv->data[l];
			}
		}
	} else {
//		for (k = 0; k < o; k++) {
//		for (d = 0; d < v; d++) {
//		for (l = 0; l < x; l++) {
//			MOV(k, d) += OVX(k, a, l) * VVX(d, c, l) -
//				     OVX(k, c, l) * VVX(d, a, l);
//		}}}

		for (k = 0; k < o; k++) {
		for (l = 0; l < x; l++) {
			MOX(k, l) = OVX(k, a, l);
		}}
		for (l = 0; l < x; l++) {
		for (d = 0; d < v; d++) {
			MXV(l, d) = VVX(d, c, l);
		}}
		gemm(v, o, x, 1.0, mxv, mox, 0.0, mov);
		for (k = 0; k < o; k++) {
		for (l = 0; l < x; l++) {
			MOX(k, l) = OVX(k, c, l);
		}}
		for (l = 0; l < x; l++) {
		for (d = 0; d < v; d++) {
			MXV(l, d) = VVX(d, a, l);
		}}
		gemm(v, o, x, -1.0, mxv, mox, 1.0, mov);
	}
	gemm(o*o, o, v, 1.0, mvoo, mov, 0.0, &(T3BBAC(0,0,0)));
}

static void
ccsd_t3b(size_t o, size_t v, size_t a, size_t b, size_t c,
    double *t3b, const double *f_ov, const double *t1,
    const struct st4 *t2, const struct st4 *i_oovv, double *work)
{
	double *moo1, *moo2;
	size_t i, j, k, l;

	/* t3b(i,j,k,a,b,c) = t1(i,a)*i_oovv(j,k,b,c)+f_ov(i,a)*t2(j,k,b,c) */

	moo1 = work;
	moo2 = moo1 + o*o;

	memset(moo1, 0, o*o*sizeof(double));
	for (l = i_oovv->offset[c]; l < i_oovv->offset[c+1]; l++) {
		if (i_oovv->idx[l].c == b) {
			j = i_oovv->idx[l].a;
			k = i_oovv->idx[l].b;
			MOO1(j, k) = i_oovv->data[l];
		}
	}
	memset(moo2, 0, o*o*sizeof(double));
	for (l = t2->offset[c]; l < t2->offset[c+1]; l++) {
		if (t2->idx[l].c == b) {
			j = t2->idx[l].a;
			k = t2->idx[l].b;
			MOO2(j, k) = t2->data[l];
		}
	}
	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		T3BABC(i, j, k) = T1(i, a) * MOO1(j, k) +
		    F_OV(i, a) * MOO2(j, k);
	}}}

	memset(moo1, 0, o*o*sizeof(double));
	for (l = i_oovv->offset[c]; l < i_oovv->offset[c+1]; l++) {
		if (i_oovv->idx[l].c == a) {
			j = i_oovv->idx[l].a;
			k = i_oovv->idx[l].b;
			MOO1(j, k) = i_oovv->data[l];
		}
	}
	memset(moo2, 0, o*o*sizeof(double));
	for (l = t2->offset[c]; l < t2->offset[c+1]; l++) {
		if (t2->idx[l].c == a) {
			j = t2->idx[l].a;
			k = t2->idx[l].b;
			MOO2(j, k) = t2->data[l];
		}
	}
	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		T3BBAC(i, j, k) = T1(i, b) * MOO1(j, k) +
		    F_OV(i, b) * MOO2(j, k);
	}}}

	memset(moo1, 0, o*o*sizeof(double));
	for (l = i_oovv->offset[a]; l < i_oovv->offset[a+1]; l++) {
		if (i_oovv->idx[l].c == b) {
			j = i_oovv->idx[l].a;
			k = i_oovv->idx[l].b;
			MOO1(j, k) = i_oovv->data[l];
		}
	}
	memset(moo2, 0, o*o*sizeof(double));
	for (l = t2->offset[a]; l < t2->offset[a+1]; l++) {
		if (t2->idx[l].c == b) {
			j = t2->idx[l].a;
			k = t2->idx[l].b;
			MOO2(j, k) = t2->data[l];
		}
	}
	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
		T3BCBA(i, j, k) = T1(i, c) * MOO1(j, k) +
		    F_OV(i, c) * MOO2(j, k);
	}}}
}

static double
ccsd_pt_energy(size_t o, size_t v, size_t a, size_t b, size_t c,
    const double *t3a, const double *t3b, const double *d_ov)
{
	double dn, e_pt = 0.0;
	size_t i, j, k;

	for (i = 0; i < o/2; i++) {
	for (j = i+1; j < o/2; j++) {
	for (k = j+1; k < o; k++) {
		dn = D_OV(i, a) + D_OV(j, b) + D_OV(k, c);
		e_pt += T3AABC(i, j, k) * T3BABC(i, j, k) / dn;
	}}}

	return (e_pt);
}

static double
ccsd_pt_worker(int id, int nid, size_t o, size_t v, size_t x,
    const double *d_ov, const double *f_ov, const double *t1,
    const struct st4 *t2, const struct st4 *i_ooov, const struct st4 *i_oovv,
    const struct st4 *i_ovvv, const double *ovx, const double *vvx)
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
	work = malloc((o*n + o*o*n + o*x + v*x) * sizeof(double));
	if (work == NULL)
		err(1, "malloc");

	for (a = 0, iter = 0; a < v/2; a++) {
	for (b = a+1; b < v/2; b++) {
	for (c = b+1; c < v; c++, iter++) {
		if (iter % nid != id)
			continue;
		ccsd_t3a(o, v, x, a, b, c, t3a, t3b, t2, i_ooov, i_ovvv,
		    ovx, vvx, work);
		ccsd_asymm_t3a(o, t3a);
		ccsd_asymm_t3b(o, t3b);

		for (n = 0; n < o*o*o; n++) t3a[n] = t3b[n]-t3a[n];

		ccsd_t3b(o, v, a, b, c, t3b, f_ov, t1, t2, i_oovv, work);
		ccsd_asymm_t3b(o, t3b);

		for (n = 0; n < o*o*o; n++) t3b[n] += t3a[n];

		e_pt += ccsd_pt_energy(o, v, a, b, c, t3a, t3b, d_ov);
	}}}
	e_pt *= 2.0;

	free(t3a);
	free(t3b);
	free(work);
	return (e_pt);
}

static double
do_ccsd_pt(size_t o, size_t v, size_t x, const double *d_ov, const double *f_ov,
    const double *t1, const struct st4 *t2, const struct st4 *i_ooov,
    const struct st4 *i_oovv, const struct st4 *i_ovvv,
    const double *ovx, const double *vvx)
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

		e_pt = ccsd_pt_worker(id, nid, o, v, x, d_ov, f_ov, t1, t2,
		    i_ooov, i_oovv, i_ovvv, ovx, vvx);
	}

	MPI_Allreduce(MPI_IN_PLACE, &e_pt, 1, MPI_DOUBLE, MPI_SUM,
	    MPI_COMM_WORLD);
	return (e_pt);
}

double
ccsd_pt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *t1, const struct st4 *t2, const struct st4 *i_ooov,
    const struct st4 *i_oovv, const struct st4 *i_ovvv)
{
	return (do_ccsd_pt(o, v, 0, d_ov, f_ov, t1, t2,
	    i_ooov, i_oovv, i_ovvv, NULL, NULL));
}

double
ccsd_ri_pt(size_t o, size_t v, size_t x, const double *d_ov,
    const double *f_ov, const double *t1, const struct st4 *t2,
    const struct st4 *i_ooov, const struct st4 *i_oovv,
    const double *ovx, const double *vvx)
{
	return (do_ccsd_pt(o, v, x, d_ov, f_ov, t1, t2,
	    i_ooov, i_oovv, NULL, ovx, vvx));
}
