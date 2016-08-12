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
#include <stdio.h> /*XXX*/

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "pt.h"

#define D_OV(i, a) d_ov[i*v+a]
#define F_OV(i, a) f_ov[i*v+a]
#define I_OOOV(i, j, k, a) i_ooov[i*o*o*v+j*o*v+k*v+a]
#define I_OOVO(i, j, a, k) i_oovo[i*o*o*v+j*o*v+a*o+k]
#define I_OOVV(i, j, a, b) i_oovv[i*o*v*v+j*v*v+a*v+b]
#define I_OVVV(i, a, b, c) i_ovvv[i*v*v*v+a*v*v+b*v+c]
#define I_VVOV(b, c, i, a) i_vvov[b*v*o*v+c*o*v+i*v+a]
#define T1(i, a) t1[i*v+a]
#define T2(i, j, a, b) t2[i*o*v*v+j*v*v+a*v+b]
#define T2T(a, b, i, j) t2t[a*v*o*o+b*o*o+i*o+j]
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

void dgemm_(char *, char *, int *, int *, int *, double *, double *,
    int *, double *, int *, double *, double *, int *);

static void
gemm(char transa, char transb, int m, int n, int k, double alpha,
    const double *a, int lda, const double *b, int ldb, double beta,
    double *c, int ldc)
{
	dgemm_(&transa, &transb, &m, &n, &k, &alpha, (double *)a, &lda,
	    (double *)b, &ldb, &beta, c, &ldc);
}

static double
comp_t3a_ijkabc_11(size_t o, size_t v, size_t i, size_t j, size_t k,
    size_t a, size_t b, size_t c, const double *t2, const double *i_vvov)
{
	double s = 0.0;
	const double *t2_p = &T2(i,j,a,0);
	const double *i_vvov_p = &I_VVOV(b,c,k,0);
	size_t l;

	/* t3b(i,j,k,a,b,c) = contract(d, t2(i,j,d,a), i_ovvv(k,d,b,c)) */

	for (l = 0; l < v/2; l++)
		s += t2_p[l] * i_vvov_p[l];

	return (s);
}

static double
comp_t3a_ijkabc_11h(size_t o, size_t v, size_t i, size_t j, size_t k,
    size_t a, size_t b, size_t c, const double *t2, const double *i_vvov)
{
	double s = 0.0;
	const double *t2_p = &T2(i,j,a,0);
	const double *i_vvov_p = &I_VVOV(b,c,k,0);
	size_t l;

	/* t3b(i,j,k,a,b,c) = contract(d, t2(i,j,d,a), i_ovvv(k,d,b,c)) */

	for (l = 0; l < v; l++)
		s += t2_p[l] * i_vvov_p[l];

	return (s);
}

static void
comp_t3a_ijk_1(size_t o, size_t v, size_t a, size_t b, size_t c,
    double *ijk, const double *t2, const double *i_vvov)
{
	double s = 0.0;
	const double *t2_p = &T2(0,0,a,0);
	const double *i_vvov_p = &I_VVOV(b,c,0,0);
//	size_t l;

	/* t3b(i,j,k,a,b,c) = contract(d, t2(i,j,d,a), i_ovvv(k,d,b,c)) */

	int lda = v*v;
	int ldb = v;
	gemm('T', 'N', o*o, o, v, 1.0, t2_p, lda, i_vvov_p, ldb, 0.0, ijk, o*o);

//	for (l = 0; l < v; l++)
//		s += t2_p[l] * i_vvov_p[l];

//	return (s);
}

static void
comp_t3a_ijk_2(size_t o, size_t v, size_t a, size_t b, size_t c,
    double *ijk, const double *t2t, const double *i_oovo)
{
	double s = 0.0;
	const double *t2t_p = &T2T(a,b,0,0);
	const double *i_oovo_p = &I_OOVO(0,0,c,0);
//	size_t l;

	/* t3a(i,j,k,a,b,c) = contract(l, t2(i,l,a,b), i_ooov(j,k,l,c)) */

	int lda = o;
	int ldb = o*v;
	gemm('T', 'N', o, o*o, o, 1.0, t2t_p, lda, i_oovo_p, ldb, 0.0, ijk, o);

//	for (l = 0; l < o; l++)
//		s += t2t_p[l] * i_oovo_p[l];

//	return (s);
}

static double
comp_t3a_ijkabc_12(size_t o, size_t v, size_t i, size_t j, size_t k,
    size_t a, size_t b, size_t c, const double *t2, const double *i_vvov)
{
	double s = 0.0;
	const double *t2_p = &T2(i,j,a,0);
	const double *i_vvov_p = &I_VVOV(b,c,k,0);
	size_t l;

	/* t3b(i,j,k,a,b,c) = contract(d, t2(i,j,d,a), i_ovvv(k,d,b,c)) */

	for (l = v/2; l < v; l++)
		s += t2_p[l] * i_vvov_p[l];

	return (s);
}

static double
comp_t3a_ijkabc_21(size_t o, size_t v, size_t i, size_t j, size_t k,
    size_t a, size_t b, size_t c, const double *t2t, const double *i_oovo)
{
	double s = 0.0;
	const double *t2t_p = &T2T(a,b,i,0);
	const double *i_oovo_p = &I_OOVO(j,k,c,0);
	size_t l;

	/* t3a(i,j,k,a,b,c) = contract(l, t2(i,l,a,b), i_ooov(j,k,l,c)) */

	for (l = 0; l < o/2; l++)
		s += t2t_p[l] * i_oovo_p[l];

	return (s);
}

static double
comp_t3a_ijkabc_21h(size_t o, size_t v, size_t i, size_t j, size_t k,
    size_t a, size_t b, size_t c, const double *t2t, const double *i_oovo)
{
	double s = 0.0;
	const double *t2t_p = &T2T(a,b,i,0);
	const double *i_oovo_p = &I_OOVO(j,k,c,0);
	size_t l;

	/* t3a(i,j,k,a,b,c) = contract(l, t2(i,l,a,b), i_ooov(j,k,l,c)) */

	for (l = 0; l < o; l++)
		s += t2t_p[l] * i_oovo_p[l];

	return (s);
}

static double
comp_t3a_ijkabc_22(size_t o, size_t v, size_t i, size_t j, size_t k,
    size_t a, size_t b, size_t c, const double *t2t, const double *i_oovo)
{
	double s = 0.0;
	const double *t2t_p = &T2T(a,b,i,0);
	const double *i_oovo_p = &I_OOVO(j,k,c,0);
	size_t l;

	/* t3a(i,j,k,a,b,c) = contract(l, t2(i,l,a,b), i_ooov(j,k,l,c)) */

	for (l = o/2; l < o; l++)
		s += t2t_p[l] * i_oovo_p[l];

	return (s);
}

static double
comp_t3b_ijkabc(size_t o, size_t v, size_t i, size_t j, size_t k,
    size_t a, size_t b, size_t c, const double *t1, const double *i_oovv,
    const double *f_ov, const double *t2)
{
	return T1(i,a)*I_OOVV(j,k,b,c) + F_OV(i,a)*T2(j,k,b,c);
}

static double
ccsd_pt_energy(size_t o, size_t v,
    const double *d_ov, const double *f_ov, const double *t1,
    const double *t2, const double *t2t, const double *i_oovo,
    const double *i_oovv, const double *i_vvov)
{
	double e_pt1 = 0.0, e_pt2 = 0.0;

#pragma omp parallel for reduction(+:e_pt1) schedule(dynamic)
	for (size_t a = 0; a < v/2; a++) {
	for (size_t b = a+1; b < v/2; b++) {
	for (size_t c = b+1; c < v/2; c++) {
	for (size_t i = 0; i < o/2; i++) {
	for (size_t j = i+1; j < o/2; j++) {
	for (size_t k = j+1; k < o/2; k++) {
		double t3ax1, t3ax2, t3bx, dn;
		double ijkabc, jikabc, kjiabc, ikjabc;
		ijkabc = comp_t3a_ijkabc_11(o,v,i,j,k,a,b,c,t2,i_vvov) -
			 comp_t3a_ijkabc_11(o,v,i,j,k,b,a,c,t2,i_vvov) -
			 comp_t3a_ijkabc_11(o,v,i,j,k,c,b,a,t2,i_vvov);
		kjiabc = comp_t3a_ijkabc_11(o,v,k,j,i,a,b,c,t2,i_vvov) -
			 comp_t3a_ijkabc_11(o,v,k,j,i,b,a,c,t2,i_vvov) -
			 comp_t3a_ijkabc_11(o,v,k,j,i,c,b,a,t2,i_vvov);
		ikjabc = comp_t3a_ijkabc_11(o,v,i,k,j,a,b,c,t2,i_vvov) -
			 comp_t3a_ijkabc_11(o,v,i,k,j,b,a,c,t2,i_vvov) -
			 comp_t3a_ijkabc_11(o,v,i,k,j,c,b,a,t2,i_vvov);
		t3ax1 = ijkabc - kjiabc - ikjabc;

		ijkabc = comp_t3a_ijkabc_21(o,v,i,j,k,a,b,c,t2t,i_oovo) -
			 comp_t3a_ijkabc_21(o,v,i,j,k,c,b,a,t2t,i_oovo) -
			 comp_t3a_ijkabc_21(o,v,i,j,k,a,c,b,t2t,i_oovo);
		jikabc = comp_t3a_ijkabc_21(o,v,j,i,k,a,b,c,t2t,i_oovo) -
			 comp_t3a_ijkabc_21(o,v,j,i,k,c,b,a,t2t,i_oovo) -
			 comp_t3a_ijkabc_21(o,v,j,i,k,a,c,b,t2t,i_oovo);
		kjiabc = comp_t3a_ijkabc_21(o,v,k,j,i,a,b,c,t2t,i_oovo) -
			 comp_t3a_ijkabc_21(o,v,k,j,i,c,b,a,t2t,i_oovo) -
			 comp_t3a_ijkabc_21(o,v,k,j,i,a,c,b,t2t,i_oovo);
		t3ax2 = ijkabc - jikabc - kjiabc;

		ijkabc = comp_t3b_ijkabc(o,v,i,j,k,a,b,c,t1,i_oovv,f_ov,t2) -
			 comp_t3b_ijkabc(o,v,i,j,k,b,a,c,t1,i_oovv,f_ov,t2) -
			 comp_t3b_ijkabc(o,v,i,j,k,c,b,a,t1,i_oovv,f_ov,t2);
		jikabc = comp_t3b_ijkabc(o,v,j,i,k,a,b,c,t1,i_oovv,f_ov,t2) -
			 comp_t3b_ijkabc(o,v,j,i,k,b,a,c,t1,i_oovv,f_ov,t2) -
			 comp_t3b_ijkabc(o,v,j,i,k,c,b,a,t1,i_oovv,f_ov,t2);
		kjiabc = comp_t3b_ijkabc(o,v,k,j,i,a,b,c,t1,i_oovv,f_ov,t2) -
			 comp_t3b_ijkabc(o,v,k,j,i,b,a,c,t1,i_oovv,f_ov,t2) -
			 comp_t3b_ijkabc(o,v,k,j,i,c,b,a,t1,i_oovv,f_ov,t2);
		t3bx = ijkabc - jikabc - kjiabc;

		dn = D_OV(i, a) + D_OV(j, b) + D_OV(k, c);
		e_pt1 += (t3ax1+t3ax2) * (t3ax1+t3ax2-t3bx) / dn;
	}}}}}}

	e_pt1 *= 2.0;
	printf("aaaaaa %g\n", e_pt1);

	double *t2_aaaa, *t2_abab;//, *t2_abba;
	double *t2t_aaaa, *t2t_abab;//, *t2t_abba;
	double *i_vvov_aaaa, *i_vvov_abab;//, *i_vvov_abba;
	double *i_oovo_aaaa, *i_oovo_abab;//, *i_oovo_abba;

	t2_aaaa = malloc(o/2*o/2*v/2*v/2*sizeof(double));
	t2_abab = malloc(o/2*o/2*v/2*v/2*sizeof(double));
//	t2_abba = malloc(o/2*o/2*v/2*v/2*sizeof(double));
	t2t_aaaa = malloc(o/2*o/2*v/2*v/2*sizeof(double));
	t2t_abab = malloc(o/2*o/2*v/2*v/2*sizeof(double));
//	t2t_abba = malloc(o/2*o/2*v/2*v/2*sizeof(double));
	i_vvov_aaaa = malloc(v/2*v/2*o/2*v/2*sizeof(double));
	i_vvov_abab = malloc(v/2*v/2*o/2*v/2*sizeof(double));
//	i_vvov_abba = malloc(v/2*v/2*o/2*v/2*sizeof(double));
	i_oovo_aaaa = malloc(o/2*o/2*v/2*o/2*sizeof(double));
	i_oovo_abab = malloc(o/2*o/2*v/2*o/2*sizeof(double));
//	i_oovo_abba = malloc(o/2*o/2*v/2*o/2*sizeof(double));

	for (size_t i = 0; i < o/2; i++) {
	for (size_t j = 0; j < o/2; j++) {
	for (size_t a = 0; a < v/2; a++) {
	for (size_t b = 0; b < v/2; b++) {
		t2_aaaa[i*o/2*v/2*v/2+j*v/2*v/2+a*v/2+b] =
		    T2(i,j,a,b);
		t2_abab[i*o/2*v/2*v/2+j*v/2*v/2+a*v/2+b] =
		    T2(i,(j+o/2),a,(b+v/2));
//		t2_abba[i*o/2*v/2*v/2+j*v/2*v/2+a*v/2+b] =
//		    T2(i,(j+o/2),(a+v/2),b);
		t2t_aaaa[a*v/2*o/2*o/2+b*o/2*o/2+i*o/2+j] =
		    T2T(a,b,i,j);
		t2t_abab[a*v/2*o/2*o/2+b*o/2*o/2+i*o/2+j] =
		    T2T(a,(b+v/2),i,(j+o/2));
//		t2t_abba[a*v/2*o/2*o/2+b*o/2*o/2+i*o/2+j] =
//		    T2T(a,(b+v/2),(i+o/2),j);
	}}}}
	for (size_t i = 0; i < o/2; i++) {
	for (size_t a = 0; a < v/2; a++) {
	for (size_t b = 0; b < v/2; b++) {
	for (size_t c = 0; c < v/2; c++) {
		i_vvov_aaaa[b*v/2*o/2*v/2+c*o/2*v/2+i*v/2+a] =
		    I_VVOV(b,c,i,a);
		i_vvov_abab[b*v/2*o/2*v/2+c*o/2*v/2+i*v/2+a] =
		    I_VVOV(b,(c+v/2),i,(a+v/2));
//		i_vvov_abba[b*v/2*o/2*v/2+c*o/2*v/2+i*v/2+a] =
//		    I_VVOV(b,(c+v/2),(i+o/2),a);
	}}}}
	for (size_t i = 0; i < o/2; i++) {
	for (size_t j = 0; j < o/2; j++) {
	for (size_t k = 0; k < o/2; k++) {
	for (size_t a = 0; a < v/2; a++) {
		i_oovo_aaaa[i*o/2*o/2*v/2+j*o/2*v/2+a*o/2+k] =
		    I_OOVO(i,j,a,k);
		i_oovo_abab[i*o/2*o/2*v/2+j*o/2*v/2+a*o/2+k] =
		    I_OOVO(i,(j+o/2),a,(k+o/2));
//		i_oovo_abba[i*o/2*o/2*v/2+j*o/2*v/2+a*o/2+k] =
//		    I_OOVO(i,(j+o/2),(a+v/2),k);
	}}}}

#pragma omp parallel
	{
	double *ijk11 = malloc(o/2*o/2*o/2*sizeof(double));
	double *ijk12 = malloc(o/2*o/2*o/2*sizeof(double));
	double *ijk13 = malloc(o/2*o/2*o/2*sizeof(double));
//	double *ijk14 = malloc(o/2*o/2*o/2*sizeof(double));
	double *ijk15 = malloc(o/2*o/2*o/2*sizeof(double));
//	double *ijk16 = malloc(o/2*o/2*o/2*sizeof(double));
	double *ijk17 = malloc(o/2*o/2*o/2*sizeof(double));
//	double *ijk18 = malloc(o/2*o/2*o/2*sizeof(double));
	double *ijk21 = malloc(o/2*o/2*o/2*sizeof(double));
//	double *ijk22 = malloc(o/2*o/2*o/2*sizeof(double));
	double *ijk23 = malloc(o/2*o/2*o/2*sizeof(double));
	double *ijk24 = malloc(o/2*o/2*o/2*sizeof(double));
//	double *ijk25 = malloc(o/2*o/2*o/2*sizeof(double));
//	double *ijk26 = malloc(o/2*o/2*o/2*sizeof(double));
	double *ijk27 = malloc(o/2*o/2*o/2*sizeof(double));
	double *ijk28 = malloc(o/2*o/2*o/2*sizeof(double));

#pragma omp for reduction(+:e_pt2) schedule(dynamic)
	for (size_t a = 0; a < v/2; a++) {
	for (size_t b = a+1; b < v/2; b++) {
	for (size_t c = 0; c < v/2; c++) {

	comp_t3a_ijk_1(o/2,v/2,a,c,b,ijk11,t2_aaaa,i_vvov_abab);
	comp_t3a_ijk_1(o/2,v/2,b,c,a,ijk12,t2_aaaa,i_vvov_abab);
	comp_t3a_ijk_1(o/2,v/2,a,b,c,ijk13,t2_abab,i_vvov_abab);
//	comp_t3a_ijk_1(o/2,v/2,a,b,c,ijk14,t2_abab,i_vvov_abab);
	comp_t3a_ijk_1(o/2,v/2,b,a,c,ijk15,t2_abab,i_vvov_abab);
//	comp_t3a_ijk_1(o/2,v/2,b,a,c,ijk16,t2_abab,i_vvov_abab);
	comp_t3a_ijk_1(o/2,v/2,c,a,b,ijk17,t2_abab,i_vvov_aaaa);
//	comp_t3a_ijk_1(o/2,v/2,c,a,b,ijk18,t2_abab,i_vvov_aaaa);

	comp_t3a_ijk_2(o/2,v/2,a,b,c,ijk21,t2t_aaaa,i_oovo_abab);
//	comp_t3a_ijk_2(o/2,v/2,a,b,c,ijk22,t2t_aaaa,i_oovo_abab);
	comp_t3a_ijk_2(o/2,v/2,a,c,b,ijk23,t2t_abab,i_oovo_abab);
	comp_t3a_ijk_2(o/2,v/2,b,c,a,ijk24,t2t_abab,i_oovo_abab);
//	comp_t3a_ijk_2(o/2,v/2,b,c,a,ijk25,t2t_abab,i_oovo_abab);
//	comp_t3a_ijk_2(o/2,v/2,a,c,b,ijk26,t2t_abab,i_oovo_abab);
	comp_t3a_ijk_2(o/2,v/2,c,a,b,ijk27,t2t_abab,i_oovo_aaaa);
	comp_t3a_ijk_2(o/2,v/2,c,b,a,ijk28,t2t_abab,i_oovo_aaaa);

	for (size_t i = 0; i < o/2; i++) {
	for (size_t j = i+1; j < o/2; j++) {
	for (size_t k = 0; k < o/2; k++) {
		double t3ax1, t3ax2, t3bx, dn;
		double ijkabc, jikabc, kjiabc, ikjabc;

		t3ax1 =
      -ijk11[j+i*o/2+k*o/2*o/2] //- comp_t3a_ijkabc_11h(o/2,v/2,i,j,k,a,c,b,t2_aaaa,i_vvov_abab)
      +ijk12[j+i*o/2+k*o/2*o/2] //+ comp_t3a_ijkabc_11h(o/2,v/2,i,j,k,b,c,a,t2_aaaa,i_vvov_abab)
      -ijk13[k+i*o/2+j*o/2*o/2] //- comp_t3a_ijkabc_11h(o/2,v/2,i,k,j,a,b,c,t2_abab,i_vvov_abab)
      +ijk13[k+j*o/2+i*o/2*o/2] //+ comp_t3a_ijkabc_11h(o/2,v/2,j,k,i,a,b,c,t2_abab,i_vvov_abab)
      +ijk15[k+i*o/2+j*o/2*o/2] //+ comp_t3a_ijkabc_11h(o/2,v/2,i,k,j,b,a,c,t2_abab,i_vvov_abab)
      -ijk15[k+j*o/2+i*o/2*o/2] //- comp_t3a_ijkabc_11h(o/2,v/2,j,k,i,b,a,c,t2_abab,i_vvov_abab)
      -ijk17[j+k*o/2+i*o/2*o/2] //- comp_t3a_ijkabc_11h(o/2,v/2,k,j,i,c,a,b,t2_abab,i_vvov_aaaa)
      +ijk17[i+k*o/2+j*o/2*o/2];//+ comp_t3a_ijkabc_11h(o/2,v/2,k,i,j,c,a,b,t2_abab,i_vvov_aaaa);

		t3ax2 =
      -ijk21[i+j*o/2+k*o/2*o/2] //- comp_t3a_ijkabc_21h(o/2,v/2,i,k,j,a,b,c,t2t_aaaa,i_oovo_abab)
      +ijk21[j+i*o/2+k*o/2*o/2] //+ comp_t3a_ijkabc_21h(o/2,v/2,j,k,i,a,b,c,t2t_aaaa,i_oovo_abab)
      -ijk23[i+k*o/2+j*o/2*o/2] //- comp_t3a_ijkabc_21h(o/2,v/2,i,j,k,a,c,b,t2t_abab,i_oovo_abab)
      +ijk24[i+k*o/2+j*o/2*o/2] //+ comp_t3a_ijkabc_21h(o/2,v/2,i,j,k,b,c,a,t2t_abab,i_oovo_abab)
      -ijk24[j+k*o/2+i*o/2*o/2] //- comp_t3a_ijkabc_21h(o/2,v/2,j,i,k,b,c,a,t2t_abab,i_oovo_abab)
      +ijk23[j+k*o/2+i*o/2*o/2] //+ comp_t3a_ijkabc_21h(o/2,v/2,j,i,k,a,c,b,t2t_abab,i_oovo_abab)
      -ijk27[k+i*o/2+j*o/2*o/2] //- comp_t3a_ijkabc_21h(o/2,v/2,k,j,i,c,a,b,t2t_abab,i_oovo_aaaa)
      +ijk28[k+i*o/2+j*o/2*o/2];//+ comp_t3a_ijkabc_21h(o/2,v/2,k,j,i,c,b,a,t2t_abab,i_oovo_aaaa);

	ijkabc = comp_t3b_ijkabc(o,v,i,j,k+o/2,a,b,c+v/2,t1,i_oovv,f_ov,t2) -
		 comp_t3b_ijkabc(o,v,i,j,k+o/2,b,a,c+v/2,t1,i_oovv,f_ov,t2) -
		 comp_t3b_ijkabc(o,v,i,j,k+o/2,c+v/2,b,a,t1,i_oovv,f_ov,t2);
	jikabc = comp_t3b_ijkabc(o,v,j,i,k+o/2,a,b,c+v/2,t1,i_oovv,f_ov,t2) -
		 comp_t3b_ijkabc(o,v,j,i,k+o/2,b,a,c+v/2,t1,i_oovv,f_ov,t2) -
		 comp_t3b_ijkabc(o,v,j,i,k+o/2,c+v/2,b,a,t1,i_oovv,f_ov,t2);
	kjiabc = comp_t3b_ijkabc(o,v,k+o/2,j,i,a,b,c+v/2,t1,i_oovv,f_ov,t2) -
		 comp_t3b_ijkabc(o,v,k+o/2,j,i,b,a,c+v/2,t1,i_oovv,f_ov,t2) -
		 comp_t3b_ijkabc(o,v,k+o/2,j,i,c+v/2,b,a,t1,i_oovv,f_ov,t2);
	t3bx = ijkabc - jikabc - kjiabc;

		dn = D_OV(i, a) + D_OV(j, b) + D_OV(k, c);
		e_pt2 += (t3ax1+t3ax2) * (t3ax1+t3ax2-t3bx) / dn;
	}}}}}}
	}

	e_pt2 *= 2.0;
	printf("aabaab %g\n", e_pt2);

	return (e_pt1+e_pt2);
}

static double
ccsd_pt_worker(int id, int nid, size_t o, size_t v, size_t x,
    const double *d_ov, const double *f_ov, const double *t1,
    const double *t2, const double *t2t, const double *i_oovo,
    const double *i_oovv, const double *i_vvov, const double *ovx,
    const double *vvx)
{
//	double e_pt = 0.0, *t3a, *t3b, *work;
//	size_t a, b, c, iter;

//	t3a = malloc(3*o*o*o*sizeof(double));
//	if (t3a == NULL)
//		err(1, "malloc");
//	t3b = malloc(3*o*o*o*sizeof(double));
//	if (t3b == NULL)
//		err(1, "malloc");
//	n = o > v ? o : v;
//	work = malloc((o*n + o*o*n + o*x + v*x) * sizeof(double));
//	if (work == NULL)
//		err(1, "malloc");

//	for (a = 0, iter = 0; a < v/2; a++) {
//	for (b = a+1; b < v/2; b++) {
//	for (c = b+1; c < v; c++, iter++) {
//		if (iter % nid != id)
//			continue;
//		ccsd_t3a(o, v, x, a, b, c, t3a, t3b, t2, i_ooov, i_ovvv,
//		    ovx, vvx, work);
//		ccsd_asymm_t3a(o, t3a);
//		ccsd_asymm_t3b(o, t3b);
//
//		for (n = 0; n < o*o*o; n++) t3a[n] = t3b[n]-t3a[n];
//
//		ccsd_t3b(o, v, a, b, c, t3b, f_ov, t1, t2, i_oovv, work);
//		ccsd_asymm_t3b(o, t3b);

//		comp_t3a_asymm(o,v,a,b,c,t2,i_ooov,i_ovvv,t3a,t3b);
//		for (n = 0; n < o*o*o; n++) t3a[n] = t3b[n]-t3a[n];
//		comp_t3b_asymm(o,v,a,b,c,t1,i_oovv,f_ov,t2,t3b);

///		e_pt += ccsd_pt_energy(o, v, a, b, c, d_ov, f_ov, t1,
///		    t2, i_ooov, i_oovv, i_ovvv);
///	}}}
///	e_pt *= 2.0;

//	free(t3a);
//	free(t3b);
//	free(work);
	return (ccsd_pt_energy(o, v, d_ov, f_ov, t1,
	    t2, t2t, i_oovo, i_oovv, i_vvov));
}

static double
do_ccsd_pt(size_t o, size_t v, size_t x, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2, const double *t2t, const double *i_oovo,
    const double *i_oovv, const double *i_vvov,
    const double *ovx, const double *vvx)
{
	double e_pt = 0.0;
	int pid, npid;

	if (o == 0 || v == 0)
		return (0.0);
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &npid);

//#ifdef _OPENMP
//#pragma omp parallel reduction(+:e_pt)
//#endif
	{
		int id, nid, tid = 0, ntid = 1;
#ifdef _OPENMP
//		tid = omp_get_thread_num();
//		ntid = omp_get_num_threads();
#endif
		id = pid * ntid + tid;
		nid = npid * ntid;

		e_pt = ccsd_pt_worker(id, nid, o, v, x, d_ov, f_ov, t1, t2, t2t,
		    i_oovo, i_oovv, i_vvov, ovx, vvx);
	}

	MPI_Allreduce(MPI_IN_PLACE, &e_pt, 1, MPI_DOUBLE, MPI_SUM,
	    MPI_COMM_WORLD);
	return (e_pt);
}

double
ccsd_pt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2, const double *t2t, const double *i_oovo,
    const double *i_oovv, const double *i_vvov)
{
	return (do_ccsd_pt(o, v, 0, d_ov, f_ov, t1, t2, t2t,
	    i_oovo, i_oovv, i_vvov, NULL, NULL));
}

//double
//ccsd_ri_pt(size_t o, size_t v, size_t x, const double *d_ov,
//    const double *f_ov, const double *t1, const double *t2,
//    const double *i_ooov, const double *i_oovv,
//    const double *ovx, const double *vvx)
//{
//	return (do_ccsd_pt(o, v, x, d_ov, f_ov, t1, t2,
//	    i_ooov, i_oovv, NULL, ovx, vvx));
//}
