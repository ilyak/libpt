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
#include <time.h>

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "pt.h"

//#define D_OV(i, a) d_ov[i*v+a]
//#define F_OV(i, a) f_ov[i*v+a]
//#define I_OOOV(i, j, k, a) i_ooov[i*o*o*v+j*o*v+k*v+a]
#define I_OOVO(i, j, a, k) i_oovo[i*o*o*v+j*o*v+a*o+k]
#define I_OOVV(i, j, a, b) i_oovv[i*o*v*v+j*v*v+a*v+b]
#define I_OVVV(i, a, b, c) i_ovvv[i*v*v*v+a*v*v+b*v+c]
//#define I_VVOV(b, c, i, a) i_vvov[b*v*o*v+c*o*v+i*v+a]
//#define T1(i, a) t1[i*v+a]
#define T2(i, j, a, b) t2[i*o*v*v+j*v*v+a*v+b]
//#define T2T(a, b, i, j) t2t[a*v*o*o+b*o*o+i*o+j]
//#define T3AABC(i, j, k) t3a[0*o*o*o+i*o*o+j*o+k]
//#define T3ACBA(i, j, k) t3a[1*o*o*o+i*o*o+j*o+k]
//#define T3AACB(i, j, k) t3a[2*o*o*o+i*o*o+j*o+k]
//#define T3BABC(i, j, k) t3b[0*o*o*o+i*o*o+j*o+k]
//#define T3BBAC(i, j, k) t3b[1*o*o*o+i*o*o+j*o+k]
//#define T3BCBA(i, j, k) t3b[2*o*o*o+i*o*o+j*o+k]
//#define OVX(i, j, k) ovx[i*v*x+j*x+k]
//#define VVX(i, j, k) vvx[i*v*x+j*x+k]
//#define MOV(i, a) mov[i*v+a]
//#define MOO1(i, j) moo1[i*o+j]
//#define MOO2(i, j) moo2[i*o+j]
//#define MOOO(i, j, k) mooo[i*o*o+j*o+k]
//#define MVOO(a, i, j) mvoo[a*o*o+i*o+j]
//#define MOX(i, j) mox[i*x+j]
//#define MXV(i, j) mxv[i*v+j]

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

static void
comp_t3a_abc_1(size_t o, size_t v, size_t i, size_t j, size_t k,
    double *abc, const double *t2, const double *i_ovvv)
{
	const double *t2_p = &T2(i,j,0,0);
	const double *i_ovvv_p = &i_ovvv[k*v*v*(v-1)/2];
	int lda = v;
	int ldb = v*(v-1)/2;

	/* t3a1(i,j,k,a,b,c) = contract(d, t2(i,j,a,d), i_ovvv(k,d,b,c)) */

	gemm('T', 'T', v, v*(v-1)/2, v, 1.0, t2_p, lda,
	    i_ovvv_p, ldb, 0, abc, v);
}

static void
comp_t3a_abc_12(size_t o, size_t v, size_t i, size_t j, size_t k,
    double *abc, const double *t2, const double *i_ovvv)
{
	const double *t2_p = &T2(i,j,0,0);
	const double *i_ovvv_p = &I_OVVV(k,0,0,0);
	int lda = v;
	int ldb = v*v;

	/* t3a1(i,j,k,a,b,c) = contract(d, t2(i,j,a,d), i_ovvv(k,d,b,c)) */

	gemm('T', 'T', v, v*v, v, 1.0, t2_p, lda, i_ovvv_p, ldb, 0, abc, v);
}

static void
comp_t3a_abc_2(size_t o, size_t v, size_t i, size_t j, size_t k,
    double *abc, const double *t2, const double *i_oovo)
{
	const double *t2_p = &T2(i,0,0,0);
	const double *i_oovo_p = &I_OOVO(j,k,0,0);
	int lda = v*v;
	int ldb = o;

	/* t3a2(i,j,k,a,b,c) = contract(l, t2(i,l,a,b), i_oovo(j,k,c,l)) */

	gemm('N', 'N', v*v, v, o, 1.0, t2_p, lda, i_oovo_p, ldb, 0, abc, v*v);
}

static double
comp_t3b_ijkabc(size_t o, size_t v, size_t i, size_t j, size_t k,
    size_t a, size_t b, size_t c, const double *t1, const double *i_oovv,
    const double *f_ov, const double *t2)
{
	double t1_ia = t1[i*v+a];
	double f_ov_ia = f_ov[i*v+a];

	return t1_ia*I_OOVV(j,k,b,c) + f_ov_ia*T2(j,k,b,c);
}

static double
ccsd_pt_aaaa(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2, const double *i_oovo,
    const double *i_oovv, const double *i_ovvv)
{
	double e_pt = 0.0;
	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

#pragma omp parallel
{
	size_t nij = 0, *ij = malloc(o*(o-1)*sizeof(size_t));
	if (ij == NULL)
		err(1, "malloc ij");
	for (size_t i = 0, n = 0; i < o; i++)
	for (size_t j = i+1; j < o; j++, n++) {
		if (n % size == rank) {
			ij[2*nij+0] = i;
			ij[2*nij+1] = j;
			nij++;
		}
	}

	const double *t2_aaaa;
	const double *i_ovvv_aaaa;
	const double *i_oovo_aaaa;
	const double *i_oovv_aaaa;

	t2_aaaa = t2;
	i_ovvv_aaaa = i_ovvv;
	i_oovo_aaaa = i_oovo;
	i_oovv_aaaa = i_oovv;

	double *t3ax1, *abc1, *abc2, *abc3;
	double *work = malloc(4*v*v*v*sizeof(double));
	if (work == NULL)
		err(1, "malloc work");
	t3ax1 = work;
	abc1 = work + 1*v*v*v;
	abc2 = work + 2*v*v*v;
	abc3 = work + 3*v*v*v;

	/* aaaa spin-block */
#pragma omp for reduction(+:e_pt) schedule(dynamic)
	for (size_t it = 0; it < nij; it++) {
		size_t i = ij[2*it+0];
		size_t j = ij[2*it+1];
	for (size_t k = j+1; k < o; k++) {

	memset(t3ax1, 0, v*v*v*sizeof(double));

	comp_t3a_abc_1(o,v,i,j,k,abc1,t2_aaaa,i_ovvv_aaaa);
	comp_t3a_abc_1(o,v,i,k,j,abc2,t2_aaaa,i_ovvv_aaaa);
	comp_t3a_abc_1(o,v,k,j,i,abc3,t2_aaaa,i_ovvv_aaaa);
	for (size_t a = 0; a < v; a++) {
	for (size_t b = 0; b < a; b++) {
	for (size_t c = 0; c < b; c++) {
		t3ax1[c+b*v+a*v*v] +=
+abc1[a+v*b*(b-1)/2+v*c] //+t3a_ijkabc_11(o,v,i,j,k,a,b,c,t2_aaaa,i_vvov_aaaa)
-abc1[b+v*a*(a-1)/2+v*c] //-t3a_ijkabc_11(o,v,k,j,i,a,b,c,t2_aaaa,i_vvov_aaaa)
+abc1[c+v*a*(a-1)/2+v*b] //-t3a_ijkabc_11(o,v,i,k,j,a,b,c,t2_aaaa,i_vvov_aaaa)
-abc2[a+v*b*(b-1)/2+v*c] //-t3a_ijkabc_11(o,v,i,j,k,b,a,c,t2_aaaa,i_vvov_aaaa)
+abc2[b+v*a*(a-1)/2+v*c] //+t3a_ijkabc_11(o,v,k,j,i,b,a,c,t2_aaaa,i_vvov_aaaa)
-abc2[c+v*a*(a-1)/2+v*b] //+t3a_ijkabc_11(o,v,i,k,j,b,a,c,t2_aaaa,i_vvov_aaaa)
-abc3[a+v*b*(b-1)/2+v*c] //-t3a_ijkabc_11(o,v,i,j,k,c,b,a,t2_aaaa,i_vvov_aaaa)
+abc3[b+v*a*(a-1)/2+v*c] //+t3a_ijkabc_11(o,v,k,j,i,c,b,a,t2_aaaa,i_vvov_aaaa)
-abc3[c+v*a*(a-1)/2+v*b];//+t3a_ijkabc_11(o,v,i,k,j,c,b,a,t2_aaaa,i_vvov_aaaa);
	}}}

	comp_t3a_abc_2(o,v,i,j,k,abc1,t2_aaaa,i_oovo_aaaa);
	comp_t3a_abc_2(o,v,j,i,k,abc2,t2_aaaa,i_oovo_aaaa);
	comp_t3a_abc_2(o,v,k,j,i,abc3,t2_aaaa,i_oovo_aaaa);
	for (size_t a = 0; a < v; a++) {
	for (size_t b = 0; b < a; b++) {
	for (size_t c = 0; c < b; c++) {
		double t3ax, t3bx, dn;

		t3ax1[c+b*v+a*v*v] +=
+abc1[b+a*v+c*v*v] //+comp_t3a_ijkabc_21(o,v,i,j,k,a,b,c,t2t_aaaa,i_oovo_aaaa)
-abc1[c+a*v+b*v*v] //-comp_t3a_ijkabc_21(o,v,j,i,k,a,b,c,t2t_aaaa,i_oovo_aaaa)
-abc1[b+c*v+a*v*v] //-comp_t3a_ijkabc_21(o,v,k,j,i,a,b,c,t2t_aaaa,i_oovo_aaaa)
-abc2[b+a*v+c*v*v] //-comp_t3a_ijkabc_21(o,v,i,j,k,a,c,b,t2t_aaaa,i_oovo_aaaa)
+abc2[c+a*v+b*v*v] //+comp_t3a_ijkabc_21(o,v,j,i,k,a,c,b,t2t_aaaa,i_oovo_aaaa)
+abc2[b+c*v+a*v*v] //+comp_t3a_ijkabc_21(o,v,k,j,i,a,c,b,t2t_aaaa,i_oovo_aaaa)
-abc3[b+a*v+c*v*v] //-comp_t3a_ijkabc_21(o,v,i,j,k,c,b,a,t2t_aaaa,i_oovo_aaaa)
+abc3[c+a*v+b*v*v] //+comp_t3a_ijkabc_21(o,v,j,i,k,c,b,a,t2t_aaaa,i_oovo_aaaa)
+abc3[b+c*v+a*v*v];//+comp_t3a_ijkabc_21(o,v,k,j,i,c,b,a,t2t_aaaa,i_oovo_aaaa);

		t3bx =
+comp_t3b_ijkabc(o,v,i,j,k,a,b,c,t1,i_oovv_aaaa,f_ov,t2_aaaa)
-comp_t3b_ijkabc(o,v,j,i,k,a,b,c,t1,i_oovv_aaaa,f_ov,t2_aaaa)
-comp_t3b_ijkabc(o,v,k,j,i,a,b,c,t1,i_oovv_aaaa,f_ov,t2_aaaa)
-comp_t3b_ijkabc(o,v,i,j,k,b,a,c,t1,i_oovv_aaaa,f_ov,t2_aaaa)
+comp_t3b_ijkabc(o,v,j,i,k,b,a,c,t1,i_oovv_aaaa,f_ov,t2_aaaa)
+comp_t3b_ijkabc(o,v,k,j,i,b,a,c,t1,i_oovv_aaaa,f_ov,t2_aaaa)
-comp_t3b_ijkabc(o,v,i,j,k,c,b,a,t1,i_oovv_aaaa,f_ov,t2_aaaa)
+comp_t3b_ijkabc(o,v,j,i,k,c,b,a,t1,i_oovv_aaaa,f_ov,t2_aaaa)
+comp_t3b_ijkabc(o,v,k,j,i,c,b,a,t1,i_oovv_aaaa,f_ov,t2_aaaa);

		dn = d_ov[i*v+a] + d_ov[j*v+b] + d_ov[k*v+c];
		t3ax = t3ax1[c+b*v+a*v*v];
		e_pt += t3ax * (t3ax-t3bx) / dn;
	}}}
	}}
	free(ij);
	free(work);
}
	MPI_Allreduce(MPI_IN_PLACE, &e_pt, 1, MPI_DOUBLE, MPI_SUM,
	    MPI_COMM_WORLD);
	return (e_pt);
}

static double
ccsd_pt_abab(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2, const double *i_oovo,
    const double *i_oovv, const double *i_ovvv)
{
	double e_pt = 0.0;
	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

#pragma omp parallel
{
	size_t nij = 0, *ij = malloc(o*(o-1)*sizeof(size_t));
	if (ij == NULL)
		err(1, "malloc ij");
	for (size_t i = 0, n = 0; i < o; i++)
	for (size_t j = i+1; j < o; j++, n++) {
		if (n % size == rank) {
			ij[2*nij+0] = i;
			ij[2*nij+1] = j;
			nij++;
		}
	}

	const double *t2_aaaa, *t2_abab;
	const double *i_ovvv_aaaa, *i_ovvv_abab;
	const double *i_oovo_aaaa, *i_oovo_abab;
	const double *i_oovv_aaaa, *i_oovv_abab;

	t2_aaaa = t2;
	t2_abab = t2 + o*o*v*v;
	i_ovvv_aaaa = i_ovvv;
	i_ovvv_abab = i_ovvv + o*v*v*(v-1)/2;
	i_oovo_aaaa = i_oovo;
	i_oovo_abab = i_oovo + o*o*o*v;
	i_oovv_aaaa = i_oovv;
	i_oovv_abab = i_oovv + o*o*v*v;

	double *abc11, *abc13, *abc14, *abc17x, *abc18x;
	double *abc21, *abc22, *abc23, *abc25, *abc27;
	double *t3ax1;

	double *work = malloc(4*v*v*v*sizeof(double));
	if (work == NULL)
		err(1, "malloc work");
	t3ax1 = work;
	abc11 = work + 1*v*v*v;
	abc13 = work + 2*v*v*v;
	abc14 = work + 1*v*v*v;
	abc17x = work + 2*v*v*v;
	abc18x = work + 3*v*v*v;
	abc21 = work + 1*v*v*v;
	abc22 = work + 2*v*v*v;
	abc27 = work + 3*v*v*v;
	abc23 = work + 1*v*v*v;
	abc25 = work + 2*v*v*v;

	/* abab spin-block */
#pragma omp for reduction(+:e_pt) schedule(dynamic)
	for (size_t it = 0; it < nij; it++) {
		size_t i = ij[2*it+0];
		size_t j = ij[2*it+1];
	for (size_t k = 0; k < o; k++) {

	memset(t3ax1, 0, v*v*v*sizeof(double));

	comp_t3a_abc_12(o,v,i,j,k,abc11,t2_aaaa,i_ovvv_abab);
//	comp_t3a_abc_12(o,v,i,j,k,abc12,t2_aaaa,i_ovvv_abab);
	comp_t3a_abc_12(o,v,i,k,j,abc13,t2_abab,i_ovvv_abab);
	for (size_t a = 0; a < v; a++) {
	for (size_t b = 0; b < a; b++) {
	for (size_t c = 0; c < v; c++) {
		t3ax1[a+b*v+c*v*v] +=
-abc11[a+b*v+c*v*v] //-t3a_ijkabc_11h(o,v,i,j,k,a,c,b,t2_aaaa,i_vvov_abab)
+abc11[b+a*v+c*v*v] //+t3a_ijkabc_11h(o,v,i,j,k,b,c,a,t2_aaaa,i_vvov_abab)
-abc13[a+c*v+b*v*v] //-t3a_ijkabc_11h(o,v,i,k,j,a,b,c,t2_abab,i_vvov_abab)
+abc13[b+c*v+a*v*v];//+t3a_ijkabc_11h(o,v,i,k,j,b,a,c,t2_abab,i_vvov_abab)
	}}}

	comp_t3a_abc_12(o,v,j,k,i,abc14,t2_abab,i_ovvv_abab);
//	comp_t3a_abc_12(o,v,i,k,j,abc15,t2_abab,i_ovvv_abab);
//	comp_t3a_abc_12(o,v,j,k,i,abc16,t2_abab,i_ovvv_abab);
	comp_t3a_abc_1(o,v,k,j,i,abc17x,t2_abab,i_ovvv_aaaa);
	comp_t3a_abc_1(o,v,k,i,j,abc18x,t2_abab,i_ovvv_aaaa);
	for (size_t a = 0; a < v; a++) {
	for (size_t b = 0; b < a; b++) {
	for (size_t c = 0; c < v; c++) {
		t3ax1[a+b*v+c*v*v] +=
+abc14[a+c*v+b*v*v] //+t3a_ijkabc_11h(o,v,j,k,i,a,b,c,t2_abab,i_vvov_abab)
-abc14[b+c*v+a*v*v] //-t3a_ijkabc_11h(o,v,j,k,i,b,a,c,t2_abab,i_vvov_abab)
-abc17x[c+v*a*(a-1)/2+v*b] //-ijkabc_11h(o,v,k,j,i,c,a,b,t2_abab,i_vvov_aaaa)
+abc18x[c+v*a*(a-1)/2+v*b];//+ijkabc_11h(o,v,k,i,j,c,a,b,t2_abab,i_vvov_aaaa)
	}}}

	comp_t3a_abc_2(o,v,i,k,j,abc21,t2_aaaa,i_oovo_abab);
	comp_t3a_abc_2(o,v,j,k,i,abc22,t2_aaaa,i_oovo_abab);
	comp_t3a_abc_2(o,v,k,j,i,abc27,t2_abab,i_oovo_aaaa);
	for (size_t a = 0; a < v; a++) {
	for (size_t b = 0; b < a; b++) {
	for (size_t c = 0; c < v; c++) {
		t3ax1[a+b*v+c*v*v] +=
-abc21[b+a*v+c*v*v] //-comp_t3a_ijkabc_21h(o,v,i,k,j,a,b,c,t2t_aaaa,i_oovo_abab)
+abc22[b+a*v+c*v*v] //+comp_t3a_ijkabc_21h(o,v,j,k,i,a,b,c,t2t_aaaa,i_oovo_abab)
-abc27[a+c*v+b*v*v] //-comp_t3a_ijkabc_21h(o,v,k,j,i,c,a,b,t2t_abab,i_oovo_aaaa)
+abc27[b+c*v+a*v*v];//+comp_t3a_ijkabc_21h(o,v,k,j,i,c,b,a,t2t_abab,i_oovo_aaaa)
	}}}

	comp_t3a_abc_2(o,v,i,j,k,abc23,t2_abab,i_oovo_abab);
//	comp_t3a_abc_2(o,v,i,j,k,abc24,t2_abab,i_oovo_abab);
	comp_t3a_abc_2(o,v,j,i,k,abc25,t2_abab,i_oovo_abab);
//	comp_t3a_abc_2(o,v,j,i,k,abc26,t2_abab,i_oovo_abab);
//	comp_t3a_abc_2(o,v,k,j,i,abc28,t2_abab,i_oovo_aaaa);

	for (size_t a = 0; a < v; a++) {
	for (size_t b = 0; b < a; b++) {
	for (size_t c = 0; c < v; c++) {
		double t3ax, t3bx, dn;

		t3ax1[a+b*v+c*v*v] +=
-abc23[c+a*v+b*v*v] //-comp_t3a_ijkabc_21h(o,v,i,j,k,a,c,b,t2t_abab,i_oovo_abab)
+abc23[c+b*v+a*v*v] //+comp_t3a_ijkabc_21h(o,v,i,j,k,b,c,a,t2t_abab,i_oovo_abab)
-abc25[c+b*v+a*v*v] //-comp_t3a_ijkabc_21h(o,v,j,i,k,b,c,a,t2t_abab,i_oovo_abab)
+abc25[c+a*v+b*v*v];//+comp_t3a_ijkabc_21h(o,v,j,i,k,a,c,b,t2t_abab,i_oovo_abab)

		t3bx =
	+comp_t3b_ijkabc(o,v,i,j,k,a,b,c,t1,i_oovv_abab,f_ov,t2_abab)
	-comp_t3b_ijkabc(o,v,i,j,k,b,a,c,t1,i_oovv_abab,f_ov,t2_abab)
//	-comp_t3b_ijkabc(o,v,i,j,k+o,c+v,b,a,t1,i_oovv,f_ov,t2)
	-comp_t3b_ijkabc(o,v,j,i,k,a,b,c,t1,i_oovv_abab,f_ov,t2_abab)
	+comp_t3b_ijkabc(o,v,j,i,k,b,a,c,t1,i_oovv_abab,f_ov,t2_abab)
//	+comp_t3b_ijkabc(o,v,j,i,k+o,c+v,b,a,t1,i_oovv,f_ov,t2)
//	-comp_t3b_ijkabc(o,v,k+o,j,i,a,b,c+v,t1,i_oovv,f_ov,t2)
//	+comp_t3b_ijkabc(o,v,k+o,j,i,b,a,c+v,t1,i_oovv,f_ov,t2)
	+comp_t3b_ijkabc(o,v,k,j,i,c,b,a,t1,i_oovv_aaaa,f_ov,t2_aaaa);

		dn = d_ov[i*v+a] + d_ov[j*v+b] + d_ov[k*v+c];
		t3ax = t3ax1[a+b*v+c*v*v];
		e_pt += t3ax * (t3ax-t3bx) / dn;
	}}}
	}}
	free(ij);
	free(work);
}
	MPI_Allreduce(MPI_IN_PLACE, &e_pt, 1, MPI_DOUBLE, MPI_SUM,
	    MPI_COMM_WORLD);
	return (e_pt);
}

double
ccsd_rpt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2, const double *i_oovo,
    const double *i_oovv, const double *i_ovvv)
{
	double e_pt1, e_pt2;
	time_t wall;

	if (o < 2 || v < 2)
		return (0.0);

	e_pt1 = ccsd_pt_aaaa(o, v, d_ov, f_ov, t1, t2, i_oovo, i_oovv, i_ovvv);
	wall = time(NULL);
	printf("aaaa %g\n", 2.0 * e_pt1);
	printf("ccsd_pt: %s", ctime(&wall));
	e_pt2 = ccsd_pt_abab(o, v, d_ov, f_ov, t1, t2, i_oovo, i_oovv, i_ovvv);
	printf("abab %g\n", 2.0 * e_pt2);

	return 2.0 * (e_pt1 + e_pt2);
}

double
ccsd_upt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2, const double *i_oovo,
    const double *i_oovv, const double *i_ovvv)
{
	if (o < 2 || v < 2)
		return (0.0);
	return ccsd_pt_aaaa(o, v, d_ov, f_ov, t1, t2, i_oovo, i_oovv, i_ovvv);
}
