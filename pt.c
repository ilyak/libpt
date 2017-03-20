/*
 * Copyright (c) 2016-2017 Ilya Kaliman
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

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "pt.h"

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
t2_i_ovvv_half(size_t o, size_t v, size_t i, size_t j, size_t k,
    double *abc, const double *t2, const double *i_ovvv)
{
	const double *t2_p = &t2[i*o*v*v+j*v*v];
	const double *i_ovvv_p = &i_ovvv[k*v*v*(v-1)/2];

	/* out(i,j,k,a,b,c) = contract(d, t2(i,j,a,d), i_ovvv(k,d,b,c)) */

	gemm('T', 'T', v, v*(v-1)/2, v, 1.0, t2_p, v,
	    i_ovvv_p, v*(v-1)/2, 0.0, abc, v);
}

static void
t2_i_ovvv(size_t o, size_t v, size_t i, size_t j, size_t k,
    double *abc, const double *t2, const double *i_ovvv)
{
	const double *t2_p = &t2[i*o*v*v+j*v*v];
	const double *i_ovvv_p = &i_ovvv[k*v*v*v];

	/* out(i,j,k,a,b,c) = contract(d, t2(i,j,a,d), i_ovvv(k,d,b,c)) */

	gemm('T', 'T', v, v*v, v, 1.0, t2_p, v,
	    i_ovvv_p, v*v, 0.0, abc, v);
}

static void
t2_i_oovo(size_t o, size_t v, size_t i, size_t j, size_t k,
    double *abc, const double *t2, const double *i_oovo)
{
	const double *t2_p = &t2[i*o*v*v];
	const double *i_oovo_p = &i_oovo[j*o*o*v+k*o*v];

	/* out(i,j,k,a,b,c) = contract(l, t2(i,l,a,b), i_oovo(j,k,c,l)) */

	gemm('N', 'N', v*v, v, o, 1.0, t2_p, v*v,
	    i_oovo_p, o, 0.0, abc, v*v);
}

static double
i_jk_a_bc_ov_oovv(size_t o, size_t v, const double *ov, const double *oovv,
    size_t i, size_t j, size_t k, size_t a, size_t b, size_t c)
{
	return +ov[i*v+a]*oovv[j*o*v*v+k*v*v+b*v+c]
	       -ov[j*v+a]*oovv[i*o*v*v+k*v*v+b*v+c]
	       -ov[k*v+a]*oovv[j*o*v*v+i*v*v+b*v+c]
	       -ov[i*v+b]*oovv[j*o*v*v+k*v*v+a*v+c]
	       +ov[j*v+b]*oovv[i*o*v*v+k*v*v+a*v+c]
	       +ov[k*v+b]*oovv[j*o*v*v+i*v*v+a*v+c]
	       -ov[i*v+c]*oovv[j*o*v*v+k*v*v+b*v+a]
	       +ov[j*v+c]*oovv[i*o*v*v+k*v*v+b*v+a]
	       +ov[k*v+c]*oovv[j*o*v*v+i*v*v+b*v+a];
}

static double
comp_t3b_ijkabc(size_t o, size_t v, size_t i, size_t j, size_t k,
    size_t a, size_t b, size_t c, const double *t1, const double *i_oovv,
    const double *f_ov, const double *t2)
{
	double t1_ia = t1[i*v+a];
	double f_ov_ia = f_ov[i*v+a];
	double t2_jkbc = t2[j*o*v*v+k*v*v+b*v+c];
	double i_oovv_jkbc = i_oovv[j*o*v*v+k*v*v+b*v+c];

	return t1_ia*i_oovv_jkbc + f_ov_ia*t2_jkbc;
}

static double
asymm_ijk_a_bc(size_t v, const double *abc1, const double *abc2,
    const double *abc3, size_t a, size_t b, size_t c)
{
	return +abc1[a*v*v+b*v+c]
	       -abc1[b*v*v+a*v+c]
	       -abc1[c*v*v+b*v+a]
	       -abc2[a*v*v+b*v+c]
	       +abc2[b*v*v+a*v+c]
	       +abc2[c*v*v+b*v+a]
	       -abc3[a*v*v+b*v+c]
	       +abc3[b*v*v+a*v+c]
	       +abc3[c*v*v+b*v+a];
}

//static double
//asymm_ijk_ab_c(size_t v, const double *abc1, const double *abc2,
//    const double *abc3, size_t a, size_t b, size_t c)
//{
//	return +abc1[a*v*v+b*v+c]
//	       -abc1[a*v*v+c*v+b]
//	       -abc1[c*v*v+b*v+a]
//	       -abc2[a*v*v+b*v+c]
//	       +abc2[a*v*v+c*v+b]
//	       +abc2[c*v*v+b*v+a]
//	       -abc3[a*v*v+b*v+c]
//	       +abc3[a*v*v+c*v+b]
//	       +abc3[c*v*v+b*v+a];
//}

static double
asymm_ijk_ab_c_half(size_t v, const double *abc1, const double *abc2,
    const double *abc3, size_t a, size_t b, size_t c)
{
	return +abc1[a*(a-1)/2*v+b*v+c]
	       -abc1[a*(a-1)/2*v+c*v+b]
	       +abc1[b*(b-1)/2*v+c*v+a]
	       -abc2[a*(a-1)/2*v+b*v+c]
	       +abc2[a*(a-1)/2*v+c*v+b]
	       -abc2[b*(b-1)/2*v+c*v+a]
	       -abc3[a*(a-1)/2*v+b*v+c]
	       +abc3[a*(a-1)/2*v+c*v+b]
	       -abc3[b*(b-1)/2*v+c*v+a];
}

static double
cc_pt_aaa(size_t oa, size_t va, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2_aaaa, const double *i_oovo_aaaa,
    const double *i_oovv_aaaa, const double *i_ovvv_aaaa)
{
	double e_pt = 0.0;
	int rank = 0, size = 1;

#ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
{
	size_t i, j, k, a, b, c, it, *ij, nij = 0;
	double *t3ax1, *abc1, *abc2, *abc3;

	if ((ij = malloc(oa*(oa-1)*sizeof(size_t))) == NULL)
		err(1, "libpt malloc ij");
	for (i = 0, it = 0; i < oa; i++) {
		for (j = i+1; j < oa; j++, it++) {
			if ((int)it % size == rank) {
				ij[2*nij+0] = i;
				ij[2*nij+1] = j;
				nij++;
			}
		}
	}

	if ((t3ax1 = malloc(4*va*va*va*sizeof(double))) == NULL)
		err(1, "libpt malloc work");
	abc1 = t3ax1 + 1*va*va*va;
	abc2 = t3ax1 + 2*va*va*va;
	abc3 = t3ax1 + 3*va*va*va;

#ifdef _OPENMP
#pragma omp for reduction(+:e_pt) schedule(dynamic)
#endif
	for (it = 0; it < nij; it++) {
		i = ij[2*it+0];
		j = ij[2*it+1];
	for (k = j+1; k < oa; k++) {

	memset(t3ax1, 0, va*va*va*sizeof(double));

	t2_i_ovvv_half(oa,va,i,j,k,abc1,t2_aaaa,i_ovvv_aaaa);
	t2_i_ovvv_half(oa,va,i,k,j,abc2,t2_aaaa,i_ovvv_aaaa);
	t2_i_ovvv_half(oa,va,k,j,i,abc3,t2_aaaa,i_ovvv_aaaa);
	for (a = 2; a < va; a++) {
	for (b = 1; b < a; b++) {
	for (c = 0; c < b; c++) {
		t3ax1[a*va*va+b*va+c] +=
		    asymm_ijk_ab_c_half(va,abc1,abc2,abc3,a,b,c);
	}}}

	t2_i_oovo(oa,va,i,j,k,abc1,t2_aaaa,i_oovo_aaaa);
	t2_i_oovo(oa,va,j,i,k,abc2,t2_aaaa,i_oovo_aaaa);
	t2_i_oovo(oa,va,k,j,i,abc3,t2_aaaa,i_oovo_aaaa);
	for (a = 2; a < va; a++) {
	for (b = 1; b < a; b++) {
	for (c = 0; c < b; c++) {
		double t3ax, t3bx, dn;

		t3ax1[a*va*va+b*va+c] +=
		    asymm_ijk_a_bc(va,abc1,abc2,abc3,a,b,c);
		dn = d_ov[i*va+a] + d_ov[j*va+b] + d_ov[k*va+c];
		t3ax = t3ax1[a*va*va+b*va+c];
		t3bx = +i_jk_a_bc_ov_oovv(oa,va,t1,i_oovv_aaaa,i,j,k,a,b,c)
		       +i_jk_a_bc_ov_oovv(oa,va,f_ov,t2_aaaa,i,j,k,a,b,c);
		e_pt += t3ax * (t3ax-t3bx) / dn;
	}}}
	}}
	free(ij);
	free(t3ax1);
}
#ifdef WITH_MPI
	MPI_Allreduce(MPI_IN_PLACE, &e_pt, 1, MPI_DOUBLE,
	    MPI_SUM, MPI_COMM_WORLD);
#endif
	return (e_pt);
}

static double
cc_pt_aab(size_t oa, size_t ob, size_t va, size_t vb,
    const double *d_ov, const double *f_ov, const double *t1,
    const double *t2_aaaa, const double *t2_abab,
    const double *i_oovo_aaaa, const double *i_oovo_abab,
    const double *i_oovv_aaaa, const double *i_oovv_abab,
    const double *i_ovvv_aaaa, const double *i_ovvv_abab)
{
	double e_pt = 0.0;
	int rank = 0, size = 1;

#ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
{
	size_t i, j, k, a, b, c, it, *ij, nij = 0;
	double *abc11, *abc13, *abc14, *abc17x, *abc18x;
	double *abc21, *abc22, *abc23, *abc25, *abc27, *t3ax1;

	if ((ij = malloc(oa*(oa-1)*sizeof(size_t))) == NULL)
		err(1, "libpt malloc ij");
	for (i = 0, it = 0; i < oa; i++) {
		for (j = i+1; j < oa; j++, it++) {
			if ((int)it % size == rank) {
				ij[2*nij+0] = i;
				ij[2*nij+1] = j;
				nij++;
			}
		}
	}

	if ((t3ax1 = malloc(4*va*va*va*sizeof(double))) == NULL)
		err(1, "libpt malloc work");
	abc11 = t3ax1 + 1*va*va*va;
	abc13 = t3ax1 + 2*va*va*va;
	abc14 = t3ax1 + 3*va*va*va;
	abc17x = t3ax1 + 1*va*va*va;
	abc18x = t3ax1 + 1*va*va*va + va*va*(va-1)/2;
	abc21 = t3ax1 + 2*va*va*va;
	abc22 = t3ax1 + 3*va*va*va;
	abc23 = t3ax1 + 1*va*va*va;
	abc25 = t3ax1 + 2*va*va*va;
	abc27 = t3ax1 + 3*va*va*va;

#ifdef _OPENMP
#pragma omp for reduction(+:e_pt) schedule(dynamic)
#endif
	for (it = 0; it < nij; it++) {
		i = ij[2*it+0];
		j = ij[2*it+1];
	for (k = 0; k < ob; k++) {

	memset(t3ax1, 0, va*va*va*sizeof(double));

	t2_i_ovvv(oa,va,i,j,k,abc11,t2_aaaa,i_ovvv_abab);
	t2_i_ovvv(oa,va,i,k,j,abc13,t2_abab,i_ovvv_abab);
	t2_i_ovvv(oa,va,j,k,i,abc14,t2_abab,i_ovvv_abab);
	for (a = 1; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		t3ax1[a*va*va+b*va+c] +=
		    -abc11[a+b*va+c*va*va]
		    +abc11[b+a*va+c*va*va]
		    -abc13[a+c*va+b*va*va]
		    +abc13[b+c*va+a*va*va]
		    +abc14[a+c*va+b*va*va]
		    -abc14[b+c*va+a*va*va];
	}}}

	t2_i_ovvv_half(oa,va,k,j,i,abc17x,t2_abab,i_ovvv_aaaa);
	t2_i_ovvv_half(oa,va,k,i,j,abc18x,t2_abab,i_ovvv_aaaa);
	t2_i_oovo(oa,va,i,k,j,abc21,t2_aaaa,i_oovo_abab);
	t2_i_oovo(oa,va,j,k,i,abc22,t2_aaaa,i_oovo_abab);
	for (a = 1; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		t3ax1[a*va*va+b*va+c] +=
		    -abc17x[c+va*a*(a-1)/2+va*b]
		    +abc18x[c+va*a*(a-1)/2+va*b]
		    -abc21[b+a*va+c*va*va]
		    +abc22[b+a*va+c*va*va];
	}}}

	t2_i_oovo(oa,va,i,j,k,abc23,t2_abab,i_oovo_abab);
	t2_i_oovo(oa,va,j,i,k,abc25,t2_abab,i_oovo_abab);
	t2_i_oovo(oa,va,k,j,i,abc27,t2_abab,i_oovo_aaaa);
	for (a = 1; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		double t3ax, t3bx, dn;

		t3ax1[a*va*va+b*va+c] += -abc23[c+a*va+b*va*va]
					 +abc23[c+b*va+a*va*va]
					 -abc25[c+b*va+a*va*va]
					 +abc25[c+a*va+b*va*va]
					 -abc27[a+c*va+b*va*va]
					 +abc27[b+c*va+a*va*va];
		t3bx = +comp_t3b_ijkabc(oa,va,i,j,k,a,b,c,
			   t1,i_oovv_abab,f_ov,t2_abab)
		       -comp_t3b_ijkabc(oa,va,i,j,k,b,a,c,
			   t1,i_oovv_abab,f_ov,t2_abab)
		       -comp_t3b_ijkabc(oa,va,j,i,k,a,b,c,
			   t1,i_oovv_abab,f_ov,t2_abab)
		       +comp_t3b_ijkabc(oa,va,j,i,k,b,a,c,
			   t1,i_oovv_abab,f_ov,t2_abab)
		       +comp_t3b_ijkabc(oa,va,k,j,i,c,b,a,
			   t1,i_oovv_aaaa,f_ov,t2_aaaa);
		dn = d_ov[i*va+a] + d_ov[j*va+b] + d_ov[k*va+c];
		t3ax = t3ax1[a*va*va+b*va+c];
		e_pt += t3ax * (t3ax-t3bx) / dn;
	}}}
	}}
	free(ij);
	free(t3ax1);
}
#ifdef WITH_MPI
	MPI_Allreduce(MPI_IN_PLACE, &e_pt, 1, MPI_DOUBLE,
	    MPI_SUM, MPI_COMM_WORLD);
#endif
	return (e_pt);
}

double
cc_rpt(size_t oa, size_t va, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2, const double *i_oovo,
    const double *i_oovv, const double *i_ovvv)
{
	double e_pt1, e_pt2;
	const double *t2_aaaa = t2;
	const double *t2_abab = t2 + oa*oa*va*va;
	const double *i_ovvv_aaaa = i_ovvv;
	const double *i_ovvv_abab = i_ovvv + oa*va*va*(va-1)/2;
	const double *i_oovo_aaaa = i_oovo;
	const double *i_oovo_abab = i_oovo + oa*oa*oa*va;
	const double *i_oovv_aaaa = i_oovv;
	const double *i_oovv_abab = i_oovv + oa*oa*va*va;

	if (oa < 2 || va < 2)
		return (0.0);

	e_pt1 = cc_pt_aaa(oa, va, d_ov, f_ov, t1, t2_aaaa,
	    i_oovo_aaaa, i_oovv_aaaa, i_ovvv_aaaa);
	e_pt2 = cc_pt_aab(oa, oa, va, va, d_ov, f_ov, t1,
	    t2_aaaa, t2_abab, i_oovo_aaaa, i_oovo_abab,
	    i_oovv_aaaa, i_oovv_abab, i_ovvv_aaaa, i_ovvv_abab);

	return 2.0 * (e_pt1 + e_pt2);
}

double
cc_upt(size_t oa, size_t ob, size_t va, size_t vb, const double *d_ov,
    const double *f_ov, const double *t1, const double *t2,
    const double *i_oovo, const double *i_oovv, const double *i_ovvv)
{
	size_t o = oa + ob;
	size_t v = va + vb;

	if (o < 2 || v < 2)
		return (0.0);
	return cc_pt_aaa(o, v, d_ov, f_ov, t1, t2, i_oovo, i_oovv, i_ovvv);

	// e1 = cc_pt_aaa(aaaa);
	// e2 = cc_pt_aaa(bbbb);
	// e3 = cc_pt_aab(aaaa, abab);
	// e4 = cc_pt_aab(bbbb, baba);
}

double
cc_gft(size_t o, size_t v, const double *d_ov, const double *f2_ov,
    const double *l1, const double *t2, const double *l2, const double *i_oovv,
    const double *i2_t2f2_oovo, const double *i3_ovvv, const double *i6_oovo,
    const double *i7_ovvv)
{
	double e_pt = 0.0;
	int rank = 0, size = 1;

#ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
{
	size_t i, j, k, a, b, c, it, *ij, nij = 0;
	double *sigvvvl, *sigvvvr, *abc1, *abc2, *abc3;

	if ((ij = malloc(o*(o-1)*sizeof(size_t))) == NULL)
		err(1, "libpt malloc ij");
	for (i = 0, it = 0; i < o; i++) {
		for (j = i+1; j < o; j++, it++) {
			if ((int)it % size == rank) {
				ij[2*nij+0] = i;
				ij[2*nij+1] = j;
				nij++;
			}
		}
	}

	if ((sigvvvl = malloc(5*v*v*v*sizeof(*sigvvvl))) == NULL)
		err(1, "libpt malloc work");
	sigvvvr = sigvvvl + v*v*v;
	abc1 = sigvvvl + 2*v*v*v;
	abc2 = sigvvvl + 3*v*v*v;
	abc3 = sigvvvl + 4*v*v*v;

#ifdef _OPENMP
#pragma omp for reduction(+:e_pt) schedule(dynamic)
#endif
	for (it = 0; it < nij; it++) {
		i = ij[2*it+0];
		j = ij[2*it+1];
	for (k = j+1; k < o; k++) {

	t2_i_ovvv_half(o,v,i,j,k,abc1,l2,i7_ovvv);
	t2_i_ovvv_half(o,v,k,j,i,abc2,l2,i7_ovvv);
	t2_i_ovvv_half(o,v,i,k,j,abc3,l2,i7_ovvv);
	for (a = 0; a < v; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < b; c++) {
		sigvvvl[a*v*v+b*v+c] =
		    asymm_ijk_ab_c_half(v,abc1,abc2,abc3,a,b,c);
	}}}

	t2_i_oovo(o,v,i,j,k,abc1,l2,i6_oovo);
	t2_i_oovo(o,v,j,i,k,abc2,l2,i6_oovo);
	t2_i_oovo(o,v,k,j,i,abc3,l2,i6_oovo);
	for (a = 0; a < v; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < b; c++) {
		sigvvvl[a*v*v+b*v+c] +=
		    asymm_ijk_a_bc(v,abc1,abc2,abc3,a,b,c);
	}}}

	t2_i_ovvv_half(o,v,i,j,k,abc1,t2,i3_ovvv);
	t2_i_ovvv_half(o,v,k,j,i,abc2,t2,i3_ovvv);
	t2_i_ovvv_half(o,v,i,k,j,abc3,t2,i3_ovvv);
	for (a = 0; a < v; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < b; c++) {
		sigvvvr[a*v*v+b*v+c] =
		    asymm_ijk_ab_c_half(v,abc1,abc2,abc3,a,b,c);
	}}}

	t2_i_oovo(o,v,i,j,k,abc1,t2,i2_t2f2_oovo);
	t2_i_oovo(o,v,j,i,k,abc2,t2,i2_t2f2_oovo);
	t2_i_oovo(o,v,k,j,i,abc3,t2,i2_t2f2_oovo);
	for (a = 0; a < v; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < b; c++) {
		double dn, l1t, sigvvvl1, sigvvvr1;

		sigvvvr[a*v*v+b*v+c] +=
		    asymm_ijk_a_bc(v,abc1,abc2,abc3,a,b,c);
		dn = d_ov[i*v+a] + d_ov[j*v+b] + d_ov[k*v+c];
		l1t = +i_jk_a_bc_ov_oovv(o,v,l1,i_oovv,i,j,k,a,b,c)
		      +i_jk_a_bc_ov_oovv(o,v,f2_ov,l2,i,j,k,a,b,c);
		sigvvvl1 = sigvvvl[a*v*v+b*v+c];
		sigvvvr1 = sigvvvr[a*v*v+b*v+c];
		e_pt += (sigvvvl1 - l1t) * sigvvvr1 / dn;
	}}}
	}}
	free(ij);
	free(sigvvvl);
}
#ifdef WITH_MPI
	MPI_Allreduce(MPI_IN_PLACE, &e_pt, 1, MPI_DOUBLE,
	    MPI_SUM, MPI_COMM_WORLD);
#endif
	return (e_pt);
}
