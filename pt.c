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
t2_baba_i_ovvv_aaaa_half(size_t oa, size_t va, size_t ob, size_t vb,
    size_t i, size_t j, size_t k, double *abc, const double *t2,
    const double *i_ovvv)
{
	const double *t2_p = &t2[i*oa*vb*va+j*vb*va];
	const double *i_ovvv_p = &i_ovvv[k*va*va*(va-1)/2];

	(void)ob; /* unused */

	/* out(i,j,k,a,b,c) = contract(d, t2(i,j,a,d), i_ovvv(k,d,b,c)) */

	gemm('T', 'T', vb, va*(va-1)/2, va, 1.0, t2_p, va,
	    i_ovvv_p, va*(va-1)/2, 0.0, abc, vb);
}

static void
t2_aaaa_i_ovvv_baba(size_t oa, size_t va, size_t ob, size_t vb,
    size_t i, size_t j, size_t k, double *abc, const double *t2,
    const double *i_ovvv)
{
	const double *t2_p = &t2[i*oa*va*va+j*va*va];
	const double *i_ovvv_p = &i_ovvv[k*va*vb*va];

	(void)ob; /* unused */

	/* out(i,j,k,a,b,c) = contract(d, t2(i,j,a,d), i_ovvv(k,d,b,c)) */

	gemm('T', 'T', va, va*vb, va, 1.0, t2_p, va,
	    i_ovvv_p, va*vb, 0.0, abc, va);
}

static void
t2_abab_i_ovvv_abab(size_t oa, size_t va, size_t ob, size_t vb,
    size_t i, size_t j, size_t k, double *abc, const double *t2,
    const double *i_ovvv)
{
	const double *t2_p = &t2[i*ob*va*vb+j*va*vb];
	const double *i_ovvv_p = &i_ovvv[k*vb*va*vb];

	(void)oa; /* unused */

	/* out(i,j,k,a,b,c) = contract(d, t2(i,j,a,d), i_ovvv(k,d,b,c)) */

	gemm('T', 'T', va, va*vb, vb, 1.0, t2_p, vb,
	    i_ovvv_p, va*vb, 0.0, abc, va);
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

static void
t2_aaaa_i_oovo_baba(size_t oa, size_t va, size_t ob, size_t vb,
    size_t i, size_t j, size_t k, double *abc, const double *t2,
    const double *i_oovo)
{
	const double *t2_p = &t2[i*oa*va*va];
	const double *i_oovo_p = &i_oovo[j*oa*vb*oa+k*vb*oa];

	(void)ob; /* unused */

	/* out(i,j,k,a,b,c) = contract(l, t2(i,l,a,b), i_oovo(j,k,c,l)) */

	gemm('N', 'N', va*va, vb, oa, 1.0, t2_p, va*va,
	    i_oovo_p, oa, 0.0, abc, va*va);
}

static void
t2_abab_i_oovo_abab(size_t oa, size_t va, size_t ob, size_t vb,
    size_t i, size_t j, size_t k, double *abc, const double *t2,
    const double *i_oovo)
{
	const double *t2_p = &t2[i*ob*va*vb];
	const double *i_oovo_p = &i_oovo[j*ob*va*ob+k*va*ob];

	(void)oa; /* unused */

	/* out(i,j,k,a,b,c) = contract(l, t2(i,l,a,b), i_oovo(j,k,c,l)) */

	gemm('N', 'N', va*vb, va, ob, 1.0, t2_p, va*vb,
	    i_oovo_p, ob, 0.0, abc, va*vb);
}

static void
t2_baba_i_oovo_aaaa(size_t oa, size_t va, size_t ob, size_t vb,
    size_t i, size_t j, size_t k, double *abc, const double *t2,
    const double *i_oovo)
{
	const double *t2_p = &t2[i*oa*vb*va];
	const double *i_oovo_p = &i_oovo[j*oa*va*oa+k*va*oa];

	(void)ob; /* unused */

	/* out(i,j,k,a,b,c) = contract(l, t2(i,l,a,b), i_oovo(j,k,c,l)) */

	gemm('N', 'N', va*vb, va, oa, 1.0, t2_p, va*vb,
	    i_oovo_p, oa, 0.0, abc, va*vb);
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
comp_t3b_ijkabc(size_t v1, size_t o2, size_t v2a, size_t v2b,
    size_t i, size_t j, size_t k, size_t a, size_t b, size_t c,
    const double *t1, const double *i_oovv, const double *f_ov,
    const double *t2)
{
	return t1[i*v1+a] * i_oovv[j*o2*v2a*v2b+k*v2a*v2b+b*v2b+c] +
	       f_ov[i*v1+a] * t2[j*o2*v2a*v2b+k*v2a*v2b+b*v2b+c];
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

	t2_i_ovvv_half(oa,va,i,j,k,abc1,t2_aaaa,i_ovvv_aaaa);
	t2_i_ovvv_half(oa,va,i,k,j,abc2,t2_aaaa,i_ovvv_aaaa);
	t2_i_ovvv_half(oa,va,k,j,i,abc3,t2_aaaa,i_ovvv_aaaa);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < b; c++) {
		t3ax1[a*va*va+b*va+c] =
		    asymm_ijk_ab_c_half(va,abc1,abc2,abc3,a,b,c);
	}}}

	t2_i_oovo(oa,va,i,j,k,abc1,t2_aaaa,i_oovo_aaaa);
	t2_i_oovo(oa,va,j,i,k,abc2,t2_aaaa,i_oovo_aaaa);
	t2_i_oovo(oa,va,k,j,i,abc3,t2_aaaa,i_oovo_aaaa);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
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
cc_pt_aab(size_t oa, size_t va, size_t ob, size_t vb,
    const double *d_ov_aa, const double *d_ov_bb,
    const double *f_ov_aa, const double *f_ov_bb,
    const double *t1_aa, const double *t1_bb,
    const double *t2_aaaa, const double *t2_abab, const double *t2_baba,
    const double *i_oovo_aaaa, const double *i_oovo_abab,
    const double *i_oovo_baba, const double *i_oovv_aaaa,
    const double *i_oovv_abab, const double *i_ovvv_aaaa,
    const double *i_ovvv_abab, const double *i_ovvv_baba)
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
	double *t3ax1, *abc1, *abc2, *abc3, *abc11, *abc12;

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

	if ((t3ax1 = malloc(4*va*va*vb*sizeof(double))) == NULL)
		err(1, "libpt malloc work");
	abc1 = t3ax1 + 1*va*va*vb;
	abc2 = t3ax1 + 2*va*va*vb;
	abc3 = t3ax1 + 3*va*va*vb;
	abc11 = t3ax1 + 1*va*va*vb;
	abc12 = t3ax1 + 1*va*va*vb + vb*va*(va-1)/2;

#ifdef _OPENMP
#pragma omp for reduction(+:e_pt) schedule(dynamic)
#endif
	for (it = 0; it < nij; it++) {
		i = ij[2*it+0];
		j = ij[2*it+1];
	for (k = 0; k < ob; k++) {

	t2_aaaa_i_ovvv_baba(oa,va,ob,vb,i,j,k,abc1,t2_aaaa,i_ovvv_baba);
	t2_abab_i_ovvv_abab(oa,va,ob,vb,i,k,j,abc2,t2_abab,i_ovvv_abab);
	t2_abab_i_ovvv_abab(oa,va,ob,vb,j,k,i,abc3,t2_abab,i_ovvv_abab);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		t3ax1[a*va*vb+b*vb+c] =
		    -abc1[a+b*va+c*va*va]
		    +abc1[b+a*va+c*va*va]
		    -abc2[a+c*va+b*va*vb]
		    +abc2[b+c*va+a*va*vb]
		    +abc3[a+c*va+b*va*vb]
		    -abc3[b+c*va+a*va*vb];
	}}}

	t2_baba_i_ovvv_aaaa_half(oa,va,ob,vb,k,j,i,abc11,t2_baba,i_ovvv_aaaa);
	t2_baba_i_ovvv_aaaa_half(oa,va,ob,vb,k,i,j,abc12,t2_baba,i_ovvv_aaaa);
	t2_aaaa_i_oovo_baba(oa,va,ob,vb,i,k,j,abc2,t2_aaaa,i_oovo_baba);
	t2_aaaa_i_oovo_baba(oa,va,ob,vb,j,k,i,abc3,t2_aaaa,i_oovo_baba);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		t3ax1[a*va*vb+b*vb+c] +=
		    -abc11[c+vb*a*(a-1)/2+vb*b]
		    +abc12[c+vb*a*(a-1)/2+vb*b]
		    -abc2[b+a*va+c*va*va]
		    +abc3[b+a*va+c*va*va];
	}}}

	t2_abab_i_oovo_abab(oa,va,ob,vb,i,j,k,abc1,t2_abab,i_oovo_abab);
	t2_abab_i_oovo_abab(oa,va,ob,vb,j,i,k,abc2,t2_abab,i_oovo_abab);
	t2_baba_i_oovo_aaaa(oa,va,ob,vb,k,j,i,abc3,t2_baba,i_oovo_aaaa);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		double t3ax, t3bx, dn;

		t3ax1[a*va*vb+b*vb+c] += -abc1[c+a*vb+b*vb*va]
					 +abc1[c+b*vb+a*vb*va]
					 -abc2[c+b*vb+a*vb*va]
					 +abc2[c+a*vb+b*vb*va]
					 -abc3[a+c*va+b*va*vb]
					 +abc3[b+c*va+a*va*vb];
		t3bx = +comp_t3b_ijkabc(va,ob,va,vb,i,j,k,a,b,c,
			   t1_aa,i_oovv_abab,f_ov_aa,t2_abab)
		       -comp_t3b_ijkabc(va,ob,va,vb,i,j,k,b,a,c,
			   t1_aa,i_oovv_abab,f_ov_aa,t2_abab)
		       -comp_t3b_ijkabc(va,ob,va,vb,j,i,k,a,b,c,
			   t1_aa,i_oovv_abab,f_ov_aa,t2_abab)
		       +comp_t3b_ijkabc(va,ob,va,vb,j,i,k,b,a,c,
			   t1_aa,i_oovv_abab,f_ov_aa,t2_abab)
		       +comp_t3b_ijkabc(vb,oa,va,va,k,j,i,c,b,a,
			   t1_bb,i_oovv_aaaa,f_ov_bb,t2_aaaa);
		dn = d_ov_aa[i*va+a] + d_ov_aa[j*va+b] + d_ov_bb[k*vb+c];
		t3ax = t3ax1[a*va*vb+b*vb+c];
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
libpt_rpt(size_t oa, size_t va, const double *d_ov, const double *f_ov,
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

	e_pt1 = cc_pt_aaa(oa, va, d_ov, f_ov, t1, t2_aaaa,
	    i_oovo_aaaa, i_oovv_aaaa, i_ovvv_aaaa);
	e_pt2 = cc_pt_aab(oa, va, oa, va, d_ov, d_ov, f_ov, f_ov, t1, t1,
	    t2_aaaa, t2_abab, t2_abab, i_oovo_aaaa, i_oovo_abab, i_oovo_abab,
	    i_oovv_aaaa, i_oovv_abab, i_ovvv_aaaa, i_ovvv_abab, i_ovvv_abab);

	return 2.0 * (e_pt1 + e_pt2);
}

double
libpt_upt(size_t oa, size_t va, size_t ob, size_t vb, const double *d_ov,
    const double *f_ov, const double *t1, const double *t2,
    const double *i_oovo, const double *i_oovv, const double *i_ovvv)
{
	double e_pt1, e_pt2, e_pt3, e_pt4;
	const double *d_ov_aa = d_ov;
	const double *d_ov_bb = d_ov_aa + oa*va;
	const double *f_ov_aa = f_ov;
	const double *f_ov_bb = f_ov_aa + oa*va;
	const double *t1_aa = t1;
	const double *t1_bb = t1_aa + oa*va;

	const double *t2_aaaa = t2;
	const double *t2_abab = t2_aaaa + oa*oa*va*va;
	const double *t2_bbbb = t2_abab + oa*ob*va*vb;
	const double *t2_baba = t2_bbbb + ob*ob*vb*vb;

	const double *i_oovo_aaaa = i_oovo;
	const double *i_oovo_abab = i_oovo_aaaa + oa*oa*va*oa;
	const double *i_oovo_bbbb = i_oovo_abab + oa*ob*va*ob;
	const double *i_oovo_baba = i_oovo_bbbb + ob*ob*vb*ob;

	const double *i_oovv_aaaa = i_oovv;
	const double *i_oovv_abab = i_oovv_aaaa + oa*oa*va*va;
	const double *i_oovv_bbbb = i_oovv_abab + oa*ob*va*vb;
	const double *i_oovv_baba = i_oovv_bbbb + ob*ob*vb*vb;

	const double *i_ovvv_aaaa = i_ovvv;
	const double *i_ovvv_abab = i_ovvv_aaaa + oa*va*va*(va-1)/2;
	const double *i_ovvv_bbbb = i_ovvv_abab + oa*vb*va*vb;
	const double *i_ovvv_baba = i_ovvv_bbbb + ob*vb*vb*(vb-1)/2;

	/* aaaaaa */
	e_pt1 = cc_pt_aaa(oa, va, d_ov_aa, f_ov_aa, t1_aa, t2_aaaa,
	    i_oovo_aaaa, i_oovv_aaaa, i_ovvv_aaaa);
	/* bbbbbb */
	e_pt2 = cc_pt_aaa(ob, vb, d_ov_bb, f_ov_bb, t1_bb, t2_bbbb,
	    i_oovo_bbbb, i_oovv_bbbb, i_ovvv_bbbb);
	/* aabaab */
	e_pt3 = cc_pt_aab(oa, va, ob, vb, d_ov_aa, d_ov_bb, f_ov_aa, f_ov_bb,
	    t1_aa, t1_bb, t2_aaaa, t2_abab, t2_baba, i_oovo_aaaa, i_oovo_abab,
	    i_oovo_baba, i_oovv_aaaa, i_oovv_abab, i_ovvv_aaaa, i_ovvv_abab,
	    i_ovvv_baba);
	/* bbabba */
	e_pt4 = cc_pt_aab(ob, vb, oa, va, d_ov_bb, d_ov_aa, f_ov_bb, f_ov_aa,
	    t1_bb, t1_aa, t2_bbbb, t2_baba, t2_abab, i_oovo_bbbb, i_oovo_baba,
	    i_oovo_abab, i_oovv_bbbb, i_oovv_baba, i_ovvv_bbbb, i_ovvv_baba,
	    i_ovvv_abab);

	return (e_pt1 + e_pt2 + e_pt3 + e_pt4);
}

static double
cc_ft_aaa(size_t oa, size_t va, const double *d_ov, const double *f2_ov,
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

	if ((sigvvvl = malloc(5*va*va*va*sizeof(*sigvvvl))) == NULL)
		err(1, "libpt malloc work");
	sigvvvr = sigvvvl + va*va*va;
	abc1 = sigvvvl + 2*va*va*va;
	abc2 = sigvvvl + 3*va*va*va;
	abc3 = sigvvvl + 4*va*va*va;

#ifdef _OPENMP
#pragma omp for reduction(+:e_pt) schedule(dynamic)
#endif
	for (it = 0; it < nij; it++) {
		i = ij[2*it+0];
		j = ij[2*it+1];
	for (k = j+1; k < oa; k++) {

	t2_i_ovvv_half(oa,va,i,j,k,abc1,l2,i7_ovvv);
	t2_i_ovvv_half(oa,va,k,j,i,abc2,l2,i7_ovvv);
	t2_i_ovvv_half(oa,va,i,k,j,abc3,l2,i7_ovvv);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < b; c++) {
		sigvvvl[a*va*va+b*va+c] =
		    asymm_ijk_ab_c_half(va,abc1,abc2,abc3,a,b,c);
	}}}

	t2_i_oovo(oa,va,i,j,k,abc1,l2,i6_oovo);
	t2_i_oovo(oa,va,j,i,k,abc2,l2,i6_oovo);
	t2_i_oovo(oa,va,k,j,i,abc3,l2,i6_oovo);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < b; c++) {
		sigvvvl[a*va*va+b*va+c] +=
		    asymm_ijk_a_bc(va,abc1,abc2,abc3,a,b,c);
	}}}

	t2_i_ovvv_half(oa,va,i,j,k,abc1,t2,i3_ovvv);
	t2_i_ovvv_half(oa,va,k,j,i,abc2,t2,i3_ovvv);
	t2_i_ovvv_half(oa,va,i,k,j,abc3,t2,i3_ovvv);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < b; c++) {
		sigvvvr[a*va*va+b*va+c] =
		    asymm_ijk_ab_c_half(va,abc1,abc2,abc3,a,b,c);
	}}}

	t2_i_oovo(oa,va,i,j,k,abc1,t2,i2_t2f2_oovo);
	t2_i_oovo(oa,va,j,i,k,abc2,t2,i2_t2f2_oovo);
	t2_i_oovo(oa,va,k,j,i,abc3,t2,i2_t2f2_oovo);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < b; c++) {
		double dn, l1t, sigvvvl1, sigvvvr1;

		sigvvvr[a*va*va+b*va+c] +=
		    asymm_ijk_a_bc(va,abc1,abc2,abc3,a,b,c);
		dn = d_ov[i*va+a] + d_ov[j*va+b] + d_ov[k*va+c];
		l1t = +i_jk_a_bc_ov_oovv(oa,va,l1,i_oovv,i,j,k,a,b,c)
		      +i_jk_a_bc_ov_oovv(oa,va,f2_ov,l2,i,j,k,a,b,c);
		sigvvvl1 = sigvvvl[a*va*va+b*va+c];
		sigvvvr1 = sigvvvr[a*va*va+b*va+c];
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

static double
cc_ft_aab(size_t oa, size_t va, size_t ob, size_t vb,
    const double *d_ov_aa, const double *d_ov_bb,
    const double *f2_ov_aa, const double *f2_ov_bb,
    const double *l1_aa, const double *l1_bb,
    const double *t2_aaaa, const double *t2_abab, const double *t2_baba,
    const double *l2_aaaa, const double *l2_abab, const double *l2_baba,
    const double *i_oovv_aaaa, const double *i_oovv_abab,
    const double *i2_t2f2_oovo_aaaa, const double *i2_t2f2_oovo_abab,
    const double *i2_t2f2_oovo_baba,
    const double *i3_ovvv_aaaa, const double *i3_ovvv_abab,
    const double *i3_ovvv_baba,
    const double *i6_oovo_aaaa, const double *i6_oovo_abab,
    const double *i6_oovo_baba,
    const double *i7_ovvv_aaaa, const double *i7_ovvv_abab,
    const double *i7_ovvv_baba)
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
	double *sigvvvl, *sigvvvr, *abc1, *abc2, *abc3, *abc11, *abc12;

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

	if ((sigvvvl = malloc(5*va*va*vb*sizeof(double))) == NULL)
		err(1, "libpt malloc work");
	sigvvvr = sigvvvl + 1*va*va*vb;
	abc1 = sigvvvl + 2*va*va*vb;
	abc2 = sigvvvl + 3*va*va*vb;
	abc3 = sigvvvl + 4*va*va*vb;
	abc11 = sigvvvl + 2*va*va*vb;
	abc12 = sigvvvl + 2*va*va*vb + vb*va*(va-1)/2;

#ifdef _OPENMP
#pragma omp for reduction(+:e_pt) schedule(dynamic)
#endif
	for (it = 0; it < nij; it++) {
		i = ij[2*it+0];
		j = ij[2*it+1];
	for (k = 0; k < ob; k++) {

	t2_aaaa_i_ovvv_baba(oa,va,ob,vb,i,j,k,abc1,l2_aaaa,i7_ovvv_baba);
	t2_abab_i_ovvv_abab(oa,va,ob,vb,i,k,j,abc2,l2_abab,i7_ovvv_abab);
	t2_abab_i_ovvv_abab(oa,va,ob,vb,j,k,i,abc3,l2_abab,i7_ovvv_abab);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		sigvvvl[a*va*vb+b*vb+c] =
		    -abc1[a+b*va+c*va*va]
		    +abc1[b+a*va+c*va*va]
		    -abc2[a+c*va+b*va*vb]
		    +abc2[b+c*va+a*va*vb]
		    +abc3[a+c*va+b*va*vb]
		    -abc3[b+c*va+a*va*vb];
	}}}

	t2_baba_i_ovvv_aaaa_half(oa,va,ob,vb,k,j,i,abc11,l2_baba,i7_ovvv_aaaa);
	t2_baba_i_ovvv_aaaa_half(oa,va,ob,vb,k,i,j,abc12,l2_baba,i7_ovvv_aaaa);
	t2_aaaa_i_oovo_baba(oa,va,ob,vb,i,k,j,abc2,l2_aaaa,i6_oovo_baba);
	t2_aaaa_i_oovo_baba(oa,va,ob,vb,j,k,i,abc3,l2_aaaa,i6_oovo_baba);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		sigvvvl[a*va*vb+b*vb+c] +=
		    -abc11[c+vb*a*(a-1)/2+vb*b]
		    +abc12[c+vb*a*(a-1)/2+vb*b]
		    -abc2[b+a*va+c*va*va]
		    +abc3[b+a*va+c*va*va];
	}}}

	t2_abab_i_oovo_abab(oa,va,ob,vb,i,j,k,abc1,l2_abab,i6_oovo_abab);
	t2_abab_i_oovo_abab(oa,va,ob,vb,j,i,k,abc2,l2_abab,i6_oovo_abab);
	t2_baba_i_oovo_aaaa(oa,va,ob,vb,k,j,i,abc3,l2_baba,i6_oovo_aaaa);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		sigvvvl[a*va*vb+b*vb+c] +=
		    -abc1[c+a*vb+b*vb*va]
		    +abc1[c+b*vb+a*vb*va]
		    -abc2[c+b*vb+a*vb*va]
		    +abc2[c+a*vb+b*vb*va]
		    -abc3[a+c*va+b*va*vb]
		    +abc3[b+c*va+a*va*vb];
	}}}

	t2_aaaa_i_ovvv_baba(oa,va,ob,vb,i,j,k,abc1,t2_aaaa,i3_ovvv_baba);
	t2_abab_i_ovvv_abab(oa,va,ob,vb,i,k,j,abc2,t2_abab,i3_ovvv_abab);
	t2_abab_i_ovvv_abab(oa,va,ob,vb,j,k,i,abc3,t2_abab,i3_ovvv_abab);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		sigvvvr[a*va*vb+b*vb+c] =
		    -abc1[a+b*va+c*va*va]
		    +abc1[b+a*va+c*va*va]
		    -abc2[a+c*va+b*va*vb]
		    +abc2[b+c*va+a*va*vb]
		    +abc3[a+c*va+b*va*vb]
		    -abc3[b+c*va+a*va*vb];
	}}}

	t2_baba_i_ovvv_aaaa_half(oa,va,ob,vb,k,j,i,abc11,t2_baba,i3_ovvv_aaaa);
	t2_baba_i_ovvv_aaaa_half(oa,va,ob,vb,k,i,j,abc12,t2_baba,i3_ovvv_aaaa);
	t2_aaaa_i_oovo_baba(oa,va,ob,vb,i,k,j,abc2,t2_aaaa,i2_t2f2_oovo_baba);
	t2_aaaa_i_oovo_baba(oa,va,ob,vb,j,k,i,abc3,t2_aaaa,i2_t2f2_oovo_baba);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		sigvvvr[a*va*vb+b*vb+c] +=
		    -abc11[c+vb*a*(a-1)/2+vb*b]
		    +abc12[c+vb*a*(a-1)/2+vb*b]
		    -abc2[b+a*va+c*va*va]
		    +abc3[b+a*va+c*va*va];
	}}}

	t2_abab_i_oovo_abab(oa,va,ob,vb,i,j,k,abc1,t2_abab,i2_t2f2_oovo_abab);
	t2_abab_i_oovo_abab(oa,va,ob,vb,j,i,k,abc2,t2_abab,i2_t2f2_oovo_abab);
	t2_baba_i_oovo_aaaa(oa,va,ob,vb,k,j,i,abc3,t2_baba,i2_t2f2_oovo_aaaa);
	for (a = 0; a < va; a++) {
	for (b = 0; b < a; b++) {
	for (c = 0; c < vb; c++) {
		double sigvvvl1, sigvvvr1, l1t, dn;

		sigvvvr[a*va*vb+b*vb+c] +=
		    -abc1[c+a*vb+b*vb*va]
		    +abc1[c+b*vb+a*vb*va]
		    -abc2[c+b*vb+a*vb*va]
		    +abc2[c+a*vb+b*vb*va]
		    -abc3[a+c*va+b*va*vb]
		    +abc3[b+c*va+a*va*vb];
		l1t = +comp_t3b_ijkabc(va,ob,va,vb,i,j,k,a,b,c,
			  l1_aa,i_oovv_abab,f2_ov_aa,l2_abab)
		      -comp_t3b_ijkabc(va,ob,va,vb,i,j,k,b,a,c,
			  l1_aa,i_oovv_abab,f2_ov_aa,l2_abab)
		      -comp_t3b_ijkabc(va,ob,va,vb,j,i,k,a,b,c,
			  l1_aa,i_oovv_abab,f2_ov_aa,l2_abab)
		      +comp_t3b_ijkabc(va,ob,va,vb,j,i,k,b,a,c,
			  l1_aa,i_oovv_abab,f2_ov_aa,l2_abab)
		      +comp_t3b_ijkabc(vb,oa,va,va,k,j,i,c,b,a,
			  l1_bb,i_oovv_aaaa,f2_ov_bb,l2_aaaa);
		dn = d_ov_aa[i*va+a] + d_ov_aa[j*va+b] + d_ov_bb[k*vb+c];
		sigvvvl1 = sigvvvl[a*va*vb+b*vb+c];
		sigvvvr1 = sigvvvr[a*va*vb+b*vb+c];
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

double
libpt_rft(size_t oa, size_t va, const double *d_ov, const double *f2_ov,
    const double *l1, const double *t2, const double *l2, const double *i_oovv,
    const double *i2_t2f2_oovo, const double *i3_ovvv, const double *i6_oovo,
    const double *i7_ovvv)
{
	double e_pt1, e_pt2;
	const double *t2_aaaa = t2;
	const double *t2_abab = t2 + oa*oa*va*va;
	const double *l2_aaaa = l2;
	const double *l2_abab = l2 + oa*oa*va*va;
	const double *i_oovv_aaaa = i_oovv;
	const double *i_oovv_abab = i_oovv + oa*oa*va*va;
	const double *i2_t2f2_oovo_aaaa = i2_t2f2_oovo;
	const double *i2_t2f2_oovo_abab = i2_t2f2_oovo + oa*oa*oa*va;
	const double *i3_ovvv_aaaa = i3_ovvv;
	const double *i3_ovvv_abab = i3_ovvv + oa*va*va*(va-1)/2;
	const double *i6_oovo_aaaa = i6_oovo;
	const double *i6_oovo_abab = i6_oovo + oa*oa*oa*va;
	const double *i7_ovvv_aaaa = i7_ovvv;
	const double *i7_ovvv_abab = i7_ovvv + oa*va*va*(va-1)/2;

	e_pt1 = cc_ft_aaa(oa, va, d_ov, f2_ov, l1, t2_aaaa, l2_aaaa,
	    i_oovv_aaaa, i2_t2f2_oovo_aaaa, i3_ovvv_aaaa, i6_oovo_aaaa,
	    i7_ovvv_aaaa);
	e_pt2 = cc_ft_aab(oa, va, oa, va, d_ov, d_ov, f2_ov, f2_ov,
	    l1, l1, t2_aaaa, t2_abab, t2_abab, l2_aaaa, l2_abab, l2_abab,
	    i_oovv_aaaa, i_oovv_abab, i2_t2f2_oovo_aaaa, i2_t2f2_oovo_abab,
	    i2_t2f2_oovo_abab, i3_ovvv_aaaa, i3_ovvv_abab, i3_ovvv_abab,
	    i6_oovo_aaaa, i6_oovo_abab, i6_oovo_abab,
	    i7_ovvv_aaaa, i7_ovvv_abab, i7_ovvv_abab);

	return 2.0 * (e_pt1 + e_pt2);
}

double
libpt_uft(size_t oa, size_t va, size_t ob, size_t vb, const double *d_ov,
    const double *f2_ov, const double *l1, const double *t2, const double *l2,
    const double *i_oovv, const double *i2_t2f2_oovo, const double *i3_ovvv,
    const double *i6_oovo, const double *i7_ovvv)
{
	double e_pt1, e_pt2, e_pt3, e_pt4;
	const double *d_ov_aa = d_ov;
	const double *d_ov_bb = d_ov_aa + oa*va;
	const double *f2_ov_aa = f2_ov;
	const double *f2_ov_bb = f2_ov_aa + oa*va;
	const double *l1_aa = l1;
	const double *l1_bb = l1_aa + oa*va;

	const double *t2_aaaa = t2;
	const double *t2_abab = t2_aaaa + oa*oa*va*va;
	const double *t2_bbbb = t2_abab + oa*ob*va*vb;
	const double *t2_baba = t2_bbbb + ob*ob*vb*vb;

	const double *l2_aaaa = l2;
	const double *l2_abab = l2_aaaa + oa*oa*va*va;
	const double *l2_bbbb = l2_abab + oa*ob*va*vb;
	const double *l2_baba = l2_bbbb + ob*ob*vb*vb;

	const double *i_oovv_aaaa = i_oovv;
	const double *i_oovv_abab = i_oovv_aaaa + oa*oa*va*va;
	const double *i_oovv_bbbb = i_oovv_abab + oa*ob*va*vb;
	const double *i_oovv_baba = i_oovv_bbbb + ob*ob*vb*vb;

	const double *i2_t2f2_oovo_aaaa = i2_t2f2_oovo;
	const double *i2_t2f2_oovo_abab = i2_t2f2_oovo_aaaa + oa*oa*va*oa;
	const double *i2_t2f2_oovo_bbbb = i2_t2f2_oovo_abab + oa*ob*va*ob;
	const double *i2_t2f2_oovo_baba = i2_t2f2_oovo_bbbb + ob*ob*vb*ob;

	const double *i3_ovvv_aaaa = i3_ovvv;
	const double *i3_ovvv_abab = i3_ovvv_aaaa + oa*va*va*(va-1)/2;
	const double *i3_ovvv_bbbb = i3_ovvv_abab + oa*vb*va*vb;
	const double *i3_ovvv_baba = i3_ovvv_bbbb + ob*vb*vb*(vb-1)/2;

	const double *i6_oovo_aaaa = i6_oovo;
	const double *i6_oovo_abab = i6_oovo_aaaa + oa*oa*va*oa;
	const double *i6_oovo_bbbb = i6_oovo_abab + oa*ob*va*ob;
	const double *i6_oovo_baba = i6_oovo_bbbb + ob*ob*vb*ob;

	const double *i7_ovvv_aaaa = i7_ovvv;
	const double *i7_ovvv_abab = i7_ovvv_aaaa + oa*va*va*(va-1)/2;
	const double *i7_ovvv_bbbb = i7_ovvv_abab + oa*vb*va*vb;
	const double *i7_ovvv_baba = i7_ovvv_bbbb + ob*vb*vb*(vb-1)/2;

	/* aaaaaa */
	e_pt1 = cc_ft_aaa(oa, va, d_ov_aa, f2_ov_aa, l1_aa, t2_aaaa, l2_aaaa,
	    i_oovv_aaaa, i2_t2f2_oovo_aaaa, i3_ovvv_aaaa, i6_oovo_aaaa,
	    i7_ovvv_aaaa);
	/* bbbbbb */
	e_pt2 = cc_ft_aaa(ob, vb, d_ov_bb, f2_ov_bb, l1_bb, t2_bbbb, l2_bbbb,
	    i_oovv_bbbb, i2_t2f2_oovo_bbbb, i3_ovvv_bbbb, i6_oovo_bbbb,
	    i7_ovvv_bbbb);
	/* aabaab */
	e_pt3 = cc_ft_aab(oa, va, ob, vb, d_ov_aa, d_ov_bb, f2_ov_aa, f2_ov_bb,
	    l1_aa, l1_bb, t2_aaaa, t2_abab, t2_baba, l2_aaaa, l2_abab, l2_baba,
	    i_oovv_aaaa, i_oovv_abab, i2_t2f2_oovo_aaaa, i2_t2f2_oovo_abab,
	    i2_t2f2_oovo_baba, i3_ovvv_aaaa, i3_ovvv_abab, i3_ovvv_baba,
	    i6_oovo_aaaa, i6_oovo_abab, i6_oovo_baba,
	    i7_ovvv_aaaa, i7_ovvv_abab, i7_ovvv_baba);
	/* bbabba */
	e_pt4 = cc_ft_aab(ob, vb, oa, va, d_ov_bb, d_ov_aa, f2_ov_bb, f2_ov_aa,
	    l1_bb, l1_aa, t2_bbbb, t2_baba, t2_abab, l2_bbbb, l2_baba, l2_abab,
	    i_oovv_bbbb, i_oovv_baba, i2_t2f2_oovo_bbbb, i2_t2f2_oovo_baba,
	    i2_t2f2_oovo_abab, i3_ovvv_bbbb, i3_ovvv_baba, i3_ovvv_abab,
	    i6_oovo_bbbb, i6_oovo_baba, i6_oovo_abab,
	    i7_ovvv_bbbb, i7_ovvv_baba, i7_ovvv_abab);

	return (e_pt1 + e_pt2 + e_pt3 + e_pt4);
}
