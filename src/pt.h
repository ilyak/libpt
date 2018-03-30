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

#ifndef LIBPT_PT_H
#define LIBPT_PT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Print libpt banner to standard output. */
void libpt_print_banner(void);

void libpt_set_malloc(void *(*fn)(size_t));

void libpt_set_free(void (*fn)(void *));

/* Compute coupled-cluster (T) energy correction for the restricted case.
 *
 * This routine is MPI/OpenMP parallel. All MPI processes must receive same
 * input data.
 *
 * Arguments:
 *   oa - size of alpha (equals beta) spin-block occupied space
 *   va - size of alpha (equals beta) spin-block virtual space
 *   d_ov - Delta matrix (size oa*va)
 *   f_ov - Fock matrix (size oa*va)
 *   t1 - CCSD T1 amplitudes (size oa*va)
 *   t2 - CCSD T2 amplitudes, aaaa/abab blocks (size 2*oa*oa*va*va)
 *   i_oovo - OOVO integrals, aaaa/abab blocks (size 2*oa*oa*va*oa)
 *   i_oovv - OOVV integrals, aaaa/abab blocks (size 2*oa*oa*va*va)
 *   i_ovvv - OVVV integrals, aaaa block with vv symmetry followed by full
 *     abab block (size oa*va*va*(va-1)/2 + oa*va*va*va)
 *
 * All arrays should be arranged contiguously in memory by last index first.
 * E.g., for d_ov the first v contiguous elements in memory are d_ov[o=0,v=0],
 * d_ov[o=0,v=1], d_ov[o=0,v=2], and so forth. The tensors are expected to
 * be properly (anti-)symmetrized.
 *
 * The function returns coupled-cluster (T) energy correction.
 *
 * References:
 *   J. Chem. Phys. 98, 8718 (1993); http://dx.doi.org/10.1063/1.464480
 */
double libpt_rpt(size_t oa, size_t va, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2, const double *i_oovo,
    const double *i_oovv, const double *i_ovvv);

/* Compute coupled-cluster (T) energy correction for the unrestricted case.
 *
 * This routine is MPI/OpenMP parallel. All MPI processes must receive same
 * input data.
 *
 * Arguments:
 *   oa - size of alpha spin-block occupied space
 *   va - size of alpha spin-block virtual space
 *   ob - size of beta spin-block occupied space
 *   vb - size of beta spin-block virtual space
 *   d_ov - Delta matrix (size oa*va + ob*vb)
 *   f_ov - Fock matrix (size oa*va + ob*vb)
 *   t1 - CCSD T1 amplitudes (size oa*va + ob*vb)
 *   t2 - CCSD T2 amplitudes aaaa/abab/bbbb/baba spin-blocks
 *     (size oa*oa*va*va + oa*ob*va*vb + ob*ob*vb*vb + ob*oa*vb*va)
 *   i_oovo - OOVO integrals aaaa/abab/bbbb/baba spin-blocks
 *     (size oa*oa*va*oa + oa*ob*va*ob + ob*ob*vb*ob + ob*oa*vb*oa)
 *   i_oovv - OOVV integrals aaaa/abab/bbbb/baba spin-blocks
 *     (size oa*oa*va*va + oa*ob*va*vb + ob*ob*vb*vb + ob*oa*vb*va)
 *   i_ovvv - OVVV integrals aaaa/abab/bbbb/baba spin-blocks
 *     (size oa*va*va*(va-1)/2 + oa*vb*va*vb + ob*vb*vb*(vb-1)/2 + ob*va*vb*va)
 *
 * All arrays should be arranged contiguously in memory by last index first.
 * E.g., for d_ov the first v contiguous elements in memory are d_ov[o=0,v=0],
 * d_ov[o=0,v=1], d_ov[o=0,v=2], and so forth. The tensors are expected to
 * be properly (anti-)symmetrized.
 *
 * The function returns coupled-cluster (T) energy correction.
 *
 * References:
 *   J. Chem. Phys. 98, 8718 (1993); http://dx.doi.org/10.1063/1.464480
 */
double libpt_upt(size_t oa, size_t va, size_t ob, size_t vb, const double *d_ov,
    const double *f_ov, const double *t1, const double *t2,
    const double *i_oovo, const double *i_oovv, const double *i_ovvv);

/* Compute coupled-cluster (fT) energy correction for the restricted case.
 *
 * This routine is MPI/OpenMP parallel. All MPI processes must receive same
 * input data.
 *
 * See the reference paper for the description of the intermediates passed to
 * this function as arguments. The memory layout of integrals is the same as
 * for restricted (T) correction.
 *
 * aaaa spin-blocks of i3_ovvv and i7_ovvv are stored with symmetry
 * (size o*v*v*(v-1)/2)
 *
 * All arrays should be arranged contiguously in memory by last index first.
 * E.g., for d_ov the first v contiguous elements in memory are d_ov[o=0,v=0],
 * d_ov[o=0,v=1], d_ov[o=0,v=2], and so forth. The tensors are expected to
 * be properly (anti-)symmetrized.
 *
 * The function returns the coupled-cluster (fT) energy correction.
 *
 * References:
 *   J. Chem. Phys. 129, 194105 (2008); http://dx.doi.org/10.1063/1.3013087
 */
double libpt_rft(size_t oa, size_t va, const double *d_ov, const double *f2_ov,
    const double *l1, const double *t2, const double *l2, const double *i_oovv,
    const double *i2_t2f2_oovo, const double *i3_ovvv, const double *i6_oovo,
    const double *i7_ovvv);

/* Compute coupled-cluster (fT) energy correction for the unrestricted case.
 *
 * This routine is MPI/OpenMP parallel. All MPI processes must receive same
 * input data.
 *
 * See the reference paper for the description of the intermediates passed to
 * this function as arguments. The memory layout of integrals is the same as
 * for unrestricted (T) correction.
 *
 * aaaa and bbbb spin-blocks of i3_ovvv and i7_ovvv are stored with symmetry
 * (size o*v*v*(v-1)/2)
 *
 * All arrays should be arranged contiguously in memory by last index first.
 * E.g., for d_ov the first v contiguous elements in memory are d_ov[o=0,v=0],
 * d_ov[o=0,v=1], d_ov[o=0,v=2], and so forth. The tensors are expected to
 * be properly (anti-)symmetrized.
 *
 * The function returns the coupled-cluster (fT) energy correction.
 *
 * References:
 *   J. Chem. Phys. 129, 194105 (2008); http://dx.doi.org/10.1063/1.3013087
 */
double libpt_uft(size_t oa, size_t va, size_t ob, size_t vb, const double *d_ov,
    const double *f2_ov, const double *l1, const double *t2, const double *l2,
    const double *i_oovv, const double *i2_t2f2_oovo, const double *i3_ovvv,
    const double *i6_oovo, const double *i7_ovvv);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBPT_PT_H */
