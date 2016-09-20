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

#ifndef LIBPT_PT_H
#define LIBPT_PT_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Compute CCSD(T) energy correction in parallel.
 *
 * The routine is MPI/OpenMP parallel. All MPI processes must receive same
 * input data.
 *
 * Arguments:
 *   o - size of spin-block occupied space
 *   v - size of spin-block virtual space
 *   d_ov - Delta matrix (size o*v)
 *   f_ov - Fock matrix (size o*v)
 *   t1 - CCSD T1 amplitudes (size o*v)
 *   t2 - CCSD T2 amplitudes, aaaa and abab blocks (size 2*o*o*v*v)
 *   i_oovo - OOOV integrals transposed, aaaa and abab blocks (size 2*o*o*o*v)
 *   i_oovv - OOVV integrals, aaaa and abab blocks (size 2*o*o*v*v)
 *   i_ovvv - OVVV integrals, aaaa and abab blocks (size 2*o*v*v*v)
 *
 * All arrays should be arranged contiguously in memory by last index first.
 * E.g., for d_ov the first v contiguous elements in memory are d_ov[o=0,v=0],
 * d_ov[o=0,v=1], d_ov[o=0,v=2], and so forth. The tensors are expected to
 * be properly (anti-)symmetrized.
 *
 * The function returns CCSD(T) energy correction.
 *
 * References:
 *   J. Chem. Phys. 98, 8718 (1993); http://dx.doi.org/10.1063/1.464480
 */
double ccsd_rpt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2, const double *i_oovo,
    const double *i_oovv, const double *i_ovvv);

/* Compute CCSD(T) energy correction in parallel for the unrestricted case.
 *
 * The routine is MPI/OpenMP parallel. All MPI processes must receive same
 * input data.
 *
 * Arguments:
 *   o - full size of occupied space
 *   v - full size of virtual space
 *   d_ov - Delta matrix (size o*v)
 *   f_ov - Fock matrix (size o*v)
 *   t1 - CCSD T1 amplitudes (size o*v)
 *   t2 - CCSD T2 amplitudes (size o*o*v*v)
 *   i_oovo - OOOV integrals transposed (size o*o*o*v)
 *   i_oovv - OOVV integrals (size o*o*v*v)
 *   i_ovvv - OVVV integrals (size o*v*v*v)
 *
 * All arrays should be arranged contiguously in memory by last index first.
 * E.g., for d_ov the first v contiguous elements in memory are d_ov[o=0,v=0],
 * d_ov[o=0,v=1], d_ov[o=0,v=2], and so forth. The tensors are expected to
 * be properly (anti-)symmetrized.
 *
 * The function returns CCSD(T) energy correction.
 *
 * References:
 *   J. Chem. Phys. 98, 8718 (1993); http://dx.doi.org/10.1063/1.464480
 */
double ccsd_upt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *t1, const double *t2, const double *i_oovo,
    const double *i_oovv, const double *i_ovvv);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBPT_PT_H */
