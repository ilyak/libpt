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

struct i4 {
	uint16_t a, b, c, d;
};

struct st4 {
	size_t len;
	struct i4 *idx;
	double *data;
};

/* Compute CCSD(T) energy correction in parallel.
 *
 * Arguments:
 *   o - size of occupied space
 *   v - size of virtual space
 *   d_ov - Delta matrix (size o*v)
 *   f_ov - Fock matrix (size o*v)
 *   t1 - CCSD T1 amplitudes (size o*v)
 *   t2 - CCSD T2 amplitudes (sparse tensor, full size o*o*v*v)
 *   i_ooov - OOOV integrals (sparse tensor, full size o*o*o*v)
 *   i_oovv - OOVV integrals (sparse tensor, full size o*o*v*v)
 *   i_ovvv - OVVV integrals (sparse tensor, full size o*v*v*v)
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
double ccsd_pt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *t1, const struct st4 *t2, const struct st4 *i_ooov,
    const struct st4 *i_oovv, const struct st4 *i_ovvv);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBPT_PT_H */
