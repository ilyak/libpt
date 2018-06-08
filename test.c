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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef LIBPT_USE_MPI
#include <mpi.h>
#endif

#include "pt.h"

#define EPSILON 1.0e-8

static void *
xmalloc(size_t size)
{
	void *ptr;

	if ((ptr = malloc(size)) == NULL)
		err(1, "malloc");
	return ptr;
}

static void
read_test_header(FILE *fp, int unrestricted, size_t *oa, size_t *va,
    size_t *ob, size_t *vb, double *e_ref)
{
	int rc;

	if (unrestricted) {
		rc = fscanf(fp, "%zu %zu %zu %zu %lf", oa, va, ob, vb, e_ref);
		if (rc != 5)
			errx(1, "unable to read test header");
	} else {
		rc = fscanf(fp, "%zu %zu %lf", oa, va, e_ref);
		if (rc != 3)
			errx(1, "unable to read test header");
		*ob = *oa;
		*vb = *va;
	}
}

static void
skip_line(FILE *fp)
{
	int ch;

	while ((ch = fgetc(fp)) != '\n')
		if (ch == EOF)
			errx(1, "unexpected end of file");
}

static double
read_next_double(FILE *fp)
{
	double el;

	if (fscanf(fp, "%lf", &el) != 1)
		errx(1, "error parsing test file data");
	return (el);
}

static void
read_test_data_rpt(FILE *fp, size_t oa, size_t va, double *d_ov, double *f_ov,
    double *t1, double *t2, double *i_oovo, double *i_oovv, double *i_ovvv)
{
	size_t i;

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*va; i++)
		d_ov[i] = read_next_double(fp);
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*va; i++)
		f_ov[i] = read_next_double(fp);
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*va; i++)
		t1[i] = read_next_double(fp);
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*va*va; i++)
		t2[i] = read_next_double(fp);
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*va*va; i++)
		t2[oa*oa*va*va+i] = read_next_double(fp);
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*oa*va; i++)
		i_oovo[i] = read_next_double(fp);
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*oa*va; i++)
		i_oovo[oa*oa*oa*va+i] = read_next_double(fp);
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*va*va; i++)
		i_oovv[i] = read_next_double(fp);
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*va*va; i++)
		i_oovv[oa*oa*va*va+i] = read_next_double(fp);
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*va*va*(va-1)/2; i++)
		i_ovvv[i] = read_next_double(fp);
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*va*va*va; i++)
		i_ovvv[oa*va*va*(va-1)/2+i] = read_next_double(fp);
}

static void
read_test_data_upt(FILE *fp, size_t oa, size_t va, size_t ob, size_t vb,
    double *d_ov, double *f_ov, double *t1, double *t2, double *i_oovo,
    double *i_oovv, double *i_ovvv)
{
	size_t i, j, k, a, b, c;
	size_t o = oa+ob;
	size_t v = va+vb;
	size_t x = o > v ? o : v;
	double *tmp;

	tmp = xmalloc(x*x*x*x*sizeof(*tmp));

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa; i++) {
		for (j = 0; j < va; j++)
			*d_ov++ = read_next_double(fp);
		for (j = 0; j < vb; j++)
			read_next_double(fp);
	}
	for (i = 0; i < ob; i++) {
		for (j = 0; j < va; j++)
			read_next_double(fp);
		for (j = 0; j < vb; j++)
			*d_ov++ = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa; i++) {
		for (j = 0; j < va; j++)
			*f_ov++ = read_next_double(fp);
		for (j = 0; j < vb; j++)
			read_next_double(fp);
	}
	for (i = 0; i < ob; i++) {
		for (j = 0; j < va; j++)
			read_next_double(fp);
		for (j = 0; j < vb; j++)
			*f_ov++ = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa; i++) {
		for (j = 0; j < va; j++)
			*t1++ = read_next_double(fp);
		for (j = 0; j < vb; j++)
			read_next_double(fp);
	}
	for (i = 0; i < ob; i++) {
		for (j = 0; j < va; j++)
			read_next_double(fp);
		for (j = 0; j < vb; j++)
			*t1++ = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
		*t2++ = tmp[i*o*v*v+j*v*v+a*v+b];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
		*t2++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < vb; b++)
		*t2++ = tmp[(i+oa)*o*v*v+(j+oa)*v*v+(a+va)*v+(b+va)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < va; b++)
		*t2++ = tmp[(i+oa)*o*v*v+j*v*v+(a+va)*v+b];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*o; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (k = 0; k < oa; k++)
		*i_oovo++ = tmp[i*o*v*o+j*v*o+a*o+k];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (k = 0; k < ob; k++)
		*i_oovo++ = tmp[i*o*v*o+(j+oa)*v*o+a*o+(k+oa)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < vb; a++)
	for (k = 0; k < ob; k++)
		*i_oovo++ = tmp[(i+oa)*o*v*o+(j+oa)*v*o+(a+va)*o+(k+oa)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < vb; a++)
	for (k = 0; k < oa; k++)
		*i_oovo++ = tmp[(i+oa)*o*v*o+j*v*o+(a+va)*o+k];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
		*i_oovv++ = tmp[i*o*v*v+j*v*v+a*v+b];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
		*i_oovv++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < vb; b++)
		*i_oovv++ = tmp[(i+oa)*o*v*v+(j+oa)*v*v+(a+va)*v+(b+va)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < va; b++)
		*i_oovv++ = tmp[(i+oa)*o*v*v+j*v*v+(a+va)*v+b];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
	for (c = 0; c < b; c++) {
		double t = read_next_double(fp);
		tmp[i*v*v*v+a*v*v+b*v+c] = t;
		tmp[i*v*v*v+a*v*v+c*v+b] = -t;
	}
	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
		tmp[i*v*v*v+a*v*v+b*v+b] = 0;
	for (i = 0; i < oa; i++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
	for (c = 0; c < b; c++)
		*i_ovvv++ = tmp[i*v*v*v+a*v*v+b*v+c];
	for (i = 0; i < oa; i++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < va; b++)
	for (c = 0; c < vb; c++)
		*i_ovvv++ = tmp[i*v*v*v+(a+va)*v*v+b*v+(c+va)];
	for (i = 0; i < ob; i++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < vb; b++)
	for (c = 0; c < b; c++)
		*i_ovvv++ = tmp[(i+oa)*v*v*v+(a+va)*v*v+(b+va)*v+(c+va)];
	for (i = 0; i < ob; i++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
	for (c = 0; c < va; c++)
		*i_ovvv++ = tmp[(i+oa)*v*v*v+a*v*v+(b+va)*v+c];
	free(tmp);
}

static void
read_test_data_rft(FILE *fp, size_t oa, size_t va, double *d_ov, double *f2_ov,
    double *l1, double *t2, double *l2, double *i_oovv, double *i2_t2f2_oovo,
    double *i3_ovvv, double *i6_oovo, double *i7_ovvv)
{
	size_t i, j, k, a, b, c;
	size_t ob = oa, vb = va, o = 2*oa, v = 2*va;
	size_t x = o > v ? o : v;
	double *tmp;

	tmp = xmalloc(x*x*x*x*sizeof(*tmp));

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa; i++) {
		for (j = 0; j < va; j++)
			*d_ov++ = read_next_double(fp);
		for (j = 0; j < vb; j++)
			read_next_double(fp);
	}
	for (i = 0; i < ob; i++) {
		for (j = 0; j < v; j++)
			read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa; i++) {
		for (j = 0; j < va; j++)
			*f2_ov++ = read_next_double(fp);
		for (j = 0; j < vb; j++)
			read_next_double(fp);
	}
	for (i = 0; i < ob; i++) {
		for (j = 0; j < v; j++)
			read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa; i++) {
		for (j = 0; j < va; j++)
			*l1++ = read_next_double(fp);
		for (j = 0; j < vb; j++)
			read_next_double(fp);
	}
	for (i = 0; i < ob; i++) {
		for (j = 0; j < v; j++)
			read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
		*t2++ = tmp[i*o*v*v+j*v*v+a*v+b];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
		*t2++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
		*l2++ = tmp[i*o*v*v+j*v*v+a*v+b];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
		*l2++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
		*i_oovv++ = tmp[i*o*v*v+j*v*v+a*v+b];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
		*i_oovv++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*o; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (k = 0; k < oa; k++)
		*i2_t2f2_oovo++ = tmp[i*o*v*o+j*v*o+a*o+k];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (k = 0; k < ob; k++)
		*i2_t2f2_oovo++ = tmp[i*o*v*o+(j+oa)*v*o+a*o+(k+oa)];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
	for (c = 0; c < b; c++) {
		double t = read_next_double(fp);
		tmp[i*v*v*v+a*v*v+b*v+c] = t;
		tmp[i*v*v*v+a*v*v+c*v+b] = -t;
	}
	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
		tmp[i*v*v*v+a*v*v+b*v+b] = 0;
	for (i = 0; i < oa; i++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
	for (c = 0; c < b; c++)
		*i3_ovvv++ = tmp[i*v*v*v+a*v*v+b*v+c];
	for (i = 0; i < oa; i++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < va; b++)
	for (c = 0; c < vb; c++)
		*i3_ovvv++ = tmp[i*v*v*v+(a+va)*v*v+b*v+(c+va)];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*o; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (k = 0; k < oa; k++)
		*i6_oovo++ = tmp[i*o*v*o+j*v*o+a*o+k];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (k = 0; k < ob; k++)
		*i6_oovo++ = tmp[i*o*v*o+(j+oa)*v*o+a*o+(k+oa)];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
	for (c = 0; c < b; c++) {
		double t = read_next_double(fp);
		tmp[i*v*v*v+a*v*v+b*v+c] = t;
		tmp[i*v*v*v+a*v*v+c*v+b] = -t;
	}
	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
		tmp[i*v*v*v+a*v*v+b*v+b] = 0;
	for (i = 0; i < oa; i++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
	for (c = 0; c < b; c++)
		*i7_ovvv++ = tmp[i*v*v*v+a*v*v+b*v+c];
	for (i = 0; i < oa; i++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < va; b++)
	for (c = 0; c < vb; c++)
		*i7_ovvv++ = tmp[i*v*v*v+(a+va)*v*v+b*v+(c+va)];
	free(tmp);
}

static void
read_test_data_uft(FILE *fp, size_t oa, size_t va, size_t ob, size_t vb,
    double *d_ov, double *f2_ov, double *l1, double *t2, double *l2,
    double *i_oovv, double *i2_t2f2_oovo, double *i3_ovvv, double *i6_oovo,
    double *i7_ovvv)
{
	size_t i, j, k, a, b, c;
	size_t o = oa+ob, v = va+vb;
	size_t x = o > v ? o : v;
	double *tmp;

	tmp = xmalloc(x*x*x*x*sizeof(*tmp));

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa; i++) {
		for (j = 0; j < va; j++)
			*d_ov++ = read_next_double(fp);
		for (j = 0; j < vb; j++)
			read_next_double(fp);
	}
	for (i = 0; i < ob; i++) {
		for (j = 0; j < va; j++)
			read_next_double(fp);
		for (j = 0; j < vb; j++)
			*d_ov++ = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa; i++) {
		for (j = 0; j < va; j++)
			*f2_ov++ = read_next_double(fp);
		for (j = 0; j < vb; j++)
			read_next_double(fp);
	}
	for (i = 0; i < ob; i++) {
		for (j = 0; j < va; j++)
			read_next_double(fp);
		for (j = 0; j < vb; j++)
			*f2_ov++ = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa; i++) {
		for (j = 0; j < va; j++)
			*l1++ = read_next_double(fp);
		for (j = 0; j < vb; j++)
			read_next_double(fp);
	}
	for (i = 0; i < ob; i++) {
		for (j = 0; j < va; j++)
			read_next_double(fp);
		for (j = 0; j < vb; j++)
			*l1++ = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
		*t2++ = tmp[i*o*v*v+j*v*v+a*v+b];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
		*t2++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < vb; b++)
		*t2++ = tmp[(i+oa)*o*v*v+(j+oa)*v*v+(a+va)*v+(b+va)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < va; b++)
		*t2++ = tmp[(i+oa)*o*v*v+j*v*v+(a+va)*v+b];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
		*l2++ = tmp[i*o*v*v+j*v*v+a*v+b];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
		*l2++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < vb; b++)
		*l2++ = tmp[(i+oa)*o*v*v+(j+oa)*v*v+(a+va)*v+(b+va)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < va; b++)
		*l2++ = tmp[(i+oa)*o*v*v+j*v*v+(a+va)*v+b];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
		*i_oovv++ = tmp[i*o*v*v+j*v*v+a*v+b];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
		*i_oovv++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < vb; b++)
		*i_oovv++ = tmp[(i+oa)*o*v*v+(j+oa)*v*v+(a+va)*v+(b+va)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < va; b++)
		*i_oovv++ = tmp[(i+oa)*o*v*v+j*v*v+(a+va)*v+b];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*o; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (k = 0; k < oa; k++)
		*i2_t2f2_oovo++ = tmp[i*o*v*o+j*v*o+a*o+k];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (k = 0; k < ob; k++)
		*i2_t2f2_oovo++ = tmp[i*o*v*o+(j+oa)*v*o+a*o+(k+oa)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < vb; a++)
	for (k = 0; k < ob; k++)
		*i2_t2f2_oovo++ = tmp[(i+oa)*o*v*o+(j+oa)*v*o+(a+va)*o+(k+oa)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < vb; a++)
	for (k = 0; k < oa; k++)
		*i2_t2f2_oovo++ = tmp[(i+oa)*o*v*o+j*v*o+(a+va)*o+k];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
	for (c = 0; c < b; c++) {
		double t = read_next_double(fp);
		tmp[i*v*v*v+a*v*v+b*v+c] = t;
		tmp[i*v*v*v+a*v*v+c*v+b] = -t;
	}
	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
		tmp[i*v*v*v+a*v*v+b*v+b] = 0;
	for (i = 0; i < oa; i++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
	for (c = 0; c < b; c++)
		*i3_ovvv++ = tmp[i*v*v*v+a*v*v+b*v+c];
	for (i = 0; i < oa; i++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < va; b++)
	for (c = 0; c < vb; c++)
		*i3_ovvv++ = tmp[i*v*v*v+(a+va)*v*v+b*v+(c+va)];
	for (i = 0; i < ob; i++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < vb; b++)
	for (c = 0; c < b; c++)
		*i3_ovvv++ = tmp[(i+oa)*v*v*v+(a+va)*v*v+(b+va)*v+(c+va)];
	for (i = 0; i < ob; i++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
	for (c = 0; c < va; c++)
		*i3_ovvv++ = tmp[(i+oa)*v*v*v+a*v*v+(b+va)*v+c];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*o; i++)
		tmp[i] = read_next_double(fp);
	for (i = 0; i < oa; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < va; a++)
	for (k = 0; k < oa; k++)
		*i6_oovo++ = tmp[i*o*v*o+j*v*o+a*o+k];
	for (i = 0; i < oa; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < va; a++)
	for (k = 0; k < ob; k++)
		*i6_oovo++ = tmp[i*o*v*o+(j+oa)*v*o+a*o+(k+oa)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < ob; j++)
	for (a = 0; a < vb; a++)
	for (k = 0; k < ob; k++)
		*i6_oovo++ = tmp[(i+oa)*o*v*o+(j+oa)*v*o+(a+va)*o+(k+oa)];
	for (i = 0; i < ob; i++)
	for (j = 0; j < oa; j++)
	for (a = 0; a < vb; a++)
	for (k = 0; k < oa; k++)
		*i6_oovo++ = tmp[(i+oa)*o*v*o+j*v*o+(a+va)*o+k];
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
	for (c = 0; c < b; c++) {
		double t = read_next_double(fp);
		tmp[i*v*v*v+a*v*v+b*v+c] = t;
		tmp[i*v*v*v+a*v*v+c*v+b] = -t;
	}
	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
		tmp[i*v*v*v+a*v*v+b*v+b] = 0;
	for (i = 0; i < oa; i++)
	for (a = 0; a < va; a++)
	for (b = 0; b < va; b++)
	for (c = 0; c < b; c++)
		*i7_ovvv++ = tmp[i*v*v*v+a*v*v+b*v+c];
	for (i = 0; i < oa; i++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < va; b++)
	for (c = 0; c < vb; c++)
		*i7_ovvv++ = tmp[i*v*v*v+(a+va)*v*v+b*v+(c+va)];
	for (i = 0; i < ob; i++)
	for (a = 0; a < vb; a++)
	for (b = 0; b < vb; b++)
	for (c = 0; c < b; c++)
		*i7_ovvv++ = tmp[(i+oa)*v*v*v+(a+va)*v*v+(b+va)*v+(c+va)];
	for (i = 0; i < ob; i++)
	for (a = 0; a < va; a++)
	for (b = 0; b < vb; b++)
	for (c = 0; c < va; c++)
		*i7_ovvv++ = tmp[(i+oa)*v*v*v+a*v*v+(b+va)*v+c];
	free(tmp);
}

static void
test_rpt(FILE *fp, double *e_ref, double *e_cmp)
{
	double *d_ov, *f_ov, *t1, *t2;
	double *i_oovo, *i_oovv, *i_ovvv;
	size_t oa, va, ob, vb, d_ov_sz, f_ov_sz, t1_sz, t2_sz;
	size_t i_oovo_sz, i_oovv_sz, i_ovvv_sz;

	read_test_header(fp, 0, &oa, &va, &ob, &vb, e_ref);

	d_ov_sz = oa*va;
	f_ov_sz = oa*va;
	t1_sz = oa*va;
	t2_sz = oa*oa*va*va + oa*ob*va*vb;
	i_oovo_sz = oa*oa*va*oa + oa*ob*va*ob;
	i_oovv_sz = oa*oa*va*va + oa*ob*va*vb;
	i_ovvv_sz = oa*va*va*(va-1)/2 + oa*vb*va*vb;

	d_ov = xmalloc(d_ov_sz * sizeof(*d_ov));
	f_ov = xmalloc(f_ov_sz * sizeof(*f_ov));
	t1 = xmalloc(t1_sz * sizeof(*t1));
	t2 = xmalloc(t2_sz * sizeof(*t2));
	i_oovo = xmalloc(i_oovo_sz * sizeof(*i_oovo));
	i_oovv = xmalloc(i_oovv_sz * sizeof(*i_oovv));
	i_ovvv = xmalloc(i_ovvv_sz * sizeof(*i_ovvv));

	read_test_data_rpt(fp, oa, va, d_ov, f_ov, t1, t2,
	    i_oovo, i_oovv, i_ovvv);
	*e_cmp = libpt_rpt(oa, va, d_ov, f_ov, t1, t2,
	    i_oovo, i_oovv, i_ovvv);

	free(d_ov);
	free(f_ov);
	free(t1);
	free(t2);
	free(i_oovo);
	free(i_oovv);
	free(i_ovvv);
}

static void
test_upt(FILE *fp, double *e_ref, double *e_cmp)
{
	double *d_ov, *f_ov, *t1, *t2;
	double *i_oovo, *i_oovv, *i_ovvv;
	size_t oa, va, ob, vb, d_ov_sz, f_ov_sz, t1_sz, t2_sz;
	size_t i_oovo_sz, i_oovv_sz, i_ovvv_sz;

	read_test_header(fp, 1, &oa, &va, &ob, &vb, e_ref);

	d_ov_sz = oa*va + ob*vb;
	f_ov_sz = oa*va + ob*vb;
	t1_sz = oa*va + ob*vb;
	t2_sz = oa*oa*va*va + 2*oa*ob*va*vb + ob*ob*vb*vb;
	i_oovo_sz = oa*oa*va*oa + oa*ob*va*ob + ob*oa*vb*oa +
	    ob*ob*vb*ob;
	i_oovv_sz = oa*oa*va*va + 2*oa*ob*va*vb + ob*ob*vb*vb;
	i_ovvv_sz = oa*va*va*(va-1)/2 + oa*vb*va*vb + ob*va*vb*va +
	    ob*vb*vb*(vb-1)/2;

	d_ov = xmalloc(d_ov_sz * sizeof(*d_ov));
	f_ov = xmalloc(f_ov_sz * sizeof(*f_ov));
	t1 = xmalloc(t1_sz * sizeof(*t1));
	t2 = xmalloc(t2_sz * sizeof(*t2));
	i_oovo = xmalloc(i_oovo_sz * sizeof(*i_oovo));
	i_oovv = xmalloc(i_oovv_sz * sizeof(*i_oovv));
	i_ovvv = xmalloc(i_ovvv_sz * sizeof(*i_ovvv));

	read_test_data_upt(fp, oa, va, ob, vb, d_ov, f_ov, t1, t2,
	    i_oovo, i_oovv, i_ovvv);
	*e_cmp = libpt_upt(oa, va, ob, vb, d_ov, f_ov, t1, t2,
	    i_oovo, i_oovv, i_ovvv);

	free(d_ov);
	free(f_ov);
	free(t1);
	free(t2);
	free(i_oovo);
	free(i_oovv);
	free(i_ovvv);
}

static void
test_rft(FILE *fp, double *e_ref, double *e_cmp)
{
	double *d_ov, *f2_ov, *l1, *t2, *l2;
	double *i_oovv, *i2_t2f2_oovo, *i3_ovvv, *i6_oovo, *i7_ovvv;
	size_t oa, va, ob, vb, d_ov_sz, f2_ov_sz, l1_sz, t2_sz, l2_sz;
	size_t i_oovv_sz, i2_t2f2_oovo_sz, i3_ovvv_sz, i6_oovo_sz, i7_ovvv_sz;

	read_test_header(fp, 0, &oa, &va, &ob, &vb, e_ref);

	d_ov_sz = oa*va;
	f2_ov_sz = oa*va;
	l1_sz = oa*va;
	t2_sz = oa*oa*va*va + oa*ob*va*vb;
	l2_sz = oa*oa*va*va + oa*ob*va*vb;
	i_oovv_sz = oa*oa*va*va + oa*ob*va*vb;
	i2_t2f2_oovo_sz = oa*oa*va*oa + oa*ob*va*ob;
	i3_ovvv_sz = oa*va*va*(va-1)/2 + oa*vb*va*vb;
	i6_oovo_sz = oa*oa*va*oa + oa*ob*va*ob;
	i7_ovvv_sz = oa*va*va*(va-1)/2 + oa*vb*va*vb;

	d_ov = xmalloc(d_ov_sz * sizeof(*d_ov));
	f2_ov = xmalloc(f2_ov_sz * sizeof(*f2_ov));
	l1 = xmalloc(l1_sz * sizeof(*l1));
	t2 = xmalloc(t2_sz * sizeof(*t2));
	l2 = xmalloc(l2_sz * sizeof(*l2));
	i_oovv = xmalloc(i_oovv_sz * sizeof(*i_oovv));
	i2_t2f2_oovo = xmalloc(i2_t2f2_oovo_sz * sizeof(*i2_t2f2_oovo));
	i3_ovvv = xmalloc(i3_ovvv_sz * sizeof(*i3_ovvv));
	i6_oovo = xmalloc(i6_oovo_sz * sizeof(*i6_oovo));
	i7_ovvv = xmalloc(i7_ovvv_sz * sizeof(*i7_ovvv));

	read_test_data_rft(fp, oa, va, d_ov, f2_ov, l1, t2, l2,
	    i_oovv, i2_t2f2_oovo, i3_ovvv, i6_oovo, i7_ovvv);
	*e_cmp = libpt_rft(oa, va, d_ov, f2_ov, l1, t2, l2,
	    i_oovv, i2_t2f2_oovo, i3_ovvv, i6_oovo, i7_ovvv);

	free(d_ov);
	free(f2_ov);
	free(l1);
	free(t2);
	free(l2);
	free(i_oovv);
	free(i2_t2f2_oovo);
	free(i3_ovvv);
	free(i6_oovo);
	free(i7_ovvv);
}

static void
test_uft(FILE *fp, double *e_ref, double *e_cmp)
{
	double *d_ov, *f2_ov, *l1, *t2, *l2;
	double *i_oovv, *i2_t2f2_oovo, *i3_ovvv, *i6_oovo, *i7_ovvv;
	size_t oa, va, ob, vb, d_ov_sz, f2_ov_sz, l1_sz, t2_sz, l2_sz;
	size_t i_oovv_sz, i2_t2f2_oovo_sz, i3_ovvv_sz, i6_oovo_sz, i7_ovvv_sz;

	read_test_header(fp, 1, &oa, &va, &ob, &vb, e_ref);

	d_ov_sz = oa*va + ob*vb;
	f2_ov_sz = oa*va + ob*vb;
	l1_sz = oa*va + ob*vb;
	t2_sz = oa*oa*va*va + 2*oa*ob*va*vb + ob*ob*vb*vb;
	l2_sz = oa*oa*va*va + 2*oa*ob*va*vb + ob*ob*vb*vb;
	i_oovv_sz = oa*oa*va*va + 2*oa*ob*va*vb + ob*ob*vb*vb;
	i2_t2f2_oovo_sz = oa*oa*va*oa + oa*ob*va*ob + ob*oa*vb*oa +
	    ob*ob*vb*ob;
	i3_ovvv_sz = oa*va*va*(va-1)/2 + oa*vb*va*vb + ob*va*vb*va +
	    ob*vb*vb*(vb-1)/2;
	i6_oovo_sz = oa*oa*va*oa + oa*ob*va*ob + ob*oa*vb*oa +
	    ob*ob*vb*ob;
	i7_ovvv_sz = oa*va*va*(va-1)/2 + oa*vb*va*vb + ob*va*vb*va +
	    ob*vb*vb*(vb-1)/2;

	d_ov = xmalloc(d_ov_sz * sizeof(*d_ov));
	f2_ov = xmalloc(f2_ov_sz * sizeof(*f2_ov));
	l1 = xmalloc(l1_sz * sizeof(*l1));
	t2 = xmalloc(t2_sz * sizeof(*t2));
	l2 = xmalloc(l2_sz * sizeof(*l2));
	i_oovv = xmalloc(i_oovv_sz * sizeof(*i_oovv));
	i2_t2f2_oovo = xmalloc(i2_t2f2_oovo_sz * sizeof(*i2_t2f2_oovo));
	i3_ovvv = xmalloc(i3_ovvv_sz * sizeof(*i3_ovvv));
	i6_oovo = xmalloc(i6_oovo_sz * sizeof(*i6_oovo));
	i7_ovvv = xmalloc(i7_ovvv_sz * sizeof(*i7_ovvv));

	read_test_data_uft(fp, oa, va, ob, vb, d_ov, f2_ov, l1, t2, l2,
	    i_oovv, i2_t2f2_oovo, i3_ovvv, i6_oovo, i7_ovvv);
	*e_cmp = libpt_uft(oa, va, ob, vb, d_ov, f2_ov, l1, t2, l2,
	    i_oovv, i2_t2f2_oovo, i3_ovvv, i6_oovo, i7_ovvv);

	free(d_ov);
	free(f2_ov);
	free(l1);
	free(t2);
	free(l2);
	free(i_oovv);
	free(i2_t2f2_oovo);
	free(i3_ovvv);
	free(i6_oovo);
	free(i7_ovvv);
}

static const struct {
	const char *type;
	void (*fn)(FILE *, double *, double *);
} tests[] = {
	{ "rpt", test_rpt },
	{ "upt", test_upt },
	{ "rft", test_rft },
	{ "uft", test_uft },
};
static const size_t ntests = sizeof tests / sizeof *tests;

int
main(int argc, char **argv)
{
	FILE *fp;
	double e_ref = 0, e_cmp = 0;
	size_t i;
	int rank = 0;
	char path[256];

#ifdef LIBPT_USE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (argc < 2)
		errx(1, "usage: test rpt|upt|rft|uft id");
	if (rank == 0)
		printf("%s%s\n", argv[1], argv[2]);
	snprintf(path, sizeof path, "tests/%s%s.dat", argv[1], argv[2]);
	if ((fp = fopen(path, "r")) == NULL)
		err(1, "fopen");
	for (i = 0; i < ntests; i++)
		if (strcmp(tests[i].type, argv[1]) == 0) {
			tests[i].fn(fp, &e_ref, &e_cmp);
			break;
		}
	if (i == ntests)
		errx(1, "bad test type");
	if (rank == 0) {
		printf("%s energy: % .8lf\n", argv[1], e_cmp);
		printf("%s ref:    % .8lf\n", argv[1], e_ref);
	}
	fclose(fp);
#ifdef LIBPT_USE_MPI
	MPI_Finalize();
#endif
	return (fabs(e_cmp - e_ref) < EPSILON ? 0 : 1);
}
