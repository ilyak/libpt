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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <err.h>

#ifdef LIBPT_USE_MPI
#include <mpi.h>
#endif

#include "pt.h"

#define EPSILON 1.0e-8

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

	if ((tmp = malloc(x*x*x*x*sizeof(*tmp))) == NULL)
		err(1, "malloc");

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

int
main(int argc, char **argv)
{
	FILE *fp;
	double e_cmp, e_ref, *d_ov, *f_ov, *t1, *t2;
	double *i_oovo, *i_oovv, *i_ovvv;
	size_t oa, va, ob, vb, d_ov_sz, f_ov_sz, t1_sz, t2_sz;
	size_t i_oovo_sz, i_oovv_sz, i_ovvv_sz;
	int rank = 0, unrestricted;

#ifdef LIBPT_USE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (argc < 2)
		errx(1, "specify test file");
	if ((fp = fopen(argv[1], "r")) == NULL)
		err(1, "fopen");
	unrestricted = strstr(argv[0], "upt") != NULL;
	read_test_header(fp, unrestricted, &oa, &va, &ob, &vb, &e_ref);

	if (unrestricted) {
		d_ov_sz = oa*va + ob*vb;
		f_ov_sz = oa*va + ob*vb;
		t1_sz = oa*va + ob*vb;
		t2_sz = oa*oa*va*va + 2*oa*ob*va*vb + ob*ob*vb*vb;
		i_oovo_sz = oa*oa*va*oa + oa*ob*va*ob + ob*oa*vb*oa +
		    ob*ob*vb*ob;
		i_oovv_sz = oa*oa*va*va + 2*oa*ob*va*vb + ob*ob*vb*vb;
		i_ovvv_sz = oa*va*va*(va-1)/2 + oa*vb*va*vb + ob*va*vb*va +
		    ob*vb*vb*(vb-1)/2;
	} else {
		d_ov_sz = oa*va;
		f_ov_sz = oa*va;
		t1_sz = oa*va;
		t2_sz = oa*oa*va*va + oa*ob*va*vb;
		i_oovo_sz = oa*oa*va*oa + oa*ob*va*ob;
		i_oovv_sz = oa*oa*va*va + oa*ob*va*vb;
		i_ovvv_sz = oa*va*va*(va-1)/2 + oa*vb*va*vb;
	}
	if ((d_ov = malloc(d_ov_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((f_ov = malloc(f_ov_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((t1 = malloc(t1_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((t2 = malloc(t2_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((i_oovo = malloc(i_oovo_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((i_oovv = malloc(i_oovv_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((i_ovvv = malloc(i_ovvv_sz * sizeof(double))) == NULL)
		err(1, "malloc");

	if (unrestricted) {
		read_test_data_upt(fp, oa, va, ob, vb, d_ov, f_ov, t1, t2,
		    i_oovo, i_oovv, i_ovvv);
	} else {
		read_test_data_rpt(fp, oa, va, d_ov, f_ov, t1, t2,
		    i_oovo, i_oovv, i_ovvv);
	}
	fclose(fp);

	if (unrestricted) {
		e_cmp = libpt_upt(oa, va, ob, vb, d_ov, f_ov, t1, t2,
		    i_oovo, i_oovv, i_ovvv);
	} else {
		e_cmp = libpt_rpt(oa, va, d_ov, f_ov, t1, t2,
		    i_oovo, i_oovv, i_ovvv);
	}
	if (rank == 0) {
		const char *s = unrestricted ? "upt" : "rpt";
		printf("%s energy: % .8lf\n", s, e_cmp);
		printf("%s ref:    % .8lf\n", s, e_ref);
	}
	free(d_ov);
	free(f_ov);
	free(t1);
	free(t2);
	free(i_oovo);
	free(i_oovv);
	free(i_ovvv);
#ifdef LIBPT_USE_MPI
	MPI_Finalize();
#endif
	return (fabs(e_cmp - e_ref) < EPSILON ? 0 : 1);
}
