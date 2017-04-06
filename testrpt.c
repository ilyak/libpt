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

#include <err.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "pt.h"

#define EPSILON 1.0e-8

static void
read_test_header(FILE *fp, size_t *oa, size_t *va, double *e_ref)
{
	if (fscanf(fp, "%zu %zu %lf", oa, va, e_ref) != 3)
		errx(1, "unable to read test header");
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
read_test_data(FILE *fp, size_t oa, size_t va, double *d_ov, double *f_ov,
    double *t1, double *t2, double *i_oovo, double *i_oovv, double *i_ovvv)
{
	size_t i;

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*va; i++) {
		d_ov[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*va; i++) {
		f_ov[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*va; i++) {
		t1[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*va*va; i++) {
		t2[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*va*va; i++) {
		t2[oa*oa*va*va+i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*oa*va; i++) {
		i_oovo[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*oa*va; i++) {
		i_oovo[oa*oa*oa*va+i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*va*va; i++) {
		i_oovv[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*oa*va*va; i++) {
		i_oovv[oa*oa*va*va+i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*va*va*(va-1)/2; i++) {
		i_ovvv[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < oa*va*va*va; i++) {
		i_ovvv[oa*va*va*(va-1)/2+i] = read_next_double(fp);
	}
}

int
main(int argc, char **argv)
{
	FILE *fp;
	double e_cmp, e_ref, *d_ov, *f_ov, *t1, *t2;
	double *i_oovo, *i_oovv, *i_ovvv;
	size_t oa, va, size;
	int rank = 0;

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (argc < 2)
		errx(1, "specify test file");
	if ((fp = fopen(argv[1], "r")) == NULL)
		err(1, "fopen");
	read_test_header(fp, &oa, &va, &e_ref);

	size = 3*oa*va + 2*oa*oa*oa*va + 4*oa*oa*va*va + oa*va*va*va +
	    oa*va*va*(va-1)/2;
	if ((d_ov = malloc(size * sizeof(double))) == NULL)
		err(1, "malloc");
	f_ov = d_ov + oa*va;
	t1 = f_ov + oa*va;
	t2 = t1 + oa*va;
	i_oovo = t2 + 2*oa*oa*va*va;
	i_oovv = i_oovo + 2*oa*oa*va*oa;
	i_ovvv = i_oovv + 2*oa*oa*va*va;

	read_test_data(fp, oa, va, d_ov, f_ov, t1, t2, i_oovo, i_oovv, i_ovvv);
	fclose(fp);

	e_cmp = libpt_rpt(oa, va, d_ov, f_ov, t1, t2, i_oovo, i_oovv, i_ovvv);
	if (rank == 0) {
		printf("rpt energy: % .8lf\n", e_cmp);
		printf("rpt ref:    % .8lf\n", e_ref);
	}
	free(d_ov);
#ifdef WITH_MPI
	MPI_Finalize();
#endif
	return (fabs(e_cmp - e_ref) < EPSILON ? 0 : 1);
}
