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

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <err.h>
#include <getopt.h>
#include <unistd.h>

#include <mpi.h>

#include "pt.h"

#define EPSILON 1.0e-8

long long strtonum(const char *, long long, long long, const char **);

static void
usage(void)
{
	fprintf(stderr, "usage: pt [-o nocc] [-v nvirt] [-t test]\n");
	exit(1);
}

static void *
xmalloc(size_t sz)
{
	void *p;

	if ((p = malloc(sz)) == NULL)
		err(1, "malloc");
	return (p);
}

static void
load_test_header(const char *testpath, size_t *o, size_t *v, double *e_ref)
{
	FILE *fp;

	if ((fp = fopen(testpath, "r")) == NULL)
		err(1, "unable to open %s", testpath);
	if (fscanf(fp, "%zu %zu %lf", o, v, e_ref) != 3)
		errx(1, "error parsing test file header");
	fclose(fp);
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

	if ((fscanf(fp, "%lf", &el)) != 1)
		errx(1, "error parsing test file data");
	return (el);
}

static void
load_test_data(const char *testpath, size_t o, size_t v, double *d_ov,
    double *f_ov, double *i_ooov, double *i_oovv, double *i_ovvv,
    double *t1, double *t2)
{
	FILE *fp;
	size_t i;

	if ((fp = fopen(testpath, "r")) == NULL)
		err(1, "unable to open %s", testpath);

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v; i++) {
		d_ov[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v; i++) {
		f_ov[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*o*v; i++) {
		i_ooov[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		i_oovv[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v*v*v; i++) {
		i_ovvv[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v; i++) {
		t1[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		t2[i] = read_next_double(fp);
	}
	fclose(fp);
}

static double
random_double(void)
{
	return (drand48() / 100.0);
}

static void
load_random_data(size_t o, size_t v, double *d_ov,
    double *f_ov, double *i_ooov, double *i_oovv, double *i_ovvv,
    double *t1, double *t2)
{
	size_t i;

	for (i = 0; i < o*v; i++) {
		d_ov[i] = random_double();
	}
	for (i = 0; i < o*v; i++) {
		f_ov[i] = random_double();
	}
	for (i = 0; i < o*o*o*v; i++) {
		i_ooov[i] = random_double();
	}
	for (i = 0; i < o*o*v*v; i++) {
		i_oovv[i] = random_double();
	}
	for (i = 0; i < o*v*v*v; i++) {
		i_ovvv[i] = random_double();
	}
	for (i = 0; i < o*v; i++) {
		t1[i] = random_double();
	}
	for (i = 0; i < o*o*v*v; i++) {
		t2[i] = random_double();
	}
}

int
main(int argc, char **argv)
{
	size_t o = 4, v = 20;
	double *d_ov, *f_ov;
	double *i_ooov, *i_oovv, *i_ovvv;
	double *t1, *t2;
	double e_pt, e_ref = 0.0;
	int rank;
	const char *errstr, *testpath = NULL;
	char ch;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	while ((ch = getopt(argc, argv, "o:t:v:")) != -1) {
		switch (ch) {
		case 'o':
			o = strtonum(optarg, 1, INT_MAX, &errstr);
			if (errstr)
				errx(1, "bad o value: %s", errstr);
			break;
		case 't':
			testpath = optarg;
			break;
		case 'v':
			v = strtonum(optarg, 1, INT_MAX, &errstr);
			if (errstr)
				errx(1, "bad v value: %s", errstr);
			break;
		default:
			usage();
			break;
		}
	}
	argv += optind;
	argc -= optind;

	if (testpath)
		load_test_header(testpath, &o, &v, &e_ref);
	d_ov = xmalloc(o*v * sizeof(double));
	f_ov = xmalloc(o*v * sizeof(double));
	i_ooov = xmalloc(o*o*o*v * sizeof(double));
	i_oovv = xmalloc(o*o*v*v * sizeof(double));
	i_ovvv = xmalloc(o*v*v*v * sizeof(double));
	t1 = xmalloc(o*v * sizeof(double));
	t2 = xmalloc(o*o*v*v * sizeof(double));

	if (rank == 0) {
		if (testpath) {
			load_test_data(testpath, o, v, d_ov, f_ov, i_ooov,
			    i_oovv, i_ovvv, t1, t2);
		} else {
			load_random_data(o, v, d_ov, f_ov, i_ooov,
			    i_oovv, i_ovvv, t1, t2);
		}
	}

	MPI_Bcast(d_ov, o*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(f_ov, o*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(i_ooov, o*o*o*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(i_oovv, o*o*v*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(i_ovvv, o*v*v*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(t1, o*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(t2, o*o*v*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	e_pt = ccsd_pt(o, v, d_ov, f_ov, i_ooov, i_oovv, i_ovvv, t1, t2);
	if (rank == 0)
		printf("ccsd(t) energy: % .8lf\n", e_pt);
	if (testpath) {
		if (rank == 0)
			printf("ccsd(t) ref:    % .8lf\n", e_ref);
	} else
		e_ref = e_pt;

	free(d_ov);
	free(f_ov);
	free(i_ooov);
	free(i_oovv);
	free(i_ovvv);
	free(t1);
	free(t2);
	MPI_Finalize();
	return (fabs(e_pt - e_ref) < EPSILON ? 0 : 1);
}
