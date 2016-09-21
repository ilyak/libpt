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
#include <time.h>

#include <err.h>
#include <getopt.h>
#include <unistd.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif
#include "pt.h"

#define EPSILON 1.0e-8

double drand48(void);
long long strtonum(const char *, long long, long long, const char **);

static void
usage(void)
{
	fprintf(stderr, "usage: pt [-o nocc] [-v nvirt] [-u] [-t test]\n");
	exit(1);
}

static void *
xmalloc(size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
		err(1, "malloc");
	return (p);
}

static void
load_test_header(const char *testpath, size_t *o, size_t *v, int *is_rpt,
    double *e_ref)
{
	FILE *fp;
	char buf[16];

	if ((fp = fopen(testpath, "r")) == NULL)
		err(1, "unable to open %s", testpath);
	if (fscanf(fp, "%12s", buf) != 1)
		errx(1, "error parsing test file header");
	*is_rpt = strcmp(buf, "unrestricted") != 0;
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
load_test_data(const char *testpath, size_t o, size_t v, int is_rpt,
    double *d_ov, double *f_ov, double *t1, double *t2, double *i_oovo,
    double *i_oovv, double *i_ovvv)
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
	for (i = 0; i < o*v; i++) {
		t1[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		t2[i] = read_next_double(fp);
	}
	if (is_rpt) {
		skip_line(fp);
		skip_line(fp);
		for (i = 0; i < o*o*v*v; i++) {
			t2[o*o*v*v+i] = read_next_double(fp);
		}
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*o*v; i++) {
		i_oovo[i] = read_next_double(fp);
	}
	if (is_rpt) {
		skip_line(fp);
		skip_line(fp);
		for (i = 0; i < o*o*o*v; i++) {
			i_oovo[o*o*o*v+i] = read_next_double(fp);
		}
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		i_oovv[i] = read_next_double(fp);
	}
	if (is_rpt) {
		skip_line(fp);
		skip_line(fp);
		for (i = 0; i < o*o*v*v; i++) {
			i_oovv[o*o*v*v+i] = read_next_double(fp);
		}
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v*v*v; i++) {
		i_ovvv[i] = read_next_double(fp);
	}
	if (is_rpt) {
		skip_line(fp);
		skip_line(fp);
		for (i = 0; i < o*v*v*v; i++) {
			i_ovvv[o*v*v*v+i] = read_next_double(fp);
		}
	}
	fclose(fp);
}

static double
random_double(void)
{
	return (drand48() / 100000.0);
}

static void
load_random_data(size_t o, size_t v, int is_rpt, double *d_ov,
    double *f_ov, double *t1, double *t2, double *i_oovo,
    double *i_oovv, double *i_ovvv)
{
	size_t i, nsp = is_rpt ? 2 : 1;

	for (i = 0; i < o*v; i++) {
		d_ov[i] = random_double();
	}
	for (i = 0; i < o*v; i++) {
		f_ov[i] = random_double();
	}
	for (i = 0; i < o*v; i++) {
		t1[i] = random_double();
	}
	for (i = 0; i < nsp*o*o*v*v; i++) {
		t2[i] = random_double();
	}
	for (i = 0; i < nsp*o*o*o*v; i++) {
		i_oovo[i] = random_double();
	}
	for (i = 0; i < nsp*o*o*v*v; i++) {
		i_oovv[i] = random_double();
	}
	for (i = 0; i < o*v*v*(v-1)/2; i++) {
		i_ovvv[i] = random_double();
	}
	if (is_rpt) {
		for (i = 0; i < o*v*v*v; i++) {
			i_ovvv[o*v*v*(v-1)/2+i] = random_double();
		}
	}
}

int
main(int argc, char **argv)
{
	size_t o = 0, v = 0, nsp;
	double e_pt = 0.0, e_ref = 0.0;
	double *d_ov, *f_ov, *t1, *t2, *i_oovv, *i_oovo, *i_ovvv;
	const char *errstr, *testpath = NULL;
	int is_rpt = 1, rank = 0;
	char ch;

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	while ((ch = getopt(argc, argv, "o:t:uv:")) != -1) {
		switch (ch) {
		case 'o':
			o = strtonum(optarg, 1, INT_MAX, &errstr);
			if (errstr)
				errx(1, "bad o value: %s", errstr);
			break;
		case 't':
			testpath = optarg;
			break;
		case 'u':
			is_rpt = 0;
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
		load_test_header(testpath, &o, &v, &is_rpt, &e_ref);
	nsp = is_rpt ? 2 : 1;

	d_ov = xmalloc(o*v*sizeof(double));
	f_ov = xmalloc(o*v*sizeof(double));
	t1 = xmalloc(o*v*sizeof(double));
	t2 = xmalloc(nsp*o*o*v*v*sizeof(double));
	i_oovo = xmalloc(nsp*o*o*o*v*sizeof(double));
	i_oovv = xmalloc(nsp*o*o*v*v*sizeof(double));
	i_ovvv = xmalloc(o*v*(v*(v-1)/2+(nsp-1)*v*v)*sizeof(double));

	if (testpath) {
		double *i_ovvv2 = xmalloc(nsp*o*v*v*v*sizeof(double));
		load_test_data(testpath, o, v, is_rpt, d_ov, f_ov,
		    t1, t2, i_oovo, i_oovv, i_ovvv2);
		for (size_t i = 0; i < o; i++) {
		for (size_t a = 0; a < v; a++) {
		for (size_t b = 0; b < v; b++) {
		for (size_t c = 0; c < b; c++) {
			i_ovvv[i*v*v*(v-1)/2+a*v*(v-1)/2+b*(b-1)/2+c] =
			    i_ovvv2[i*v*v*v+a*v*v+b*v+c];
		}}}}
		if (is_rpt)
			memcpy(i_ovvv+o*v*v*(v-1)/2, i_ovvv2+o*v*v*v,
			    o*v*v*v*sizeof(double));
		free(i_ovvv2);
	} else {
		load_random_data(o, v, is_rpt, d_ov, f_ov,
		    t1, t2, i_oovo, i_oovv, i_ovvv);
	}

	if (rank == 0) {
		time_t t = time(NULL);
		printf("ccsd_pt: %s", ctime(&t));
	}
	if (is_rpt) {
		e_pt = ccsd_rpt(o, v, d_ov, f_ov, t1, t2,
		    i_oovo, i_oovv, i_ovvv);
	} else {
		e_pt = ccsd_upt(o, v, d_ov, f_ov, t1, t2,
		    i_oovo, i_oovv, i_ovvv);
	}
	if (rank == 0) {
		time_t t = time(NULL);
		printf("ccsd_pt: %s", ctime(&t));
		printf("ccsd(t) energy: % .8lf\n", e_pt);
	}
	if (testpath) {
		if (rank == 0)
			printf("ccsd(t) ref:    % .8lf\n", e_ref);
	} else
		e_ref = e_pt;

	free(d_ov);
	free(f_ov);
	free(t1);
	free(t2);
	free(i_oovo);
	free(i_oovv);
	free(i_ovvv);
#ifdef WITH_MPI
	MPI_Finalize();
#endif
	return (fabs(e_pt - e_ref) < EPSILON ? 0 : 1);
}
