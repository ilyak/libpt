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
load_test_header(const char *testpath, size_t *oa, size_t *ob,
    size_t *va, size_t *vb, int *is_rpt, double *e_ref)
{
	FILE *fp;
	char buf[16];

	if ((fp = fopen(testpath, "r")) == NULL)
		err(1, "unable to open %s", testpath);
	if (fscanf(fp, "%12s", buf) != 1)
		errx(1, "error parsing test file header");
	*is_rpt = strcmp(buf, "unrestricted") != 0;
	if (fscanf(fp, "%zu %zu %lf", oa, va, e_ref) != 3)
		errx(1, "error parsing test file header");
	if (*is_rpt) {
		*ob = *oa;
		*vb = *va;
	} else {
		*ob = 0;
		*vb = 0;
	}
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
load_test_data(const char *testpath, size_t oa, size_t ob,
    size_t va, size_t vb, int is_rpt, double *d_ov, double *f_ov,
    double *t1, double *t2, double *i_oovo, double *i_oovv, double *i_ovvv)
{
	FILE *fp;
	size_t i, o, v;

	if ((fp = fopen(testpath, "r")) == NULL)
		err(1, "unable to open %s", testpath);

	if (is_rpt) {
		o = oa;
		v = va;
	} else {
		o = oa + ob;
		v = va + vb;
	}

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
	for (i = 0; i < o*v*v*(v-1)/2; i++) {
		i_ovvv[i] = read_next_double(fp);
	}
	if (is_rpt) {
		skip_line(fp);
		skip_line(fp);
		for (i = 0; i < o*v*v*v; i++) {
			i_ovvv[o*v*v*(v-1)/2+i] = read_next_double(fp);
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
load_random_data(size_t oa, size_t ob, size_t va, size_t vb, int is_rpt,
    double *d_ov, double *f_ov, double *t1, double *t2, double *i_oovo,
    double *i_oovv, double *i_ovvv)
{
	size_t i, o, v, nsp = is_rpt ? 2 : 1;

	if (is_rpt) {
		o = oa;
		v = va;
	} else {
		o = oa + ob;
		v = va + vb;
	}

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
	size_t oa = 0, ob = 0, va = 0, vb = 0;
	double e_pt = 0.0, e_ref = 0.0;
	double *d_ov, *f_ov, *t1, *t2, *i_oovv, *i_oovo, *i_ovvv;
	const char *testpath = NULL;
	time_t wall;
	long num;
	int ch, is_rpt = 1, rank = 0;

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	while ((ch = getopt(argc, argv, "o:t:uv:")) != -1) {
		switch (ch) {
		case 'o':
			num = strtol(optarg, NULL, 10);
			if (num < 1)
				err(1, "bad o value");
			oa = ob = (size_t)num;
			break;
		case 't':
			testpath = optarg;
			break;
		case 'u':
			is_rpt = 0;
			break;
		case 'v':
			num = strtol(optarg, NULL, 10);
			if (num < 1)
				err(1, "bad v value");
			va = vb = (size_t)num;
			break;
		default:
			usage();
			break;
		}
	}
	argv += optind;
	argc -= optind;

	if (testpath)
		load_test_header(testpath, &oa, &ob, &va, &vb, &is_rpt, &e_ref);

	if (is_rpt) {
		d_ov = xmalloc(oa*va*sizeof(double));
		f_ov = xmalloc(oa*va*sizeof(double));
		t1 = xmalloc(oa*va*sizeof(double));
		t2 = xmalloc(2*oa*oa*va*va*sizeof(double));
		i_oovo = xmalloc(2*oa*oa*oa*va*sizeof(double));
		i_oovv = xmalloc(2*oa*oa*va*va*sizeof(double));
		i_ovvv = xmalloc(oa*va*(va*(va-1)/2+va*va)*sizeof(double));
	} else {
		size_t o = oa + ob;
		size_t v = va + vb;
		d_ov = xmalloc(o*v*sizeof(double));
		f_ov = xmalloc(o*v*sizeof(double));
		t1 = xmalloc(o*v*sizeof(double));
		t2 = xmalloc(o*o*v*v*sizeof(double));
		i_oovo = xmalloc(o*o*o*v*sizeof(double));
		i_oovv = xmalloc(o*o*v*v*sizeof(double));
		i_ovvv = xmalloc(o*v*v*(v-1)/2*sizeof(double));
	}

	if (testpath) {
		load_test_data(testpath, oa, ob, va, vb, is_rpt, d_ov, f_ov,
		    t1, t2, i_oovo, i_oovv, i_ovvv);
	} else {
		load_random_data(oa, ob, va, vb, is_rpt, d_ov, f_ov,
		    t1, t2, i_oovo, i_oovv, i_ovvv);
	}

	if (rank == 0)
		printf("starting calculation...\n");
	wall = time(NULL);
	if (is_rpt) {
		e_pt = cc_rpt(oa, va, d_ov, f_ov, t1, t2,
		    i_oovo, i_oovv, i_ovvv);
	} else {
		e_pt = cc_upt(oa, ob, va, vb, d_ov, f_ov, t1, t2,
		    i_oovo, i_oovv, i_ovvv);
	}
	wall = time(NULL) - wall;
	if (rank == 0) {
		printf("done in %d sec\n", (int)wall);
		printf("cc (t) energy: % .8lf\n", e_pt);
	}
	if (testpath) {
		if (rank == 0)
			printf("cc (t) ref:    % .8lf\n", e_ref);
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
