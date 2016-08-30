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

#include <mpi.h>

#include "pt.h"

#define EPSILON 1.0e-8

double drand48(void);
long long strtonum(const char *, long long, long long, const char **);

static void
usage(void)
{
	fprintf(stderr, "usage: pt [-o nocc] [-v nvirt] [-t test]\n");
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
load_test_header(const char *testpath, size_t *o, size_t *v, size_t *x,
    double *e_ref)
{
	FILE *fp;
	char buf[16];

	if ((fp = fopen(testpath, "r")) == NULL)
		err(1, "unable to open %s", testpath);
	if (fscanf(fp, "%4s", buf) != 1)
		errx(1, "error parsing test file header");
	if (strcmp(buf, "ri") == 0) {
		if (fscanf(fp, "%zu %zu %zu %lf", o, v, x, e_ref) != 4)
			errx(1, "error parsing test file header");
	} else if (strcmp(buf, "full") == 0) {
		if (fscanf(fp, "%zu %zu %lf", o, v, e_ref) != 3)
			errx(1, "error parsing test file header");
	} else
		errx(1, "unknown test type %s", buf);
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

//static void
//scan_next_line_st4(FILE *fp, struct st4 *st)
//{
//	double val;
//	unsigned a, b, c, d;
//	int rv;
//
//	rv = fscanf(fp, "%u %u %u %u %lf\n", &a, &b, &c, &d, &val);
//	if (rv != 5)
//		errx(1, "bad file format");
//	st->len++;
//	st->idx = xreallocarray(st->idx, st->len, sizeof(*st->idx));
//	st->data = xreallocarray(st->data, st->len, sizeof(*st->data));
//	st->idx[st->len-1].a = a;
//	st->idx[st->len-1].b = b;
//	st->idx[st->len-1].c = c;
//	st->idx[st->len-1].d = d;
//	st->data[st->len-1] = val;
//}

static void
load_test_data(const char *testpath, size_t o, size_t v, double *d_ov,
    double *f_ov, double *t1, double *t2, double *i_oovo, double *i_oovv,
    double *i_vvov)
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
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		t2[o*o*v*v+i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*o*v; i++) {
		i_oovo[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*o*v; i++) {
		i_oovo[o*o*o*v+i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		i_oovv[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		i_oovv[o*o*v*v+i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v*v*v; i++) {
		i_vvov[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v*v*v; i++) {
		i_vvov[o*v*v*v+i] = read_next_double(fp);
	}
	fclose(fp);
}

//static void
//load_test_data(const char *testpath, size_t o, size_t v, size_t x, double *d_ov,
//    double *f_ov, struct st4 *i_ooov, struct st4 *i_oovv, struct st4 *i_ovvv,
//    double *t1, struct st4 *t2, double *ovx, double *vvx)
//{
//	FILE *fp;
//	size_t i, len;
//
//	memset(t2, 0, sizeof(*t2));
//	memset(i_ooov, 0, sizeof(*i_ooov));
//	memset(i_oovv, 0, sizeof(*i_oovv));
//	memset(i_ovvv, 0, sizeof(*i_ovvv));
//
//	if ((fp = fopen(testpath, "r")) == NULL)
//		err(1, "unable to open %s", testpath);
//
//	skip_line(fp);
//	skip_line(fp);
//	for (i = 0; i < o*v; i++) {
//		d_ov[i] = read_next_double(fp);
//	}
//	skip_line(fp);
//	skip_line(fp);
//	for (i = 0; i < o*v; i++) {
//		f_ov[i] = read_next_double(fp);
//	}
//	skip_line(fp);
//	skip_line(fp);
//	for (i = 0; i < o*v; i++) {
//		t1[i] = read_next_double(fp);
//	}
//	skip_line(fp);
//	skip_line(fp);
//	if (fscanf(fp, "%zu\n", &len) != 1)
//		errx(1, "bad file format");
//	for (i = 0; i < len; i++)
//		scan_next_line_st4(fp, t2);
//	skip_line(fp);
//	if (fscanf(fp, "%zu\n", &len) != 1)
//		errx(1, "bad file format");
//	for (i = 0; i < len; i++)
//		scan_next_line_st4(fp, i_ooov);
//	skip_line(fp);
//	if (fscanf(fp, "%zu\n", &len) != 1)
//		errx(1, "bad file format");
//	for (i = 0; i < len; i++)
//		scan_next_line_st4(fp, i_oovv);
//	skip_line(fp);
//	if (x == 0) { /* full */
//		if (fscanf(fp, "%zu\n", &len) != 1)
//			errx(1, "bad file format");
//		for (i = 0; i < len; i++)
//			scan_next_line_st4(fp, i_ovvv);
//	} else { /* ri */
//		for (i = 0; i < o*v*x; i++)
//			ovx[i] = read_next_double(fp);
//		skip_line(fp);
//		skip_line(fp);
//		for (i = 0; i < v*v*x; i++)
//			vvx[i] = read_next_double(fp);
//	}
//	fclose(fp);
//}

static double
random_double(void)
{
	return (drand48() / 10000.0);
}

static void
load_random_data(size_t o, size_t v, double *d_ov,
    double *f_ov, double *t1, double *t2, double *i_oovo,
    double *i_oovv, double *i_vvov)
{
	size_t i;

	for (i = 0; i < o*v; i++) {
		d_ov[i] = random_double();
	}
	for (i = 0; i < o*v; i++) {
		f_ov[i] = random_double();
	}
	for (i = 0; i < o*v; i++) {
		t1[i] = random_double();
	}
	for (i = 0; i < 2*o*o*v*v; i++) {
		t2[i] = random_double();
	}
	for (i = 0; i < 2*o*o*o*v; i++) {
		i_oovo[i] = random_double();
	}
	for (i = 0; i < 2*o*o*v*v; i++) {
		i_oovv[i] = random_double();
	}
	for (i = 0; i < 2*o*v*v*v; i++) {
		i_vvov[i] = random_double();
	}
}

#if 0
#define T2(i, j, a, b) t2[i*o*v*v+j*v*v+a*v+b]
#define I_OOOV(i, j, k, a) i_ooov[i*o*o*v+j*o*v+k*v+a]
#define I_OOVV(i, j, a, b) i_oovv[i*o*v*v+j*v*v+a*v+b]
#define I_OVVV(i, a, b, c) i_ovvv[i*v*v*v+a*v*v+b*v+c]

static void
convert_t2(size_t o, size_t v, const double *t2, struct st4 *st)
{
	size_t i, j, b, c;

	memset(st, 0, sizeof(*st));

	for (i = 0; i < o; i++)
	for (j = 0; j < o; j++)
	for (b = 0; b < v; b++)
	for (c = 0; c < v; c++)
		if (T2(i, j, b, c) != 0.0) {
			st->len++;
			st->idx = xreallocarray(st->idx, st->len,
			    sizeof(*st->idx));
			st->data = xreallocarray(st->data, st->len,
			    sizeof(*st->data));
			st->idx[st->len-1].a = i;
			st->idx[st->len-1].b = j;
			st->idx[st->len-1].c = b;
			st->idx[st->len-1].d = c;
			st->data[st->len-1] = T2(i, j, b, c);
		}

	//printf("t2 %zu %zu\n", o*o*v*v, st->len);
}

static void
convert_ooov(size_t o, size_t v, const double *i_ooov, struct st4 *st)
{
	size_t i, j, k, c;

	memset(st, 0, sizeof(*st));

	for (i = 0; i < o; i++)
	for (j = 0; j < o; j++)
	for (k = 0; k < o; k++)
	for (c = 0; c < v; c++)
		if (I_OOOV(i, j, k, c) != 0.0) {
			st->len++;
			st->idx = xreallocarray(st->idx, st->len,
			    sizeof(*st->idx));
			st->data = xreallocarray(st->data, st->len,
			    sizeof(*st->data));
			st->idx[st->len-1].a = i;
			st->idx[st->len-1].b = j;
			st->idx[st->len-1].c = k;
			st->idx[st->len-1].d = c;
			st->data[st->len-1] = I_OOOV(i, j, k, c);
		}

	//printf("i_ooov %zu %zu\n", o*o*o*v, st->len);
}

static void
convert_oovv(size_t o, size_t v, const double *i_oovv, struct st4 *st)
{
	size_t i, j, b, c;

	memset(st, 0, sizeof(*st));

	for (i = 0; i < o; i++)
	for (j = 0; j < o; j++)
	for (b = 0; b < v; b++)
	for (c = 0; c < v; c++)
		if (I_OOVV(i, j, b, c) != 0.0) {
			st->len++;
			st->idx = xreallocarray(st->idx, st->len,
			    sizeof(*st->idx));
			st->data = xreallocarray(st->data, st->len,
			    sizeof(*st->data));
			st->idx[st->len-1].a = i;
			st->idx[st->len-1].b = j;
			st->idx[st->len-1].c = b;
			st->idx[st->len-1].d = c;
			st->data[st->len-1] = I_OOVV(i, j, b, c);
		}

	//printf("i_oovv %zu %zu\n", o*o*v*v, st->len);
}

static void
convert_ovvv(size_t o, size_t v, const double *i_ovvv, struct st4 *st)
{
	size_t i, a, b, c;

	memset(st, 0, sizeof(*st));

	for (i = 0; i < o; i++)
	for (a = 0; a < v; a++)
	for (b = 0; b < v; b++)
	for (c = 0; c < v; c++)
		if (I_OVVV(i, a, b, c) != 0.0) {
			st->len++;
			st->idx = xreallocarray(st->idx, st->len,
			    sizeof(*st->idx));
			st->data = xreallocarray(st->data, st->len,
			    sizeof(*st->data));
			st->idx[st->len-1].a = i;
			st->idx[st->len-1].b = a;
			st->idx[st->len-1].c = b;
			st->idx[st->len-1].d = c;
			st->data[st->len-1] = I_OVVV(i, a, b, c);
		}

	//printf("i_ovvv %zu %zu\n", o*v*v*v, st->len);
}

static void
print_st(const struct st4 *st)
{
	size_t i;

	printf("%zu\n", st->len);
	for (i = 0; i < st->len; i++) {
		printf("%3u %3u %3u %3u %24.16e\n", st->idx[i].a, st->idx[i].b,
		    st->idx[i].c, st->idx[i].d, st->data[i]);
	}
}

static void
setup_offsets(size_t ldim, struct st4 *st)
{
	size_t cur, i;
	int diff, j;

	st->offset = xmalloc((ldim+1)*sizeof(*st->offset));

	cur = 0;
	st->offset[cur] = 0;
	for (i = 0; i < st->len; i++) {
		diff = st->idx[i].d - st->idx[st->offset[cur]].d;
		if (diff > 0) {
			for (j = 0; j < diff; j++) {
				cur++;
				if (cur >= ldim)
					errx(1, "index out of range");
				st->offset[cur] = i;
			}
		} else if (diff < 0)
			errx(1, "last sparse index not sorted");
	}
	while (cur < ldim)
		st->offset[++cur] = st->len;
}
#endif

int
main(int argc, char **argv)
{
	size_t o = 0, v = 0, x = 0;
	double e_pt = 0.0, e_ref = 0.0;
	double *d_ov, *f_ov, *t1, *t2, *i_oovv, *i_oovo, *i_vvov;
	const char *errstr, *testpath = NULL;
	int rank;
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
		load_test_header(testpath, &o, &v, &x, &e_ref);
	d_ov = xmalloc(o*v*sizeof(double));
	f_ov = xmalloc(o*v*sizeof(double));
	t1 = xmalloc(o*v*sizeof(double));
	t2 = xmalloc(2*o*o*v*v*sizeof(double));
	i_oovo = xmalloc(2*o*o*o*v*sizeof(double));
	i_oovv = xmalloc(2*o*o*v*v*sizeof(double));
	i_vvov = xmalloc(2*o*v*v*v*sizeof(double));

	if (rank == 0) {
		if (testpath) {
			load_test_data(testpath, o, v, d_ov, f_ov, t1, t2,
			    i_oovo, i_oovv, i_vvov);
		} else {
			load_random_data(o, v, d_ov, f_ov, t1, t2,
			    i_oovo, i_oovv, i_vvov);
		}
	}

//	for (i = 0; i < o; i++) {
//	for (j = 0; j < o; j++) {
//	for (a = 0; a < v; a++) {
//	for (b = 0; b < v; b++) {
//		T2T(a,b,i,j) = T2(i,j,a,b);
//	}}}}
//	for (i = 0; i < o; i++) {
//	for (j = 0; j < o; j++) {
//	for (k = 0; k < o; k++) {
//	for (a = 0; a < v; a++) {
//		I_OOVO(i,j,a,k) = I_OOOV(i,j,k,a);
//	}}}}
//	for (i = 0; i < o; i++) {
//	for (a = 0; a < v; a++) {
//	for (b = 0; b < v; b++) {
//	for (c = 0; c < v; c++) {
//		I_VVOV(b,c,i,a) = I_OVVV(i,a,b,c);
//	}}}}
//	t2 += o*o*v*v;
//	t2t += o*o*v*v;
//	i_ooov += o*o*o*v;
//	i_oovo += o*o*o*v;
//	i_ovvv += o*v*v*v;
//	i_vvov += o*v*v*v;
//	for (i = 0; i < o; i++) {
//	for (j = 0; j < o; j++) {
//	for (a = 0; a < v; a++) {
//	for (b = 0; b < v; b++) {
//		T2T(a,b,i,j) = T2(i,j,a,b);
//	}}}}
//	for (i = 0; i < o; i++) {
//	for (j = 0; j < o; j++) {
//	for (k = 0; k < o; k++) {
//	for (a = 0; a < v; a++) {
//		I_OOVO(i,j,a,k) = I_OOOV(i,j,k,a);
//	}}}}
//	for (i = 0; i < o; i++) {
//	for (a = 0; a < v; a++) {
//	for (b = 0; b < v; b++) {
//	for (c = 0; c < v; c++) {
//		I_VVOV(b,c,i,a) = I_OVVV(i,a,b,c);
//	}}}}
//	t2 -= o*o*v*v;
//	t2t -= o*o*v*v;
//	i_ooov -= o*o*o*v;
//	i_oovo -= o*o*o*v;
//	i_ovvv -= o*v*v*v;
//	i_vvov -= o*v*v*v;

//	setup_offsets(v, &tt2);
//	setup_offsets(v, &it_ooov);
//	setup_offsets(v, &it_oovv);
//	setup_offsets(v, &it_ovvv);

	MPI_Bcast(d_ov, o*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(f_ov, o*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(t1, o*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(t2, 2*o*o*v*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(i_oovo, 2*o*o*o*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(i_oovv, 2*o*o*v*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(i_vvov, 2*o*v*v*v, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//	convert_t2(o, v, t2, &tt2);
//	convert_ooov(o, v, i_ooov, &it_ooov);
//	convert_oovv(o, v, i_oovv, &it_oovv);
//	convert_ovvv(o, v, i_ovvv, &it_ovvv);
//	printf("t2\n");
//	print_st(&tt2);
//	printf("ooov\n");
//	print_st(&it_ooov);
//	printf("oovv\n");
//	print_st(&it_oovv);
//	printf("ovvv\n");
//	print_st(&it_ovvv);

	if (rank == 0) {
		time_t t = time(NULL);
		printf("ccsd_pt: %s", ctime(&t));
	}

	e_pt = ccsd_pt(o, v, d_ov, f_ov, t1, t2, i_oovo, i_oovv, i_vvov);

	if (rank == 0) {
		time_t t = time(NULL);
		printf("ccsd_pt: %s", ctime(&t));
		printf("ccsd(t) energy: % .8lf\n", e_pt);
	}

//		size_t o11 = 0, o12 = o, v11 = 0, v12 = v;
//		size_t o1 = o12 - o11;
//		size_t v1 = v12 - v11;
//		double *d_ov1 = xmalloc(o1*v1*sizeof(double));
//		double *f_ov1 = xmalloc(o1*v1*sizeof(double));
//		double *t11 = xmalloc(o1*v1*sizeof(double));
//		double *t21 = xmalloc(2*o1*o1*v1*v1*sizeof(double));
//		double *i_oovo1 = xmalloc(2*o1*o1*o1*v1*sizeof(double));
//		double *i_oovv1 = xmalloc(2*o1*o1*v1*v1*sizeof(double));
//		double *i_vvov1 = xmalloc(2*o1*v1*v1*v1*sizeof(double));
//		for (size_t i = 0; i < o1; i++)
//		for (size_t j = 0; j < v1; j++) {
//			d_ov1[i*v1+j] = d_ov[(i+o11)*v+(j+v11)];
//			f_ov1[i*v1+j] = f_ov[(i+o11)*v+(j+v11)];
//			t11[i*v1+j] = t1[(i+o11)*v+(j+v11)];
//		}
//		for (size_t i = 0; i < o1; i++)
//		for (size_t j = 0; j < o1; j++)
//		for (size_t k = 0; k < v1; k++)
//		for (size_t l = 0; l < o1; l++) {
//			i_oovo1[i*o1*v1*o1+j*o1*v1+k*o1+l] =
//			    i_oovo[(i+o11)*o*v*o+(j+o11)*o*v+(k+v11)*o+(l+o11)];
//			i_oovo1[o1*o1*v1*o1+
//				i*o1*v1*o1+j*o1*v1+k*o1+l] =
//			    i_oovo[o*o*v*o+
//				(i+o11)*o*v*o+(j+o11)*o*v+(k+v11)*o+(l+o11)];
//		}
//		for (size_t i = 0; i < o1; i++)
//		for (size_t j = 0; j < o1; j++)
//		for (size_t k = 0; k < v1; k++)
//		for (size_t l = 0; l < v1; l++) {
//			i_oovv1[i*o1*v1*v1+j*v1*v1+k*v1+l] =
//			    i_oovv[(i+o11)*o*v*v+(j+o11)*v*v+(k+v11)*v+(l+v11)];
//			t21[i*o1*v1*v1+j*v1*v1+k*v1+l] =
//			    t2[(i+o11)*o*v*v+(j+o11)*v*v+(k+v11)*v+(l+v11)];
//			i_oovv1[o1*o1*v1*v1+
//				i*o1*v1*v1+j*v1*v1+k*v1+l] =
//			    i_oovv[o*o*v*v+
//				(i+o11)*o*v*v+(j+o11)*v*v+(k+v11)*v+(l+v11)];
//			t21[o1*o1*v1*v1+
//				i*o1*v1*v1+j*v1*v1+k*v1+l] =
//			    t2[o*o*v*v+
//				(i+o11)*o*v*v+(j+o11)*v*v+(k+v11)*v+(l+v11)];
//		}
//		for (size_t i = 0; i < v1; i++)
//		for (size_t j = 0; j < v1; j++)
//		for (size_t k = 0; k < o1; k++)
//		for (size_t l = 0; l < v1; l++) {
//			i_vvov1[i*v1*o1*v1+j*v1*o1+k*v1+l] =
//			    i_vvov[(i+v11)*v*o*v+(j+v11)*v*o+(k+o11)*v+(l+v11)];
//			i_vvov1[v1*v1*o1*v1+
//				i*v1*o1*v1+j*v1*o1+k*v1+l] =
//			    i_vvov[v*v*o*v+
//				(i+v11)*v*o*v+(j+v11)*v*o+(k+o11)*v+(l+v11)];
//		}
//		e_pt = ccsd_pt(o1, v1, d_ov1, f_ov1, t11, t21, i_oovo1,
//		    i_oovv1, i_vvov1);

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
	free(i_vvov);
	MPI_Finalize();
	return (fabs(e_pt - e_ref) < EPSILON ? 0 : 1);
}
