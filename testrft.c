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
read_test_header(FILE *fp, size_t *o, size_t *v, double *e_ref)
{
	if (fscanf(fp, "%zu %zu %lf", o, v, e_ref) != 3)
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
read_test_data(FILE *fp, size_t oa, size_t va, double *d_ov, double *f2_ov,
    double *l1, double *t2, double *l2, double *i_oovv, double *i2_t2f2_oovo,
    double *i3_ovvv, double *i6_oovo, double *i7_ovvv)
{
	size_t i, j, k, a, b, c;
	size_t ob = oa, vb = va, o = 2*oa, v = 2*va;
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
			/* *d_ov++ = */read_next_double(fp);
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
			/* *f2_ov++ = */read_next_double(fp);
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
			/* *l1++ = */read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		tmp[i] = read_next_double(fp);
	}
	for (i = 0; i < oa; i++) {
	for (j = 0; j < oa; j++) {
	for (a = 0; a < va; a++) {
	for (b = 0; b < va; b++) {
		*t2++ = tmp[i*o*v*v+j*v*v+a*v+b];
	}}}}
	for (i = 0; i < oa; i++) {
	for (j = 0; j < ob; j++) {
	for (a = 0; a < va; a++) {
	for (b = 0; b < vb; b++) {
		*t2++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	}}}}
//	for (i = 0; i < ob; i++) {
//	for (j = 0; j < ob; j++) {
//	for (a = 0; a < vb; a++) {
//	for (b = 0; b < vb; b++) {
//		*t2++ = tmp[(i+oa)*o*v*v+(j+oa)*v*v+(a+va)*v+(b+va)];
//	}}}}
//	for (i = 0; i < ob; i++) {
//	for (j = 0; j < oa; j++) {
//	for (a = 0; a < vb; a++) {
//	for (b = 0; b < va; b++) {
//		*t2++ = tmp[(i+oa)*o*v*v+j*v*v+(a+va)*v+b];
//	}}}}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		tmp[i] = read_next_double(fp);
	}
	for (i = 0; i < oa; i++) {
	for (j = 0; j < oa; j++) {
	for (a = 0; a < va; a++) {
	for (b = 0; b < va; b++) {
		*l2++ = tmp[i*o*v*v+j*v*v+a*v+b];
	}}}}
	for (i = 0; i < oa; i++) {
	for (j = 0; j < ob; j++) {
	for (a = 0; a < va; a++) {
	for (b = 0; b < vb; b++) {
		*l2++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	}}}}
//	for (i = 0; i < ob; i++) {
//	for (j = 0; j < ob; j++) {
//	for (a = 0; a < vb; a++) {
//	for (b = 0; b < vb; b++) {
//		*l2++ = tmp[(i+oa)*o*v*v+(j+oa)*v*v+(a+va)*v+(b+va)];
//	}}}}
//	for (i = 0; i < ob; i++) {
//	for (j = 0; j < oa; j++) {
//	for (a = 0; a < vb; a++) {
//	for (b = 0; b < va; b++) {
//		*l2++ = tmp[(i+oa)*o*v*v+j*v*v+(a+va)*v+b];
//	}}}}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		tmp[i] = read_next_double(fp);
	}
	for (i = 0; i < oa; i++) {
	for (j = 0; j < oa; j++) {
	for (a = 0; a < va; a++) {
	for (b = 0; b < va; b++) {
		*i_oovv++ = tmp[i*o*v*v+j*v*v+a*v+b];
	}}}}
	for (i = 0; i < oa; i++) {
	for (j = 0; j < ob; j++) {
	for (a = 0; a < va; a++) {
	for (b = 0; b < vb; b++) {
		*i_oovv++ = tmp[i*o*v*v+(j+oa)*v*v+a*v+(b+va)];
	}}}}
//	for (i = 0; i < ob; i++) {
//	for (j = 0; j < ob; j++) {
//	for (a = 0; a < vb; a++) {
//	for (b = 0; b < vb; b++) {
//		*i_oovv++ = tmp[(i+oa)*o*v*v+(j+oa)*v*v+(a+va)*v+(b+va)];
//	}}}}
//	for (i = 0; i < ob; i++) {
//	for (j = 0; j < oa; j++) {
//	for (a = 0; a < vb; a++) {
//	for (b = 0; b < va; b++) {
//		*i_oovv++ = tmp[(i+oa)*o*v*v+j*v*v+(a+va)*v+b];
//	}}}}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*o; i++) {
		tmp[i] = read_next_double(fp);
	}
	for (i = 0; i < oa; i++) {
	for (j = 0; j < oa; j++) {
	for (a = 0; a < va; a++) {
	for (k = 0; k < oa; k++) {
		*i2_t2f2_oovo++ = tmp[i*o*v*o+j*v*o+a*o+k];
	}}}}
	for (i = 0; i < oa; i++) {
	for (j = 0; j < ob; j++) {
	for (a = 0; a < va; a++) {
	for (k = 0; k < ob; k++) {
		*i2_t2f2_oovo++ = tmp[i*o*v*o+(j+oa)*v*o+a*o+(k+oa)];
	}}}}
//	for (i = 0; i < ob; i++) {
//	for (j = 0; j < ob; j++) {
//	for (a = 0; a < vb; a++) {
//	for (k = 0; k < ob; k++) {
//		*i2_t2f2_oovo++ = tmp[(i+oa)*o*v*o+(j+oa)*v*o+(a+va)*o+(k+oa)];
//	}}}}
//	for (i = 0; i < ob; i++) {
//	for (j = 0; j < oa; j++) {
//	for (a = 0; a < vb; a++) {
//	for (k = 0; k < oa; k++) {
//		*i2_t2f2_oovo++ = tmp[(i+oa)*o*v*o+j*v*o+(a+va)*o+k];
//	}}}}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < b; c++) {
		double t = read_next_double(fp);
		tmp[i*v*v*v+a*v*v+b*v+c] = t;
		tmp[i*v*v*v+a*v*v+c*v+b] = -t;
	}}}}
	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
		tmp[i*v*v*v+a*v*v+b*v+b] = 0;
	}}}
	for (i = 0; i < oa; i++) {
	for (a = 0; a < va; a++) {
	for (b = 0; b < va; b++) {
	for (c = 0; c < b; c++) {
		*i3_ovvv++ = tmp[i*v*v*v+a*v*v+b*v+c];
	}}}}
	for (i = 0; i < oa; i++) {
	for (a = 0; a < vb; a++) {
	for (b = 0; b < va; b++) {
	for (c = 0; c < vb; c++) {
		*i3_ovvv++ = tmp[i*v*v*v+(a+va)*v*v+b*v+(c+va)];
	}}}}
//	for (i = 0; i < ob; i++) {
//	for (a = 0; a < vb; a++) {
//	for (b = 0; b < vb; b++) {
//	for (c = 0; c < b; c++) {
//		*i3_ovvv++ = tmp[(i+oa)*v*v*v+(a+va)*v*v+(b+va)*v+(c+va)];
//	}}}}
//	for (i = 0; i < ob; i++) {
//	for (a = 0; a < va; a++) {
//	for (b = 0; b < vb; b++) {
//	for (c = 0; c < va; c++) {
//		*i3_ovvv++ = tmp[(i+oa)*v*v*v+a*v*v+(b+va)*v+c];
//	}}}}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*o; i++) {
		tmp[i] = read_next_double(fp);
	}
	for (i = 0; i < oa; i++) {
	for (j = 0; j < oa; j++) {
	for (a = 0; a < va; a++) {
	for (k = 0; k < oa; k++) {
		*i6_oovo++ = tmp[i*o*v*o+j*v*o+a*o+k];
	}}}}
	for (i = 0; i < oa; i++) {
	for (j = 0; j < ob; j++) {
	for (a = 0; a < va; a++) {
	for (k = 0; k < ob; k++) {
		*i6_oovo++ = tmp[i*o*v*o+(j+oa)*v*o+a*o+(k+oa)];
	}}}}
//	for (i = 0; i < ob; i++) {
//	for (j = 0; j < ob; j++) {
//	for (a = 0; a < vb; a++) {
//	for (k = 0; k < ob; k++) {
//		*i6_oovo++ = tmp[(i+oa)*o*v*o+(j+oa)*v*o+(a+va)*o+(k+oa)];
//	}}}}
//	for (i = 0; i < ob; i++) {
//	for (j = 0; j < oa; j++) {
//	for (a = 0; a < vb; a++) {
//	for (k = 0; k < oa; k++) {
//		*i6_oovo++ = tmp[(i+oa)*o*v*o+j*v*o+(a+va)*o+k];
//	}}}}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < b; c++) {
		double t = read_next_double(fp);
		tmp[i*v*v*v+a*v*v+b*v+c] = t;
		tmp[i*v*v*v+a*v*v+c*v+b] = -t;
	}}}}
	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
		tmp[i*v*v*v+a*v*v+b*v+b] = 0;
	}}}
	for (i = 0; i < oa; i++) {
	for (a = 0; a < va; a++) {
	for (b = 0; b < va; b++) {
	for (c = 0; c < b; c++) {
		*i7_ovvv++ = tmp[i*v*v*v+a*v*v+b*v+c];
	}}}}
	for (i = 0; i < oa; i++) {
	for (a = 0; a < vb; a++) {
	for (b = 0; b < va; b++) {
	for (c = 0; c < vb; c++) {
		*i7_ovvv++ = tmp[i*v*v*v+(a+va)*v*v+b*v+(c+va)];
	}}}}
//	for (i = 0; i < ob; i++) {
//	for (a = 0; a < vb; a++) {
//	for (b = 0; b < vb; b++) {
//	for (c = 0; c < b; c++) {
//		*i7_ovvv++ = tmp[(i+oa)*v*v*v+(a+va)*v*v+(b+va)*v+(c+va)];
//	}}}}
//	for (i = 0; i < ob; i++) {
//	for (a = 0; a < va; a++) {
//	for (b = 0; b < vb; b++) {
//	for (c = 0; c < va; c++) {
//		*i7_ovvv++ = tmp[(i+oa)*v*v*v+a*v*v+(b+va)*v+c];
//	}}}}
	free(tmp);
/*
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v; i++) {
		d_ov[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v; i++) {
		f2_ov[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v; i++) {
		l1[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		t2[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		l2[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*v; i++) {
		i_oovv[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*o; i++) {
		i2_t2f2_oovo[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v*v*(v-1)/2; i++) {
		i3_ovvv[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*o*v*o; i++) {
		i6_oovo[i] = read_next_double(fp);
	}
	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o*v*v*(v-1)/2; i++) {
		i7_ovvv[i] = read_next_double(fp);
	}
*/
}

int
main(int argc, char **argv)
{
	FILE *fp;
	double e_cmp, e_ref, *d_ov, *f2_ov, *l1, *t2, *l2;
	double *i_oovv, *i2_t2f2_oovo, *i3_ovvv, *i6_oovo, *i7_ovvv;
	size_t o, v, size;
	int rank = 0;

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (argc < 2)
		errx(1, "specify test file");
	if ((fp = fopen(argv[1], "r")) == NULL)
		err(1, "fopen");
	read_test_header(fp, &o, &v, &e_ref);

	size = 3*o*v + 6*o*o*v*v + 4*o*o*o*v + 2*o*v*v*v + 2*o*v*v*(v-1)/2;
	if ((d_ov = malloc(size * sizeof(double))) == NULL)
		err(1, "malloc");
	f2_ov = d_ov + o*v;
	l1 = f2_ov + o*v;
	t2 = l1 + o*v;
	l2 = t2 + 2*o*o*v*v;
	i_oovv = l2 + 2*o*o*v*v;
	i2_t2f2_oovo = i_oovv + 2*o*o*v*v;
	i3_ovvv = i2_t2f2_oovo + 2*o*o*o*v;
	i6_oovo = i3_ovvv + o*v*v*(v-1)/2 + o*v*v*v;
	i7_ovvv = i6_oovo + 2*o*o*o*v;

	read_test_data(fp, o, v, d_ov, f2_ov, l1, t2, l2, i_oovv,
	    i2_t2f2_oovo, i3_ovvv, i6_oovo, i7_ovvv);
	fclose(fp);

	e_cmp = libpt_rft(o, v, d_ov, f2_ov, l1, t2, l2, i_oovv,
	    i2_t2f2_oovo, i3_ovvv, i6_oovo, i7_ovvv);
	if (rank == 0) {
		printf("rft energy: % .8lf\n", e_cmp);
		printf("rft ref:    % .8lf\n", e_ref);
	}
	free(d_ov);
#ifdef WITH_MPI
	MPI_Finalize();
#endif
	return (fabs(e_cmp - e_ref) < EPSILON ? 0 : 1);
}
