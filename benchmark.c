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
#include <time.h>

#ifdef LIBPT_USE_MPI
#include <mpi.h>
#endif

#include "pt.h"

static void *
xmalloc(size_t size)
{
	void *ptr;

	if ((ptr = malloc(size)) == NULL)
		err(1, "malloc");
	return ptr;
}

static void
randomfill(double *arr, size_t arrlen)
{
	size_t i;

	for (i = 0; i < arrlen; i++)
		arr[i] = drand48() / 10000.0;
}

static void
benchmark_rpt(size_t oa, size_t va, size_t ob, size_t vb)
{
	double *d_ov, *f_ov, *t1, *t2;
	double *i_oovo, *i_oovv, *i_ovvv;
	size_t d_ov_sz, f_ov_sz, t1_sz, t2_sz;
	size_t i_oovo_sz, i_oovv_sz, i_ovvv_sz;

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

	randomfill(d_ov, d_ov_sz);
	randomfill(f_ov, f_ov_sz);
	randomfill(t1, t1_sz);
	randomfill(t2, t2_sz);
	randomfill(i_oovo, i_oovo_sz);
	randomfill(i_oovv, i_oovv_sz);
	randomfill(i_ovvv, i_ovvv_sz);

	libpt_rpt(oa, va, d_ov, f_ov, t1, t2, i_oovo, i_oovv, i_ovvv);

	free(d_ov);
	free(f_ov);
	free(t1);
	free(t2);
	free(i_oovo);
	free(i_oovv);
	free(i_ovvv);
}

static void
benchmark_upt(size_t oa, size_t va, size_t ob, size_t vb)
{
	double *d_ov, *f_ov, *t1, *t2;
	double *i_oovo, *i_oovv, *i_ovvv;
	size_t d_ov_sz, f_ov_sz, t1_sz, t2_sz;
	size_t i_oovo_sz, i_oovv_sz, i_ovvv_sz;

	d_ov_sz = oa*va + ob*vb;
	f_ov_sz = oa*va + ob*vb;
	t1_sz = oa*va + ob*vb;
	t2_sz = oa*oa*va*va + 2*oa*ob*va*vb + ob*ob*vb*vb;
	i_oovo_sz = oa*oa*va*oa + oa*ob*va*ob + ob*oa*vb*oa + ob*ob*vb*ob;
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

	randomfill(d_ov, d_ov_sz);
	randomfill(f_ov, f_ov_sz);
	randomfill(t1, t1_sz);
	randomfill(t2, t2_sz);
	randomfill(i_oovo, i_oovo_sz);
	randomfill(i_oovv, i_oovv_sz);
	randomfill(i_ovvv, i_ovvv_sz);

	libpt_upt(oa, va, ob, vb, d_ov, f_ov, t1, t2, i_oovo, i_oovv, i_ovvv);

	free(d_ov);
	free(f_ov);
	free(t1);
	free(t2);
	free(i_oovo);
	free(i_oovv);
	free(i_ovvv);
}

static void
benchmark_rft(size_t oa, size_t va, size_t ob, size_t vb)
{
	double *d_ov, *f2_ov, *l1, *t2, *l2;
	double *i_oovv, *i2_t2f2_oovo, *i3_ovvv, *i6_oovo, *i7_ovvv;
	size_t d_ov_sz, f2_ov_sz, l1_sz, t2_sz, l2_sz;
	size_t i_oovv_sz, i2_t2f2_oovo_sz, i3_ovvv_sz, i6_oovo_sz, i7_ovvv_sz;

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

	randomfill(d_ov, d_ov_sz);
	randomfill(f2_ov, f2_ov_sz);
	randomfill(l1, l1_sz);
	randomfill(t2, t2_sz);
	randomfill(l2, l2_sz);
	randomfill(i_oovv, i_oovv_sz);
	randomfill(i2_t2f2_oovo, i2_t2f2_oovo_sz);
	randomfill(i3_ovvv, i3_ovvv_sz);
	randomfill(i6_oovo, i6_oovo_sz);
	randomfill(i7_ovvv, i7_ovvv_sz);

	libpt_rft(oa, va, d_ov, f2_ov, l1, t2, l2, i_oovv,
	    i2_t2f2_oovo, i3_ovvv, i6_oovo, i7_ovvv);

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
benchmark_uft(size_t oa, size_t va, size_t ob, size_t vb)
{
	double *d_ov, *f2_ov, *l1, *t2, *l2;
	double *i_oovv, *i2_t2f2_oovo, *i3_ovvv, *i6_oovo, *i7_ovvv;
	size_t d_ov_sz, f2_ov_sz, l1_sz, t2_sz, l2_sz;
	size_t i_oovv_sz, i2_t2f2_oovo_sz, i3_ovvv_sz, i6_oovo_sz, i7_ovvv_sz;

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
	i6_oovo_sz = oa*oa*va*oa + oa*ob*va*ob + ob*oa*vb*oa + ob*ob*vb*ob;
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

	randomfill(d_ov, d_ov_sz);
	randomfill(f2_ov, f2_ov_sz);
	randomfill(l1, l1_sz);
	randomfill(t2, t2_sz);
	randomfill(l2, l2_sz);
	randomfill(i_oovv, i_oovv_sz);
	randomfill(i2_t2f2_oovo, i2_t2f2_oovo_sz);
	randomfill(i3_ovvv, i3_ovvv_sz);
	randomfill(i6_oovo, i6_oovo_sz);
	randomfill(i7_ovvv, i7_ovvv_sz);

	libpt_uft(oa, va, ob, vb, d_ov, f2_ov, l1, t2, l2, i_oovv,
	    i2_t2f2_oovo, i3_ovvv, i6_oovo, i7_ovvv);

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
	void (*fn)(size_t, size_t, size_t, size_t);
} benchmarks[] = {
	{ "rpt", benchmark_rpt },
	{ "upt", benchmark_upt },
	{ "rft", benchmark_rft },
	{ "uft", benchmark_uft },
};
static const size_t nbenchmarks = sizeof benchmarks / sizeof *benchmarks;

int
main(int argc, char **argv)
{
	size_t i, oa, va, ob, vb;
	time_t wall;
	int rank = 0;

#ifdef LIBPT_USE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (argc < 4)
		errx(1, "usage: benchmark rpt|upt|rft|uft o v");
	oa = ob = strtol(argv[2], NULL, 10);
	va = vb = strtol(argv[3], NULL, 10);
	if (rank == 0)
		printf("%s, oa = ob = %zu, va = vb = %zu\n", argv[1], oa, va);
	wall = time(NULL);
	for (i = 0; i < nbenchmarks; i++)
		if (strcmp(benchmarks[i].type, argv[1]) == 0) {
			benchmarks[i].fn(oa, va, ob, vb);
			break;
		}
	if (i == nbenchmarks)
		errx(1, "bad benchmark type");
	if (rank == 0)
		printf("done in %d sec\n", (int)(time(NULL)-wall));
#ifdef LIBPT_USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
