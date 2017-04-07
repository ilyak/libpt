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
#include <time.h>

#include <err.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "pt.h"

static void
randomfill(double *arr, size_t arrlen)
{
	size_t i;

	for (i = 0; i < arrlen; i++)
		*arr++ = drand48() / 10000.0;
}

int
main(int argc, char **argv)
{
	double *d_ov, *f2_ov, *l1, *t2, *l2;
	double *i_oovv, *i2_t2f2_oovo, *i3_ovvv, *i6_oovo, *i7_ovvv;
	size_t d_ov_sz, f2_ov_sz, l1_sz, t2_sz, l2_sz;
	size_t i_oovv_sz, i2_t2f2_oovo_sz, i3_ovvv_sz, i6_oovo_sz, i7_ovvv_sz;
	size_t oa, va, ob, vb;
	time_t wall = 0;
	int rank = 0;

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (argc < 3)
		errx(1, "specify oa and va");
	oa = ob = strtol(argv[1], NULL, 10);
	va = vb = strtol(argv[2], NULL, 10);

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
	if ((d_ov = malloc(d_ov_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((f2_ov = malloc(f2_ov_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((l1 = malloc(l1_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((t2 = malloc(t2_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((l2 = malloc(l2_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((i_oovv = malloc(i_oovv_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((i2_t2f2_oovo = malloc(i2_t2f2_oovo_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((i3_ovvv = malloc(i3_ovvv_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((i6_oovo = malloc(i6_oovo_sz * sizeof(double))) == NULL)
		err(1, "malloc");
	if ((i7_ovvv = malloc(i7_ovvv_sz * sizeof(double))) == NULL)
		err(1, "malloc");
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

	if (rank == 0) {
		printf("rft with oa = %zu, va = %zu\n", oa, va);
		wall = time(NULL);
	}
	libpt_rft(oa, va, d_ov, f2_ov, l1, t2, l2, i_oovv,
	    i2_t2f2_oovo, i3_ovvv, i6_oovo, i7_ovvv);
	if (rank == 0) {
		wall = time(NULL) - wall;
		printf("rft done in %d sec\n", (int)wall);
	}
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
#ifdef WITH_MPI
	MPI_Finalize();
#endif
	return (0);
}
