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
	double *d_ov, *f_ov, *t1, *t2, *i_oovo, *i_oovv, *i_ovvv;
	size_t d_ov_sz, f_ov_sz, t1_sz, t2_sz, i_oovo_sz, i_oovv_sz, i_ovvv_sz;
	size_t oa, va, ob, vb;
	time_t wall = 0;
	int rank = 0, unrestricted;

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
	if (argc < 3)
		errx(1, "specify oa and va");
	oa = ob = strtol(argv[1], NULL, 10);
	va = vb = strtol(argv[2], NULL, 10);
	unrestricted = strstr(argv[0], "upt") != NULL;

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
	randomfill(d_ov, d_ov_sz);
	randomfill(f_ov, f_ov_sz);
	randomfill(t1, t1_sz);
	randomfill(t2, t2_sz);
	randomfill(i_oovo, i_oovo_sz);
	randomfill(i_oovv, i_oovv_sz);
	randomfill(i_ovvv, i_ovvv_sz);

	if (rank == 0) {
		printf("%s with oa = %zu, va = %zu, ob = %zu, vb = %zu\n",
		    unrestricted ? "upt" : "rpt", oa, va, ob, vb);
		wall = time(NULL);
	}
	if (unrestricted) {
		libpt_upt(oa, va, ob, vb, d_ov, f_ov, t1, t2,
		    i_oovo, i_oovv, i_ovvv);
	} else {
		libpt_rpt(oa, va, d_ov, f_ov, t1, t2,
		    i_oovo, i_oovv, i_ovvv);
	}
	if (rank == 0) {
		wall = time(NULL) - wall;
		printf("done in %d sec\n", (int)wall);
	}
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
	return (0);
}
