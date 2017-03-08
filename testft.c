#include <err.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "pt.h"

#define EPSILON 1.0e-8

static void
read_test_header(FILE *fp, size_t *o, size_t *v, double *e_ref)
{
	if (fscanf(fp, "%zu %zu %lf", o, v, e_ref) != 3)
		errx(1, "unable to read test header");
}

static void
read_test_data(FILE *fp)
{
}

int
main(int argc, char **argv)
{
	FILE *fp;
	double e_ft, e_ref, *f_ov, *d_ov, *l1, *t2, *l2, *i_oovv;
	double *i2_oovo, *i3_ovvv, *i6_oovo, *i7_ovvv;
	size_t o, v, size;

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
#endif
	if (argc < 2)
		errx(1, "specify test file");
	if ((fp = fopen(argv[1], "r")) == NULL)
		err(1, "fopen");
	read_test_header(fp, &o, &v, &e_ref);

	size = 3*o*v+3*o*o*v*v+2*o*o*o*v+2*o*v*v*v;
	if ((f_ov = malloc(size * sizeof(double))) == NULL)
		err(1, "malloc");
	d_ov = f_ov + o*v;
	l1 = d_ov + o*v;
	t2 = l1 + o*v;
	l2 = t2 + o*o*v*v;
	i_oovv = l2 + o*o*v*v;
	i2_oovo = i_oovv + o*o*v*v;
	i3_ovvv = i_oovv + o*o*o*v;
	i6_oovo = i_oovv + o*v*v*v;
	i7_ovvv = i_oovv + o*o*o*v;

	read_test_data(fp);
	fclose(fp);

	e_ft = cc_ft(o, v, f_ov, d_ov, l1, t2, l2, i_oovv,
	    i2_oovo, i3_ovvv, i6_oovo, i7_ovvv);

	free(f_ov);
#ifdef WITH_MPI
	MPI_Finalize();
#endif
	return (fabs(e_ft - e_ref) < EPSILON ? 0 : 1);
}
