#include <err.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "pt.h"

#define EPSILON 1.0e-8

int
main(int argc, char **argv)
{
	double e_ft, e_ref, *f_ov, *d_ov, *l1, *t2, *l2, *i_oovv;
	double *i2_oovo, *i3_ovvv, *i6_oovo, *i7_ovvv;
	size_t o, v;

#ifdef WITH_MPI
	MPI_Init(&argc, &argv);
#endif
	if (argc < 2)
		errx(1, "specify test file");

	e_ft = cc_ft(o, v, f_ov, d_ov, l1, t2, l2, i_oovv, i2_oovo,
	    i3_ovvv, i6_oovo, i7_ovvv);

#ifdef WITH_MPI
	MPI_Finalize();
#endif
	return (fabs(e_ft - e_ref) < EPSILON ? 0 : 1);
}
