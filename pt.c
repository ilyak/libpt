#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <err.h>
#include <getopt.h>
#include "xutil.h"

static double
ccsd_pt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *i_ooov, const double *i_oovv, const double *i_ovvv,
    const double *t1, const double *t2)
{
	double e_pt = 0.0;

	return (e_pt);
}

int
main(int argc, char **argv)
{
	size_t o = 2, v = 8;
	double *d_ov, *f_ov;
	double *i_ooov, *i_oovv, *i_ovvv;
	double *t1, *t2;
	double e_pt;
	const char *errstr;
	char ch;

	while ((ch = getopt(argc, argv, "lo:v:")) != -1) {
		switch (ch) {
		case 'l':
			log_add_level();
			log_open("pt");
			break;
		case 'o':
			o = strtonum(optarg, 1, INT_MAX, &errstr);
			if (errstr)
				errx(1, "bad o value: %s", errstr);
			break;
		case 'v':
			v = strtonum(optarg, 1, INT_MAX, &errstr);
			if (errstr)
				errx(1, "bad v value: %s", errstr);
			break;
		}
	}
	argv += optind;
	argc -= optind;

	d_ov = xmalloc(o * v * sizeof(double));
	f_ov = xmalloc(o * v * sizeof(double));
	i_ooov = xmalloc(o * o * o * v * sizeof(double));
	i_oovv = xmalloc(o * o * v * v * sizeof(double));
	i_ovvv = xmalloc(o * v * v * v * sizeof(double));
	t1 = xmalloc(o * v * sizeof(double));
	t2 = xmalloc(o * o * v * v * sizeof(double));

	e_pt = ccsd_pt(o, v, d_ov, f_ov, i_ooov, i_oovv, i_ovvv, t1, t2);
	printf("ccsd(t) energy: %lf\n", e_pt);

	free(d_ov);
	free(f_ov);
	free(i_ooov);
	free(i_oovv);
	free(i_ovvv);
	free(t1);
	free(t2);
	log_close();
	return (0);
}
