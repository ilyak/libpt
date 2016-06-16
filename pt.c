#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <err.h>
#include <getopt.h>

#include <xutil.h>

#define EPSILON 1.0e-8

static void __dead
usage(void)
{
	fprintf(stderr, "usage: pt [-l] [-o nocc] [-v nvirt] [-t test]\n");
	exit(1);
}

#define D_OV(i, a) \
    d_ov[i*v+a]
#define F_OV(i, a) \
    f_ov[i*v+a]
#define I_OOOV(i, j, k, a) \
    i_ooov[i*o*o*v+j*o*v+k*v+a]
#define I_OOVV(i, j, a, b) \
    i_oovv[i*o*v*v+j*v*v+a*v+b]
#define I_OVVV(i, a, b, c) \
    i_ovvv[i*v*v*v+a*v*v+b*v+c]
#define T1(i, a) \
    t1[i*v+a]
#define T2(i, j, a, b) \
    t2[i*o*v*v+j*v*v+a*v+b]
#define T3A(i, j, k, a, b, c) \
    t3a[i*o*o*v*v*v+j*o*v*v*v+k*v*v*v+a*v*v+b*v+c]
#define T3B(i, j, k, a, b, c) \
    t3b[i*o*o*v*v*v+j*o*v*v*v+k*v*v*v+a*v*v+b*v+c]

static void
ccsd_asymm_t3(size_t o, size_t v, double *t3a)
{
	size_t i, j, k, a, b, c;

	for (i = 0; i < o; i++) {
	for (j = i; j < o; j++) {
	for (k = j; k < o; k++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		double x;
		x = T3A(i, j, k, a, b, c) -
		    T3A(j, i, k, a, b, c) -
		    T3A(k, j, i, a, b, c) -
		    T3A(i, k, j, a, b, c) +
		    T3A(j, k, i, a, b, c) +
		    T3A(k, i, j, a, b, c);
		T3A(i, j, k, a, b, c) = x;
		T3A(j, i, k, a, b, c) = x;
		T3A(k, j, i, a, b, c) = x;
		T3A(i, k, j, a, b, c) = x;
		T3A(j, k, i, a, b, c) = x;
		T3A(k, i, j, a, b, c) = x;
	}}}}}}

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
	for (a = 0; a < v; a++) {
	for (b = a; b < v; b++) {
	for (c = b; c < v; c++) {
		double x;
		x = T3A(i, j, k, a, b, c) -
		    T3A(i, j, k, b, a, c) -
		    T3A(i, j, k, c, b, a) -
		    T3A(i, j, k, a, c, b) +
		    T3A(i, j, k, b, c, a) +
		    T3A(i, j, k, c, a, b);
		T3A(i, j, k, a, b, c) = x;
		T3A(i, j, k, b, a, c) = x;
		T3A(i, j, k, c, b, a) = x;
		T3A(i, j, k, a, c, b) = x;
		T3A(i, j, k, b, c, a) = x;
		T3A(i, j, k, c, a, b) = x;
	}}}}}}
}

static void
ccsd_t3a(size_t o, size_t v, double *t3a, const double *t2,
    const double *i_ooov, const double *i_ovvv)
{
	size_t i, j, k, l, a, b, c, d;
	size_t ooovvv = o * o * o * v * v * v;

	memset(t3a, 0, ooovvv * sizeof(double));

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
	for (d = 0; d < v; d++) {
		T3A(i, j, k, a, b, c) +=
		    T2(i, j, a, d) * I_OVVV(k, d, c, b);
	}}}}}}}

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
	for (l = 0; l < o; l++) {
		T3A(i, j, k, a, b, c) -=
		    T2(i, l, a, b) * I_OOOV(j, k, l, c);
	}}}}}}}

	ccsd_asymm_t3(o, v, t3a);
}

static void
ccsd_t3b(size_t o, size_t v, double *t3b, const double *t3a, const double *t1,
    const double *t2, const double *i_oovv, const double *f_ov)
{
	size_t i, j, k, a, b, c;
	size_t ooovvv = o * o * o * v * v * v;

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		T3B(i, j, k, a, b, c) =
		    T1(i, a) * I_OOVV(j, k, b, c) +
		    F_OV(i, a) * T2(j, k, b, c);
	}}}}}}

	ccsd_asymm_t3(o, v, t3b);

	for (i = 0; i < ooovvv; i++) {
		t3b[i] += t3a[i];
	}
}

static double
ccsd_pt_energy(size_t o, size_t v, const double *t3a, const double *t3b,
    const double *d_ov)
{
	double e_pt = 0.0;
	size_t i, j, k, a, b, c;

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		double dn = D_OV(i, a) + D_OV(i, b) + D_OV(i, c) +
			    D_OV(j, a) + D_OV(j, b) + D_OV(j, c) +
			    D_OV(k, a) + D_OV(k, b) + D_OV(k, c);
		e_pt += T3A(i, j, k, a, b, c) * T3B(i, j, k, a, b, c) / dn;
	}}}}}}

	return (e_pt);
}

static double
ccsd_pt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *i_ooov, const double *i_oovv, const double *i_ovvv,
    const double *t1, const double *t2)
{
	double e_pt = 0.0;
	double *t3a, *t3b;
	size_t ooovvv = o * o * o * v * v * v;

	t3a = xmalloc(ooovvv * sizeof(double));
	ccsd_t3a(o, v, t3a, t2, i_ooov, i_ovvv);

	t3b = xmalloc(ooovvv * sizeof(double));
	ccsd_t3b(o, v, t3b, t3a, t1, t2, i_oovv, f_ov);

	e_pt = ccsd_pt_energy(o, v, t3a, t3b, d_ov);
	e_pt *= (1.0 / 12.0 / 16.0);

	free(t3a);
	free(t3b);
	return (e_pt);
}

static void
load_test_header(const char *testpath, size_t *o, size_t *v, double *e_ref)
{
	FILE *fp;

	if ((fp = fopen(testpath, "r")) == NULL)
		err(1, "unable to open %s", testpath);
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
load_test_data(const char *testpath, size_t o, size_t v, double *d_ov,
    double *f_ov, double *i_ooov, double *i_oovv, double *i_ovvv,
    double *t1, double *t2)
{
	FILE *fp;
	size_t i, j, k, a, b, c;

	if ((fp = fopen(testpath, "r")) == NULL)
		err(1, "unable to open %s", testpath);

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
		D_OV(i, a) = read_next_double(fp);
	}}

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
		F_OV(i, a) = read_next_double(fp);
	}}

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
	for (a = 0; a < v; a++) {
		I_OOOV(i, j, k, a) = read_next_double(fp);
	}}}}

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
		I_OOVV(i, j, a, b) = read_next_double(fp);
	}}}}

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		I_OVVV(i, a, b, c) = read_next_double(fp);
	}}}}

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
		T1(i, a) = read_next_double(fp);
	}}

	skip_line(fp);
	skip_line(fp);
	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
		T2(i, j, a, b) = read_next_double(fp);
	}}}}

	fclose(fp);
}

static double
random_double(void)
{
	return (drand48());
}

static void
load_random_data(size_t o, size_t v, double *d_ov,
    double *f_ov, double *i_ooov, double *i_oovv, double *i_ovvv,
    double *t1, double *t2)
{
	size_t i, j, k, a, b, c;

	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
		D_OV(i, a) = random_double();
	}}

	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
		F_OV(i, a) = random_double();
	}}

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (k = 0; k < o; k++) {
	for (a = 0; a < v; a++) {
		I_OOOV(i, j, k, a) = random_double();
	}}}}

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
		I_OOVV(i, j, a, b) = random_double();
	}}}}

	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		I_OVVV(i, a, b, c) = random_double();
	}}}}

	for (i = 0; i < o; i++) {
	for (a = 0; a < v; a++) {
		T1(i, a) = random_double();
	}}

	for (i = 0; i < o; i++) {
	for (j = 0; j < o; j++) {
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
		T2(i, j, a, b) = random_double();
	}}}}
}

int
main(int argc, char **argv)
{
	size_t o = 4, v = 20;
	double *d_ov, *f_ov;
	double *i_ooov, *i_oovv, *i_ovvv;
	double *t1, *t2;
	double e_pt, e_ref = 0.0;
	const char *errstr, *testpath = NULL;
	char ch;

	while ((ch = getopt(argc, argv, "lo:t:v:")) != -1) {
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
		load_test_header(testpath, &o, &v, &e_ref);

	d_ov = xmalloc(o * v * sizeof(double));
	f_ov = xmalloc(o * v * sizeof(double));
	i_ooov = xmalloc(o * o * o * v * sizeof(double));
	i_oovv = xmalloc(o * o * v * v * sizeof(double));
	i_ovvv = xmalloc(o * v * v * v * sizeof(double));
	t1 = xmalloc(o * v * sizeof(double));
	t2 = xmalloc(o * o * v * v * sizeof(double));

	if (testpath) {
		load_test_data(testpath, o, v, d_ov, f_ov, i_ooov,
		    i_oovv, i_ovvv, t1, t2);
	} else {
		load_random_data(o, v, d_ov, f_ov, i_ooov,
		    i_oovv, i_ovvv, t1, t2);
	}

	e_pt = ccsd_pt(o, v, d_ov, f_ov, i_ooov, i_oovv, i_ovvv, t1, t2);
	if (testpath) {
		printf("ccsd(t) corr energy: % .6lf\n", e_pt);
		printf("ccsd(t) corr ref:    % .6lf\n", e_ref);
	} else {
		e_ref = e_pt;
	}

	free(d_ov);
	free(f_ov);
	free(i_ooov);
	free(i_oovv);
	free(i_ovvv);
	free(t1);
	free(t2);
	log_close();
	return (fabs(e_pt - e_ref) < EPSILON ? 0 : 1);
}
