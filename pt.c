#include <limits.h>
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
#define T3A(a, b, c) \
    t3a[a*v*v+b*v+c]
#define T3B(a, b, c) \
    t3b[a*v*v+b*v+c]
#define T3AIJK(a, b, c) t3a[0*v*v*v+a*v*v+b*v+c]
#define T3AJIK(a, b, c) t3a[1*v*v*v+a*v*v+b*v+c]
#define T3AKJI(a, b, c) t3a[2*v*v*v+a*v*v+b*v+c]
#define T3AIKJ(a, b, c) t3a[3*v*v*v+a*v*v+b*v+c]
#define T3AJKI(a, b, c) t3a[4*v*v*v+a*v*v+b*v+c]
#define T3AKIJ(a, b, c) t3a[5*v*v*v+a*v*v+b*v+c]
#define T3BIJK(a, b, c) t3b[0*v*v*v+a*v*v+b*v+c]
#define T3BJIK(a, b, c) t3b[1*v*v*v+a*v*v+b*v+c]
#define T3BKJI(a, b, c) t3b[2*v*v*v+a*v*v+b*v+c]
#define T3BIKJ(a, b, c) t3b[3*v*v*v+a*v*v+b*v+c]
#define T3BJKI(a, b, c) t3b[4*v*v*v+a*v*v+b*v+c]
#define T3BKIJ(a, b, c) t3b[5*v*v*v+a*v*v+b*v+c]

static void
ccsd_asymm_t3(size_t v, double *t3a)
{
	size_t a, b, c, n;

#pragma omp parallel for private(a,b,c)
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		double x;
		x = T3AIJK(a, b, c) -
		    T3AJIK(a, b, c) -
		    T3AKJI(a, b, c) -
		    T3AIKJ(a, b, c) +
		    T3AJKI(a, b, c) +
		    T3AKIJ(a, b, c);
		T3AIJK(a, b, c) = x;
		T3AJIK(a, b, c) = x;
		T3AKJI(a, b, c) = x;
		T3AIKJ(a, b, c) = x;
		T3AJKI(a, b, c) = x;
		T3AKIJ(a, b, c) = x;
	}}}

	for (n = 0; n < 6; n++) {
#pragma omp parallel for private(a,b,c)
	for (a = 0; a < v; a++) {
	for (b = a; b < v; b++) {
	for (c = b; c < v; c++) {
		double x;
		x = T3A(a, b, c) -
		    T3A(b, a, c) -
		    T3A(c, b, a) -
		    T3A(a, c, b) +
		    T3A(b, c, a) +
		    T3A(c, a, b);
		T3A(a, b, c) = x;
		T3A(b, a, c) = x;
		T3A(c, b, a) = x;
		T3A(a, c, b) = x;
		T3A(b, c, a) = x;
		T3A(c, a, b) = x;
	}}}
		t3a += v * v * v;
	}
}

static void
ccsd_t3a(size_t o, size_t v, size_t i, size_t j, size_t k, double *t3a,
    const double *t2, const double *i_ooov, const double *i_ovvv)
{
	size_t l, a, b, c, d;
	size_t vvv = v * v * v;

	memset(t3a, 0, 6 * vvv * sizeof(double));

#pragma omp parallel for private(a,b,c,d)
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
	for (d = 0; d < v; d++) {
		T3AIJK(a, b, c) += T2(i, j, a, d) * I_OVVV(k, d, c, b);
		T3AJIK(a, b, c) += T2(j, i, a, d) * I_OVVV(k, d, c, b);
		T3AKJI(a, b, c) += T2(k, j, a, d) * I_OVVV(i, d, c, b);
		T3AIKJ(a, b, c) += T2(i, k, a, d) * I_OVVV(j, d, c, b);
		T3AJKI(a, b, c) += T2(j, k, a, d) * I_OVVV(i, d, c, b);
		T3AKIJ(a, b, c) += T2(k, i, a, d) * I_OVVV(j, d, c, b);
	}}}}

#pragma omp parallel for private(a,b,c,l)
	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
	for (l = 0; l < o; l++) {
		T3AIJK(a, b, c) -= T2(i, l, a, b) * I_OOOV(j, k, l, c);
		T3AJIK(a, b, c) -= T2(j, l, a, b) * I_OOOV(i, k, l, c);
		T3AKJI(a, b, c) -= T2(k, l, a, b) * I_OOOV(j, i, l, c);
		T3AIKJ(a, b, c) -= T2(i, l, a, b) * I_OOOV(k, j, l, c);
		T3AJKI(a, b, c) -= T2(j, l, a, b) * I_OOOV(k, i, l, c);
		T3AKIJ(a, b, c) -= T2(k, l, a, b) * I_OOOV(i, j, l, c);
	}}}}
}

static void
ccsd_t3b(size_t o, size_t v, size_t i, size_t j, size_t k,
    double *t3b, const double *t1, const double *t2,
    const double *i_oovv, const double *f_ov)
{
	size_t a, b, c;

	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		T3BIJK(a, b, c) = T1(i, a) * I_OOVV(j, k, b, c) +
		    F_OV(i, a) * T2(j, k, b, c);
		T3BJIK(a, b, c) = T1(j, a) * I_OOVV(i, k, b, c) +
		    F_OV(i, a) * T2(j, k, b, c);
		T3BKJI(a, b, c) = T1(k, a) * I_OOVV(j, i, b, c) +
		    F_OV(i, a) * T2(j, k, b, c);
		T3BIKJ(a, b, c) = T1(i, a) * I_OOVV(k, j, b, c) +
		    F_OV(i, a) * T2(j, k, b, c);
		T3BJKI(a, b, c) = T1(j, a) * I_OOVV(k, i, b, c) +
		    F_OV(i, a) * T2(j, k, b, c);
		T3BKIJ(a, b, c) = T1(k, a) * I_OOVV(i, j, b, c) +
		    F_OV(i, a) * T2(j, k, b, c);
	}}}
}

static double
ccsd_pt_energy(size_t v, size_t i, size_t j, size_t k,
    const double *t3a, const double *t3b, const double *d_ov)
{
	double e_pt = 0.0;
	size_t a, b, c;

	for (a = 0; a < v; a++) {
	for (b = 0; b < v; b++) {
	for (c = 0; c < v; c++) {
		double dn = D_OV(i, a) + D_OV(i, b) + D_OV(i, c) +
			    D_OV(j, a) + D_OV(j, b) + D_OV(j, c) +
			    D_OV(k, a) + D_OV(k, b) + D_OV(k, c);
		e_pt += T3A(a, b, c) * T3B(a, b, c) / dn;
	}}}

	return (e_pt);
}

static double
ccsd_pt(size_t o, size_t v, const double *d_ov, const double *f_ov,
    const double *i_ooov, const double *i_oovv, const double *i_ovvv,
    const double *t1, const double *t2)
{
	double e_pt = 0.0;
	double *t3a, *t3b;
	size_t i, j, k, n, vvv = v * v * v;

	t3a = xmalloc(6 * vvv * sizeof(double));
	t3b = xmalloc(6 * vvv * sizeof(double));

	for (i = 0; i < o; i++) {
	for (j = i+1; j < o; j++) {
	for (k = j+1; k < o; k++) {
		ccsd_t3a(o, v, i, j, k, t3a, t2, i_ooov, i_ovvv);
		ccsd_asymm_t3(v, t3a);
		ccsd_t3b(o, v, i, j, k, t3b, t1, t2, i_oovv, f_ov);
		ccsd_asymm_t3(v, t3b);

		for (n = 0; n < vvv; n++) {
			t3b[n] += t3a[n];
		}

		e_pt += ccsd_pt_energy(v, i, j, k, t3a, t3b, d_ov);
		e_pt += ccsd_pt_energy(v, j, i, k, t3a, t3b, d_ov);
		e_pt += ccsd_pt_energy(v, k, j, i, t3a, t3b, d_ov);
		e_pt += ccsd_pt_energy(v, i, k, j, t3a, t3b, d_ov);
		e_pt += ccsd_pt_energy(v, j, k, i, t3a, t3b, d_ov);
		e_pt += ccsd_pt_energy(v, k, i, j, t3a, t3b, d_ov);
	}}}
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
		printf("ccsd(t) corr energy: % .8lf\n", e_pt);
		printf("ccsd(t) corr ref:    % .8lf\n", e_ref);
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
