#ifndef LIBPT_PT_H
#define LIBPT_PT_H

#include <sys/cdefs.h>

__BEGIN_DECLS

/* Computes CCSD(T) energy correction. */
double ccsd_pt(size_t o, size_t v, const double *d_ov,
    const double *f_ov, const double *i_ooov, const double *i_oovv,
    const double *i_ovvv, const double *t1, const double *t2);

__END_DECLS

#endif /* LIBPT_PT_H */
