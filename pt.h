#ifndef LIBPT_PT_H
#define LIBPT_PT_H

#include <sys/cdefs.h>
#include <sys/types.h>

__BEGIN_DECLS

/* Compute CCSD(T) energy correction in parallel.
 *
 * o/v - sizes of occupied/virtual spaces
 * d_ov - Delta matrix (size o*v)
 * f_ov - Fock matrix (size o*v)
 * i_ooov - OOOV integrals (size o*o*o*v)
 * i_oovv - OOVV integrals (size o*o*v*v)
 * i_ovvv - OVVV integrals (size o*v*v*v)
 * t1 - CCSD T1 amplitudes (size o*v)
 * t2 - CCSD T2 amplitudes (size o*o*v*v)
 *
 * The function returns CCSD(T) energy correction.
 */
double ccsd_pt(size_t o, size_t v, const double *d_ov,
    const double *f_ov, const double *i_ooov, const double *i_oovv,
    const double *i_ovvv, const double *t1, const double *t2);

__END_DECLS

#endif /* LIBPT_PT_H */
