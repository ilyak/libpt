#ifndef UT_H_INCLUDED
#define UT_H_INCLUDED

#ifndef __GNUC__
#define __attribute__(a)
#endif
#ifndef __dead
#define __dead __attribute__ ((__noreturn__))
#endif

long long strtonum(const char *, long long, long long, const char **);

#endif /* UT_H_INCLUDED */
