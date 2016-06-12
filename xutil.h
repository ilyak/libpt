#ifndef XUTIL_H_INCLUDED
#define XUTIL_H_INCLUDED

#include <stdarg.h>
#include <stdlib.h>

#include "vis.h"

#ifndef __GNUC__
#define __attribute__(a)
#endif

#ifndef __unused
#define __unused __attribute__ ((__unused__))
#endif
#ifndef __dead
#define __dead __attribute__ ((__noreturn__))
#endif
#ifndef __packed
#define __packed __attribute__ ((__packed__))
#endif

#if !defined(__bounded__)
# define __bounded__(x, y, z)
#endif

/* Attribute to make gcc check printf-like arguments. */
#define printflike(a, b) __attribute__ ((format (printf, a, b)))

/* Number of items in array. */
#ifndef nitems
#define nitems(_a) (sizeof((_a)) / sizeof((_a)[0]))
#endif

/* asprintf.c */
int		 asprintf(char **, const char *, ...);
int		 vasprintf(char **, const char *, va_list);

/* reallocarray.c */
void		*reallocarray(void *, size_t, size_t size);

/* strcasestr.c */
char		*strcasestr(const char *, const char *);

/* strlcpy.c */
size_t	 	 strlcpy(char *, const char *, size_t);

/* strlcat.c */
size_t	 	 strlcat(char *, const char *, size_t);

/* strsep.c */
char		*strsep(char **, const char *);

/* strtonum.c */
long long	 strtonum(const char *, long long, long long, const char **);

/* xmalloc.c */
void	*xmalloc(size_t);
void	*xcalloc(size_t, size_t);
void	*xrealloc(void *, size_t);
void	*xreallocarray(void *, size_t, size_t);
char	*xstrdup(const char *);
int	 xasprintf(char **, const char *, ...)
		__attribute__((__format__ (printf, 2, 3)))
		__attribute__((__nonnull__ (2)));
int	 xvasprintf(char **, const char *, va_list)
		__attribute__((__nonnull__ (2)));
int	 xsnprintf(char *, size_t, const char *, ...)
		__attribute__((__format__ (printf, 3, 4)))
		__attribute__((__nonnull__ (3)))
		__attribute__((__bounded__ (__string__, 1, 2)));
int	 xvsnprintf(char *, size_t, const char *, va_list)
		__attribute__((__nonnull__ (3)))
		__attribute__((__bounded__ (__string__, 1, 2)));

/* log.c */
void	log_add_level(void);
int	log_get_level(void);
void	log_open(const char *);
void	log_close(void);
void printflike(1, 2) log_debug(const char *, ...);
__dead void printflike(1, 2) fatal(const char *, ...);
__dead void printflike(1, 2) fatalx(const char *, ...);

#endif /* XUTIL_H_INCLUDED */
