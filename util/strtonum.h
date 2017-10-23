#ifndef STRTONUM_H
#define STRTONUM_H

#include <stdio.h>      /* printf */
#include <stdlib.h>     /* strtol */
#include <ctype.h>
#include <errno.h>
#include <limits.h>

int strtoint(const char* str, long *ret, const char **ret_ep);
int strtodouble(const char* str, double *ret, const char **ret_ep);
const char *strtonum_error(int i);

#endif // STRTONUM_H