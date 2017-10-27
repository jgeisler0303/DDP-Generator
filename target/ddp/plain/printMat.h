#ifndef PRINTMAT_H
#define PRINTMAT_H

#include <string.h>
#include <stdio.h>

#ifndef PRNT
#define PRNT printf
#endif

void printVec(const double *A, const int n, const char *nm);
void printTri(const double *A, const int n, const char *nm);
void printMat(const double *A, const int n, const int m, const char *nm);

#endif // PRINTMAT_H
