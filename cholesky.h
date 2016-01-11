#ifndef CHOLESKY_H
#define CHOLESKY_H

int cholesky_tri(const double *A, int n, double *L);
void cholesky_tri_solve(const double *L_, const double *b, double *x, int n);
void cholesky_tri_inv(const double *L_, double *invA, const int n, double *x);

#endif // CHOLESKY_H
