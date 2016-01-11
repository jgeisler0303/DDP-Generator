#include <math.h>
#include "cholesky.h"
#include "printMat.h"
#include "matMult.h"

int cholesky_tri(const double *A, int n, double *L) {
    int i, j, k;
    double s;
    
    for(i = 0; i < n; i++)
        for(j = 0; j < (i+1); j++) {
            s= 0;
            for(k = 0; k < j; k++)
                s += L[UTRI_MAT_IDX(k, i)] * L[UTRI_MAT_IDX(k, j)];

            s= A[UTRI_MAT_IDX(j, i)] - s;
            if(i == j) {
                if(s<=0.0)
                    return 0;
                L[UTRI_MAT_IDX(j, i)]= sqrt(s);
            } else
//                L[i * n + j]= s / L[j * n + j];
                L[UTRI_MAT_IDX(j, i)]= 1.0 / L[UTRI_MAT_IDX(j, j)] * s;
        }
 
    return 1;
}

void cholesky_solve_tri(const double *L_, const double *b, double *x, int n) {
    int i, k;
    double sum;

    // Solve L*y= b;
    for(k= 0; k<n; k++) {
        sum= b[k];
        for(i= 0; i<k; i++)
            sum-= x[i]*L_[UTRI_MAT_IDX(i, k)];
        
        x[k]= sum/L_[UTRI_MAT_IDX(k, k)];
    }

    // Solve L'*X= Y;
    for(k= n-1; k >= 0; k--) {
       sum= x[k];
       for(i= k+1; i<n; i++)
             sum-= x[i]*L_[UTRI_MAT_IDX(k, i)];
       x[k]= sum/L_[UTRI_MAT_IDX(k, k)];
    }
}

void cholesky_tri_inv(const double *L_, double *invA, const int n, double *x) {
    int i, k, l;

    for(l= 0; l<n; l++) {
        x[l]= 1.0;
        for(k= l+1; k<n; k++) x[k]= 0.0;

        // Solve L*y= b;
        for(k= l; k<n; k++) {
            for(i= l; i<k; i++)
                x[k]-= x[i]*L_[UTRI_MAT_IDX(i, k)];
            x[k]/= L_[UTRI_MAT_IDX(k, k)];
        }

        // Solve L'*X= Y;
        for(k= n-1; k >= l; k--) {
            for(i= k+1; i<n; i++)
                 x[k]-= x[i]*L_[UTRI_MAT_IDX(k, i)];
            x[k]/= L_[UTRI_MAT_IDX(k, k)];
            
            invA[UTRI_MAT_IDX(l, k)]= x[k];
        }
    }
}