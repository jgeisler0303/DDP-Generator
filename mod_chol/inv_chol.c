#include "mex.h"
#include "mod_chol_inplace.h"


void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    int i, m, n, *P;
    double *A, *E, *invA, d;
    
    /* Check for proper number of input and output arguments */
    if (nrhs!=1) {
        mexErrMsgTxt("1 input argument required.");
    }
    if(nlhs<1 | nlhs>2){
        mexErrMsgTxt("1 to 2 output argument required.");
    }
    
    m= mxGetM(prhs[0]);
    n= mxGetN(prhs[0]);
    if(m!=n) {
        mexErrMsgTxt("Input must be a square matrix.\n");
    }

    A= mxMalloc(n*n*sizeof(double));
    for(i= 0; i<n*n; i++) A[i]= mxGetPr(prhs[0])[i];
	plhs[0]= mxCreateDoubleMatrix(n, n, mxREAL);
    invA= mxGetPr(plhs[0]);
    P= mxMalloc(n*sizeof(int));
    E= mxMalloc(n*sizeof(double));

    d= mod_chol(A, n, E, P);
    chol_inv(A, P, invA, n);
    if(nlhs>1) {
        plhs[1]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[1])[0]= d;
    }

    mxFree(P);
    mxFree(E);
    mxFree(A);
}
