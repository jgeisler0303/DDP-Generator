#include "mex.h"
#include "mod_chol_inplace.h"


void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    int i, m, n, *P;
    double *A, *L, *E;
    
    /* Check for proper number of input and output arguments */
    if (nrhs!=1) {
        mexErrMsgTxt("1 input argument required.");
    }
    if(nlhs<1 | nlhs>3){
        mexErrMsgTxt("1 to 3 output arguments required.");
    }
    
    m= mxGetM(prhs[0]);
    n= mxGetN(prhs[0]);
    if(m!=n || m<2) {
        mexErrMsgTxt("Input must be a square matrix.\n");
    }

    plhs[0]= mxCreateDoubleMatrix(n, n, mxREAL);
    A= mxGetPr(plhs[0]);
    for(i= 0; i<n*n; i++) A[i]= mxGetPr(prhs[0])[i];
    P= mxMalloc(n*sizeof(int));
    E= mxMalloc(n*sizeof(double));

    mod_chol(A, n, E, P);

    if(nlhs>1) {
        plhs[1]= mxCreateDoubleMatrix(1, n, mxREAL);
        for(i= 0; i<n; i++) mxGetPr(plhs[1])[i]= E[i];
    }
    if(nlhs>2) {
        plhs[2]= mxCreateDoubleMatrix(1, n, mxREAL);
        for(i= 0; i<n; i++) mxGetPr(plhs[2])[i]= P[i]+1;
    }

    mxFree(P);
    mxFree(E);
}
