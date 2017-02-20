// mkoctfile --mex -I../common -DPRNT=mexPrintf -DN_U=4 testBoxQP_mex.cpp boxQP.cpp
// [x, is_clamped, lltHfree, status]= testBoxQP_mex(eye(4), ones(4, 1), ones(4, 1)*-inf, ones(4, 1)*inf, zeros(4, 1))

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <time.h>

#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

#include "boxQP.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    if(nrhs!=5) { mexErrMsgIdAndTxt("MATLAB:minrhs", "wrong number of arguments (expected: H, g, lower, upper, x0)"); return; }
    if(nlhs!=4) { mexErrMsgIdAndTxt("MATLAB:minlhs", "wrong number of return values (expected: x, is_clamped, lltHfree, status)"); return; }

    int m= mxGetM(prhs[0]);
    int n= mxGetN(prhs[0]);
    
    if(n!=N_U || m!=N_U) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong dimensions of H (%d x %d expected)", N_U, N_U); return; }

    m= mxGetNumberOfElements(prhs[1]);
    if(m!=N_U) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong dimensions of g (%d expected)", N_U); return; }
    
    m= mxGetNumberOfElements(prhs[2]);
    if(m!=N_U) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong dimensions of lower (%d expected)", N_U); return; }

    m= mxGetNumberOfElements(prhs[3]);
    if(m!=N_U) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong dimensions of upper (%d expected)", N_U); return; }

    m= mxGetNumberOfElements(prhs[4]);
    if(m!=N_U) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong dimensions of x0 (%d expected)", N_U); return; }
    
    MatrixUU H(mxGetPr(prhs[0]));
    VectorU g(mxGetPr(prhs[1]));
    VectorU lower(mxGetPr(prhs[2]));
    VectorU upper(mxGetPr(prhs[3]));

    VectorU x(mxGetPr(prhs[4]));
    VectorU_int is_clamped;
    LLT<MatrixUU_dyn, Upper> llt;
    
    plhs[3]= mxCreateDoubleMatrix(1, 1, mxREAL);
    
    int n_free= 0;
    
    mxGetPr(plhs[3])[0]= boxQP(H, g, lower, upper, x, is_clamped, n_free, llt);
    
    plhs[0]= mxCreateDoubleMatrix(N_U, 1, mxREAL);
    for(int i= 0; i<N_U; i++) mxGetPr(plhs[0])[i]= x(i);

    plhs[1]= mxCreateNumericMatrix(N_U, 1, mxINT32_CLASS, mxREAL);
    for(int i= 0; i<N_U; i++) ((int32_t *)mxGetData(plhs[1]))[i]= is_clamped(i);

    plhs[2]= mxCreateDoubleMatrix(n_free, n_free, mxREAL);
    for(int j= 0; j<llt.cols(); j++) for(int i= 0; i<=j; i++) mxGetPr(plhs[2])[i + j*n_free]= llt.matrixU()(i, j);
}

 