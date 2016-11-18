// Copyright (c) 2016 Jens Geisler


#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

#include "iLQG.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // dims
    int N, n, m, m_, n_, si, i, j, k, l;
    int dims[4];
    // inputs
    double alpha;
    // outputs
    double *fx;
    double *fu;
#if FULL_DDP
    double *fxx;
    double *fxu;
    double *fuu;
#endif
    double *cx;
    double *cu;
    double *cxx;
    double *cxu;
    double *cuu;
    int out_idx;
    // aux
    tOptSet o= INIT_OPTSET;
    const mxArray *mxParams, *mxParam;

    if(nrhs!=6) { mexErrMsgIdAndTxt("MATLAB:minrhs", "wrong number of arguments (expected: x, u, l, L, alpha, params)"); return; }
#if FULL_DDP
    if(nlhs!=2 && nlhs!=12) { mexErrMsgIdAndTxt("MATLAB:minlhs", "wrong number of return values (expected: x_new,c[,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu])"); return; }
#else
    if(nlhs!=2 && nlhs!=9) { mexErrMsgIdAndTxt("MATLAB:minlhs", "wrong number of return values (expected: x_new,c[,fx,fu,cx,cu,cxx,cxu,cuu])"); return; }
#endif

    n= mxGetM(prhs[0]);
    m= mxGetM(prhs[1]);
    N= mxGetN(prhs[1]);
    
    if(n!=N_X) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of states (%d expected)", N_X); return; }
    if(m!=N_U) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of inputs (%d expected)", N_U); return; }
    
    if(mxGetNumberOfElements(prhs[0])!=n) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of elements in x (%dx1 expected)", n, N); return; }
    if(mxGetNumberOfElements(prhs[1])!=m*N) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of elements in u (%dx%d expected)", m, N); return; }
    
    if(!mxIsEmpty(prhs[2]) && mxGetNumberOfElements(prhs[2])!=m*N) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of elements in l (%dx%d or empty expected)", m, N); return; }
    if(mxIsEmpty(prhs[2]))
        if(!mxIsEmpty(prhs[3])) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of elements in L (l was empty, so L must be too)", m, n, N); return; }
    else        
        if(mxGetNumberOfElements(prhs[3])!=m*n*N) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of elements in L (%dx%dx%d or empty expected)", m, n, N); return; }
        
    if(mxIsEmpty(prhs[2]))
        if(!mxIsEmpty(prhs[4])) { mexErrMsgIdAndTxt("MATLAB:dimagree", "l was empty, so alpha must be too"); return; }
    else
        if(mxGetNumberOfElements(prhs[4])!=1) { mexErrMsgIdAndTxt("MATLAB:dimagree", "alpha must be a scalar"); return; }

    mxParams= prhs[5];
    if(!mxIsStruct(mxParams) || !mxGetNumberOfElements(mxParams)==1) {
        mexErrMsgIdAndTxt("MATLAB:dimagree", "Input 6 must be a scalar struct.\n");
    }

        
    o.n_hor= N;
    o.p= mxMalloc(n_params*sizeof(double *));
    for(i=0; i<n_params; i++) {
        si= (paramdesc[i]->size==-1)? o.n_hor+1: paramdesc[i]->size;
        if((mxParam= mxGetField(mxParams, 0, paramdesc[i]->name))==NULL) {
            mxFree(o.p);
            mexErrMsgIdAndTxt("MATLAB:dimagree", "Parameter name '%s' is not member of parameters struct.\n", paramdesc[i]->name);
        }
        m_= mxGetM(mxParam);
        n_= mxGetN(mxParam);
        if(mxIsSparse(mxParam) || !mxIsDouble(mxParam) || (m_!=1 && n_!=1) || (m_*n_!=si)) {
            mxFree(o.p);
            mexErrMsgIdAndTxt("MATLAB:dimagree", "Parameter name '%s' must be a vector length %d.\n", paramdesc[i]->name, si);
        }
        o.p[i]= mxGetPr(mxParam);
    }

    // aux
    o.trajectory= mxMalloc(sizeof(trajEl_t)*(N+1));
    o.candidates[0]= mxMalloc(sizeof(trajEl_t)*(N+1));

    for(k= 0; k<N; k++)
        for(i= 0; i<N_U; i++)
            o.trajectory[k].u[i]= mxGetPr(prhs[1])[MAT_IDX(i, k, N_U)];
    
    
    if(init_trajectory(o.trajectory, &o) && init_trajectory(o.candidates[0], &o)) {
        if(!mxIsEmpty(prhs[2])) {
            for(k= 0; k<N; k++) {
                for(i= 0; i<N_U; i++)
                    o.trajectory[k].l[i]= mxGetPr(prhs[2])[MAT_IDX(i, k, N_U)];

                for(i= 0; i<N_U*N_X; i++)
                    o.trajectory[k].L[i]= mxGetPr(prhs[2])[MAT_IDX(i, k, N_U*N_X)];
            }
            alpha= mxGetScalar(prhs[4]);
        } else
            alpha= 0.0;
        for(i= 0; i<N_X; i++)
            o.x0[i]= mxGetPr(prhs[0])[i];
        
        if(forward_pass(o.candidates[0], &o, alpha, &o.cost)) {
            makeCandidateNominal(&o, 0);
            plhs[0]= mxCreateDoubleMatrix(n, N+1, mxREAL);
            plhs[1]= mxCreateDoubleMatrix(1, N+1, mxREAL);
            for(k= 0; k<N+1; k++) {
                for(i= 0; i<N_X; i++)
                    mxGetPr(plhs[0])[MAT_IDX(i, k, N_X)]= o.trajectory[k].x[i];
                
                mxGetPr(plhs[1])[k]= o.trajectory[k].c;                    
            }
            if(nlhs>2) {
                if(calc_derivs(o)) {
                    out_idx= 2;
                    dims[0]= n; dims[1]= n; dims[2]= N;
                    fx= mxGetPr(plhs[out_idx]= mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL)); out_idx++;
                    dims[0]= n; dims[1]= m; dims[2]= N;
                    fu= mxGetPr(plhs[out_idx]= mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL)); out_idx++;
#if FULL_DDP
                    dims[0]= n; dims[1]= n; dims[2]= n; dims[3]= N;
                    fxx= mxGetPr(plhs[out_idx]= mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL)); out_idx++;
                    dims[0]= n; dims[1]= n; dims[2]= m; dims[3]= N;
                    fxu= mxGetPr(plhs[out_idx]= mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL)); out_idx++;
                    dims[0]= n; dims[1]= m; dims[2]= m; dims[3]= N;
                    fuu= mxGetPr(plhs[out_idx]= mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL)); out_idx++;
#endif
                    cx= mxGetPr(plhs[out_idx]= mxCreateDoubleMatrix(n, N+1, mxREAL)); out_idx++;
                    cu= mxGetPr(plhs[out_idx]= mxCreateDoubleMatrix(m, N, mxREAL)); out_idx++;
                    dims[0]= n; dims[1]= n; dims[2]= N+1;
                    cxx= mxGetPr(plhs[out_idx]= mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL)); out_idx++;
                    dims[0]= n; dims[1]= m; dims[2]= N;
                    cxu= mxGetPr(plhs[out_idx]= mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL)); out_idx++;
                    dims[0]= m; dims[1]= m; dims[2]= N;
                    cuu= mxGetPr(plhs[out_idx]= mxCreateDoubleMatrix(n*m, N, mxREAL)); out_idx++;
                    
                    for(k= 0; k<N+1; k++) {
                        for(i= 0; i<N_X; i++)
                            cx[MAT_IDX(i, k, N_X)]= o.trajectory[k].cx[i];
                        for(j= 0; j<N_X; j++)
                            for(i= 0; i<=j; i++) {
                                cxx[MAT_IDX3(i, j, k, N_X, N_X)]= o.trajectory[k].cxx[UTRI_MAT_IDX(i, j)];
                                if(i!=j)
                                    cxx[MAT_IDX3(j, i, k, N_X, N_X)]= o.trajectory[k].cxx[UTRI_MAT_IDX(i, j)];
                            }
                        if(k==N) break;
                        
                        for(i= 0; i<N_U; i++)
                            cu[MAT_IDX(i, k, N_U)]= o.trajectory[k].cu[i];
                        for(j= 0; j<N_U; j++)
                            for(i= 0; i<=j; i++) {
                                cuu[MAT_IDX3(i, j, k, N_U, N_U)]= o.trajectory[k].cuu[UTRI_MAT_IDX(i, j)];
                                if(i!=j)
                                    cuu[MAT_IDX3(j, i, k, N_U, N_U)]= o.trajectory[k].cuu[UTRI_MAT_IDX(i, j)];
                            }
                        for(i= 0; i<N_X; i++)
                            for(j= 0; j<N_U; j++)
                                cxu[MAT_IDX3(i, j, k, N_X, N_U)]= o.trajectory[k].cxu[MAT_IDX(i, j, N_X)];
                        
                        for(i= 0; i<N_X; i++)
                            for(j= 0; j<N_X; j++)
                                fx[MAT_IDX3(i, j, k, N_X, N_X)]= o.trajectory[k].fx[MAT_IDX(i, j, N_X)];
                        for(i= 0; i<N_X; i++)
                            for(j= 0; j<N_U; j++)
                                fu[MAT_IDX3(i, j, k, N_X, N_U)]= o.trajectory[k].fu[MAT_IDX(i, j, N_X)];
#if FULL_DDP
                        for(i= 0; i<N_X; i++)
                            for(l= 0; l<N_X; l++)
                                for(j= 0; j<=l; j++) {
                                    fxx[MAT_IDX4(i, j, l, k, N_X, N_X, N_X)]= o.trajectory[k].fxx[UTRI_MAT_IDX(j, l)+sizeofQxx*i];
                                    if(l!=j)
                                        fxx[MAT_IDX4(i, l, j, k, N_X, N_X, N_X)]= o.trajectory[k].fxx[UTRI_MAT_IDX(j, l)+sizeofQxx*i];
                                }
                        for(i= 0; i<N_X; i++)
                            for(l= 0; l<N_U; l++)
                                for(j= 0; j<=l; j++) {
                                    fuu[MAT_IDX4(i, j, l, k, N_X, N_U, N_U)]= o.trajectory[k].fuu[UTRI_MAT_IDX(j, l)+sizeofQuu*i];
                                    if(l!=j)
                                        fuu[MAT_IDX4(i, l, j, k, N_X, N_U, N_U)]= o.trajectory[k].fuu[UTRI_MAT_IDX(j, l)+sizeofQuu*i];
                                }
                        for(i= 0; i<N_X; i++)
                            for(j= 0; j<N_X; j++)
                                for(l= 0; l<N_U; l++)
                                    fxu[MAT_IDX4(i, j, l, k, N_X, N_X, N_U)]= o.trajectory[k].fxu[MAT_IDX(j, l, N_X)+sizeofQxu*i];
#endif
                    }
                }
            }
        }
        
    }
            
    mxFree(o.p);
    mxFree(o.trajectory);
    mxFree(o.candidates[0]);
}
