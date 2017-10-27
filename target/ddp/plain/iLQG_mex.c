// MATLAB Mex function wrapper for iLQG algorithm
// Copyright (c) 2016 Jens Geisler


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

#include "iLQG.h"
#include "printMat.h"
#include "matMult.h"

#define OUTPUT_FUNC_VEC(e, n, ef, Nf) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [2]; \
                                dims[0]= n; dims[1]= N+Nf; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++) \
                                    for(int i_n= 0; i_n<n; i_n++) \
                                        a[MAT_IDX(i_n, i_N, n)] = t->t[i_N].e[i_n];\
                                \
                                if(Nf) { \
                                    for(int i_n= 0; i_n<n; i_n++) \
                                        a[MAT_IDX(i_n, N, n)] = t->f.ef[i_n];\
                                } \
                            }
            
#define OUTPUT_FUNC_TRI(e, n, ef, Nf) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [3]; \
                                dims[0]= n; dims[1]= n; dims[2]= N+Nf; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++) \
                                    for(int i_m= 0; i_m<n; i_m++) \
                                        for(int i_n= 0; i_n<n; i_n++) \
                                            a[MAT_IDX3(i_n, i_m, i_N, n, n)] = t->t[i_N].e[SYMTRI_MAT_IDX(i_n, i_m)];\
                                \
                                if(Nf) { \
                                    for(int i_m= 0; i_m<n; i_m++) \
                                        for(int i_n= 0; i_n<n; i_n++) \
                                            a[MAT_IDX3(i_n, i_m, N, n, n)] = t->f.ef[SYMTRI_MAT_IDX(i_n, i_m)];\
                                } \
                            }
            
#define OUTPUT_FUNC_MAT(e, n, m, ef, Nf) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [3]; \
                                dims[0]= n; dims[1]= m; dims[2]= N+Nf; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++) \
                                    for(int i_m= 0; i_m<m; i_m++) \
                                        for(int i_n= 0; i_n<n; i_n++) \
                                            a[MAT_IDX3(i_n, i_m, i_N, n, m)] = t->t[i_N].e[MAT_IDX(i_n, i_m, n)];\
                                \
                                if(Nf+Nf) { \
                                    for(int i_m= 0; i_m<m; i_m++) \
                                        for(int i_n= 0; i_n<n; i_n++) \
                                            a[MAT_IDX3(i_n, i_m, N, n, m)] = t->f.ef[MAT_IDX(i_n, i_m, n)];\
                                } \
                            }
                            
#define OUTPUT_FUNC_MAT3(e, n1, n2, n3) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [4]; \
                                dims[0]= n1; dims[1]= n2; dims[2]= n3; dims[3]= N; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++) \
                                    for(int i3= 0; i3<n3; i3++) \
                                        for(int i2= 0; i2<n2; i2++) \
                                            for(int i1= 0; i1<n1; i1++) \
                                                a[MAT_IDX4(i1, i2, i3, i_N, n1, n2, n3)] = t->t[i_N].e[MAT_IDX3(i1, i2, i3, n1, n2)];\
                            }

#define OUTPUT_FUNC_TRI3(e, n1, n12, n3) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [4]; \
                                dims[0]= n1; dims[1]= n1; dims[2]= n3; dims[3]= N; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++) \
                                    for(int i3= 0; i3<n3; i3++) \
                                        for(int i2= 0; i2<n1; i2++) \
                                            for(int i1= 0; i1<n1; i1++) \
                                                a[MAT_IDX4(i1, i2, i3, i_N, n1, n1, n3)] = t->t[i_N].e[SYMTRI_MAT_IDX(i1, i2) + i3*n12];\
                            }

OUTPUT_FUNC_VEC(cx, N_X, cx, 1)
OUTPUT_FUNC_VEC(cu, N_U, cxx, 0)
OUTPUT_FUNC_TRI(cxx, N_X, cxx, 1)
OUTPUT_FUNC_TRI(cuu, N_U, cxx, 0)
OUTPUT_FUNC_MAT(cxu, N_X, N_U, cxx, 0)
OUTPUT_FUNC_MAT(fx, N_X, N_X, cxx, 0)
OUTPUT_FUNC_MAT(fu, N_X, N_U, cxx, 0)
#if FULL_DDP
OUTPUT_FUNC_TRI3(fxx, N_X, sizeofQxx, N_X)
OUTPUT_FUNC_TRI3(fuu, N_U, sizeofQuu, N_X)
OUTPUT_FUNC_MAT3(fxu, N_X, N_U, N_X)
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // dims
    int N, n, m, m_, n_, si, i, k;
    // inputs
    tOptSet o= INIT_OPTSET;
    // outputs
    double *success, *new_cost, *u_nom, *x_new, *u_new;
    // aux
    const mxArray *mxParams, *mxOptParam, *mxParam;
    const char *err_msg, *fname;
    clock_t begin, end;
    int calc_system_only= 0;

    if(nrhs!=4) { mexErrMsgIdAndTxt("MATLAB:minrhs", "wrong number of arguments (expected: x0, u_nom, params, opt_params)"); return; }
    if(nlhs<4 || nlhs>14) { mexErrMsgIdAndTxt("MATLAB:minlhs", "wrong number of return values (expected: success, x_new, u_new, new_cost, [cx, cu, cxx, cuu, cxu, fx, fu, fxx, fuu, fxu])"); return; }

    n= mxGetNumberOfElements(prhs[0]);
    m= mxGetM(prhs[1]);
    N= mxGetN(prhs[1])+1;
    
    if(n!=N_X) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of states (%d expected)", N_X); return; }
    if(m!=N_U) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of inputs (%d expected)", N_U); return; }
    
    if(mxGetNumberOfElements(prhs[0])!=n) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of elements in x0 (%d expected)", n); return; }
    if(mxGetNumberOfElements(prhs[1])!=m*(N-1)) { mexErrMsgIdAndTxt("MATLAB:dimagree", "wrong number of elements in u_nom (%dx%d expected)", m, N-1); return; }

    mxParams= prhs[2];
    if(!mxIsStruct(mxParams) || !mxGetNumberOfElements(mxParams)==1) {
        mexErrMsgIdAndTxt("MATLAB:dimagree", "Input 3 must be a scalar struct.\n");
    }
    mxOptParam= prhs[3];
    if(!mxIsStruct(mxOptParam) || !mxGetNumberOfElements(mxOptParam)==1) {
        mexErrMsgIdAndTxt("MATLAB:dimagree", "Input 4 must be a scalar struct of optimization parameters.\n");
    }

    // inputs
    for(i= 0; i<N_X; i++)
        o.x0[i]= mxGetPr(prhs[0])[i];
    u_nom= mxGetPr(prhs[1]);
    o.n_hor= N-1;
    
    standard_parameters(&o);
    for(i= 0; i<mxGetNumberOfFields(mxOptParam); i++) {
        mxParam= mxGetFieldByNumber(mxOptParam, 0, i);
        fname= mxGetFieldNameByNumber(mxOptParam, i);
        if(strcmp(fname, "max_iter")==0 && mxGetPr(mxParam)[0]<0.0)
            calc_system_only= 1;
        else {
            err_msg= setOptParam(&o, fname, mxGetPr(mxParam), mxGetNumberOfElements(mxParam));
            if(err_msg) {
                mexErrMsgIdAndTxt("MATLAB:dimagree", "Error setting optimization parameter '%s': %s.\n", fname, err_msg);
            }
        }
    }

    
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
    // printParams(o.p, o.n_hor);

    // outputs
    plhs[0]= mxCreateDoubleMatrix(1, 1, mxREAL);
    success= mxGetPr(plhs[0]);

    plhs[1]= mxCreateDoubleMatrix(n, N, mxREAL);
    x_new= mxGetPr(plhs[1]);

    plhs[2]= mxCreateDoubleMatrix(m, N-1, mxREAL);
    u_new= mxGetPr(plhs[2]);

    plhs[3]= mxCreateDoubleMatrix(1, 1, mxREAL);
    new_cost= mxGetPr(plhs[3]);

    // aux
    for(i= 0; i<NUMBER_OF_THREADS+1; i++) {
        o.trajectories[i].t= mxMalloc(sizeof(trajEl_t)*(N-1));
        // memset(o.trajectories[i].t, 1, sizeof(trajEl_t)*(N-1));
    }
    o.multipliers.t= mxMalloc(sizeof(multipliersEl_t)*N);
//     mexPrintf("sizeof trajEl_t: %d\n", sizeof(trajEl_t));

    
    mexPrintf("Set const vars\n");
    if(!init_opt(&o)) {
        success[0]= 0;
        new_cost[0]= o.cost;
    } else {
        mexPrintf("Initializing trajectory\n");
        for(k= 0; k<N-1; k++)
            for(i= 0; i<N_U; i++)
                o.nominal->t[k].u[i]= u_nom[MAT_IDX(i, k, N_U)];
        if(!forward_pass(o.candidates[0], &o, 0.0, &o.cost, 0)) {
            success[0]= 0;
            new_cost[0]= o.cost;
        } else {  
            makeCandidateNominal(&o, 0);
            if(calc_system_only) {
                mexPrintf("System calculation only\n");
                if(!calc_derivs(&o)) {
                    mexPrintf("Calculating derivatives failed.\n");
                }
            } else {
                mexPrintf("Starting iLQG\n");
                begin = clock();
                success[0]= iLQG(&o);
                end = clock();
                
                printLog(&o);
                mexPrintf("Time for iLQG: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
            }
            
            for(k= 0; k<N-1; k++)
                for(i= 0; i<N_X; i++)
                    x_new[MAT_IDX(i, k, N_X)]= o.nominal->t[k].x[i];
            
            for(i= 0; i<N_X; i++)
                x_new[MAT_IDX(i, N-1, N_X)]= o.nominal->f.x[i];
                
            for(k= 0; k<N-1; k++)
                for(i= 0; i<N_U; i++)
                    u_new[MAT_IDX(i, k, N_U)]= o.nominal->t[k].u[i];
                    
            new_cost[0]= o.cost;

            // TODO: add output of stat log
            if(nlhs>4) output_cx(plhs, 4, o.nominal, N-1);
            if(nlhs>5) output_cu(plhs, 5, o.nominal, N-1);
            if(nlhs>6) output_cxx(plhs, 6, o.nominal, N-1);
            if(nlhs>7) output_cuu(plhs, 7, o.nominal, N-1);
            if(nlhs>8) output_cxu(plhs, 8, o.nominal, N-1);
            if(nlhs>9) output_fx(plhs, 9, o.nominal, N-1);
            if(nlhs>10) output_fu(plhs, 10, o.nominal, N-1);
#if FULL_DDP
            if(nlhs>11) output_fxx(plhs, 11, o.nominal, N-1);
            if(nlhs>12) output_fuu(plhs, 12, o.nominal, N-1);
            if(nlhs>13) output_fxu(plhs, 13, o.nominal, N-1);
#endif
        }
    }
    
    mxFree(o.p);
    for(i= 0; i<NUMBER_OF_THREADS+1; i++)
        mxFree(o.trajectories[i].t);
    
    mxFree(o.multipliers.t);
}

 