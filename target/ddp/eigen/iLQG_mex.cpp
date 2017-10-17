// MATLAB Mex function wrapper for iLQG algorithm
// Copyright (c) 2016 Jens Geisler


#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <time.h>

#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

#include "iLQG.hpp"

using namespace Eigen;

#define OUTPUT_FUNC_VEC_F(T, e, n) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [2]; \
                                dims[0]= n; dims[1]= N+1; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++, a+= n) {\
                                    Map<T> d(a); \
                                    d= t->t[i_N].e;\
                                } \
                                \
                                Map<T> d(a); \
                                d= t->f.e;\
                            }
            
#define OUTPUT_FUNC_VEC(T, e, n) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [2]; \
                                dims[0]= n; dims[1]= N; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++, a+= n) {\
                                    Map<T> d(a); \
                                    d= t->t[i_N].e;\
                                } \
                            }
            
#define OUTPUT_FUNC_TRI_F(T, e, n) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [3]; \
                                dims[0]= n; dims[1]= n; dims[2]= N+1; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++, a+= n*n) {\
                                    Map<T> d(a); \
                                    d= t->t[i_N].e.selfadjointView<Upper>();\
                                } \
                                \
                                Map<T> d(a); \
                                d= t->f.e.selfadjointView<Upper>();\
                            }
            
#define OUTPUT_FUNC_TRI(T, e, n) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [3]; \
                                dims[0]= n; dims[1]= n; dims[2]= N; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++, a+= n*n) {\
                                    Map<T> d(a); \
                                    d= t->t[i_N].e.selfadjointView<Upper>();\
                                } \
                            }
            
#define OUTPUT_FUNC_MAT(T, e, n, m) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [3]; \
                                dims[0]= n; dims[1]= m; dims[2]= N; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++, a+= n*m) { \
                                    Map<T> d(a); \
                                    d= t->t[i_N].e;\
                                } \
                            }
                            
#define OUTPUT_FUNC_MAT3(T, e, n1, n2, n3) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [4]; \
                                dims[0]= n1; dims[1]= n2; dims[2]= n3; dims[3]= N; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++) \
                                    for(int i3= 0; i3<n3; i3++, a+= n1*n2) { \
                                        Map<T> d(a); \
                                        d= t->t[i_N].e[i3];\
                                    } \
                            }

#define OUTPUT_FUNC_TRI3(T, e, n1, n2, n3) void output_ ## e(mxArray *plhs[], int out_idx, traj_t *t, int N) { \
                                double *a; \
                                int dims [4]; \
                                dims[0]= n1; dims[1]= n2; dims[2]= n3; dims[3]= N; \
                                a= mxGetPr(plhs[out_idx]= mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL)); \
                                for(int i_N= 0; i_N<N; i_N++) \
                                    for(int i3= 0; i3<n3; i3++, a+= n1*n2) { \
                                        Map<T> d(a); \
                                        d= t->t[i_N].e[i3].selfadjointView<Upper>(); \
                                     } \
                            }

OUTPUT_FUNC_VEC_F(VectorX, cx, N_X)
OUTPUT_FUNC_VEC(VectorU, cu, N_U)
OUTPUT_FUNC_TRI_F(MatrixXX, cxx, N_X)
OUTPUT_FUNC_TRI(MatrixUU, cuu, N_U)
OUTPUT_FUNC_MAT(MatrixXU, cxu, N_X, N_U)
OUTPUT_FUNC_MAT(MatrixXX, fx, N_X, N_X)
OUTPUT_FUNC_MAT(MatrixXU, fu, N_X, N_U)
#if FULL_DDP
OUTPUT_FUNC_TRI3(MatrixXX, fxx, N_X, N_X, N_X)
OUTPUT_FUNC_TRI3(MatrixUU, fuu, N_U, N_U, N_X)
OUTPUT_FUNC_MAT3(MatrixXU, fxu, N_X, N_U, N_X)
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int N, n, m, m_, n_, si;
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

    tOptSet *o= new optSet();
    standard_parameters(o);
    for(int i= 0; i<mxGetNumberOfFields(mxOptParam); i++) {
        mxParam= mxGetFieldByNumber(mxOptParam, 0, i);
        fname= mxGetFieldNameByNumber(mxOptParam, i);
        if(strcmp(fname, "max_iter")==0 && mxGetPr(mxParam)[0]<0.0)
            calc_system_only= 1;
        else {
            err_msg= setOptParam(o, fname, mxGetPr(mxParam), mxGetNumberOfElements(mxParam));
            if(err_msg) {
                mexErrMsgIdAndTxt("MATLAB:dimagree", "Error setting optimization parameter '%s': %s.\n", fname, err_msg);
            }
        }
    }

    
    o->n_hor= N-1;
    
    o->p= (double **)mxMalloc(n_params*sizeof(double *));
    for(int i=0; i<n_params; i++) {
        si= (paramdesc[i]->size==-1)? o->n_hor+1: paramdesc[i]->size;
        if((mxParam= mxGetField(mxParams, 0, paramdesc[i]->name))==NULL) {
            mxFree(o->p);
            mexErrMsgIdAndTxt("MATLAB:dimagree", "Parameter name '%s' is not member of parameters struct.\n", paramdesc[i]->name);
        }
        m_= mxGetM(mxParam);
        n_= mxGetN(mxParam);
        if(mxIsSparse(mxParam) || !mxIsDouble(mxParam) || (m_!=1 && n_!=1) || (m_*n_!=si)) {
            mxFree(o->p);
            mexErrMsgIdAndTxt("MATLAB:dimagree", "Parameter name '%s' must be a vector length %d.\n", paramdesc[i]->name, si);
        }
        o->p[i]= mxGetPr(mxParam);
    }
    // printParams(o.p, o.n_hor);

    // outputs
    plhs[0]= mxCreateDoubleMatrix(1, 1, mxREAL);
    double *success= mxGetPr(plhs[0]);

    plhs[1]= mxCreateDoubleMatrix(n, N, mxREAL);
    double *x_new= mxGetPr(plhs[1]);

    plhs[2]= mxCreateDoubleMatrix(m, N-1, mxREAL);
    double *u_new= mxGetPr(plhs[2]);

    plhs[3]= mxCreateDoubleMatrix(1, 1, mxREAL);
    double *new_cost= mxGetPr(plhs[3]);

    // aux
    for(int i= 0; i<NUMBER_OF_THREADS+1; i++) {
        o->trajectories[i].t= new trajEl_t[N-1];
        // memset(o.trajectories[i].t, 1, sizeof(trajEl_t)*(N-1));
    }
    o->multipliers.t= new multipliersEl_t[N];
//     mexPrintf("sizeof trajEl_t: %d\n", sizeof(trajEl_t));

    // inputs
    Map<VectorX> mapX0(mxGetPr(prhs[0]));
    o->x0= mapX0;
    // std::cout << "x0= " << o->x0 << std::endl;
    
    mexPrintf("Set const vars\n");
    if(!init_opt(o)) {
        success[0]= 0;
        new_cost[0]= o->cost;
    } else {
        mexPrintf("Initializing trajectory\n");
        double *u_nom= mxGetPr(prhs[1]);
        for(int k= 0; k<N-1; k++, u_nom+= N_U) {
            Map<VectorU> mapU(u_nom);
            o->nominal->t[k].u= mapU;
        }
        double cost;
        if(!forward_pass(o->candidates[0], o, 0.0, cost, 0)) {
            mexPrintf("makeCandidateNominal\n");
            makeCandidateNominal(o, 0);
            
            for(int k= 0; k<N-1; k++, x_new+= N_X) {
                Map<VectorX> m(x_new);
                m= o->nominal->t[k].x;
            }
            Map<VectorX> m(x_new); m= o->nominal->f.x;
                
            
            for(int k= 0; k<N-1; k++, u_new+= N_U) {
                Map<VectorU> m(u_new);
                m= o->nominal->t[k].u;
            }
            success[0]= 0;
            new_cost[0]= cost;
        } else {
            o->cost= cost;
            mexPrintf("makeCandidateNominal\n");
            makeCandidateNominal(o, 0);
            if(calc_system_only) {
                mexPrintf("System calculation only\n");
                if(!calc_derivs(o)) {
                    mexPrintf("Calculating derivatives failed.\n");
                }
            } else {
                mexPrintf("Starting iLQG\n");
                begin = clock();
                success[0]= iterate(o);
                end = clock();
                mexPrintf("Time for iLQG: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
            }
            
            for(int k= 0; k<N-1; k++, x_new+= N_X) {
                Map<VectorX> m(x_new);
                m= o->nominal->t[k].x;
            }
            Map<VectorX> m(x_new); m= o->nominal->f.x;
                
            
            for(int k= 0; k<N-1; k++, u_new+= N_U) {
                Map<VectorU> m(u_new);
                m= o->nominal->t[k].u;
            }
                    
            new_cost[0]= o->cost;

            
            // TODO: add output of stat log
            if(nlhs>4) output_cx(plhs, 4, o->nominal, N-1);
            if(nlhs>5) output_cu(plhs, 5, o->nominal, N-1);
            if(nlhs>6) output_cxx(plhs, 6, o->nominal, N-1);
            if(nlhs>7) output_cuu(plhs, 7, o->nominal, N-1);
            if(nlhs>8) output_cxu(plhs, 8, o->nominal, N-1);
            if(nlhs>9) output_fx(plhs, 9, o->nominal, N-1);
            if(nlhs>10) output_fu(plhs, 10, o->nominal, N-1);
#if FULL_DDP
            if(nlhs>11) output_fxx(plhs, 11, o->nominal, N-1);
            if(nlhs>12) output_fuu(plhs, 12, o->nominal, N-1);
            if(nlhs>13) output_fxu(plhs, 13, o->nominal, N-1);
#endif
        }
    }
    
    
    mxFree(o->p);
    
    for(int i= 0; i<NUMBER_OF_THREADS+1; i++)
        delete[] o->trajectories[i].t;
    
    delete[] o->multipliers.t;
    
    delete o;
}

 
