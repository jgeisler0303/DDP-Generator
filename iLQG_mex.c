// MATLAB Mex function wrapper for iLQG algorithm
// Copyright (c) 2016 Jens Geisler


#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mex.h"

#include "iLQG.h"
#include "printMat.h"
#include "matMult.h"
#include "iLQG_problem.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // dims
    int N, n, m, m_, n_, si, i, k;
    // inputs
    tOptSet o= INIT_OPTSET;
    int dims[3];
    // outputs
    double *success, *new_cost, *x_nom, *u_nom, *x_new, *u_new, *l, *L;
    // aux
    const mxArray *mxParams, *mxOptParam, *mxParam;

    const char *err_msg, *fname;

    if(nrhs!=4) { mexErrMsgIdAndTxt("MATLAB:minrhs", "wrong number of arguments (expected: x0, u_nom, params, opt_params)"); return; }
    if(nlhs!=4) { mexErrMsgIdAndTxt("MATLAB:minlhs", "wrong number of return values (expected: success, x_new, u_new, new_cost)"); return; }

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
    o.x0= mxGetPr(prhs[0]);
    u_nom= mxGetPr(prhs[1]);
    o.n_hor= N-1;
    
    standard_parameters(&o);
    for(i= 0; i<mxGetNumberOfFields(mxOptParam); i++) {
        mxParam= mxGetFieldByNumber(mxOptParam, 0, i);
        fname= mxGetFieldNameByNumber(mxOptParam, i);
        err_msg= setOptParam(&o, fname, mxGetPr(mxParam), mxGetNumberOfElements(mxParam));
        if(err_msg) {
            mexErrMsgIdAndTxt("MATLAB:dimagree", "Error setting optimization parameter '%s': %s.\n", fname, err_msg);
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
    o.x_nom= mxMalloc(sizeof(double)*n*N);
    o.u_nom= mxMalloc(sizeof(double)*m*(N-1));
    o.x_new= mxMalloc(sizeof(double)*n*N);
    o.u_new= mxMalloc(sizeof(double)*m*(N-1));
    o.l= mxMalloc(sizeof(double)*m*(N-1));
    o.L= mxMalloc(sizeof(double)*m*n*(N-1));

    memcpy(o.u_nom, u_nom, sizeof(double)*m*(N-1));
    
    if(!initialize_iLQG(&o)) {
        success[0]= 0;
        new_cost[0]= o.cost;
    } else {
        success[0]= iLQG(&o);

        if(!success[0]) {
            memcpy(x_new, o.x_nom, sizeof(double)*n*N);
            memcpy(u_new, o.u_nom, sizeof(double)*m*(N-1));
            new_cost[0]= o.cost;
        } else {
            memcpy(x_new, o.x_new, sizeof(double)*n*N);
            memcpy(u_new, o.u_new, sizeof(double)*m*(N-1));
            new_cost[0]= o.new_cost;
        }
    }
    
    mxFree(o.p);
    mxFree(o.x_nom);
    mxFree(o.u_nom);
    mxFree(o.x_new);
    mxFree(o.u_new);
    mxFree(o.l);
    mxFree(o.L);
}
