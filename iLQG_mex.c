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
    clock_t begin, end;

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
    for(i= 0; i<NUMBER_OF_THREADS+1; i++)
        o.trajectories[i].t= mxMalloc(sizeof(trajEl_t)*(N-1));
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
            
            mexPrintf("Starting iLQG\n");
            begin = clock();
            success[0]= iLQG(&o);
            end = clock();
            mexPrintf("Time for iLQG: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
            for(k= 0; k<N-1; k++)
                for(i= 0; i<N_X; i++)
                    x_new[MAT_IDX(i, k, N_X)]= o.nominal->t[k].x[i];
            for(i= 0; i<N_X; i++)
                x_new[MAT_IDX(i, N-1, N_X)]= o.nominal->f.x[i];
                
            for(k= 0; k<N-1; k++)
                for(i= 0; i<N_U; i++)
                    u_new[MAT_IDX(i, k, N_U)]= o.nominal->t[k].u[i];
                    
            new_cost[0]= o.cost;
        }
    }
    
    mxFree(o.p);
    for(i= 0; i<NUMBER_OF_THREADS+1; i++)
        mxFree(o.trajectories[i].t);
}
