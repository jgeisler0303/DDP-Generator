#include "mex.h"
#include "matrix.h"
#include "ddp_wrapper.h"

#define freeme(x) do { if(x!=NULL) mxFree(x); x= NULL; } while(0)

void cleanUp(tOptSet *o) {
    freeme(o->u_nom);
    freeme(o->x_nom);
    freeme(o->x_new);
    freeme(o->u_new);
    freeme(o->alpha);
    freeme(o->beta);
    freeme(o->p);

    freeme(o->p_opt.e);
    freeme(o->p_opt.J);
    freeme(o->p_opt.a0);
    freeme(o->p_opt.linesearch);
    freeme(o->p_opt.regu_sum);
    freeme(o->p_opt.regu_max);
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    int i, j, m, n, success, si;
    mxArray *mxParams, *mxParam;
    char s[256];
    int n_opt_res= 8;
    const char *opt_res_names[]= {"iterations", "linesearch", "success", "e", "J", "a0", "regu_max", "regu_sum"};

    tOptSet o= {0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0.0, {0, 0, 0.0, 0, NULL, 0, NULL, NULL, NULL, NULL, NULL}, NULL};

    /* Check for proper number of input and output arguments */
    if (nrhs<3 || nrhs>6) {
        mexErrMsgTxt("3 to 6 input arguments required.");
    }
    if(nlhs<2 || nlhs>4) {
        mexErrMsgTxt("2 to 4 output arguments required.");
    }

    m= mxGetM(prhs[0]);
    n= mxGetN(prhs[0]);
    if(mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) || (m!=1 && n!=1) || (m!=N_X && n!=N_X)) {
        snprintf(s, 255, "Input 1 must be a vector length %d for initial conditions.\n", N_X);
        mexErrMsgTxt(s);
    }
    o.x_ic= mxGetPr(prhs[0]);

    if(mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]) || mxGetNumberOfElements(prhs[1])!=1 || floor(mxGetScalar(prhs[1]))!=mxGetScalar(prhs[1]))
        mexErrMsgTxt("Input 2 must be a scalar integer for optimization horizont.\n");
    o.n_hor= mxGetScalar(prhs[1]);

    mxParams= prhs[2];
    if(!mxIsStruct(mxParams) || mxGetNumberOfFields(mxParams)!=n_params) {
        snprintf(s, 255, "Input 3 must be a struct of %d parameters.\n", n_params);
        mexErrMsgTxt(s);
    }

    o.p= mxMalloc(n_params*sizeof(double *));
    for(i=0; i<n_params; i++) {
        si= (paramdesc[i]->size==-1)? o.n_hor+1: paramdesc[i]->size;
        if((mxParam= mxGetField(mxParams, 0, paramdesc[i]->name))==NULL) {
            cleanUp(&o);
            snprintf(s, 255, "Parameter name '%s' is not member of parameters struct.\n", paramdesc[i]->name);
            mexErrMsgTxt(s);
        }
        m= mxGetM(mxParam);
        n= mxGetN(mxParam);
        if(mxIsSparse(mxParam) || !mxIsDouble(mxParam) || (m!=1 && n!=1) || (m!=si && n!=si)) {
            cleanUp(&o);
            snprintf(s, 255, "Parameter name '%s' must be a vector length %d.\n", paramdesc[i]->name, si);
            mexErrMsgTxt(s);
        }
        o.p[i]= mxGetPr(mxParam);
    }

    if(nrhs>3 && !mxIsEmpty(prhs[3])) {
        m= mxGetM(prhs[3]);
        n= mxGetN(prhs[3]);
        if(mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]) || (m!=N_U && n!=o.n_hor)) {
            cleanUp(&o);
            snprintf(s, 255, "Input 4 must be a matrix(%d x n_hor= %d) of initial control moves.\n", N_U, o.n_hor);
            mexErrMsgTxt(s);
        }
        o.u_nom= mxMalloc(N_U*o.n_hor*sizeof(double));
        for(i= 0; i<N_U; i++) for(j= 0; j<o.n_hor; j++) o.u_nom[i + N_U*j]= mxGetPr(prhs[3])[i + N_U*j];
    } else {
        o.u_nom= mxMalloc(N_U*o.n_hor*sizeof(double));
        for(i= 0; i<N_U; i++) for(j= 0; j<o.n_hor; j++) o.u_nom[i + N_U*j]= 0.0;
    }

    o.p_opt.max_linesearch= 50;
    o.p_opt.max_iter= 100;
    o.p_opt.eta= 1e-6;
    o.p_opt.c= 0.5;
    if(nrhs>4 && !mxIsEmpty(prhs[4])) {
        if(!mxIsStruct(prhs[4])) {
            cleanUp(&o);
            mexErrMsgTxt("Input 5 must be a struct of optimization parameters.\n");
        }
        if((mxParam= mxGetField(prhs[4], 0, "max_linesearch"))!=NULL) {
            if(!mxIsDouble(mxParam) || mxGetNumberOfElements(mxParam)!=1 || floor(mxGetScalar(mxParam))!=mxGetScalar(mxParam)) {
                cleanUp(&o);
                mexErrMsgTxt("Optimization parameter 'max_linesearch' must be scalar integer.\n");
            }
            o.p_opt.max_linesearch= (int)mxGetScalar(mxParam);
        }
        if((mxParam= mxGetField(prhs[4], 0, "max_iter"))!=NULL) {
            if(!mxIsDouble(mxParam) || mxGetNumberOfElements(mxParam)!=1 || floor(mxGetScalar(mxParam))!=mxGetScalar(mxParam)) {
                cleanUp(&o);
                mexErrMsgTxt("Optimization parameter 'max_iter' must be scalar integer.\n");
            }
            o.p_opt.max_iter= (int)mxGetScalar(mxParam);
        }
        if((mxParam= mxGetField(prhs[4], 0, "eta"))!=NULL) {
            if(!mxIsDouble(mxParam) || mxGetNumberOfElements(mxParam)!=1) {
                cleanUp(&o);
                mexErrMsgTxt("Optimization parameter 'eta' must be scalar.\n");
            }
            o.p_opt.eta= mxGetScalar(mxParam);
        }
        if((mxParam= mxGetField(prhs[4], 0, "c"))!=NULL) {
            if(!mxIsDouble(mxParam) || mxGetNumberOfElements(mxParam)!=1 || mxGetScalar(mxParam)<=0 || mxGetScalar(mxParam)>=1) {
                cleanUp(&o);
                mexErrMsgTxt("Optimization parameter 'c' must be scalar (0..1).\n");
            }
            o.p_opt.c= mxGetScalar(mxParam);
        }
    }

    if(nrhs>5 && !mxIsEmpty(prhs[5])) {
        if(mxIsSparse(prhs[5]) || !mxIsDouble(prhs[5]) || mxGetNumberOfElements(prhs[5])!=1 || floor(mxGetScalar(prhs[5]))!=mxGetScalar(prhs[5]))
            mexErrMsgTxt("Input 6 must be a scalar integer debug level.\n");
        o.debug_level= mxGetScalar(prhs[5]);
    }

    o.x_nom= mxMalloc(N_X*(o.n_hor+1)*sizeof(double));
    o.x_new= mxMalloc(N_X*(o.n_hor+1)*sizeof(double));
    o.u_new= mxMalloc(N_U*o.n_hor*sizeof(double));
    o.alpha= mxMalloc(N_U*o.n_hor*sizeof(double));
    o.beta= mxMalloc(N_U*N_X*o.n_hor*sizeof(double));
    o.p_opt.linesearch= mxMalloc(o.p_opt.max_iter*sizeof(int));
    o.p_opt.e= mxMalloc(o.p_opt.max_iter*sizeof(double));
    o.p_opt.J= mxMalloc(o.p_opt.max_iter*sizeof(double));
    o.p_opt.a0= mxMalloc(o.p_opt.max_iter*sizeof(double));
    o.p_opt.regu_max= mxMalloc(o.p_opt.max_iter*sizeof(double));
    o.p_opt.regu_sum= mxMalloc(o.p_opt.max_iter*sizeof(double));

    success= iterateDDP(&o);

    plhs[0]= mxCreateDoubleMatrix(N_X, o.n_hor+1, mxREAL);
    for(i= 0; i<N_X; i++) for(j= 0; j<(o.n_hor+1); j++) mxGetPr(plhs[0])[i + N_X*j]= o.x_nom[i + N_X*j];

    plhs[1]= mxCreateDoubleMatrix(N_U, o.n_hor, mxREAL);
    for(i= 0; i<N_U; i++) for(j= 0; j<o.n_hor; j++) mxGetPr(plhs[1])[i + N_U*j]= o.u_nom[i + N_U*j];

    if(nlhs>2) {
        plhs[2]= mxCreateStructMatrix(1, 1, n_opt_res, opt_res_names);
        mxParam= mxCreateDoubleMatrix(1, 1, mxREAL); mxGetPr(mxParam)[0]= o.p_opt.iterations; mxSetField(plhs[2], 0, opt_res_names[0], mxParam);
        mxParam= mxCreateDoubleMatrix(1, o.p_opt.iterations, mxREAL); for(i= 0; i<o.p_opt.iterations; i++) mxGetPr(mxParam)[i]= o.p_opt.linesearch[i]; mxSetField(plhs[2], 0, opt_res_names[1], mxParam);
        mxParam= mxCreateDoubleMatrix(1, 1, mxREAL); mxGetPr(mxParam)[0]= o.p_opt.success; mxSetField(plhs[2], 0, opt_res_names[2], mxParam);
        mxParam= mxCreateDoubleMatrix(1, o.p_opt.iterations, mxREAL); for(i= 0; i<o.p_opt.iterations; i++) mxGetPr(mxParam)[i]= o.p_opt.e[i]; mxSetField(plhs[2], 0, opt_res_names[3], mxParam);
        mxParam= mxCreateDoubleMatrix(1, o.p_opt.iterations, mxREAL); for(i= 0; i<o.p_opt.iterations; i++) mxGetPr(mxParam)[i]= o.p_opt.J[i]; mxSetField(plhs[2], 0, opt_res_names[4], mxParam);
        mxParam= mxCreateDoubleMatrix(1, o.p_opt.iterations, mxREAL); for(i= 0; i<o.p_opt.iterations; i++) mxGetPr(mxParam)[i]= o.p_opt.a0[i]; mxSetField(plhs[2], 0, opt_res_names[5], mxParam);
        mxParam= mxCreateDoubleMatrix(1, o.p_opt.iterations, mxREAL); for(i= 0; i<o.p_opt.iterations; i++) mxGetPr(mxParam)[i]= o.p_opt.regu_max[i]; mxSetField(plhs[2], 0, opt_res_names[6], mxParam);
        mxParam= mxCreateDoubleMatrix(1, o.p_opt.iterations, mxREAL); for(i= 0; i<o.p_opt.iterations; i++) mxGetPr(mxParam)[i]= o.p_opt.regu_sum[i]; mxSetField(plhs[2], 0, opt_res_names[7], mxParam);
    }

    if(nlhs>3) {
        plhs[3]= mxCreateDoubleMatrix(1, 1, mxREAL);
        mxGetPr(plhs[3])[0]= o.V0;
    }
        
    cleanUp(&o);
    if(!success)
         mexWarnMsgTxt("The optimization was not successful\n");
}
