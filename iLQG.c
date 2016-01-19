// C implementation of iLQG algorithm from http://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization by Yuval Tassa
// Copyright (c) 2016 Jens Geisler
//
// BIBTeX:
// @INPROCEEDINGS{
// author={Tassa, Y. and Mansard, N. and Todorov, E.},
// booktitle={Robotics and Automation (ICRA), 2014 IEEE International Conference on},
// title={Control-Limited Differential Dynamic Programming},
// year={2014}, month={May}, doi={10.1109/ICRA.2014.6907001}}

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mex.h"
#include "iLQG.h"
#include "line_search.h"
#include "back_pass.h"
#include "iLQG_problem.h"

#ifndef DEBUG_ILQG
#define DEBUG_ILQG 1
#else
    #if PREFIX1(DEBUG_ILQG)==1
    #define DEBUG_ILQG 1
    #endif
#endif

#define TRACE(x) do { if (DEBUG_ILQG) PRNT x; } while (0)


double default_alpha[]= {1.0, 0.3727594, 0.1389495, 0.0517947, 0.0193070, 0.0071969, 0.0026827, 0.0010000};


void printParams(double **p, int k) {
    int i, k_;
    for(i=0; i<n_params; i++) {
        if(paramdesc[i]->size==-1)
            PRNT("%s[k]= %g\n", paramdesc[i]->name, p[i][k]);
        else if(paramdesc[i]->size==1)
            PRNT("%s= %g\n", paramdesc[i]->name, p[i][0]);
        else
            printVec(p[i], paramdesc[i]->size, paramdesc[i]->name);
    }
}

void standard_parameters(tOptSet *o) {
    o->alpha= default_alpha;
    o->n_alpha= 8;
    o->tolFun= 1e-7;
    o->tolGrad= 1e-5;
    o->max_iter= 20;
    o->lambdaInit= 1;
    o->dlambdaInit= 1;
    o->lambdaFactor= 1.6;
    o->lambdaMax= 1e10;
    o->lambdaMin= 1e-6;
    o->regType= 1;
    o->zMin= 0;
    o->debug_level= 2;
}

char setOptParamErr_not_scalar[]= "parameter must be scalar";
char setOptParamErr_alpha_range[]= "all alpha must be in the range [1.0..0.0)";
char setOptParamErr_alpha_monotonic[]= "all alpha must be monotonically decreasing";
char setOptParamErr_not_pos[]= "parameter must be positive";
char setOptParamErr_lt_one[]= "parameter must be > 1";
char setOptParamErr_range_one_two[]= "parameter must be in range [1..2]";
char setOptParamErr_range_zero_one[]= "parameter must be in range [0..1)";
char setOptParamErr_debug_level_range[]= "parameter must be in range [0..6]";
char setOptParamErr_no_such_parameter[]= "no such parameter";

char *setOptParam(tOptSet *o, const char *name, const double *value, const int n) {
    int i;
    
    if(strcmp(name, "alpha")==0) {
        for(i= 0; i<n; i++) {
            if(value[i]<0.0 || value[i]>1.0)
                return setOptParamErr_alpha_range;
            if(i>0 && value[i]>=value[i-1])
                return setOptParamErr_alpha_monotonic;
        }
        o->alpha= value;
        o->n_alpha= n;
    } else if(strcmp(name, "tolFun")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<=0.0)
            return setOptParamErr_not_pos;
        o->tolFun= value[0];
    } else if(strcmp(name, "tolGrad")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<=0.0)
            return setOptParamErr_not_pos;
        o->tolGrad= value[0];
    } else if(strcmp(name, "max_iter")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->max_iter= value[0];
    } else if(strcmp(name, "lambdaInit")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->lambdaInit= value[0];
    } else if(strcmp(name, "dlambdaInit")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->dlambdaInit= value[0];
    } else if(strcmp(name, "lambdaFactor")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<1.0)
            return setOptParamErr_lt_one;
        o->lambdaFactor= value[0];
    } else if(strcmp(name, "lambdaMax")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->lambdaMax= value[0];
    } else if(strcmp(name, "lambdaMin")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->lambdaMin= value[0];
    } else if(strcmp(name, "regType")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<1.0 || value[0]>2.0)
            return setOptParamErr_range_one_two;
        o->regType= value[0];
    } else if(strcmp(name, "zMin")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0 || value[0]>=1.0)
            return setOptParamErr_range_zero_one;
        o->zMin= value[0];
    } else if(strcmp(name, "debug_level")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0 || value[0]>6.0)
            return setOptParamErr_debug_level_range;
        o->debug_level= value[0];
    } else {
        return setOptParamErr_no_such_parameter;
    }
    
    return NULL;
}

int initialize_iLQG(tOptSet *o) {
    int i, j;
    double *x= o->x_nom;
    double *u= o->u_nom;
    double **p= o->p;    
    double lower[N_U];
    double upper[N_U];
    aux_t aux;
    
    for(i= 0; i<N_X; i++) x[i]= o->x0[i];
    
    if(!calcStaticAux(&aux, p)) {
        if(o->debug_level>=1) PRNT("initializeDDP aux failed @step %d\n", -1);
        return 0;
    }
    
    o->cost= 0.0;
    
    for(i= 0; i<o->n_hor; i++, x+= N_X, u+= N_U) {
        if(!calcXVariableAux(x, i, &aux, p)) {
            if(o->debug_level>=1) PRNT("initializeDDP aux failed @step %d\n", i);
            return 0;
        }
    
        clampU(x, u, i, &aux, p, o->n_hor);

        if(!calcXUVariableAux(x, u, i, &aux, p)) {
            if(o->debug_level>=1) PRNT("initializeDDP aux failed @step %d\n", i);
            return 0;
        }
        
        if(!ddpf(x+N_X, x, u, i, &aux, p, o->n_hor)) {
            if(o->debug_level>=1) PRNT("initializeDDP failed @step %d\n", i);
            return 0;
        }
        
        o->cost+= ddpJ(x, u, i, &aux, p, o->n_hor);        
        if(isnan(o->cost) || !finite(o->cost))  {
            if(o->debug_level>=1) PRNT("Objective is inf @step %d\n", i);
            return 0;
        }
    }
    
    if(!calcXVariableAux(x, o->n_hor, &aux, p)) {
        if(o->debug_level>=1) PRNT("initializeDDP aux failed @step %d\n", i);
        return 0;
    }
    
    o->cost+= ddpJ(x, NULL, o->n_hor, &aux, p, o->n_hor);
    if(isnan(o->cost) || !finite(o->cost)) {
        if(o->debug_level>=1) PRNT("Objective is inf @step %d\n", o->n_hor);
        return 0;
    }
    
    return 1;
}
    
int iLQG(tOptSet *o) {
    int iter, diverge, backPassDone, fwdPassDone;
    double dlambda= o->dlambdaInit;
    double *temp;
    
    o->lambda= o->lambdaInit;
    
    for(iter= 0; iter < o->max_iter; iter++) {
        // ====== STEP 1: differentiate dynamics and cost along new trajectory: integrated in back_pass
        // ====== STEP 2: backward pass, compute optimal control law and cost-to-go
        
        backPassDone= 0;
        while(!backPassDone) {
            diverge= back_pass(o);
        
            if(diverge) {
                if(o->debug_level>=1)
                    TRACE(("Back pass failed at timestep %d.\n",diverge));

                dlambda= max(dlambda * o->lambdaFactor, o->lambdaFactor);
                o->lambda= max(o->lambda * dlambda, o->lambdaMin);
                if(o->lambda > o->lambdaMax)
                    break;
            } else
                backPassDone= 1;
        }
        
        // check for termination due to small gradient
        if(o->g_norm < o->tolGrad && o->lambda < 1e-5) {
            dlambda= min(dlambda / o->lambdaFactor, 1.0/o->lambdaFactor);
            o->lambda= o->lambda * dlambda * (o->lambda > o->lambdaMin);
            if(o->debug_level>=1)
                TRACE(("\nSUCCESS: gradient norm < tolGrad\n"));
            break;
        }
    
        // ====== STEP 3: line-search to find new control sequence, trajectory, cost
        if(backPassDone)
            fwdPassDone= line_search(o, iter);

        // ====== STEP 4: accept (or not), draw graphics
        if(fwdPassDone) {
            if(o->debug_level>=1)
                TRACE(("iter: %-3d  cost: %-9.6g  reduction: %-9.3g  gradient: %-9.3g  log10lam: %3.1f\n", iter+1, o->cost, o->dcost, o->g_norm, log10(o->lambda)));
            
            // decrease lambda
            dlambda= min(dlambda / o->lambdaFactor, 1.0/o->lambdaFactor);
            o->lambda= o->lambda * dlambda * (o->lambda > o->lambdaMin);
            
            // accept changes
            temp= o->x_nom;
            o->x_nom= o->x_new;
            o->x_new= temp;
            
            temp= o->u_nom;
            o->u_nom= o->u_new;
            o->u_new= temp;

            o->cost= o->new_cost;
            
            // terminate ?
            if(o->dcost < o->tolFun) {
                if(o->debug_level>=1)
                    TRACE(("\nSUCCESS: cost change < tolFun\n"));
            
                break;
            }
        } else { // no cost improvement
            // increase lambda
            dlambda= max(dlambda * o->lambdaFactor, o->lambdaFactor);
            o->lambda= max(o->lambda * dlambda, o->lambdaMin);
            
            // print status
            if(o->debug_level>=1)
                TRACE(("iter: %-3d  REJECTED    expected: %-11.3g    actual: %-11.3g    log10lam: %3.1f\n", iter+1, o->expected , o->dcost, log10(o->lambda)));
            
            // terminate ?
            if(o->lambda > o->lambdaMax) {
                if(o->debug_level>=1)
                    TRACE(("\nEXIT: lambda > lambdaMax\n"));
                break;
            }
        }
    }
    
    
    o->iterations= iter;
    
    if(iter>=o->max_iter) {
        if(o->debug_level>=1)
            TRACE(("\nEXIT: Maximum iterations reached.\n"));
        
        return 0;
    }
    return 1;
}
