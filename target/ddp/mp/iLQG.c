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
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

#include "printMat.h"
#include "iLQG.h"
#include "line_search.h"
#include "back_pass.h"

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
    int i;
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
    o->tolConstraint= 1e-7;
    o->tolGrad= 1e-5;
    o->max_iter= 20;
    o->lambdaInit= 1;
    o->dlambdaInit= 1;
    o->lambdaFactor= 1.6;
    o->lambdaMax= 1e10;
    o->lambdaMin= 1e-6;
    o->regType= 1;
    o->zMin= 0.0;
    o->debug_level= 2;
    o->w_pen_init_l= 1.0;
    o->w_pen_init_f= 1.0;
    o->w_pen_max_l= INF;
    o->w_pen_max_f= INF;
    o->w_pen_fact1= 4.0; // 4...10 Bertsekas p. 123
    o->w_pen_fact2= 1.0;
}

char setOptParamErr_not_scalar[]= "parameter must be scalar";
char setOptParamErr_alpha_range[]= "all alpha must be in the range [1.0..0.0)";
char setOptParamErr_alpha_monotonic[]= "all alpha must be monotonically decreasing";
char setOptParamErr_not_pos[]= "parameter must be positive";
char setOptParamErr_lt_one[]= "parameter must be > 1";
char setOptParamErr_gt_one[]= "parameter must be < 1";
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
    } else if(strcmp(name, "tolConstraint")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<=0.0)
            return setOptParamErr_not_pos;
        o->tolConstraint= value[0];
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
    } else if(strcmp(name, "w_pen_init_l")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->w_pen_init_l= value[0];
    } else if(strcmp(name, "w_pen_init_f")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->w_pen_init_f= value[0];
    } else if(strcmp(name, "w_pen_max_l")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->w_pen_max_l= value[0];
    } else if(strcmp(name, "w_pen_max_f")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->w_pen_max_f= value[0];
    } else if(strcmp(name, "w_pen_fact1")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<1.0)
            return setOptParamErr_lt_one;
        o->w_pen_fact1= value[0];
    } else if(strcmp(name, "w_pen_fact2")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<1.0)
            return setOptParamErr_lt_one;
        o->w_pen_fact2= value[0];
    } else {
        return setOptParamErr_no_such_parameter;
    }

    return NULL;
}



int iLQG(tOptSet *o) {
    int iter, diverge, backPassDone, fwdPassDone, successCandidate;
    int newDeriv;
    double dlambda= o->dlambdaInit;
    o->lambda= o->lambdaInit;
    o->w_pen_l= o->w_pen_init_l;
    o->w_pen_f= o->w_pen_init_f;
    newDeriv= 1;
    
    #ifdef _OPENMP
    if(omp_get_max_threads()>NUMBER_OF_THREADS)
        omp_set_num_threads(NUMBER_OF_THREADS);
    #endif

    
    update_multipliers(o, 1);
    
    for(iter= 0; iter < o->max_iter; iter++) {
        // ====== STEP 1: differentiate dynamics and cost along new trajectory: integrated in back_pass
        // moved to back_pass
            
        // ====== STEP 2: backward pass, compute optimal control law and cost-to-go
        backPassDone= 0;
//         TRACE(("Back pass:\n"));
        while(!backPassDone) {
            if(back_pass(o, newDeriv)) {
                // TODO check if backpass or derivs failed!
                if(o->debug_level>=1)
                    TRACE(("Back pass or derivatives failed.\n"));
                
                dlambda= max(dlambda * o->lambdaFactor, o->lambdaFactor);
                o->lambda= max(o->lambda * dlambda, o->lambdaMin);
                if(o->lambda > o->lambdaMax)
                    break;
            } else {
                backPassDone= 1;
//                 TRACE(("...done\n"));
            }
            newDeriv= 0;
        }
        
        // check for termination due to small gradient
        // TODO: add constraint tolerance check
        if(o->g_norm < o->tolGrad && o->lambda < 1e-5) {
            dlambda= min(dlambda / o->lambdaFactor, 1.0/o->lambdaFactor);
            o->lambda= o->lambda * dlambda * (o->lambda > o->lambdaMin);
            if(o->debug_level>=1)
                TRACE(("\nSUCCESS: gradient norm < tolGrad\n"));
            break;
        }
    
        // ====== STEP 3: line-search to find new control sequence, trajectory, cost
        if(backPassDone) {
            successCandidate= line_search(o, iter);
            fwdPassDone= successCandidate>-1;
        } else
            break;
        
        // ====== STEP 4: accept (or not), draw graphics
        if(fwdPassDone) {
            if(o->debug_level>=1)
                TRACE(("iter: %-3d  cost: %-9.6g  reduction: %-9.3g  gradient: %-9.3g  z: %-5.3g log10(lam): %3.1f w_pen_l: %-9.3g w_pen_f: %-9.3g\n", iter+1, o->cost, o->dcost, o->g_norm, o->dcost/o->expected, log10(o->lambda), o->w_pen_l, o->w_pen_f));
            
            // decrease lambda
            dlambda= min(dlambda / o->lambdaFactor, 1.0/o->lambdaFactor);
            o->lambda= o->lambda * dlambda * (o->lambda > o->lambdaMin);

            
            // accept changes
            makeCandidateNominal(o, successCandidate);

            o->cost= o->new_cost;
            newDeriv= 1;
            
            // terminate ?
            // TODO: add constraint tolerance check
            if(o->dcost < o->tolFun) {
                if(o->debug_level>=1)
                    TRACE(("\nSUCCESS: cost change < tolFun\n"));
            
                break;
            }
            // adapt w_pen
            // TODO: add check for sufficient decrease of gradient
            update_multipliers(o, 0);
            forward_pass(o->nominal, o, 0.0, &o->cost, 1);

        } else { // no cost improvement
            // increase lambda
            dlambda= max(dlambda * o->lambdaFactor, o->lambdaFactor);
            o->lambda= max(o->lambda * dlambda, o->lambdaMin);

            if(o->w_pen_fact2>1.0) {
                o->w_pen_l= min(o->w_pen_max_l, o->w_pen_l*o->w_pen_fact2);
                o->w_pen_f= min(o->w_pen_max_f, o->w_pen_f*o->w_pen_fact2);
                forward_pass(o->nominal, o, 0.0, &o->cost, 1);
            }
            
            // print status
            if(o->debug_level>=1)
                TRACE(("iter: %-3d  REJECTED    expected: %-11.3g    actual: %-11.3g    log10lam: %3.1f w_pen_l: %-9.3g w_pen_l: %-9.3g\n", iter+1, o->expected , o->dcost, log10(o->lambda), o->w_pen_l, o->w_pen_f));
            
            // terminate ?
            if(o->lambda > o->lambdaMax) {
                if(o->debug_level>=1)
                    TRACE(("\nEXIT: lambda > lambdaMax\n"));
                break;
            }
        }
    }
    
    
    o->iterations= iter;
    
    if(!backPassDone) {
        if(o->debug_level>=1)
            TRACE(("\nEXIT: no descent direction found.\n"));
        
        return 0;    
    } else if(iter>=o->max_iter) {
        if(o->debug_level>=1)
            TRACE(("\nEXIT: Maximum iterations reached.\n"));
        
        return 0;
    }
    return 1;
}

void makeCandidateNominal(tOptSet *o, int idx) {
    traj_t *temp;
    temp= o->nominal;
    o->nominal= o->candidates[idx];
    o->candidates[idx]= temp;
}
