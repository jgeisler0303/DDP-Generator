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

#include "iLQG.hpp"
#include "line_search.h"
#include "back_pass.h"
#include "printMat.h"

double default_alpha[]= {1.0, 0.3727594, 0.1389495, 0.0517947, 0.0193070, 0.0071969, 0.0026827, 0.0010000};

void printParams(double **p, int k) {
    for(int i=0; i<n_params; i++) {
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
    o->w_pen_init_l= 1.0;
    o->w_pen_init_f= 1.0;
    o->w_pen_max_l= INF;
    o->w_pen_max_f= INF;
    o->w_pen_fact1= 4.0; // 4...10 Bertsekas p. 123
    o->w_pen_fact2= 1.0;
    o->h_fd= 7.6294e-06;
    o->log= NULL;
    o->log_line= NULL;
    o->iterations= 0;
}

const char setOptParamErr_not_scalar[]= "parameter must be scalar";
const char setOptParamErr_alpha_range[]= "all alpha must be in the range [1.0..0.0)";
const char setOptParamErr_alpha_monotonic[]= "all alpha must be monotonically decreasing";
const char setOptParamErr_not_pos[]= "parameter must be positive";
const char setOptParamErr_lt_one[]= "parameter must be > 1";
const char setOptParamErr_gt_one[]= "parameter must be < 1";
const char setOptParamErr_range_one_two[]= "parameter must be in range [1..2]";
const char setOptParamErr_range_zero_one[]= "parameter must be in range [0..1)";
const char setOptParamErr_no_such_parameter[]= "no such parameter";

const char *setOptParam(tOptSet *o, const char *name, const double *value, const int n) {
    if(strcmp(name, "alpha")==0) {
        for(int i= 0; i<n; i++) {
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
    } else if(strcmp(name, "h_fd")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->h_fd= value[0];
    } else {
        return setOptParamErr_no_such_parameter;
    }

    return NULL;
}

const char *qpErrorStr(int e) {
    static const char llt_error[]= "llt error";
    static const char no_descent[]= "no descent";
    static const char max_iterations[]= "max iterations";
    static const char max_ls[]= "max line search";
    static const char no_bounds[]= "no bounds";
    static const char dcost_le_tol[]= "dcost<tol";
    static const char grad_le_tol[]= "grad<tol";
    static const char all_clamped[]= "all clamped";
    static const char unknown[]= "unknown";
    
    if(e==-1) return llt_error;
    if(e==-2) return no_descent;
    if(e==1) return max_iterations;
    if(e==2) return max_ls;
    if(e==3) return no_bounds;
    if(e==4) return dcost_le_tol;
    if(e==5) return grad_le_tol;
    if(e==6) return all_clamped;
    return unknown;
}

void printBackPassInfo(tLogLine *l) {
    printf("new deriv= %d, n back passes= %d, g_norm= %9.3f, lambda= %3.1f, w_pen= %9.3f / %9.3f", l->new_deriv, l->n_back_pass, l->g_norm, log10(l->lambda), l->w_pen_l, l->w_pen_f);
}

void printLineSearchInfo(tLogLine *l) {
    printf("searches= %2d (#neg= %2d), z= %12.3g/%12.3g = %6.3f", l->n_line_searches, l->neg_exp_red, l->dcost, l->expected_red, l->z);    
}

void printLogLine(int i, tLogLine *l) {
    printf("%3d: ", i);
    switch(l->res) {
        case 0:
            if(l->line_search_res>0) {
                printf("improvement: ");
                printf("cost= %12.6g; ", l->cost);
                printLineSearchInfo(l);
            } else {
                printf("no improvement: ");
                printf("z= %12.3g/%12.3g (#neg= %2d))", l->dcost, l->expected_red, l->neg_exp_red);
            }
            printf("; ");
            printBackPassInfo(l);
            break;
        case -1:
            printf("ERROR nan or inf in derivatives at k= %d", l->derivs_fail);
            break;
        case -2:
            printf("ERROR max lambda reached after %d back passes at k= %d with qp error \"%s\" (lambda= %f, w_pen_l= %f, w_pen_f= %f)", l->n_back_pass, l->back_pass_failed, qpErrorStr(l->qp_res), log10(l->lambda), l->w_pen_l, l->w_pen_f);
            break;
        case -3:
            printf("ERROR nan or inf in forward pass at k= %d after %d searches (", l->forward_pass_fail, l->n_line_searches);
            printBackPassInfo(l);
            printf(")");
            break;
        case -4:
            printf("ERROR max lambda reached after line search (");
            printBackPassInfo(l);
            printf("; z= %12.3g/%12.3g (#neg= %2d))", l->dcost, l->expected_red, l->neg_exp_red);
            break;
        case -5:
            printf("ERROR max iterations reached (");
            printLineSearchInfo(l);
            printf("; ");
            printBackPassInfo(l);
            printf(")");
            break;
        case 1:
            printf("grad < tol:  ");
            printf("cost= %12.6g; ", l->cost);
            printLineSearchInfo(l);
            printf("; ");
            printBackPassInfo(l);
            break;
        case 2:
            printf("dcost < tol: ");
            printf("cost= %12.6g; ", l->cost);
            printLineSearchInfo(l);
            printf("; ");
            printBackPassInfo(l);
            break;
        default:
            printf("ERROR unknown result code %d", l->res);
            break;
    }
    printf("\n");
}

void printLog(tOptSet *o) {
    if(o->log) {
        tLogLine *l= o->log;
        for(int i= 0; i <= o->iterations; l++, i++)
            printLogLine(i, l);
    } else
        printf("No log recorded\n");
}


int iterate(tOptSet *o) {
    int iter;
    int fwdPass;
    int newDeriv;
    double dlambda= o->dlambdaInit;
    int res= 0;
    
    o->lambda= o->lambdaInit;
    o->w_pen_l= o->w_pen_init_l;
    o->w_pen_f= o->w_pen_init_f;
    newDeriv= 1;
    
    update_multipliers(o, 1);
    
    for(iter= 0; iter < o->max_iter; iter++) {
        if(o->log) o->log_line= o->log+iter;

        // ====== STEP 1: differentiate dynamics and cost along new trajectory: integrated in back_pass
        if(newDeriv) {
            if(o->log_line) o->log_line->new_deriv= 1;

            if(!calc_derivs(o)) {
                res= -1;
                break;
            }
            
            newDeriv= 0;
        }
            
        // ====== STEP 2: backward pass, compute optimal control law and cost-to-go
        if(o->log_line) o->log_line->w_pen_l= o->w_pen_l;
        if(o->log_line) o->log_line->w_pen_f= o->w_pen_f;
        while(o->lambda < o->lambdaMax) {
            if(o->log_line) o->log_line->n_back_pass++;
            if(o->log_line) o->log_line->lambda= o->lambda;
            
            if(back_pass(o)) {
                // this doesn't make sense: if dlambda==1/o->lambdaFactor then lambda will not change for one pass
                // dlambda= max(dlambda * o->lambdaFactor, o->lambdaFactor);
                dlambda= o->lambdaFactor;
                o->lambda= max(o->lambda * dlambda, o->lambdaMin);
            } else {
                break;
            }
        }
        if(o->lambda >= o->lambdaMax) {
            res= -2;
            break;
        }

        if(o->log_line) o->log_line->g_norm= o->g_norm;
        
        // check for termination due to small gradient
        // TODO: add constraint tolerance check
        // TODO: make lambda _term a parameter
        if(o->g_norm < o->tolGrad && o->lambda < 1e-5) {
            res= 1;
            break;
        }
    
        // ====== STEP 3: line-search to find new control sequence, trajectory, cost
        fwdPass= line_search(o);
        if(o->log_line) o->log_line->line_search_res= fwdPass;

        if(fwdPass==-2) {
            res= -3;
            break;
        }
        
        // ====== STEP 4: accept (or not), draw graphics
        if(fwdPass>0) {
            // decrease lambda
            dlambda= min(dlambda / o->lambdaFactor, 1.0/o->lambdaFactor);
            o->lambda= o->lambda * dlambda * (o->lambda > o->lambdaMin);

            
            // accept changes
            makeCandidateNominal(o, 0);

            o->cost= o->new_cost;
            newDeriv= 1;
            
            // terminate ?
            // TODO: add constraint tolerance check
            if(o->dcost < o->tolFun) {
                res= 2;
                break;
            }
            // adapt w_pen
            // TODO: add check for sufficient decrease of gradient
            update_multipliers(o, 0);
            forward_pass(o->nominal, o, 0.0, o->cost, 1);

        } else { // no cost improvement
            // increase lambda
            // dlambda= max(dlambda * o->lambdaFactor, o->lambdaFactor);
            dlambda= o->lambdaFactor;
            o->lambda= max(o->lambda * dlambda, o->lambdaMin);

            if(o->w_pen_fact2>1.0) {
                o->w_pen_l= min(o->w_pen_max_l, o->w_pen_l*o->w_pen_fact2);
                o->w_pen_f= min(o->w_pen_max_f, o->w_pen_f*o->w_pen_fact2);
                forward_pass(o->nominal, o, 0.0, o->cost, 1);
            }
            
            if(o->lambda > o->lambdaMax) {
                res= -4;
                break;
            }
        }
    }
    if(iter>=o->max_iter)
        res= -5;
    
    o->iterations= iter;
    
    if(o->log_line) o->log_line->res= res;
        
    if(res>0)
        return 1;
    
    return 0;
}

void makeCandidateNominal(tOptSet *o, int idx) {
    traj_t *temp;
    temp= o->nominal;
    o->nominal= o->candidates[idx];
    o->candidates[idx]= temp;
}
