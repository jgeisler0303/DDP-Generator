#include <stdio.h>
#include <string.h>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <fenv.h>

#include "ddp.h"
#include "line_search.h"
#include "back_pass.h"
#include "printMat.h"

using namespace std;

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
    o->n_alpha= sizeof(default_alpha)/sizeof(default_alpha[0]);
    for(int i= 0; i<o->n_alpha; i++)
        o->alpha[i]= default_alpha[i];
    o->tolFun= 1e-7;
    o->tolConstraint= 1e-5;
    o->tolGrad= 1e-5;
    o->max_iter= 20;
    o->lambdaInit= 1;
    o->lambdaFactor= 1.6;
    o->lambdaMax= 1e10;
    o->lambdaMin= 1e-6;
    o->regType= 1;
    o->zMin= 0.0;
    o->w_pen_init_l= 1.0;
    o->w_pen_init_f= 1.0;
    o->w_pen_max_l= INF;
    o->w_pen_max_f= INF;
    o->w_pen_fact= 4.0; // 4...10 Bertsekas p. 123
    o->h_fd= 7.6294e-06;
    o->log= NULL;
    o->log_line= NULL;
    o->iterations= 0;
    o->contractGradInit= 0.0;   // auto init
    o->contractConstrInit= 0.0; // auto init
    o->contractGradFactor= 0.7;
    o->contractConstrFactor= 0.5; // suggested by "ON AUGMENTED LAGRANGIAN METHODS WITH GENERAL LOWER-LEVEL CONSTRAINTS", https://doi.org/10.1137/060654797
    o->lambdaFactorUpdateP= 1.0;
    o->lambdaFactorUpdateM= 1.0;
}

const char setOptParamErr_not_scalar[]= "parameter must be scalar";
const char setOptParamErr_alpha_range[]= "all alpha must be in the range [1.0..0.0)";
#define str(s) #s
#define stringify(s) str(s)
#define MAX_ALPHA_STR stringify(MAX_ALPHA)
const char setOptParamErr_alpha_n[]= "currently no more than " MAX_ALPHA_STR "values can be set for alpha";
const char setOptParamErr_alpha_monotonic[]= "all alpha must be monotonically decreasing";
const char setOptParamErr_not_pos[]= "parameter must be positive";
const char setOptParamErr_lt_one[]= "parameter must be > 1";
const char setOptParamErr_gt_one[]= "parameter must be < 1";
const char setOptParamErr_range_one_two[]= "parameter must be in range [1..2]";
const char setOptParamErr_range_zero_one[]= "parameter must be in range [0..1)";
const char setOptParamErr_no_such_parameter[]= "no such parameter";

const char *setOptParam(tOptSet *o, const char *name, const double *value, const int n) {
    if(strcmp(name, "alpha")==0) {
        if(n>MAX_ALPHA)
            return setOptParamErr_alpha_n;
        for(int i= 0; i<n; i++) {
            if(value[i]<0.0 || value[i]>1.0)
                return setOptParamErr_alpha_range;
            if(i>0 && value[i]>=value[i-1])
                return setOptParamErr_alpha_monotonic;
        }
        for(int i= 0; i<n; i++)
            o->alpha[i]= value[i];
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
    } else if(strcmp(name, "w_pen_fact")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<1.0)
            return setOptParamErr_lt_one;
        o->w_pen_fact= value[0];
    } else if(strcmp(name, "h_fd")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->h_fd= value[0];
    }  else if(strcmp(name, "lambdaFactorUpdateP")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<1.0)
            return setOptParamErr_lt_one;
        o->lambdaFactorUpdateP= value[0];
    }  else if(strcmp(name, "lambdaFactorUpdateM")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<1.0)
            return setOptParamErr_lt_one;
        o->lambdaFactorUpdateM= value[0];
    } else if(strcmp(name, "contractGradFactor")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0 || value[0]>=1.0)
            return setOptParamErr_range_zero_one;
        o->contractGradFactor= value[0];
    } else if(strcmp(name, "contractConstrFactor")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0 || value[0]>=1.0)
            return setOptParamErr_range_zero_one;
        o->contractConstrFactor= value[0];
    } else if(strcmp(name, "contractGradInit")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->contractGradInit= value[0];
    } else if(strcmp(name, "contractConstrInit")==0) {
        if(n!=1)
            return setOptParamErr_not_scalar;
        if(value[0]<0.0)
            return setOptParamErr_not_pos;
        o->contractConstrInit= value[0];
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
//    printf("new deriv= %d, n= %d, g_norm= %8.3g, lambda= %3.1f", l->new_deriv, l->n_back_pass, l->g_norm, log10(l->lambda));
    printf("g_norm= %8.3g, lambda= %3.1f", l->g_norm, log10(l->lambda));
}

void printLineSearchInfo(tLogLine *l) {
    printf("searches= %2d, z= %9.3g", l->n_line_searches, l->z);    
}

void printLogLine(int i, tLogLine *l) {
    printf("%3d: ", i);
    switch(l->res) {
        case 0:
            switch(l->multiplier_action) {
                case 0:
                    if(l->line_search_res>=0) {
                        printf("improvement: ");
                        printf("cost= %6.3f, constrs= %8.3g; ", l->cost, l->maxConstraint);
                        printLineSearchInfo(l);
                    } else {
                        printf("no improvement: ");
                        printf("cost= %6.3f, constrs= %8.3g; z= %8.3g/%-8.3g = %9.3g", l->cost, l->maxConstraint, l->dcost, l->expected_red, l->z);
                    }
                    printf("; ");
                    printBackPassInfo(l);
                    break;
                case 1:
                    printf("pen. update: ");
                    printf("cost= %6.3f, constrs= %8.3g, next grad= %8.3g, next constrs= %8.3g, g_norm= %8.3g, w_pen= %5.3g/%5.3g", l->cost, l->maxConstraint, l->contractGrad, l->contractConstr, l->g_norm, l->w_pen_l, l->w_pen_f);
                    break;
                case 2:
                    printf("mult.update: ");
                    printf("cost= %6.3f, constrs= %8.3g, next grad= %8.3g, next constrs= %8.3g, g_norm= %8.3g, w_pen= %5.3g/%5.3g", l->cost, l->maxConstraint, l->contractGrad, l->contractConstr, l->g_norm, l->w_pen_l, l->w_pen_f);
                    break;
            }
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
            printf("; z= %9.3g)", l->z);
            break;
        case -5:
            printf("ERROR max iterations reached cost= %6.3f, constrs= %6.3f (", l->cost, l->maxConstraint);
            printLineSearchInfo(l);
            printf("; ");
            printBackPassInfo(l);
            printf(")");
            break;
        case -6:
            printf("ERROR nan or inf in forward in multiplier update (");
            printBackPassInfo(l);
            printf(")");
            break;
        case -7:
            printf("ERROR nan or inf in forward pass after multiplier update (");
            printBackPassInfo(l);
            printf(")");
            break;
        case -8:
            printf("ERROR nan or inf in initial forward pass");
            break;
        case 1:
            printf("grad < tol:  ");
            printf("cost= %6.3f, constrs= %8.3g; ", l->cost, l->maxConstraint);
            printLineSearchInfo(l);
            printf("; ");
            printBackPassInfo(l);
            break;
        case 2:
            printf("dcost < tol: ");
            printf("cost= %6.3f, constrs= %8.3g; ", l->cost, l->maxConstraint);
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
        for(int i= 0; i < o->iterations; l++, i++)
            printLogLine(i+1, l);
    } else
        printf("No log recorded\n");
}


int iterate(tOptSet *o) {
    int iter;
    int fwdPass;
    int newDeriv;
    int res= 0;
    
//     fenv_t curr_excepts;
//     feholdexcept(&curr_excepts);
    
    o->lambda= o->lambdaInit;
    o->w_pen_l= o->w_pen_init_l;
    o->w_pen_f= o->w_pen_init_f;
    // default init: never update multipliers
    o->contractGrad= 0.0;
    o->contractConstr= 0.0;
    newDeriv= 1;
    
    if(!forward_pass(o->candidates[0], o, 0.0, o->cost, 0)) {
        o->iterations= 1;
        if(o->log_line) o->log_line->res= -8;
        return 0;
    } else
        makeCandidateNominal(o, 0);
    
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
                o->lambda= max(o->lambda * o->lambdaFactor, o->lambdaMin);
            } else {
                break;
            }
        }
        if(o->lambda >= o->lambdaMax) {
            res= -2;
            break;
        }

        if(o->log_line) o->log_line->g_norm= o->g_norm;
        
        calc_constraint_violation(o);
        if(o->log_line) o->log_line->maxConstraint= o->maxConstraint;
        
        // check for termination due to small gradient
        // TODO: make lambda _term a parameter, lambda should not influence the gradient criterium
        // if(o->g_norm < o->tolGrad && o->lambda < 1e-5) {
        if(o->g_norm < o->tolGrad && o->maxConstraint < o->tolConstraint) {
            if(o->log_line) o->log_line->cost= o->cost;
            res= 1;
            break;
        }
        
        if(iter==0) {
            // o->contractGradInit= o->g_norm * min(1.0/o->w_pen_l, o->contractGradMin);
            if(o->contractGradInit==0.0)
                o->contractGrad= max(o->g_norm * o->contractGradFactor, o->tolGrad);
            else
                o->contractGrad= min(o->contractGradInit, o->g_norm);
            
            // o->contractConstrInit= o->maxConstraint * pow(min(1.0/o->w_pen_l, o->contractConstrMin), 0.1);
            if(o->contractConstrInit==0.0)
                o->contractConstr= o->maxConstraint * o->contractConstrFactor;
            else
                o->contractConstr= min(o->maxConstraint, o->contractConstrInit);
        }
        
        // ====== STEP 3: line-search to find new control sequence, trajectory, cost
        fwdPass= line_search(o);
        if(o->log_line) o->log_line->line_search_res= fwdPass;

        if(fwdPass==-2) {
            res= -3;
            break;
        }
        
        // ====== STEP 4: accept (or not), draw graphics
        if(fwdPass>=0) {
            // accept changes
            makeCandidateNominal(o, fwdPass);
            
            // TODO: can this be an alternative criterium to the gradient?
            // if(o->dcost < o->tolFun) {
            //     res= 2;
            //    break;
            // }
            
            // check for sufficient convergence to update multipliers or weights
#ifdef HAS_AUGMENTED_LAGRANGIAN
            if(o->g_norm < o->contractGrad) {
                if(!update_multipliers(o)) {
                    res= -6;
                    break;
                }
                // o->contractConstr*= pow(min(1.0/o->w_pen_l, o->contractConstrMin), 0.9);
                // o->contractGrad= max(o->contractGrad*min(1.0/o->w_pen_l, o->contractGradMin), o->tolGrad);
                o->lambda= min(o->lambda * o->lambdaFactorUpdateM, o->lambdaInit);
                
                if(o->log_line) o->log_line->multiplier_action= 2;
                
                if(o->maxConstraint > o->contractConstr) {
                    bool penaltyUpdated= false;
                    if(o->w_pen_l < o->w_pen_max_l) {
                        o->w_pen_l= min(o->w_pen_max_l, o->w_pen_l*o->w_pen_fact);
                        penaltyUpdated= true;
                    }
                    if(o->w_pen_f < o->w_pen_max_f) {
                        o->w_pen_f= min(o->w_pen_max_f, o->w_pen_f*o->w_pen_fact);
                        penaltyUpdated= true;
                    }
                    if(penaltyUpdated) {
                        // according to Trust Region Augmented Lagrangian Methods for Sequential Response Surface Approximation and Optimization, doi:10.1115/1.2826677
                        // o->contractConstr= o->contractConstrInit * pow(min(1.0/o->w_pen_l, o->contractConstrMin), 0.1);
                        // o->contractGrad= max(o->contractGradInit * min(1.0/o->w_pen_l, o->contractGradMin), o->tolGrad);
                        o->lambda= min(o->lambda * o->lambdaFactorUpdateP, o->lambdaInit);
                        if(o->log_line) o->log_line->multiplier_action= 1;
                        if(o->log_line) o->log_line->w_pen_l= o->w_pen_l;
                        if(o->log_line) o->log_line->w_pen_f= o->w_pen_f;
                    }
                }
                
                if(!forward_pass(o->nominal, o, 0.0, o->cost, 1)) {
                    res= -7;
                    break;
                }
                o->contractGrad= max(o->contractGrad * o->contractGradFactor, o->tolGrad);
                o->contractConstr= max(o->maxConstraint * o->contractConstrFactor, o->tolConstraint);
                
                if(o->log_line) o->log_line->cost= o->cost;
                if(o->log_line) o->log_line->contractGrad= o->contractGrad;
                if(o->log_line) o->log_line->contractConstr= o->contractConstr;
            } else 
#endif
            {
                o->cost= o->new_cost;
                if(o->log_line) o->log_line->multiplier_action= 0;
            }
                
            // decrease lambda
            // maybe make these thresholds parameters
            if(o->n_ls==1 && o->last_z > 0.75 && o->last_z < 1.5) {
                o->lambda= o->lambda / o->lambdaFactor;
                if(o->lambda < o->lambdaMin) o->lambda= 0.0;
            }
                

            newDeriv= 1;
        } else { // no cost improvement
            // increase lambda
            o->lambda= max(o->lambda * o->lambdaFactor, o->lambdaMin);

            if(o->lambda > o->lambdaMax) {
                res= -4;
                break;
            }
        }
    }
    if(iter>=o->max_iter)
        res= -5;
    else
        iter++;
    
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
