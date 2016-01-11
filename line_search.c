// C implementation of line search from http://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization by Yuval Tassa
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

#include "line_search.h"
#include "matMult.h"
#include "printMat.h"
#include "iLQG_problem.h"
 
#ifndef DEBUG_FORWARDPASS
#define DEBUG_FORWARDPASS 1
#else
    #if PREFIX1(DEBUG_FORWARDPASS)==1
    #define DEBUG_FORWARDPASS 1
    #endif
#endif


#define TRACE(x) do { if (DEBUG_FORWARDPASS) PRNT x; } while (0)
#define printVec_(x) do { if (DEBUG_FORWARDPASS) printVec x; } while (0)
#define printMat_(x) do { if (DEBUG_FORWARDPASS) printMat x; } while (0)
   

int line_search(tOptSet *o, int iter) {
    double expected, z, alpha, dcost, cnew;
    int i, success;
    
    for(i= 0; i < o->n_alpha; i++) {
        alpha= o->alpha[i];
        
        success= forward_pass(o->x0, o->x_nom, o->u_nom, o->x_new, o->u_new, &cnew, o->l, o->L, alpha, o->n_hor, o->p);
        if(success) {
            dcost= o->cost - cnew;
            expected= -alpha*(o->dV[0] + alpha*o->dV[1]);
            if(expected > 0)
                z = dcost/expected;
            else {
                z= 0;
                TRACE(("non-positive expected reduction: should not occur\n"));
            }

            if(z > o->zMin)
                break;
            else
                success= 0;
        } else {
            if(o->debug_level>=2) {
                TRACE(("line search: %-3d: prediction or objective failed with inf or nan\n", i+1));
            }
        }
    }
    
    if(o->debug_level>=2) {
        if(!success) {
            TRACE(("max number of line searches reached\n"));
        } else {
//             TRACE(("iter: %-3d  alpha: %-9.6g cost: %-9.6g  reduction: %-9.3g  z: %-9.3g\n", iter, alpha, o->cost, dcost, z));
        }
    }
    
    if(o->log_linesearch!=NULL) o->log_linesearch[iter]= i+1;
    if(o->log_z!=NULL) o->log_z[iter]= z;
    if(o->log_cost!=NULL) o->log_cost[iter]= cnew;
    o->new_cost= cnew;
    o->dcost= dcost;
    o->expected= expected;

    return success;
}


int forward_pass(const double *x0, const double *x_nom, const double *u_nom, double *x_new, double *u_new, double *cnew, const double *l, const double *L, double alpha, int N, double **params) {
    int i, k, j;
    double dx;
    cnew[0]= 0.0;
    aux_t aux;

    for(i= 0; i<N_X; i++) x_new[i]= x0[i]; // ic

    if(!calcStaticAux(&aux, params))
        return 0;
    
    for(k= 0; k<N; k++) {
        for(j= 0; j<N_U; j++)
            u_new[MAT_IDX(j, k, N_U)]= u_nom[MAT_IDX(j, k, N_U)] + l[MAT_IDX(j, k, N_U)]*alpha;
        
        for(i= 0; i<N_X; i++) {
            dx= x_new[MAT_IDX(i, k, N_X)] - x_nom[MAT_IDX(i, k, N_X)];
            
            for(j= 0; j<N_U; j++) {
                u_new[MAT_IDX(j, k, N_U)]+= L[MAT_IDX3(j, i, k, N_U, N_X)]*dx;
            }
        }
        
        if(!calcXVariableAux(&x_new[MAT_IDX(0, k, N_X)], k, &aux, params))
            return 0;
        
        clampU(&x_new[MAT_IDX(0, k, N_X)], &u_new[MAT_IDX(0, k, N_U)], k, &aux, params, N);
        if(!calcXUVariableAux(&x_new[MAT_IDX(0, k, N_X)], &u_new[MAT_IDX(0, k, N_U)], k, &aux, params))
            return 0;
        
        
        if(!ddpf(&x_new[MAT_IDX(0, k+1, N_X)], &x_new[MAT_IDX(0, k, N_X)], &u_new[MAT_IDX(0, k, N_U)], k, &aux, params, N))
            return 0;
        
        cnew[0]+= ddpJ(&x_new[MAT_IDX(0, k, N_X)], &u_new[MAT_IDX(0, k, N_U)], k, &aux, params, N);        
        if(isnan(cnew[0]) || !finite(cnew[0]))
            return 0;
    }
    
    if(!calcXVariableAux(&x_new[MAT_IDX(0, N, N_X)], k, &aux, params))
        return 0;
        
    cnew[0]+= ddpJ(&x_new[MAT_IDX(0, N, N_X)], 0, N, &aux, params, N);
    if(isnan(cnew[0]) || !finite(cnew[0]))
        return 0;
        
    return 1;
}

