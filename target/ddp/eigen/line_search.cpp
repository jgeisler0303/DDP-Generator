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

#include "iLQG.hpp"
#include "line_search.h"
 
#ifndef DEBUG_FORWARDPASS
#define DEBUG_FORWARDPASS 0
#else
    #if PREFIX1(DEBUG_FORWARDPASS)==1
    #define DEBUG_FORWARDPASS 1
    #endif
#endif


#define TRACE(x) do { if (DEBUG_FORWARDPASS) PRNT x; } while (0)
#define printVec_(x) do { if (DEBUG_FORWARDPASS) printVec x; } while (0)
#define printMat_(x) do { if (DEBUG_FORWARDPASS) printMat x; } while (0)
   
// return number of line searches if successfull, -1 if max number of line searches reached, -2 if inf or nan encountered
int line_search(tOptSet *o, int iter) {
    double expected, z, alpha, dcost, cnew;
    int i, success, ret= 0;
    
    for(i= 0; i < o->n_alpha; i++) {
        alpha= o->alpha[i];
        
        success= forward_pass(o->candidates[0], o, alpha, cnew, 0);
        if(success) {
            dcost= o->cost - cnew;
            expected= -alpha*(o->dV[0] + alpha*o->dV[1]);
            if(expected > 0)
                z = dcost/expected;
            else {
                z= 0;
                TRACE(("non-positive expected reduction: should not occur (dV[0]= %g, dV[1]= %g)\n", o->dV[0], o->dV[1]));
            }

            if(z > o->zMin)
                break;
            else
                success= 0;
        } else {
            if(o->debug_level>=2) {
                TRACE(("line search: %-3d: prediction or objective failed with inf or nan\n", i+1));
                ret= -2;
                break;
            }
        }
    }
    
    if(i>=o->n_alpha) {
        if(o->debug_level>=2)
            TRACE(("max number of line searches reached\n"));
        ret= -1;
    } else if(ret!=-2) {
//             TRACE(("iter: %-3d  alpha: %-9.6g cost: %-9.6g  reduction: %-9.3g  z: %-9.3g\n", iter, alpha, o->cost, dcost, z));
        ret= i+1;
    }
    
    if(o->log_linesearch!=NULL) o->log_linesearch[iter]= ret;
    if(o->log_z!=NULL) o->log_z[iter]= z;
    if(o->log_cost!=NULL) o->log_cost[iter]= cnew;
    o->new_cost= cnew;
    o->dcost= dcost;
    o->expected= expected;

    return ret;
}
