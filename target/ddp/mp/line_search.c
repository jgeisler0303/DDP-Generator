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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "line_search.h"
#include "matMult.h"
#include "printMat.h"
 
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
    int successCandidate= -1, successIter= o->n_alpha;
    double bestZ, bestCost, lastExpected;
    
    #pragma omp parallel for ordered schedule(dynamic)
    for(int i= 0; i < o->n_alpha; i++) {
        double expected, z, alpha, dcost, cnew;
        int useCandidate= 0;
        int success;

        #pragma omp flush(successCandidate)
        if(successCandidate>-1) continue;

        #ifdef _OPENMP
        useCandidate= omp_get_thread_num();
        #endif
        
        alpha= o->alpha[i];
        
        success= forward_pass(o->candidates[useCandidate], o, alpha, &cnew, 0);
        if(success) {
            dcost= o->cost - cnew;
            expected= -alpha*(o->dV[0] + alpha*o->dV[1]);
            if(expected > 0)
                z = dcost/expected;
            else {
                z= 0;
                TRACE(("non-positive expected reduction: should not occur (dV[0]= %g, dV[1]= %g)\n", o->dV[0], o->dV[1]));
            }

            if(z < o->zMin)
                success= 0;
        } else {
            if(o->debug_level>=2) {
                TRACE(("line search: %-3d: prediction or objective failed with inf or nan\n", i+1));
            }
        }
        
        #pragma omp ordered
        {
            #pragma omp flush(successCandidate)
            if(successCandidate==-1) {
                if(success) {
                    #pragma omp atomic write
                    successCandidate= useCandidate;
                    successIter= i;
                }
                if(success || i==o->n_alpha-1) {
                    bestZ= z;
                    bestCost= cnew;
                    lastExpected= expected;
                }
            }
        }
    }
    
    if(o->debug_level>=2) {
        if(successCandidate==-1) {
            TRACE(("  line search: max number of line searches reached\n"));
        } else {
//             TRACE(("  line search: search iter: %2d  alpha: %-9.6g\n", successIter+1, o->alpha[successIter]));
        }
    }
    
    if(o->log_linesearch!=NULL) o->log_linesearch[iter]= successIter+1;
    if(o->log_z!=NULL) o->log_z[iter]= bestZ;
    if(o->log_cost!=NULL) o->log_cost[iter]= bestCost;
    o->new_cost= bestCost;
    o->dcost= o->cost-bestCost;
    o->expected= lastExpected;

    return successCandidate;
}
