// C implementation of line search from http://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization by Yuval Tassa
// Copyright (c) 2016 Jens Geisler
//
// BIBTeX:
// @INPROCEEDINGS{
// author={Tassa, Y. and Mansard, N. and Todorov, E.},
// booktitle={Robotics and Automation (ICRA), 2014 IEEE International Conference on},
// title={Control-Limited Differential Dynamic Programming},
// year={2014}, month={May}, doi={10.1109/ICRA.2014.6907001}}

#include <math.h>

#include "line_search.h"
#include "matMult.h"
 

   
// return number of line searches if successfull, -1 if max number of line searches reached, -2 if inf or nan encountered
int line_search(tOptSet *o) {
    double expected, z, alpha, dcost= INF, cnew= INF;
    int i, ret= -1;
    
    for(i= 0; i < o->n_alpha; i++) {
        if(o->log_line) o->log_line->n_line_searches= i+1;
        
        alpha= o->alpha[i];
        
        if(forward_pass(o->candidates[0], o, alpha, &cnew, 0)) {
            dcost= o->cost - cnew;
            expected= -alpha*(o->dV[0] + alpha*o->dV[1]);
            if(expected > 0)
                z = dcost/expected;
            else {
                z= 0;
                if(o->log_line) o->log_line->neg_exp_red++;
            }

            if(z > o->zMin) {
                ret= 1;
                break;
            }
        } else {
            ret= -2;
            break;
        }
    }
    
    if(ret==1 && o->log_line) {
        o->log_line->alpha= alpha;
        o->log_line->z= z;
        o->log_line->cost= cnew;
        o->log_line->dcost= dcost;
        o->log_line->expected_red= expected;
    }
    
    o->new_cost= cnew;
    o->dcost= dcost;

    return ret;
}
