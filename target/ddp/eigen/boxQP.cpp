// Box constrained quadratic optimizer c implementation of http://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization by Yuval Tassa
// Copyright (c) 2016 Jens Geisler
//
// BIBTeX:
// @INPROCEEDINGS{
// author={Tassa, Y. and Mansard, N. and Todorov, E.},
// booktitle={Robotics and Automation (ICRA), 2014 IEEE International Conference on},
// title={Control-Limited Differential Dynamic Programming},
// year={2014}, month={May}, doi={10.1109/ICRA.2014.6907001}}

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>

#include "boxQP.h"
#define DO_PREFIX1(VAL)  1 ## VAL
#define PREFIX1(VAL)     DO_PREFIX1(VAL)


#ifndef DEBUG_BOXQP
#define DEBUG_BOXQP 0
#else
    #if PREFIX1(DEBUG_BOXQP)==1
    #define DEBUG_BOXQP 1
    #endif
#endif

#define TRACE(x) do { if (DEBUG_BOXQP) PRNT x; } while (0)


int boxQP(const Ref<const MatrixUU> &H, const Ref<const VectorU> &g, const Ref<const VectorU> &lower, const Ref<const VectorU> &upper, Ref<VectorU> x, Ref<VectorU_int> is_clamped, int &n_free, LLT<MatrixUU_dyn, Upper> &llt) {
// int boxQP(MatrixUU &H, const Ref<const VectorU> &g, const Ref<const VectorU> &lower, const Ref<const VectorU> &upper, Ref<VectorU> x, Ref<VectorU_int> is_clamped, int &n_free, Ref<MatrixUU_dyn> lltHfree) {
    double value, oldvalue= 0.0;
    VectorU grad, search, xc;
    MatrixUU_dyn Hfree;
    VectorU_dyn grad_clamped, search_free;
    MatrixUU L;
    double gnorm= 0.0;
    int nfactor= 0;
    double trace= 0.0;
    double sdotg;
    int all_clamped, clamps_changed, was_clamped;
    double step, vc;
    int nstep;
    
    // TODO make settable parameters
    const int maxIter           = 100;    // maximum number of iterations
    const double minGrad        = 1e-8;   // minimum norm of non-fixed gradient
    const double minRelImprove  = 1e-8;   // minimum relative improvement
    const double stepDec        = 0.6;    // factor for decreasing stepsize
    const double minStep        = 1e-22;  // minimal stepsize for linesearch
    const double Armijo         = 0.1;    // Armijo parameter (fraction of linear improvement required)
    
    for(int i= 0; i<nQP; i++) {
        // clamp to limits
        if(x(i)>upper(i)) x(i)= upper(i);
        if(x(i)<lower(i)) x(i)= lower(i);
    }
    is_clamped= VectorU_int::Zero();

    // initial objective value
    grad.noalias()= H.selfadjointView<Upper>()*x;
    value= x.dot(0.5*grad+g);
    
    for(int iter= 0; iter<maxIter; iter++) {
        oldvalue= value;
        
        grad= g;
        grad.noalias()+= H.selfadjointView<Upper>()*x;

        all_clamped= 1;
        clamps_changed= 0;
        n_free= 0;
        gnorm= 0.0;
        for(int i= 0; i<nQP; i++) {
            was_clamped= is_clamped(i);
            if(x(i)<=lower(i) && grad(i)>0)
                is_clamped(i)= 1;
            else if(x(i)>=upper(i) && grad(i)<0)
                is_clamped(i)= 2;
            else {
                is_clamped(i)= 0;
                all_clamped= 0;
                gnorm+= grad(i)*grad(i);
                n_free++;
            }
            if((!was_clamped) != (!is_clamped(i)))
                clamps_changed= 1;
        }
        
        // TRACE(("gnorm= %g\n", sqrt(gnorm)));
        
        if(all_clamped)
            return 6;
        
        // [Hfree, indef]  = chol(H(free,free));
        Hfree.resize(n_free, n_free);
        grad_clamped.resize(n_free);
        if(iter==0 || clamps_changed) {
            for(int j= 0, j_free= 0; j<nQP; j++) { // cols
                if(!is_clamped(j)) {
                    for(int i= 0, i_free= 0; i<=j; i++) { // rows
                        if(!is_clamped(i)) {
                            Hfree(i_free, j_free)= H(i, j);
                            i_free++;
                        }
                    }
                    grad_clamped(j_free)= g(j);
                    j_free++;
                }
            }

            llt.compute(Hfree);

            if(llt.info()!=Eigen::Success)
                return -1;
            
            nfactor++;
            // std::cout << "lltHfree: " << llt.matrixLLT() << std::endl;
        }
        
        // gnorm= sqrt(gnorm);
        if(gnorm<minGrad*minGrad)
            return 5;
        

        // get search direction
        // grad_clamped   = g  + H*(x.*clamped);
        VectorU x_clamped= x;
        for(int j= 0; j<nQP; j++)
            if(!is_clamped(j))
                x_clamped(j)= 0.0;

        VectorU Hx_clamped= H.selfadjointView<Upper>()*x_clamped;
        
        for(int j= 0, j_free= 0; j<nQP; j++) {
            if(!is_clamped(j)) {
                grad_clamped(j_free)+= Hx_clamped(j);
                j_free++;
            }
        }

        // search(free)   = -Hfree\(Hfree'\grad_clamped(free)) - x(free);
        search_free= -llt.solve(grad_clamped);
        for(int i= 0, i_free= 0; i<nQP; i++) {
            if(!is_clamped(i)) {
                search(i)= search_free(i_free) - x(i);
                i_free++;
            } else
                search(i)= 0.0;
        }
        
        // std::cout << "search: " << search << std::endl;
        
        // check for descent direction
        // sdotg          = sum(search.*grad);
        sdotg= search.dot(grad);
        // TRACE(("sdotg= %g\n", sdotg));
        if(sdotg>=0.0)
            return -2;
        
        // armijo linesearch
        step= 1.0;
        nstep= 0;
        while(1) {
            // xc    = clamp(x+step*search);
            xc= x + step*search;
            for(int i= 0; i<nQP; i++) {
                if(xc(i)>upper(i)) xc(i)= upper(i);
                if(xc(i)<lower(i)) xc(i)= lower(i);
            }
            // vc    = xc'*g + 0.5*xc'*H*xc;
            grad.noalias()= H.selfadjointView<Upper>()*xc;
            vc= xc.dot(0.5*grad + g);
            
            if(((vc-oldvalue)/(step*sdotg)) >= Armijo)
                break;
            
            step= step*stepDec;
            if(step<minStep) {
                if(vc<oldvalue)
                    x= xc;
                
                return 2;
            }
            
            nstep++;
        }
        
        // accept candidate
        x= xc;
        value= vc;
        
        if((oldvalue-value) < minRelImprove*fabs(oldvalue))
            return 4;
    }
    
    return 1;
}
