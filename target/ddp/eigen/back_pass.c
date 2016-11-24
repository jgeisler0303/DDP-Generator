// iLQG back pass c implementation of http://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization by Yuval Tassa
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
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

#include "back_pass.h"
#include "iLQG.h"
#include "boxQP.h"
EIGEN_NO_DEBUG // no range checking!
#ifndef DEBUG_BACKPASS
#define DEBUG_BACKPASS 0
#else
    #if PREFIX1(DEBUG_BACKPASS)==1
    #define DEBUG_BACKPASS 1
    #endif
#endif

#define TRACE(x) do { if (DEBUG_BACKPASS) PRNT x; } while (0)
#define printVec_(x) do { if (DEBUG_BACKPASS) printVec x; } while (0)
#define printMat_(x) do { if (DEBUG_BACKPASS) printMat x; } while (0)
   
int back_pass(tOptSet *o) {
    double d1, g_norm_i, g_norm_max, g_norm_sum;
    int k, i_, j_, k_, l_, i_free, j_free, k_free, m_free, qpRes;
    int N= o->n_hor;
    VectorU grad, grad_clamped, search, Qu;
    VectorU_int is_clamped;
    VectorX Vx, Qx
    MatrixXX Vxx, Qxx;
    MatrixUU invHfree, L, R, QuuF, Quu;
    MatrixXU Qxu_reg, Qxu;
    trajEl_t *t= o->nominal->t + N - 1;
    trajFin_t *f= &o->nominal->f;
    
    g_norm_sum= 0.0;

    o->dV[0]= 0.0;
    o->dV[1]= 0.0;
    
    Vx= f->cx;
    Vxx= f->cxx;

    for(k= N-1; k>=0; k--, t--) {
        Qu= t->cu + Vx * t->fu;
        Qx= t->cx + Vx * t->fx;
        Qxu= t->cxu + t->fx.transpose() * Vxx.selfadjointView<Upper>() * t->fu;
#if FULL_DDP
        addVecTensXU(Qxu, Vx, t);
#endif

        Quu.triangularView<Upper>()= t->cuu + t->fu.transpose() * Vxx.selfadjointView<Upper>() * t->fu;
#if FULL_DDP
        addVecTensUU(Quu, Vx, t);
#endif
        
        Qxx.triangularView<Upper>()= t->cxx + t->fx.transpose() * Vxx.selfadjointView<Upper>() * t->fx;
#if FULL_DDP
        addVecTensXX(Qxx, Vx, t);
#endif
        
//         TRACE(("regularization\n"));
        QuuF= Quu;
        Qxu_reg= Qxu;
        if(o->regType==2) {
//             TRACE(("type 2\n"));
            QuuF.triangularView<Upper>()+= o->lambda * (t->fu.transpose() * t->fu);
            Qxu_reg+= o->lambda * (t->fx.transpose() * t->fu);
        }
        if(o->regType==1) {
//             TRACE(("type 1\n"));
            QuuF.diagonal().array()+= o->lambda;          
        }

        // solve Quadratic Program
//         TRACE(("boxQP\n"));
        if(k==o->n_hor-1)
            t->l.setZero();
        else
            t->l= (t+1)->l;

        if((qpRes= boxQP(QuuF, Qu, t->lower, t->upper, t->l, R, L, grad, grad_clamped, search, is_clamped, &m_free, invHfree, N_U))<1) {
            TRACE(("@k= %d: qpRes= %d \n", k, qpRes));
            return 1;
        }
        
        
        for(i_= N_U-1, i_free= m_free; i_>=0; i_--) {
            if(!is_clamped[i_]) {
                for(j_= N_U-1, j_free= m_free; j_>=i_; j_--) {
                    if(!is_clamped[j_]) {
                        invHfree(i_, j_)= invHfree(i_free, j_free);
                        j_free--;
                    } else {
                        invHfree(i_, j_)= 0.0;
                    }
                }
                i_free--;
            } else {
                if(is_clamped[i_]==1) {
                    invHfree.row(i_).setZero();
                    invHfree(i_, i_)= t->lower_sign[i_];
                    
                    Qxu_reg.row(i_)= t->lower_hx.row(i_);
                } else {
                    invHfree.row(i_).setZero();
                    invHfree(i_, i_)= t->upper_sign[i_];
                    
                    Qxu_reg.row(i_)= t->upper_hx.row(i_);
                }
            }
        }
        
        for(i_= N_U-1, i_free= m_free; i_>=0; i_--) {
            if(!is_clamped[i_]) {
                for(j_= N_U-1, j_free= m_free; j_>=i_; j_--) {
                    if(is_clamped[j_]) {
                        invHfree(i_, j_)= -1.0 * (invHfree.selfadjointView<Upper>().row(i_) * QuuF.selfadjointView<Upper>().col(j_));
                        if(is_clamped[j_]==1)
                            invHfree(i_, j_)*= t->lower_sign[j_];
                        else
                            invHfree(i_, j_)*=  t->upper_sign[j_];
                    }
                }
            }
        }

        t->L= invHfree.selfadjointView<Upper>() * Qxu_reg;
        
        
        // dV          = dV + [k_i'*Qu  .5*k_i'*Quu*k_i];
        o->dV[0]= t->l.transpose() * Qu;
        o->dV[1]= 0.5 * t->l.transpose() * Quu.selfadjointView<Upper>() * t->l
//         TRACE(("dV= %g, %g, %g, %g, %g\n", o->dV[0], o->dV[1], t->l[0], Qu[0], Qu[1]));

        // Vx(:,i)     = Qx  + K_i'*Quu*k_i + K_i'*Qu  + Qux'*k_i;
        Vx= Qx + t->L.transpose()*Quu.selfadjointView<Upper>()*t->l + t->L.transpose()*Qu + Qxu*t->l;
    
        // Vxx(:,:,i)  = Qxx + K_i'*Quu*K_i + K_i'*Qux + Qux'*K_i;
        Vxx.triangularView<Upper>()= Qxx + t->L.transpose()*Quu.selfadjointView<Upper>()*t->L;
        Qxx= t->L.transpose() * Qux;
        Vxx.triangularView<Upper>()+= Qxx + Qxx.transpose();
        
        // g_norm= mean(max(abs(l) ./ (abs(u)+1),[],1));
        g_norm_max= 0.0;
        for(i_= 0; i_<N_U; i_++) {
            g_norm_i= fabs(t->l[i_]) / (fabs(t->u[i_])+1.0);
            if(g_norm_i>g_norm_max) g_norm_max= g_norm_i;
        }
        g_norm_sum+= g_norm_max;
    }
    
    o->g_norm= g_norm_sum/((double)(o->n_hor-1));
    
    return 0;
}
