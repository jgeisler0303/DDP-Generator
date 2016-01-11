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
// uncomment this for use in MATLAB
// #include "matrix.h"
#include "mex.h"

#include "back_pass.h"
#include "matMult.h"
#include "boxQP.h"
#include "printMat.h"

#include "iLQG_problem.h"
#include "iLQG_bp.h"


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
    int k, i_, j_, k_, i_free, k_free, m_free, qpRes;
    double *x, *u;
    double invHfree[sizeofQuu], L[sizeofQuu];
    double lower[N_U], upper[N_U];
    double R[sizeofQuu], grad[N_U], grad_clamped[N_U], search[N_U];
    int is_clamped[N_U];
    back_pass_t derivs;
    aux_t aux;
    
    double **p= o->p;
    
    g_norm_sum= 0.0;

    o->dV[0]= 0.0;
    o->dV[1]= 0.0;

    //Vx=, Vxx=
    k= o->n_hor;
    x= &(o->x_nom[N_X*k]);
    u= &(o->u_nom[N_X*(k-1)]);
    
    if(!calcStaticAux(&aux, p))
        return k+1;
    
    if(!calcXVariableAux(x, k, &aux, p))
        return k+1;
    
    if(!calcStaticDeriv(&aux, p))
        return k+1;
    
    if(!calcVariableDeriv(x, u, k, &aux, p))
        return k+1;
    
    if(!bp_derivsF(&derivs, &aux, x, u, p, k))
        return k+1;

    for(k= o->n_hor-1; k>=0; k--) {
        // Lx=, Lu=, Lxx=, Lxu=, Luu= is written to Q..
        x= &(o->x_nom[N_X*k]);
        u= &(o->u_nom[N_U*k]);
        
        if(!calcXVariableAux(x, k, &aux, p) || !calcXUVariableAux(x, u, k, &aux, p))
            return k+1;
        
        if(!calcVariableDeriv(x, u, k, &aux, p))
            return k+1;
        
        if(!bp_derivsL(&derivs, &aux, x, u, p, k))
            return k+1;
        
        if(o->regType==2) {
            memcpy(derivs.Qxu_reg, derivs.Qxu, sizeof(double)*sizeofQxu);
            memcpy(derivs.QuuF, derivs.Quu, sizeof(double)*sizeofQuu);
        }
        
        // Qu  = cu(:,i)      + fu(:,:,i)'*Vx(:,i+1);
        addMulVec(derivs.Qu, derivs.Vx, derivs.fu, N_X, N_U);

        // Qx  = cx(:,i)      + fx(:,:,i)'*Vx(:,i+1);
        addMulVec(derivs.Qx, derivs.Vx, derivs.fx, N_X, N_X);

        // Qux = cxu(:,:,i)'  + fu(:,:,i)'*Vxx(:,:,i+1)*fx(:,:,i);
        addMul2Tri(derivs.Qxu, derivs.Vxx, derivs.fx, N_X, N_X, derivs.fu, N_X, N_U, derivs.dummy);
        // fxuVx = vectens(Vx(:,i+1),fxu(:,:,:,i));
        // Qux   = Qux + fxuVx;
        for(j_= 0; j_<N_X*N_U; j_++) { // x, u
            d1= 0.0;
            for(i_= 0, k_= 0; i_<N_X; i_++, k_+= N_X*N_U) // f
                d1+= derivs.Vx[i_]*derivs.fxu[j_+k_];
            derivs.Qxu[j_]+= d1;
        }

        // Quu = cuu(:,:,i)   + fu(:,:,i)'*Vxx(:,:,i+1)*fu(:,:,i);
        addSquareTri(derivs.Quu, derivs.Vxx, derivs.fu, N_X, N_U, derivs.dummy);
        // fuuVx = vectens(Vx(:,i+1),fuu(:,:,:,i));
        // Quu   = Quu + fuuVx;
        for(j_= 0; j_<sizeofQuu; j_++) { // u, u
            d1= 0.0;
            for(i_= 0, k_= 0; i_<N_X; i_++, k_+= sizeofQuu) // f
                d1+= derivs.Vx[i_]*derivs.fuu[j_+k_];
            derivs.Quu[j_]+= d1;
        }                
        
        // Qxx = cxx(:,:,i)   + fx(:,:,i)'*Vxx(:,:,i+1)*fx(:,:,i);
        addSquareTri(derivs.Qxx, derivs.Vxx, derivs.fx, N_X, N_X, derivs.dummy);

        // Qxx = Qxx + vectens(Vx(:,i+1),fxx(:,:,:,i));
        for(j_= 0; j_<sizeofQxx; j_++) {// x, x
            d1= 0.0;
            for(i_= 0, k_= 0; i_<N_X; i_++, k_+= sizeofQxx) // f
                d1+= derivs.Vx[i_]*derivs.fxx[j_+k_];
            derivs.Qxx[j_]+= d1;
        }
        
        if(o->regType==2) {
            memcpy(derivs.Vxx_reg, derivs.Vxx, sizeof(double)*sizeofQxx);
            for(i_= 0; i_<N_X; i_++) derivs.Vxx_reg[UTRI_MAT_IDX(i_, i_)]+= o->lambda;

            addMul2Tri(derivs.Qxu_reg, derivs.Vxx_reg, derivs.fx, N_X, N_X, derivs.fu, N_X, N_U, derivs.dummy);
            for(j_= 0; j_<N_X*N_U; j_++) {// x, u
                d1= 0.0;
                for(i_= 0, k_= 0; i_<N_X; i_++, k_+=N_X*N_U) // f
                    d1+= derivs.Vx[i_]*derivs.fxu[j_+k_];
                derivs.Qxu_reg[j_]+= d1;
            }
            addSquareTri(derivs.QuuF, derivs.Vxx_reg, derivs.fu, N_X, N_U, derivs.dummy);
            
            for(j_= 0; j_<sizeofQuu; j_++) { // u, u
                d1= 0.0;
                for(i_= 0, k_= 0; i_<N_X; i_++, k_+= sizeofQuu) // f
                    d1+= derivs.Vx[i_]*derivs.fuu[j_+k_];
                derivs.QuuF[j_]+= d1;
            }
        } else {
            memcpy(derivs.QuuF, derivs.Quu, sizeof(double)*sizeofQuu);
            memcpy(derivs.Qxu_reg, derivs.Qxu, sizeof(double)*N_X*N_U);
        }

        if(o->regType==1) {
            for(i_= 0; i_<N_U; i_++) derivs.QuuF[UTRI_MAT_IDX(i_, i_)]+= o->lambda;          
        }

        // solve Quadratic Program
        for(i_= 0; i_<N_U; i_++) lower[i_]= -mxGetInf();
        clampU(x, lower, k, &aux, p, o->n_hor);
        for(i_= 0; i_<N_U; i_++) upper[i_]= mxGetInf();
        clampU(x, upper, k, &aux, p, o->n_hor);
        for(i_= 0; i_<N_U; i_++) {
            lower[i_]-= u[i_];
            upper[i_]-= u[i_]; 
        }
        
        if(k==o->n_hor-1)
            memset(&(o->l[MAT_IDX(0, k, N_U)]), 0, sizeof(double)*N_U);
        else
            memcpy(&(o->l[MAT_IDX(0, k, N_U)]), &(o->l[MAT_IDX(0, k+1, N_U)]), sizeof(double)*N_U);

        if((qpRes= boxQP(derivs.QuuF, derivs.Qu, lower, upper, &(o->l[MAT_IDX(0, k, N_U)]), R, L, grad, grad_clamped, search, is_clamped, &m_free, invHfree, N_U))<1) {
            return k+1;
        }
        TRACE(("qpRes= %d\n", qpRes));
        // Lfree        = -R\(R'\Qux_reg(free,:));
        // K_i(free,:)   = Lfree;
        memset(&o->L[MAT_IDX3(0, 0, k, N_U, N_X)], 0, sizeof(double)*N_U*N_X);
        for(i_= 0, i_free= 0; i_<N_U; i_++) {
            if(!is_clamped[i_]) {
                for(j_= 0; j_<N_X; j_++) {
                    for(k_= 0, k_free= 0; k_<N_U; k_++) {
                        if(!is_clamped[k_]) {
                          o->L[MAT_IDX3(i_, j_, k, N_U, N_X)]-= invHfree[SYMTRI_MAT_IDX(i_free, k_free)] * derivs.Qxu_reg[MAT_IDX(j_, k_, N_X)];
                          k_free++;
                        }
                    }
                }
                i_free++;
            }
        }

        // dV          = dV + [k_i'*Qu  .5*k_i'*Quu*k_i];
        for(i_= 0; i_<N_U; i_++) {
            o->dV[0]+= derivs.Qu[i_]*o->l[MAT_IDX(i_, k, N_U)];
        }
        for(i_= 0; i_<N_U; i_++) {
            d1= 0.0;
            for(j_= 0; j_<N_U; j_++) {
                d1+= o->l[MAT_IDX(j_, k, N_U)]*derivs.Quu[SYMTRI_MAT_IDX(j_, i_)];
            }
            o->dV[1]+= 0.5*o->l[MAT_IDX(i_, k,  N_U)]*d1;
        }

        // Vx(:,i)     = Qx  + K_i'*Quu*k_i + K_i'*Qu  + Qux'*k_i;
        memcpy(derivs.Vx, derivs.Qx, sizeof(double)*N_X);
        addMul2Tri(derivs.Vx, derivs.Quu, &o->L[MAT_IDX3(0, 0, k, N_U, N_X)], N_U, N_X, &o->l[MAT_IDX(0, k, N_U)], N_U, 1, derivs.dummy);
        for(i_= 0; i_<N_X; i_++) // row Vx
            for(j_= 0; j_<N_U; j_++) // col o->L' & row Qu
                derivs.Vx[i_]+= o->L[MAT_IDX3(j_, i_, k, N_U, N_X)]*derivs.Qu[j_];

        for(i_= 0; i_<N_X; i_++) // row Vx
            for(j_= 0; j_<N_U; j_++) // col Qxu & row o->l
                derivs.Vx[i_]+= derivs.Qxu[MAT_IDX(i_, j_, N_X)]*o->l[MAT_IDX(j_, k, N_U)];
        
    
        // Vxx(:,:,i)  = Qxx + K_i'*Quu*K_i + K_i'*Qux + Qux'*K_i;
        memcpy(derivs.Vxx, derivs.Qxx, sizeof(double)*sizeofQxx);
        addSquareTri(derivs.Vxx, derivs.Quu, &o->L[MAT_IDX3(0, 0, k, N_U, N_X)], N_U, N_X, derivs.dummy);
        for(i_= 0; i_<N_X; i_++)
            for(j_= 0; j_<N_X; j_++)
                for(k_= 0; k_<N_U; k_++) {
                    d1= o->L[MAT_IDX3(k_, i_, k, N_U, N_X)]*derivs.Qxu[MAT_IDX(j_, k_, N_X)];
                    if(i_==j_) d1*= 2.0;
                    
                    derivs.Vxx[SYMTRI_MAT_IDX(i_, j_)]+= d1;
                }
        

        // g_norm= mean(max(abs(l) ./ (abs(u)+1),[],1));
        g_norm_max= 0.0;
        for(i_= 0; i_<N_U; i_++) {
            g_norm_i= fabs(o->l[MAT_IDX(i_, k, N_U)]) / (fabs(u[i_])+1.0);
            if(g_norm_i>g_norm_max) g_norm_max= g_norm_i;
        }
        g_norm_sum+= g_norm_max;
    }
    
    o->g_norm= g_norm_sum/((double)(o->n_hor-1));
    
    return 0;
}
