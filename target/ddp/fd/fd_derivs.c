#include "iLQG.h"

int calc_derivs(tOptSet *o) {
    int k, i_;
    int N= o->n_hor;
    trajEl_t *t= o->trajectory + N;

    if(!calc_deriv(t, k, o->p, N)) return 0;
    

    t--;
    for(k= N-1; k>=0; k--, t--) {
        if(!calc_deriv(t, k, o->p, N)) return 0;
        
        limitsU(t, k, o->p, N);
    }
    return 1;
}

#define EPS 1e-16;
// https://en.wikipedia.org/wiki/Finite_difference#Finite_difference_in_several_variables
int calc_deriv(trajEl_t *t, int k, double **p, int N) {
    double xu[N_X+N_U];
    double J_p[N_X+N_U];
    double J_m[N_X+N_U];
    double f_p[N_X*(N_X+N_U)];
    double f_m[N_X*(N_X+N_U)];
    double h[N_X+N_U];
    double *u= xu+N_X;
    double xu_, xu_i, xu_j, J_pp, J_mm, f_pp[N_X], f_mm[NX], dJ1, dJ2, *df1, *df2;
    int i, j, i_u, j_u, i_f, j_f, k_f, n, n_f;
    
    for(i= 0; i<N_X; i++) xu[i]= t->x[i];
    for(i= 0; i<N_U; i++) u[i]= t->u[i];
    
    n= (k<N)? N_X+N_U: N_X;
    
    for(i= 0, i_u= 0, i_f= 0; i<n; i++, i_f+=N_X) {
        xu_= xu[i];
        h[i]= EPS*fabs(xu_);
        if(h[i]<EPS) h[i]= EPS;
        xu[i]= xu_ + h[i];
        h[i]= xu[i] - xu_; 
        
        if(!calcXVariableAux(xu, t, k, p))
            return 0;
        if(!calcXUVariableAux(xu, u, t, k, p))
            return 0;
        J_p[i]= ddpJ(xu, u, t, k, p, N);
        if(isNANorINF(J_p[i]))
            return 0;
        if(k<N)
            if(!ddpf(&f_p[i_f], xu, u, t, k, p, N)) return 0;
        
        xu[i]= xu_ - h[i];
        if(!calcXVariableAux(xu, t, k, p))
            return 0;
        if(!calcXUVariableAux(xu, u, t, k, p))
            return 0;
        J_m[i]= ddpJ(xu, u, t, k, p, N);
        if(isNANorINF(J_m[i]))
            return 0;
        if(k<N)
            if(!ddpf(&f_m[i_f], xu, u, t, k, p, N)) return 0;

        xu[i]= xu_;
        
        
        dJ1= (J_p[i] - J_m[i]) / (2.0*h[i]);
        dJ2= (J_p[i] - 2.0*t->c[k] + J_m[i]) / (h[i]*h[i]);
        if(i>=N_X) {
            t->cu[i_u]= dJ1;
            t->cuu[UTRI_MAT_IDX(i_u, i_u)]= dJ2;
            df1= t->fu + i_u*N_U;
#if FULL_DDP
            n_f= sizeofQuu;
            df2= &t->fuu[UTRI_MAT_IDX(i_u, i_u)];
#endif            
            i_u++;
        } else {
            t->cx[i]= dJ1;
            t->cxx[UTRI_MAT_IDX(i, i)]= dJ2;
            df1= t->fx + i*N_X;
#if FULL_DDP
            n_f= sizeofQuu;
            df2= &t->fxx[UTRI_MAT_IDX(i, i)];
#endif            
        }
        if(k<N) {
            for(j_f= 0, k_f= 0; j_f<N_X; j_f++, k_f+= n_f) {
                df1[k_f]= (f_p[i_f + j_f] - f_m[i_f + j_f]) / (2.0*h[i]);
#if FULL_DDP
                df2[k_f]= (f_p[i_f + j_f] - 2.0*(t+1)->x[j_f] + f_m[i_f + j_f]) / (h[i]*h[i]);
#endif            
            }
        }
    }
    
    for(j= 1, j_u= 0; j<n; j++) {
        for(i= 0, i_u= 0; i<j; i++) {
            xu_i= xu[i];
            xu_j= xu[j];
            xu[i]= xu_i + h[i];
            xu[j]= xu_j + h[j];
            if(!calcXVariableAux(xu, t, k, p))
                return 0;
            if(!calcXUVariableAux(xu, u, t, k, p))
                return 0;
            J_pp= ddpJ(xu, u, t, k, p, N);
            if(isNANorINF(J_pp))
                return 0;
#if FULL_DDP
            if(k<N)
                if(!ddpf(f_pp, xu, u, t, k, p, N)) return 0;
#endif    
            xu[i]= xu_i - h[i];
            xu[j]= xu_j - h[j];
            if(!calcXVariableAux(xu, t, k, p))
                return 0;
            if(!calcXUVariableAux(xu, u, t, k, p))
                return 0;
            J_mm= ddpJ(xu, u, t, k, p, N);
            if(isNANorINF(J_mm))
                return 0;
#if FULL_DDP
            if(k<N)
                if(!ddpf(f_mm, xu, u, t, k, p, N)) return 0;
#endif    
            xu[i]= xu_i;
            xu[j]= xu_j;
            
            dJ2= (J_pp - J_p[i] - J_p[j] + 2.0*t->c[k] - J_m[i] - J_m[j] + J_mm) / (2.0*h[i]*h[j]);
            if(j>=N_X) {
                if(i>=N_X) {
                    t->cuu[UTRI_MAT_IDX(i_u, j_u)]= dJ2;
#if FULL_DDP
                    n_f= sizeofQuu;
                    df2= &t->fuu[UTRI_MAT_IDX(i_u, j_u)];
#endif    
                    i_u++;
                } else {
                    t->cxu[MAT_IDX(i, j_u, N_X)]= dJ2;
#if FULL_DDP
                    n_f= sizeofQxu;
                    df2= &t->fxu[MAT_IDX(i, j_u, N_X)];
#endif    
                }
                j_u++;
            } else {
                t->cxx[UTRI_MAT_IDX(i, j)]= dJ2;
#if FULL_DDP
                    n_f= sizeofQxx;
                    df2= &t->fxx[UTRI_MAT_IDX(i, j)];
#endif    
            }
#if FULL_DDP
            for(j_f= 0, k_f= 0; j_f<N_X; j_f++, k_f+= n_f) {
                df2[k_f]= (f_pp[j_f] - f_p[i*N_X + j_f] - J_p[j*N_X + j_f] + 2.0*(t+1)->x[j_f] - f_m[i*N_X + j_f] - f_m[j*N_X + j_f] + f_mm[j_f]) / (2.0*h[i]*h[j]);
            }
#endif            
        }
    }
    
    return 1;
}
