#include "iLQG.h"
#include "matMult.h"
#include "printMat.h"


typedef int(fJfunc_t)(double* f, double *J, void *t_, const void *m_, int k, const tOptSet *o);

int fJfuncL(double* f, double* J, void *t_, const void *m_, int k, const tOptSet *o);
int fJfuncF(double* f, double* J, void *t_, const void *m_, int k, const tOptSet *o);


double calc_dxu(double xu, const tOptSet *o) {
    double dxu;
    double h= o->h_fd;
    
    dxu= h*fabs(xu);
    if(dxu<h) dxu= h;
    
    return dxu;
}

// https://en.wikipedia.org/wiki/Finite_difference#Finite_difference_in_several_variables
#if FULL_DDP
int calc_deriv(double *cx, double *cu, double *cxx, double *cuu, double *cxu, double *fx, double *fu, double *fxx, double *fuu, double *fxu, fJfunc_t *func, void *t, const void *m, const double *f_nom, int k, int N, const tOptSet *o)
#else
int calc_deriv(double *cx, double *cu, double *cxx, double *cuu, double *cxu, double *fx, double *fu, fJfunc_t *func, void *t, const void *m, const double *f_nom, int k, int N, const tOptSet *o)
#endif
{
    double J_p[N_X+N_U];
    double J_m[N_X+N_U];
    double J_pp, J_mm;
    double f_p[N_X*(N_X+N_U)];
    double f_m[N_X*(N_X+N_U)];
    double f_pp[N_X], f_mm[N_X];
    double dxu[N_X+N_U];
    double save_xu_i, save_xu_j;
    double dJ1, dJ2;
    double *df1, *df2; // pointers to destination for derivatives of f
    int n, n_f= 0;
    trajFin_t *tf= (trajFin_t*)t;
    trajEl_t  *tr= (trajEl_t *)t;
    
    double J_nom;
    if(k<N) {
        n= N_X+N_U;
        J_nom= tr->c;
    } else {
        n= N_X;
        J_nom= tf->c;
    }

    // calc differences, first derivatives and symmetric second derivatives
    for(int i= 0, i_u= 0, i_f= 0; i<n; i++, i_f+=N_X) {
        // calc positive difference
        if(k<N) {
            if(i<N_X) {
                dxu[i]= calc_dxu(tr->x[i], o);
                save_xu_i= tr->x[i];
                tr->x[i]+= dxu[i];
                dxu[i]= tr->x[i] - save_xu_i;
            } else {
                dxu[i]= calc_dxu(tr->u[i_u], o);
                save_xu_i= tr->u[i_u];
                tr->u[i_u]+= dxu[i];
                dxu[i]= tr->u[i_u] - save_xu_i; 
            }
        } else {
            dxu[i]= calc_dxu(tf->x[i], o);
            save_xu_i= tf->x[i];
            tf->x[i]+= dxu[i];
            dxu[i]= tf->x[i] - save_xu_i;
        }
        if(!func(&f_p[i_f], &J_p[i], t, m, k, o)) return 0;

        // calc neagtive difference
        if(k<N) {
            if(i<N_X)
                tr->x[i]= save_xu_i - dxu[i];
            else
                tr->u[i_u]= save_xu_i - dxu[i];
        } else {
            tf->x[i]= save_xu_i - dxu[i];
        }
        if(!func(&f_m[i_f], &J_m[i], t, m, k, o)) return 0;

        // restore nominal x/u
        if(k<N) {
            if(i<N_X)
                tr->x[i]= save_xu_i;
            else
                tr->u[i_u]= save_xu_i;
        } else {
            tf->x[i]= save_xu_i;
        }
        
        // first and second derivative of J
        dJ1= (J_p[i] - J_m[i]) / (2.0*dxu[i]);
        dJ2= (J_p[i] - 2.0*J_nom + J_m[i]) / (dxu[i]*dxu[i]);
        if(i<N_X) {
            cx[i]= dJ1;
            cxx[UTRI_MAT_IDX(i, i)]= dJ2;
            df1= &fx[i*N_X];
#if FULL_DDP
            n_f= sizeofQxx;
            df2= &fxx[UTRI_MAT_IDX(i, i)];
#endif            
        } else {
            cu[i_u]= dJ1;
            cuu[UTRI_MAT_IDX(i_u, i_u)]= dJ2;
            df1= &fu[i_u*N_X];
#if FULL_DDP
            n_f= sizeofQuu;
            df2= &fuu[UTRI_MAT_IDX(i_u, i_u)];
#endif            
            i_u++;
        }
        
        // calc first and second derivative of f
        if(k<N) {
            for(int j_f= 0, k_f= 0; j_f<N_X; j_f++, k_f+= n_f) {
                df1[j_f]= (f_p[i_f + j_f] - f_m[i_f + j_f]) / (2.0*dxu[i]);
#if FULL_DDP
                df2[k_f]= (f_p[i_f + j_f] - 2.0*f_nom[j_f] + f_m[i_f + j_f]) / (dxu[i]*dxu[i]);
#endif            
            }
        }
    }

    // calc asymmetric second derivatives
    for(int j= 1, j_u= 0; j<n; j++) {
        for(int i= 0, i_u= 0; i<j; i++) {
            // calc positive difference
            if(k<N) {
                if(i<N_X) {
                    save_xu_i= tr->x[i];
                    tr->x[i]= save_xu_i + dxu[i];
                } else {
                    save_xu_i= tr->u[i_u];
                    tr->u[i_u]= save_xu_i + dxu[i];
                }
                if(j<N_X) {
                    save_xu_j= tr->x[j];
                    tr->x[j]= save_xu_j + dxu[j];
                } else {
                    save_xu_j= tr->u[j_u];
                    tr->u[j_u]= save_xu_j + dxu[j];
                }
            } else {
                save_xu_i= tf->x[i];
                tf->x[i]= save_xu_i + dxu[i];
                save_xu_j= tf->x[j];
                tf->x[j]= save_xu_j + dxu[j];
            }
#if FULL_DDP
            if(!func(f_pp, &J_pp, t, m, k, o)) return 0;
#else
            if(!func(NULL, &J_pp, t, m, k, o)) return 0;
#endif
            
            // calc negative difference
            if(k<N) {
                if(i<N_X) {
                    tr->x[i]= save_xu_i - dxu[i];
                } else {
                    tr->u[i_u]= save_xu_i - dxu[i];
                }
                if(j<N_X) {
                    tr->x[j]= save_xu_j - dxu[j];
                } else {
                    tr->u[j_u]= save_xu_j - dxu[j];
                }
            } else {
                tf->x[i]= save_xu_i - dxu[i];
                tf->x[j]= save_xu_j - dxu[j];
            }
#if FULL_DDP
            if(!func(f_mm, &J_mm, t, m, k, o)) return 0;
#else
            if(!func(NULL, &J_mm, t, m, k, o)) return 0;
#endif
            
            // restore nominal x/u
            if(k<N) {
                if(i<N_X) {
                    tr->x[i]= save_xu_i;
                } else {
                    tr->u[i_u]= save_xu_i;
                }
                if(j<N_X) {
                    tr->x[j]= save_xu_j;
                } else {
                    tr->u[j_u]= save_xu_j;
                }
            } else {
                tf->x[i]= save_xu_i;
                tf->x[j]= save_xu_j;
            }
            
            
            // second derivative of J
            // mean of mixed backwards-backwards and forwards-forwards difference
            dJ2= (J_pp - J_p[i] - J_p[j] + 2.0*J_nom - J_m[i] - J_m[j] + J_mm) / (2.0*dxu[i]*dxu[j]);
            if(j<N_X) {
                cxx[UTRI_MAT_IDX(i, j)]= dJ2;
#if FULL_DDP
                n_f= sizeofQxx;
                df2= &fxx[UTRI_MAT_IDX(i, j)];
#endif    
            } else {
                if(i<N_X) {
                    cxu[MAT_IDX(i, j_u, N_X)]= dJ2;
#if FULL_DDP
                    n_f= sizeofQxu;
                    df2= &fxu[MAT_IDX(i, j_u, N_X)];
#endif    
                } else {
                    cuu[UTRI_MAT_IDX(i_u, j_u)]= dJ2;
#if FULL_DDP
                    n_f= sizeofQuu;
                    df2= &fuu[UTRI_MAT_IDX(i_u, j_u)];
#endif    
                    i_u++;
                }
            }
            
            // calc second derivative of f
#if FULL_DDP
            if(k<N)
                for(int j_f= 0, k_f= 0; j_f<N_X; j_f++, k_f+= n_f) {
                    df2[k_f]= (f_pp[j_f] - f_p[i*N_X + j_f] - f_p[j*N_X + j_f] + 2.0*f_nom[j_f] - f_m[i*N_X + j_f] - f_m[j*N_X + j_f] + f_mm[j_f]) / (2.0*dxu[i]*dxu[j]);
                    // mexPrintf("f_nom= %15.12g, nom= %g\n", f_nom[j_f], (f_pp[j_f] - f_p[i*N_X + j_f] - f_p[j*N_X + j_f] + 2.0*f_nom[j_f] - f_m[i*N_X + j_f] - f_m[j*N_X + j_f] + f_mm[j_f]));
                    // mexPrintf("k= %d, i= %d, j= %d, df2= %g, f_p= %15.12g, f_p= %15.12g, f_nom= %15.12g, f_m= %15.12g, f_m= %15.12g, f_pp= %15.12g, f_mm= %15.12g, dxu[i]= %g, dxu[j]= %g\n", k, i, j, df2[k_f], f_p[i*N_X + j_f], f_p[j*N_X + j_f] ,f_nom[j_f], f_m[i*N_X + j_f], f_m[j*N_X + j_f], f_pp[j_f], f_mm[j_f], dxu[i], dxu[j]);
                }
#endif
        }
        if(j>=N_X) j_u++;
    }
    
    return 1;
}

int calc_derivs(const tOptSet *o) {
    int k;
    int N= o->n_hor;
    double *f_nom;
    
    trajEl_t *t= o->nominal->t + N -1;
    trajFin_t *f= &o->nominal->f;
    
    multipliersEl_t *m= o->multipliers.t + N - 1;
    const multipliersFin_t *mf= &o->multipliers.f;

    
#if FULL_DDP
    if(!calc_deriv(f->cx, NULL, f->cxx, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &fJfuncF, f, mf, NULL, N, N, o)) return 0;
#else
    if(!calc_deriv(f->cx, NULL, f->cxx, NULL, NULL, NULL, NULL, &fJfuncF, f, mf, NULL, N, N, o)) return 0;
#endif

    f_nom= f->x;
    for(k= N-1; k>=0; k--, t--, m--) {
        limitsU(t, k, o->p, N);
        
#if FULL_DDP
        if(!calc_deriv(t->cx, t->cu, t->cxx, t->cuu, t->cxu, t->fx, t->fu, t->fxx, t->fuu, t->fxu, &fJfuncL, t, m, f_nom, k, N, o)) return 0;
#else
        if(!calc_deriv(t->cx, t->cu, t->cxx, t->cuu, t->cxu, t->fx, t->fu, &fJfuncL, t, m, f_nom, k, N, o)) return 0;
#endif
        f_nom= t->x;        
    }

    return 1;
}

