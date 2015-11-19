#include <math.h>
#include <tmwtypes.h>

#define PRNT(s) mexPrintf(s)
#define PRNT1(s, x1) mexPrintf(s, x1)
#define PRNT2(s, x1, x2) mexPrintf(s, x1, x2)
#define PRNT3(s, x1, x2, x3) mexPrintf(s, x1, x2, x3)
#define PRNT4(s, x1, x2, x3, x4) mexPrintf(s, x1, x2, x3, x4)

#define PR(s)
#define PR1(s, x1)
#define PR2(s, x1, x2)
#define PR3(s, x1, x2, x3)
#define PR4(s, x1, x2, x3, x4)

#define UTRI_MAT_IDX(r, c) (((c)*((c)+1))/2 + (r))
#define LTRI_MAT_IDX(r, c) (((r)*((r)+1))/2 + (c))
#define TRI_MAT_IDX(r, c) ((r>c)? UTRI_MAT_IDX(c, r): UTRI_MAT_IDX(r, c))

// int n_hor;
// int debug_level= 0;
// double *x_ic, *x_nom= NULL, *u_nom= NULL, *x_new= NULL, *u_new= NULL, *alpha= NULL, *beta= NULL, V0;

typedef struct paramDesc {
  char *name;
  int size;
} tParamDesc;

typedef struct p_Opt {
    int max_linesearch;
    int max_iter;
    double eta;
    double c;
    int16_T iterations;
    int16_T *linesearch;
    int16_T success;
    double *e;
    double *J;
    double *a0;
    double *regu_max;
    double *regu_sum;
} tPOpt;

typedef struct optSet {
    int n_hor;
    int debug_level;
    double *x_ic, *x_nom, *u_nom, *x_new, *u_new, *alpha, *beta, V0;
    tPOpt p_opt;
    double **p;
} tOptSet;

// tPOpt p_opt= {50, 100, 1e-6, 0, NULL, 0, NULL, NULL, NULL, NULL, NULL};

#define DDP_PROBLEM
    #include "ddpProblemDescription.h" // defines N_X, N_U, param_names, param_sizes
#undef DDP_PROBLEM
#include "mod_chol/mod_chol_tri.h"
#include "ddpMatMult.h"
#include "printers.h"

double backSweep(tOptSet *o) {
    int k, i, j, l;
    int Pchol[N_U];
    double regu, Echol[N_U];
    double **p= o->p, N= o->n_hor;
    double a0, Vx[N_X], Vxx[(N_X*(N_X+1))/2];
    double *x, *u;
    double Qx[N_X], Qu[N_U], Qxu[N_X*N_U];
    double Qxx[(N_X*(N_X+1))/2], Quu[(N_U*(N_U+1))/2], invQuu[(N_U*(N_U+1))/2]; // obere rechts dreiecksmatrizen
    double fx[N_X*N_X], fu[N_X*N_U];
#ifdef HAS_FXX
    double fxx[N_X*(N_X*(N_X+1))/2]; // 2. Ableitung nach x [[obere rechte dreiecksmatrix zustand 1],...,[obere rechte dreiecksmatrix zustand n]]
#endif
#ifdef HAS_FUU
    double fuu[N_X*(N_U*(N_U+1))/2]; // 2. Ableitung nach u [[obere rechte dreiecksmatrix zustand 1],...,[obere rechte dreiecksmatrix zustand n]]
#endif
#ifdef HAS_FXU
    double fxu[N_X*N_X*N_U]; // Ableitung nach x und u [[matrix in zeilen nach x, in spalten nach u zustand 1],...,[matrix in zeilen nach x, in spalten nach u zustand n]]
#endif

    double s, dummy[N_X*N_X];

    a0= 0.0;
    if(o->p_opt.regu_max!=NULL) o->p_opt.regu_max[o->p_opt.iterations]= 0.;
    if(o->p_opt.regu_sum!=NULL) o->p_opt.regu_sum[o->p_opt.iterations]= 0.;
    k= o->n_hor;

    if(o->debug_level>=3) PRNT("  backSweep:\n");
#ifndef IS_TV
    #define DDP_FXK
        #include "ddpProblemDescription.h"
    #undef DDP_FXK
    if(o->debug_level>=4) {
        printMat(fx, N_X, N_X, "        fx");
        printMat(fu, N_X, N_U, "        fu");
    }
#endif
//Vx=, Vxx=
    x= &(o->x_nom[N_X*k]);
    #define DDP_F
        #include "ddpProblemDescription.h"
    #undef DDP_F

    if(o->debug_level>=4) {
        printVec(Vx, N_X, "        Vx");
        printTri(Vxx, N_X, "        Vxx");
    }

    for(k= o->n_hor-1; k>=0; k--) {
        if(o->debug_level>=3) {
            PRNT1("    @stage %d\n", k);
        }
        // Lx=, Lu=, Lxx=, Lxu=, Luu= wird in Q.. geschrieben
        x= &(o->x_nom[N_X*k]);
        u= &(o->u_nom[N_U*k]);
        #define DDP_L
            #include "ddpProblemDescription.h"
        #undef DDP_L
        #ifdef IS_TV
            #define DDP_FXK
                #include "ddpProblemDescription.h"
            #undef DDP_FXK
        #endif
        if(o->debug_level>=4) {
            printVec(Qx, N_X, "        Lx");
            printTri(Qxx, N_X, "        Lxx");
            printVec(Qu, N_U, "        Lu");
            printTri(Quu, N_U, "        Luu");
            printMat(Qxu, N_X, N_U, "        Lxu");
            #ifdef IS_TV
                printMat(fx, N_X, N_X, "        fx");
                printMat(fu, N_X, N_U, "        fu");
            #endif
        }

        // Qx= deriv.Lx + Vx*deriv.fx;  % 1 x n
        if(o->debug_level>=6) PRNT("      Qx= deriv.Lx + Vx*deriv.fx;\n");
        addMulVec(Qx, Vx, fx, N_X, N_X);
        // Qu= deriv.Lu + Vx*deriv.fu;  % 1 x m
        if(o->debug_level>=6) PRNT("      Qu= deriv.Lu + Vx*deriv.fu;\n");
        addMulVec(Qu, Vx, fu, N_X, N_U);
        // Qxx= deriv.Lxx + deriv.fx'*Vxx*deriv.fx;  % n x n
        if(o->debug_level>=6) PRNT("      Qxx= deriv.Lxx + deriv.fx'*Vxx*deriv.fx;\n");
        addMul2Tri(Qxx, Vxx, fx, N_X, N_X, dummy);
        #ifdef HAS_FXX
            // Qxx= deriv.Lxx + deriv.fx'*Vxx*deriv.fx + hypermul(Vx, deriv.fxx);  % n x n
            for(i= 0; i<(N_X*(N_X+1))/2; i++) for(j= 0; j<N_X; j++) Qxx[i]+= Vx[j]*fxx[i + (N_X*(N_X+1))/2*j];
        #endif
        // Quu= deriv.Luu + deriv.fu'*Vxx*deriv.fu;  % m x m
        if(o->debug_level>=6) PRNT("      Quu= deriv.Luu + deriv.fu'*Vxx*deriv.fu;\n");
        addMul2Tri(Quu, Vxx, fu, N_X, N_U, dummy);
        #ifdef HAS_FUU
            // Quu= deriv.Luu + deriv.fu'*Vxx*deriv.fu + hypermul(Vx, deriv.fuu);  % m x m
            for(i= 0; i<(N_U*(N_U+1))/2; i++) for(j= 0; j<N_X; j++) Quu[i]+= Vx[j]*fuu[i + (N_U*(N_U+1))/2*j];
        #endif
        // Qxu= deriv.Lxu + deriv.fx'*Vxx*deriv.fu;  % n x m
        if(o->debug_level>=6) PRNT("      Qxu= deriv.Lxu + deriv.fx'*Vxx*deriv.fu;\n");
        addMul2Mat(Qxu, Vxx, fx, N_X, N_X, fu, N_X, N_U, dummy);
        #ifdef HAS_FXU
            // Qxu= deriv.Lxu + deriv.fx'*Vxx*deriv.fu + hypermul(Vx, deriv.fxu);  % n x m
            for(i= 0; i<N_X*N_U; i++) for(j= 0; j<N_X; j++) Qxu[i]+= Vx[j]*fxu[i + N_X*N_U*j];
        #endif

        if(o->debug_level>=4) {
            printVec(Qx, N_X, "        Qx");
            printVec(Qu, N_U, "        Qu");
            printTri(Qxx, N_X, "        Qxx");
            printTri(Quu, N_U, "        Quu");
            printMat(Qxu, N_X, N_U, "        Qxu");
        }

        if(o->debug_level>=6) PRNT("      cholesky of Quu\n");
        regu= mod_chol(Quu, N_U, Echol, Pchol);
        if(o->p_opt.regu_max!=NULL) o->p_opt.regu_max[o->p_opt.iterations]= fmax(o->p_opt.regu_max[o->p_opt.iterations], regu);
        if(o->p_opt.regu_sum!=NULL) o->p_opt.regu_sum[o->p_opt.iterations]+= regu;

        if(o->debug_level>=6) PRNT("      invert Quu\n");
        chol_inv(Quu, Pchol, invQuu, N_U);
        if(o->debug_level>=6) PRNT("      done inverting\n");

        if(o->debug_level>=4) {
            printTri(invQuu, N_U, "        invQuu");
        }

        // alpha{k}= -Quu^-1*Qu'; % m x 1
        for(i= 0; i<N_U; i++) {
            o->alpha[i + N_U*k]= 0.0;
            for(j= 0; j<N_U; j++) {
                o->alpha[i + N_U*k]-= invQuu[TRI_MAT_IDX(i, j)]*Qu[j];
            }
        }
        // beta{k}= -Quu^-1*Qxu'; % m x n
        for(i= 0; i<N_U; i++) {
            for(j= 0; j<N_X; j++) {
                o->beta[i + N_U*j + N_U*N_X*k]= 0.0;
                for(l= 0; l<N_U; l++) {
                    o->beta[i + N_U*j + N_U*N_X*k]-= invQuu[TRI_MAT_IDX(i, l)]*Qxu[j + N_X*l];
                }
            }
        }

        // a= a - Qu*Quu^-1*Qu'; % 1 x 1
        // a= a + Qu*alpha; % 1 x 1
        for(i= 0; i<N_U; i++) {
            a0+= Qu[i]*o->alpha[i + N_U*k];
        }
        // Vx= Qx + Qu*beta{k}; % 1 x n
        for(i= 0; i<N_X; i++) {
            Vx[i]= Qx[i];
            for(j= 0; j<N_U; j++) {
                Vx[i]+= Qu[j]*o->beta[j + N_U*i + N_U*N_X*k];
            }
        }
        // Vxx= Qxx + Qxu*-Quu^-1*Qxu';   % n x n
        // Vxx= Qxx + Qxu*beta{k};   % n x n
        for(j= 0; j<N_X; j++) {
            for(i= 0; i<=j; i++) {
                Vxx[UTRI_MAT_IDX(i, j)]= Qxx[UTRI_MAT_IDX(i, j)];
                s= 0.0;
                for(l= 0; l<N_U; l++) {
                    s+= Qxu[i + l*N_X]*o->beta[l + N_U*j + N_U*N_X*k];
                }
                if(i!=j) {
                    // transponiertes element zur symmetrierung addieren
                    for(l= 0; l<N_U; l++) {
                        s+= Qxu[j + N_X*l]*o->beta[l + N_U*i + N_U*N_X*k];
                    }
                    s/= 2.0;
                }
                Vxx[UTRI_MAT_IDX(i, j)]+= s;
            }
        }

        if(o->debug_level>=4) {
            printVec(&(o->alpha[N_U*k]), N_U, "        alpha");
            printMat(&(o->beta[N_U*N_X*k]), N_U, N_X, "        beta");
            printVec(Vx, N_X, "        Vx");
            printTri(Vxx, N_X, "        Vxx");
        }
        if(o->debug_level>=3) {
            PRNT2("    end stage %d: a0= %e\n", k, a0);
        }
        if(o->debug_level>=3) PRNT1("      cholesky regularisierung delta max= %g\n", regu);
    }
    return a0;
}

int predictTrajectory(tOptSet *o) {
    int i, k, l;
    double dx;

    for(i= 0; i<N_X; i++) o->x_new[i]= o->x_ic[i]; // ic
    for(k= 0; k<o->n_hor; k++) {
        for(l= 0; l<N_U; l++) o->u_new[l + N_U*k]= o->u_nom[l + N_U*k] + o->alpha[l + N_U*k];
        for(i= 0; i<N_X; i++) {
            dx= o->x_new[i + N_X*k]-o->x_nom[i + N_X*k];
            for(l= 0; l<N_U; l++) {
                o->u_new[l + N_U*k]+= o->beta[l + N_U*i + N_U*N_X*k]*dx;
            }
        }
        
        if(!ddpf(&(o->x_new[N_X*(k+1)]), &(o->x_new[N_X*k]), &(o->u_new[N_U*k]), k, o->p, o->n_hor)) return 0;
    }
    if(o->debug_level>=5) {
        printX(o->x_new, o->n_hor+1);
        printU(o->u_new, o->n_hor);
    }
    return 1;
}

double objectiveFunc(tOptSet *o, double *x, double *u) {
    int k;
    double J= 0.0;

    for(k= 0; k<o->n_hor; k++)
        J+= ddpJ(&(x[N_X*k]), &(u[N_U*k]), k, o->p, o->n_hor);

    J+= ddpJ(&(x[N_X*o->n_hor]), 0, -1, o->p, o->n_hor);

    return J;
}

int lineSearch(tOptSet *o, double a0) {
    double e= 1.0, DV0, J= 0.0, *temp;
    int i, j, k, success= 1;

    if(o->debug_level>=2) PRNT("  lineSearch:\n");

    for(i= 0; i<o->p_opt.max_linesearch; i++) {
        if(predictTrajectory(o) && !mxIsNaN(J= objectiveFunc(o, o->x_new, o->u_new)) && !mxIsInf(J)) {
            DV0= J-o->V0;
            if(o->debug_level>=2) {
                PRNT4("      @no %d: V0= %g, J= %g, DV0= %g\n", i+1, o->V0, J, DV0);
                PRNT2("      (fabs(DV0)/fabs(e*(1-e/2)*a0))= %g, e= %g\n", (fabs(DV0)/fabs(e*(1-e/2)*a0)), e);
            }
            if(DV0<=0.0 & (fabs(DV0)/fabs(e*(1-e/2)*a0))>o->p_opt.c) break;
        } else {
            if(o->debug_level>=2) {
                PRNT1("      @no %d: prediction or objective failed\n", i+1);
            }
        }
        e/= 2.0;
        for(j= 0; j<o->n_hor; j++) for(k= 0; k<N_U; k++) o->alpha[k + N_U*j]/= 2.0;
    }

    if(i>=o->p_opt.max_linesearch) {
        success= 0;
        i--;
        if(o->debug_level>=1) PRNT("-->maximale Anzahl Iterationen überschritten\n");
    } else {
        o->V0= J;
        temp= o->x_nom;
        o->x_nom= o->x_new;
        o->x_new= temp;
        temp= o->u_nom;
        o->u_nom= o->u_new;
        o->u_new= temp;
    }
    if(o->debug_level>=4) {
        printX(o->x_nom, o->n_hor+1);
        printU(o->u_nom, o->n_hor);
    }
    if(o->debug_level>=2) PRNT2("  lineSearch end s_nr= %d, e= %e\n", i+1, e);
    if(o->p_opt.linesearch!=NULL) o->p_opt.linesearch[o->p_opt.iterations]= i+1;
    if(o->p_opt.e!=NULL) o->p_opt.e[o->p_opt.iterations]= e;
    if(o->p_opt.J!=NULL) o->p_opt.J[o->p_opt.iterations]= J;

    return success;
}

int initializeDDP(tOptSet *o) {
    int i;

    for(i= 0; i<N_X; i++) o->x_nom[i]= o->x_ic[i];
    for(i= 0; i<o->n_hor; i++)
        if(!ddpf(&(o->x_nom[N_X*(i+1)]), &(o->x_nom[N_X*i]), &(o->u_nom[N_U*i]), i, o->p, o->n_hor))
            return 0;

    if(mxIsNaN(o->V0= objectiveFunc(o,o->x_nom, o->u_nom)) || mxIsInf(o->V0)) return 0;

    for(i= 0; i<o->p_opt.max_iter; i++) {
        if(o->p_opt.linesearch!=NULL) o->p_opt.linesearch[i]= 0;
        if(o->p_opt.e!=NULL) o->p_opt.e[i]= 0.;
        if(o->p_opt.J!=NULL) o->p_opt.J[i]= 0.;
        if(o->p_opt.a0!=NULL) o->p_opt.a0[i]= 0.;
    }
    if(o->debug_level>=1) PRNT1("initializeDDP, J= %g\n", o->V0);
    if(o->debug_level>=4) printX(o->x_nom, o->n_hor+1);
    if(o->debug_level>=4) printU(o->u_nom, o->n_hor);
    return 1;
}

int iterateDDP(tOptSet *o) {
    int ls_ok= 1, success= 0;
    double a0;

    o->p_opt.iterations= 0;
    if(!initializeDDP(o)) {
        o->p_opt.success= 0;
        return 0;
    }
    while((o->p_opt.iterations<o->p_opt.max_iter) && ls_ok && !success) {
        if(o->debug_level>=1) PRNT1("iteration %d\n", o->p_opt.iterations+1);
        if(mxIsNaN(a0= backSweep(o)) || mxIsInf(a0)) {
            PRNT1("Terminating because a0= %g\n", a0);
            o->p_opt.success= 0;
            return 0;
        }
        if(o->p_opt.a0!=NULL) o->p_opt.a0[o->p_opt.iterations]= a0;
        if(fabs(a0)>o->p_opt.eta)
            ls_ok= lineSearch(o, a0);
        else {
            success= 1;
            if(o->p_opt.linesearch!=NULL) o->p_opt.linesearch[o->p_opt.iterations]= 0;
            if(o->p_opt.J!=NULL) o->p_opt.J[o->p_opt.iterations]= o->V0;
        }
        o->p_opt.iterations++;
    }

    if(o->debug_level>=1) {
        PRNT2("iterateDDP end after %d iterations a0= %g\n", o->p_opt.iterations, a0);
    }
    o->p_opt.success= success && ls_ok;
    return success && ls_ok;
}
