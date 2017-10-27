#include <string.h>
#include <stdio.h>

#include "printMat.h"
#include "matMult.h"

void printVec(const double *A, const int n, const char *nm) {
    int i, nc, ns= 255-20;
    char s[256], *ps;
    
    ps= &s[0];
    nc= sprintf(ps, "%s= [%g", nm, A[0]);
    ns-= nc;
    ps+= nc;
    for(i= 1; i<n && ns>0; i++) {
        nc= sprintf(ps, ", %g", A[i]);
        ns-= nc;
        ps+= nc;
    }
    if(ns>0) sprintf(ps, "]");
    
    PRNT("%s\n",s);
}

void printTri(const double *A, const int n, const char *nm) {
    int c, r, i, sl= strlen(nm);

    PRNT("%s= [", nm);
    for(r= 0; r<n; r++) {
        for(c= 0; c<n; c++) {
            if(r>c) {
                i= (r*(r+1))/2 + c;
            } else {
                i= (c*(c+1))/2 + r;
            }
            PRNT("%g", A[i]);
            if(c<(n-1)) {
                PRNT(", ");
            } else {
                if(r<(n-1)) {
                    PRNT("\n   ");
                    for(i= 0; i<sl; i++) PRNT(" ");
                } else {
                    PRNT("]\n");
                }
            }
        }
    }
}

void printMat(const double *A, const int n, const int m, const char *nm) {
    int c, r, i, sl= strlen(nm);;

    PRNT("%s= [", nm);
    for(r= 0; r<n; r++) {
        for(c= 0; c<m; c++) {
            PRNT("%g", A[r + n*c]);
            if(c<(m-1)) {
                PRNT(", ");
            } else {
                if(r<(n-1)) {
                    PRNT("\n   ");
                    for(i= 0; i<sl; i++) PRNT(" ");
                } else {
                    PRNT("]\n");
                }
            }
        }
    }
}


