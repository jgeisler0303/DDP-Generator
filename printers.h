#include <string.h>

void printX(double *x, int n) {
    int i, k;
    for(i= 0; i<N_X; i++) {
        PRNT1("x(%d, :)= [", i+1);
        for(k= 0; k<n; k++) {
            PRNT1("%g", x[i+k*N_X]);
            if(k<(n-1)) {
                PRNT(", ");
            } else {
                PRNT("];\n");
            }
        }
    }
}

void printU(double *u, int n) {
    int i, k;
    for(i= 0; i<N_U; i++) {
        PRNT1("u(%d, :)= [", i+1);
        for(k= 0; k<n; k++) {
            PRNT1("%g", u[i+k*N_U]);
            if(k<(n-1)) {
                PRNT(", ");
            } else {
                PRNT("];\n");
            }
        }
    }
}

void printVec(double *A, int n, char *nm) {
    int i;
    PRNT2("%s= [%g", nm, A[0]);
    for(i= 1; i<n; i++) PRNT1(", %g", A[i]);
    PRNT("]\n");
}

void printTri(double *A, int n, char *nm) {
    int c, r, i, sl= strlen(nm);

    PRNT1("%s= [", nm);
    for(r= 0; r<n; r++) {
        for(c= 0; c<n; c++) {
            if(r>c) {
                i= (r*(r+1))/2 + c;
            } else {
                i= (c*(c+1))/2 + r;
            }
            PRNT1("%g", A[i]);
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

void printMat(double *A, int n, int m, char *nm) {
    int c, r, i, sl= strlen(nm);;

    PRNT1("%s= [", nm);
    for(r= 0; r<n; r++) {
        for(c= 0; c<m; c++) {
            PRNT1("%g", A[r + n*c]);
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

