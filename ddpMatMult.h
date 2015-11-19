void addMulVec(double base[], double a[], double b[], int n_r, int n_c) {
    // base= base + a*b; base= [1 x n_c]; a= [1 x n_c]; b= [n_r x n_c];
    int ci, ri;

    for(ci= 0; ci<n_c; ci++) {
        for(ri= 0; ri<n_r; ri++) {
            base[ci]+= a[ri]*b[ri + ci*n_r];
        }
    }
}

void addMul2Tri(double base[], double b[], double a[], int n_r, int n_c, double ba[]) {
    // base= base + a'*b*a; base= [n_c x n_c]; a= [n_r x n_c]; b= [n_r x n_r];
    // base: obere rechte dreiecksmatrix
    int ci, ri, si;
    double s;

    // b*a= [n_r x n_c]
    for(ri= 0; ri<n_r; ri++) {
        for(ci= 0; ci<n_c; ci++) {
            ba[ri + ci*n_r]= 0.0;
            for(si= 0; si<n_r; si++) {
                ba[ri + ci*n_r]+= b[TRI_MAT_IDX(ri, si)]*a[si + ci*n_r];
            }
        }
    }
    // a'*ba= [n_c x n_c]
    for(ci= 0; ci<n_c; ci++) {
        for(ri= 0; ri<=ci; ri++) {
            s= 0.0;
            for(si= 0; si<n_r; si++) {
                s+= a[si + ri*n_r]*ba[si + ci*n_r]; // a': ri= si, ci= ri
            }
            if(ri!=ci) {
                // transponiertes element zur symmetrierung addieren
                for(si= 0; si<n_r; si++) {
                    s+= a[si + ci*n_r]*ba[si + ri*n_r];
                }
                s/= 2.0;
            }
            base[UTRI_MAT_IDX(ri, ci)]+= s;
        }
    }
}

void addMul2Mat(double base[], double b[], double a[], int n_ra, int n_ca, double c[], int n_rc, int n_cc, double bc[]) {
    // base= base + a'*b*c; base= [n_ca x n_cc]; a= [n_ra x n_ca]; b= [n_ra x n_rc]; c= [n_rc x n_cc];
    int ci, ri, si;

    // b*c= [n_ra x n_cc]
    for(ri= 0; ri<n_ra; ri++) {
        for(ci= 0; ci<n_cc; ci++) {
            bc[ri + ci*n_ra]= 0.0;
            for(si= 0; si<n_rc; si++) {
                bc[ri + ci*n_ra]+= b[TRI_MAT_IDX(ri, si)]*c[si + ci*n_rc];
            }
        }
    }
    // a'*bc= [n_ca x n_cc]
    for(ri= 0; ri<n_ca; ri++) {
        for(ci= 0; ci<n_cc; ci++) {
            for(si= 0; si<n_ra; si++) {
                base[ri + ci*n_ca]+= a[si + ri*n_ra]*bc[si + ci*n_ra]; // a': ri= si, ci= ri
            }
        }
    }
}
