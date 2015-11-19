#define MAT_IDX(r, c, n) LTRI_MAT_IDX(r, c)

double fmax(double a, double b) {
    if(a>b)
        return a;
    else
        return b;
}

void switch_row_and_colMAT(double *A, int n, int i, int j) {
	double tmp;
	int k;

	for(k= 0; k<j; k++) {
		tmp= A[MAT_IDX(j, k, n)];
		A[MAT_IDX(j, k, n)]= A[MAT_IDX(i, k, n)];
		A[MAT_IDX(i, k, n)]= tmp;
	}
	
	for(k= j+1; k<i; k++) {
		tmp= A[MAT_IDX(k, j, n)];
		A[MAT_IDX(k, j, n)]= A[MAT_IDX(i, k, n)];
		A[MAT_IDX(i, k, n)]= tmp;
	} 

	for(k= i+1; k<n; k++) {
		tmp= A[MAT_IDX(k, j, n)];
		A[MAT_IDX(k, j, n)]= A[MAT_IDX(k, i, n)];
		A[MAT_IDX(k, i, n)]= tmp;
	}

	tmp= A[MAT_IDX(j, j, n)];
	A[MAT_IDX(j, j, n)]= A[MAT_IDX(i, i, n)];
	A[MAT_IDX(i, i, n)]= tmp;
}

void switch_row_and_colVEC(int *b, int i, int j) {
	int id;
	id= b[i];
	b[i]= b[j];
	b[j]= id;
}

void jthIteration(double *L_, int n, int j) {
    int i, k;

    L_[MAT_IDX(j, j, n)]= sqrt(L_[MAT_IDX(j, j, n)]);
    for(i=j+1; i<n; i++) {
        L_[MAT_IDX(i, j, n)]/= L_[MAT_IDX(j, j, n)];
        for(k=j+1; k<=i; k++) {
            L_[MAT_IDX(i, k, n)]-= L_[MAT_IDX(i, j, n)]*L_[MAT_IDX(k, j, n)];
        }
    }
}

double mod_chol(double *A, int n, double *E, int *P) {
	int phaseone= 1;
	int i, j, k;
	double tau= pow(mxGetEps(), 1./3.);
	double taubar= pow(mxGetEps(), 2.0/3.0);
    double mu= 0.1;
	double gamma= 0.0;
	double tmp= 0.;
	double delta= 0.;
    double deltaprev= 0.;
    double normj;
    double lambda_hi;
    double lambda_lo;

    if(n==1) {
        delta= (taubar * fabs(A[MAT_IDX(0, 0, n)])) - A[MAT_IDX(0, 0, n)];
        if(delta>0.0) E[0]= delta; else E[0]= 0.;
        if(A[MAT_IDX(0, 0, n)]==0.0) E[0]= taubar;
        A[MAT_IDX(0, 0, n)]= sqrt(A[MAT_IDX(0, 0, n)]+E[0]);
        P[0]= 0;
        return E[0];
    }

    for(i= 0; i<n; i++) {
		P[i]= i;
		E[i]= 0.0;
	}
	for(i= 0; i<n; i++) {
		tmp= fabs(A[MAT_IDX(i, i, n)]);
		if(tmp>gamma) gamma= tmp;
		if(A[MAT_IDX(i, i, n)]<0.0) phaseone= 0;
	}
    PR1("gamma= %f\n", gamma);

	/* Phase one, A potentially positive-definite */
	j= 0;
	while((j<n) && phaseone) {
		double tmp_max= A[MAT_IDX(j, j, n)];
		double tmp_min= A[MAT_IDX(j, j, n)];
		int id= j;

        // find extrema
		for(i= j+1; i<n; i++) {
			if(tmp_max<A[MAT_IDX(i, i, n)]) {
				tmp_max= A[MAT_IDX(i, i, n)];
				id= i;
			}
			if(tmp_min>A[MAT_IDX(i, i, n)])
				tmp_min= A[MAT_IDX(i, i, n)];
		}
//		if(tmp_max<taubar*gamma || tmp_min<tau*tmp_max) {
        PR3("tmp_max= %f, tmp_min= %f @ j= %d\n", tmp_max, tmp_min, j+1);
		if(tmp_max<taubar*gamma || tmp_min<-mu*tmp_max) {
            PR3("if(%f<taubar*gamma || %f<-mu*tmp_max) @ j= %d\n", tmp_max, tmp_min, j+1);
			phaseone= 0;
			break; //go to phasetwo
		} else {
			/* Pivot on maximum diagonal of remaining submatrix */
			if(id!=j) {
				//switch rows and cols of id and j of A
				switch_row_and_colMAT(A, n, id, j);
				switch_row_and_colVEC(P, id, j);
			}

            tmp_min= 0.0;
            for(i=j+1; i<n; i++) {
				tmp= A[MAT_IDX(i, i, n)] - A[MAT_IDX(i, j, n)]*A[MAT_IDX(i, j, n)]/A[MAT_IDX(j, j, n)];
				if(tmp_min>tmp) tmp_min= tmp;
            }
            PR2("tmp_min= %f @ j= %d\n", tmp_min, j+1);
//			if(tmp_min<tau*gamma) {
			if(tmp_min<-mu*gamma) {
                PR2("if(%f<-mu*gamma) @ j= %d\n", tmp_min, j+1);
				phaseone= 0;
				break;
			} else {// perform jth iteration of factorization
                jthIteration(A, n, j);
				j++;
			}
		}//end of if
	}//end of while

	/* Phase two, A not positive-definite */
	if(!phaseone && (j==(n-1))) {
        PR("if(!phaseone && (j==(n-1)))\n");

		delta= -A[MAT_IDX(n-1, n-1, n)] + fmax(tau*A[MAT_IDX(n-1, n-1, n)]/(tau-1.), taubar*gamma);
		A[MAT_IDX(n-1, n-1, n)]+= delta;
		A[MAT_IDX(n-1, n-1, n)]= sqrt(A[MAT_IDX(n-1, n-1, n)]);
		E[n-1]= delta;
        deltaprev= delta;
	}

	if(!phaseone && (j<(n-1))) {
		double *g= malloc(n*sizeof(double));
		// k= number of iterations performed in phase one
        PR1("if(!phaseone && (j<(n-1))), j= %d\n", j+1);
		k= j - 1;

		/* Caculate lower Gerschgorin bound */
		for(i=k+1; i<n; i++) {
			g[i]= A[MAT_IDX(i, i, n)];
			for(j=k+1; j<=i-1; j++) {
				g[i]-= fabs(A[MAT_IDX(i, j, n)]);
			}
			for(j=i+1; j<n; j++) {
				g[i]-= fabs(A[MAT_IDX(j, i, n)]);
			}
		}
		/* Modified Cholesky Decomposition */
		for(j=k+1; j<n-2; j++) {
			// Pivot on maximum lower Gerschgorin bound estimate
			int id= j;
			tmp= g[id];
			for(i=j+1; i<n; i++) {
				if(tmp<g[i]) {
					tmp= g[i];
					id= i;
				}
			}
			if(id!=j) {
				switch_row_and_colMAT(A, n, id, j);
				switch_row_and_colVEC(P, id, j);
				tmp= g[id];
				g[id]= g[j];
				g[j]= tmp;
			}
			//Calculate E[j, j] and add to diagonal
			normj= 0.;
			for(i=j+1; i<n; i++) {
				normj+= fabs(A[MAT_IDX(i, j, n)]);
                PR1("abs(A(i, j))= %f\n", fabs(A[MAT_IDX(i, j, n)]));
            }
			PR1("normj= %f\n", normj);
			delta= fmax(0.0, fmax(fmax(normj, taubar*gamma)-A[MAT_IDX(j, j, n)], deltaprev));
			if(delta>0) {
				A[MAT_IDX(j, j, n)]+= delta;
				deltaprev= delta;
				E[j]= delta;
			}

			//Update Gerschgorin bound estimates
			if(A[MAT_IDX(j, j, n)]!=normj) {
				tmp= 1.0 - normj/A[MAT_IDX(j, j, n)];
				for(i=j+1; i<n; i++) {
					g[i]+= fabs(A[MAT_IDX(i, j, n)])*tmp;
				}
			}
			// Perform jth iteration of factorization
            jthIteration(A, n, j);
		}
		//Final 2x2 submatrix
		tmp= sqrt((A[MAT_IDX(n-2, n-2, n)]-A[MAT_IDX(n-1, n-1, n)])*(A[MAT_IDX(n-2, n-2, n)]-A[MAT_IDX(n-1, n-1, n)]) +4.0*A[MAT_IDX(n-1, n-2, n)]*A[MAT_IDX(n-1, n-2, n)]);
		lambda_hi= ((A[MAT_IDX(n-2, n-2, n)]+A[MAT_IDX(n-1, n-1, n)]) + tmp)*0.5;
		lambda_lo= ((A[MAT_IDX(n-2, n-2, n)]+A[MAT_IDX(n-1, n-1, n)]) - tmp)*0.5;
		delta= fmax(fmax(0.0, -lambda_lo + fmax(tau*(lambda_hi-lambda_lo)/(1.0-tau), taubar*gamma)), deltaprev);
		if(delta>0) {
			A[MAT_IDX(n-2, n-2, n)]+= delta;
			A[MAT_IDX(n-1, n-1, n)]+= delta;
			deltaprev= delta;
			E[n-2]= delta;
			E[n-1]= delta;
		}
		A[MAT_IDX(n-2, n-2, n)]= sqrt(A[MAT_IDX(n-2, n-2, n)]);
		A[MAT_IDX(n-1, n-2, n)]/= A[MAT_IDX(n-2, n-2, n)];
		A[MAT_IDX(n-1, n-1, n)]= sqrt(A[MAT_IDX(n-1, n-1, n)]-A[MAT_IDX(n-1, n-2, n)]*A[MAT_IDX(n-1, n-2, n)]);

		free(g);
	}
	return deltaprev;
}

void solve(double *L_, int *P, double *b, double *x, int n) {
	int i, k;
    double *y;

    for(k= 0; k<n; k++)	x[k]= b[P[k]];

    // Solve L*y= b;
    for(k= 0; k<n; k++) {
        for(i= 0; i<k; i++)
            x[k]-= x[i]*L_[MAT_IDX(k, i, n)];
	    x[k]/= L_[MAT_IDX(k, k, n)];
 	}

    // Solve L'*X= Y;
    for(k= n-1; k >= 0; k--) {
       for(i= k+1; i<n; i++)
             x[k]-= x[i]*L_[MAT_IDX(i, k, n)];
       x[k]/= L_[MAT_IDX(k, k, n)];
    }

	//Solve P'*y= x
	y= malloc(n*sizeof(double));
    for(k= 0; k<n; k++) y[P[k]]= x[k];
    for(k= 0; k<n; k++) x[k]= y[k];
    free(y);
}

void chol_inv(double *L_, int *P, double *invA, int n) {
	int i, k, l, cp, rp, tp;
    double *x;

//    PR("invert malloc\n");
    x= malloc(n*sizeof(double));
    for(l= 0; l<n; l++) {
//        PR1("invert1: col %d\n", l+1);
        x[l]= 1.0;
        for(k= l+1; k<n; k++) x[k]= 0.0;

        // Solve L*y= b; die Einträge in x sind nach P sortiert
        for(k= l; k<n; k++) {
            for(i= l; i<k; i++)
                x[k]-= x[i]*L_[MAT_IDX(k, i, n)];
    	    x[k]/= L_[MAT_IDX(k, k, n)];
     	}

//        PR("invert2: col\n");
        // Solve L'*X= Y;
        for(k= n-1; k >= l; k--) {
//            PR1("invert2,%d\n", k);
            for(i= k+1; i<n; i++)
                 x[k]-= x[i]*L_[MAT_IDX(i, k, n)];
            x[k]/= L_[MAT_IDX(k, k, n)];

//            PR2("invert set matrix[k: %d, l:%d]\n", k, l);
//            PR2("invert set matrix[kp: %d, lp:%d]\n", P[k], P[l]);
            rp= P[k]; cp= P[l]; if(rp<cp) { tp= rp; rp= cp; cp= tp; }
//            PR3("invert set matrix[%d, %d -> %d]\n", rp, cp, MAT_IDX(rp, cp, n));
            invA[MAT_IDX(rp, cp, n)]= x[k];
        }
    }
//    PR("invert free\n");
    free(x);
}
