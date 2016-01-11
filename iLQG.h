#ifndef ILQG_H
#define ILQG_H

#define INIT_OPTSET {0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, NULL, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0, 0, NULL, NULL, NULL, {0.0, 0.0}}

#ifndef PRNT
#define PRNT printf
#endif

#define DO_PREFIX1(VAL)  1 ## VAL
#define PREFIX1(VAL)     DO_PREFIX1(VAL)


typedef struct paramDesc {
  char *name;
  int size;
  int is_var;
} tParamDesc;

typedef struct optSet {
    int n_hor;
    int debug_level;
    double *x0, *x_nom, *u_nom, *x_new, *u_new, *l, *L, new_cost, cost, dcost, lambda, g_norm, expected;
    double **p;
    const double *alpha;
    int n_alpha;
    double lambdaMax;
    double lambdaMin;
    double lambdaInit;
    double dlambdaInit;
    double lambdaFactor;
    int max_iter;
    double tolGrad;
    double tolFun;
    double zMin;
    int regType;
    int iterations;
    int *log_linesearch;
    double *log_z;
    double *log_cost;
    double dV[2];
} tOptSet;

void standard_parameters(tOptSet *o);
int iLQG(tOptSet *o);
int initialize_iLQG(tOptSet *o);
char *setOptParam(tOptSet *o, const char *name, const double *value, const int n);

static double max(double a, double b) {
    return (a > b)? a: b;
}

static double min(double a, double b) {
    return (a < b)? a: b;
}

#endif /* ILQG_H */
