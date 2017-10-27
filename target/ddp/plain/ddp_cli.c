// MATLAB Mex function wrapper for iLQG algorithm
// Copyright (c) 2016 Jens Geisler


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "iLQG.h"
#include "ini.h"
#include "strtonum.h"

int *k_counters;
double *u_nom[N_U];

static int ddp_ini_handler(void* user, const char* section, const char* name, const char* val_str) {
    long ival;
    double dval;
    int conv_res, i, si;
    const char *endptr;
    const char *err_msg;
    char u_name[10];        
    double vals[32];
    int count= 0;

    tOptSet *o= (tOptSet *)user;

    if(strcasecmp(section, "init")==0) {
        if(strcasecmp(name, "n_hor")==0) {
            if(o->n_hor > 0) {
                fprintf(stderr, "Parameter n_hor already set, ignoring subsequent definitions\n");
                return 0;
            }
            conv_res= strtoint(val_str, &ival, &endptr);
            if(conv_res) {
                fprintf(stderr, "Error reading n_hor: %s\n", strtonum_error(conv_res));
                return 0;
            }
            if(*endptr != '\0') {
                fprintf(stderr, "Error reading n_hor: number followed by other characters\n");
                return 0;
            }
            if(ival<1 || ival>10000) {
                fprintf(stderr, "Error reading n_hor: must be > 0 and < 10000\n");
                return 0;                
            }
            
            o->n_hor= ival;
            
            for(i= 0; i<N_U; i++)
                u_nom[i]= (double *)malloc(o->n_hor * sizeof(double));
            
            for(i=0; i<n_params; i++) {
                si= (paramdesc[i]->size==-1)? o->n_hor+1: paramdesc[i]->size;
                o->p[i]= (double *)malloc(si * sizeof(double));
            }
            
            return 1;
        }
        if(strcasecmp(name, "x0")==0) {
            endptr= val_str;
            while(*endptr!='\0' && k_counters[N_U]<N_X) {
                conv_res= strtodouble(val_str, &dval, &endptr);
                if(conv_res) {
                    fprintf(stderr, "Error reading x0: %s\n", strtonum_error(conv_res));
                    return 0;
                }
                o->x0[k_counters[N_U]]= dval;
                k_counters[N_U]++;
                
                if(*endptr==';' || *endptr==',') endptr++;
                val_str= endptr;
            }
            if(k_counters[N_U]>=N_X && *endptr!='\0') {
                    fprintf(stderr, "Found more values for x0 than %d, ignoring the rest\n", N_X);
                    return 0;
            }
                
            return 1;
        }
        
        for(i= 0; i<N_U; i++) {
            snprintf(u_name, 9, "u_nom%d", i+1);
            
            if(strcasecmp(name, u_name)==0) {
                if(o->n_hor == 0) {
                    fprintf(stderr, "Parameter n_hor must be set before any u_nom, ignoring %s\n", u_name);
                    return 0;
                }
                endptr= val_str;
                while(*endptr!='\0' && k_counters[i]<o->n_hor) {
                    conv_res= strtodouble(val_str, &dval, &endptr);
                    if(conv_res) {
                        fprintf(stderr, "Error reading %s: %s\n", u_name, strtonum_error(conv_res));
                        return 0;
                    }
                    u_nom[i][k_counters[i]]= dval;
                    k_counters[i]++;
                    
                    if(*endptr==';' || *endptr==',') endptr++;
                    val_str= endptr;
                }
                if(k_counters[i]>=o->n_hor && *endptr!='\0') {
                        fprintf(stderr, "Found more values for %s than %d, ignoring the rest\n", u_name, o->n_hor);
                        return 0;
                }
                    
                return 1;
            }
        }
        
        fprintf(stderr, "Unknown parameter \"%s\" in section \"init\" in ini file\n", name);
        return 0;
    }
    if(strcasecmp(section, "param")==0) {
        if(o->n_hor == 0) {
            fprintf(stderr, "Parameter n_hor must be set before any parameter, ignoring %s\n", name);
            return 0;
        }
        for(i=0; i<n_params; i++) {
            if(strcasecmp(name, paramdesc[i]->name)==0) {
                si= (paramdesc[i]->size==-1)? o->n_hor+1: paramdesc[i]->size;
                
                endptr= val_str;
                while(*endptr!='\0' && k_counters[N_U+1+i]<si) {
                    conv_res= strtodouble(val_str, &dval, &endptr);
                    if(conv_res) {
                        fprintf(stderr, "Error reading parameter \"%s\": %s\n", paramdesc[i]->name, strtonum_error(conv_res));
                        return 0;
                    }
                    o->p[i][k_counters[N_U+1+i]]= dval;
                    k_counters[N_U+1+i]++;
                    
                    if(*endptr==';' || *endptr==',') endptr++;
                    val_str= endptr;
                }
                if(k_counters[N_U+1+i]>=si && *endptr!='\0') {
                        fprintf(stderr, "Found more values for parameter \"%s\" than %d, ignoring the rest\n", paramdesc[i]->name, si);
                        return 0;
                }
                    
                return 1;
            }
        }
        
        fprintf(stderr, "Unknown parameter \"%s\" in section \"param\" in ini file\n", name);
        return 0;
    }
    if(strcasecmp(section, "opt_param")==0) {
        count= 0;
        endptr= val_str;
        while(*endptr!='\0' && count<32) {
            conv_res= strtodouble(val_str, &dval, &endptr);
            if(conv_res) {
                fprintf(stderr, "Error reading opt parameter \"%s\": %s\n", name, strtonum_error(conv_res));
                return 0;
            }
            vals[count]= dval;
            count++;
            
            if(*endptr==';' || *endptr==',') endptr++;
            val_str= endptr;
        }
        if(count>=32 && *endptr!='\0') {
                fprintf(stderr, "Found more values for parameter \"%s\" than 32, ignoring the rest\n", name);
                return 0;
        }
        
        err_msg= setOptParam(o, name, vals, count);
        if(err_msg) {
            fprintf(stderr, "Error setting optimization parameter \"%s\": %s\n", name, err_msg);
            return 0;
        }
        
        return 1;
    }
    
    fprintf(stderr, "Unknown section \"%s\" in ini file\n", section);
    return 0;
}

void writeResult(tOptSet *o, const char* file_name) {
    int i, k;
    FILE *fptr;

    fptr = fopen(file_name, "w");
    if(fptr == NULL) {
        printf("Error opening file \"%s\" for writing results", file_name);
        return;
    }
   
    for(i= 0; i<N_X; i++)
        fprintf(fptr,"x%02d              ", i);
    for(i= 0; i<N_U; i++)
        fprintf(fptr,"u%02d              ", i);
    fprintf(fptr,"\n");
    
    for(k= 0; k<o->n_hor; k++) {
        for(i= 0; i<N_X; i++)
            fprintf(fptr,"%16.8e ", o->nominal->t[k].x[i]);
        for(i= 0; i<N_U; i++)
            fprintf(fptr,"%16.8e ", o->nominal->t[k].u[i]);
        fprintf(fptr,"\n");
    }
   
    fclose(fptr);
}

void writeIter(tOptSet *o, const char* file_name) {
    int i;
    FILE *fptr;

    fptr = fopen(file_name, "w");
    if(fptr == NULL) {
        printf("Error opening file \"%s\" for writing iteration data", file_name);
        return;
    }
   
    fprintf(fptr,"cost            ");
    fprintf(fptr,"lambda          ");
    fprintf(fptr,"alpha           ");
    fprintf(fptr,"gradient        ");
    fprintf(fptr,"dcost           ");
    fprintf(fptr,"reduction_ratio ");
    fprintf(fptr,"\n");
    
    for(i= 0; i<o->iterations; i++) {
        fprintf(fptr,"%16.8e ", o->log[i].cost);
        fprintf(fptr,"%16.8e ", o->log[i].lambda);
        fprintf(fptr,"%16.8e ", o->log[i].alpha);
        fprintf(fptr,"%16.8e ", o->log[i].g_norm);
        fprintf(fptr,"%16.8e ", o->log[i].dcost);
        fprintf(fptr,"%16.8e ", o->log[i].z);
        fprintf(fptr,"\n");
    }
   
    fclose(fptr);
}

int main(int argc, char* argv[]) {
    const char *ini_file, *res_file, *iter_file;
    clock_t begin, end;
    int k, i, res, ret= 0;
    
    tOptSet o= INIT_OPTSET;
    standard_parameters(&o);
    o.p= (double **)malloc(n_params*sizeof(double *));
    k_counters= (int *)malloc((N_U+n_params+1)*sizeof(int));
    memset(k_counters, 0, (N_U+n_params+1)*sizeof(int));
    
    if(argc>=2) {
        ini_file= argv[1];
        if(argc>4)
            printf("Warning: ignoring any command line parameters past the third\n");
    } else {
        ini_file= "param.ini";
            printf("Using default parameter file \"param.ini\"\n");
    }
    if(argc>=3)
        res_file= argv[2];
    else
        res_file= "ddp_result.dat";
    if(argc>=4)
        iter_file= argv[3];
    else
        iter_file= "ddp_iter.dat";
    
    res= ini_parse(ini_file, ddp_ini_handler, &o);
    if(res < 0) {
        printf("Can't load '%s'\n", ini_file);
        ret= 1;
    } else {
        if(res > 0) {
            printf("There were errors reading '%s' in line %d\n", ini_file, res);
            ret= 1;
        }
        if(o.n_hor==0) {
            printf("Error: n_hor was not set\n");
            ret= 1;            
        } else {
            for(i= 0; i<N_U; i++) {
                if(k_counters[i]<o.n_hor) {
                    printf("Error: number of elements of nominal u[%d] is less than %d (actual= %d)\n", i, o.n_hor, k_counters[i]);
                    ret= 1;            
                }
            }
            if(k_counters[N_U]<N_X) {
                printf("Error: number of elements of x0 is less than %d (actual= %d)\n", N_X, k_counters[N_U]);
                ret= 1;            
            }
            for(i= 0; i<n_params; i++) {
                int si= (paramdesc[i]->size==-1)? o.n_hor+1: paramdesc[i]->size;
                if(k_counters[N_U+1+i]<si) {
                    printf("Error: number of elements of parameter %d is less than %d (actual= %d)\n", i, si, k_counters[N_U+1+i]);
                    ret= 1;            
                }
            }
        }
    }
    
    free(k_counters);
    if(ret) {
        if(o.n_hor>0) {
            for(i= 0; i>n_params; i++)
                free(o.p[i]);
            for(i= 0; i>N_U; i++)
                free(u_nom[i]);
        }
        free(o.p);
        return ret;
    }
    
    for(i= 0; i<NUMBER_OF_THREADS+1; i++) {
        o.trajectories[i].t= malloc(sizeof(trajEl_t)*o.n_hor);
        // memset(o.trajectories[i].t, 1, sizeof(trajEl_t)*(N-1));
    }
    o.multipliers.t= malloc(sizeof(multipliersEl_t)*(o.n_hor+1));
    
    o.log= (tLogLine *)malloc(o.max_iter*sizeof(tLogLine));
    memset(o.log, 0, o.max_iter*sizeof(tLogLine));
    
    
    
    printf("Initializing constants\n");
    if(!init_opt(&o)) {
        printf("Error initializing constants\n");
        ret= 1;
    } else {
        printf("Initializing trajectory\n");
        for(k= 0; k<o.n_hor; k++) {
            for(i= 0; i<N_U; i++)
                o.nominal->t[k].u[i]= u_nom[i][k];
        }
        
        double cost;
        if(!forward_pass(o.candidates[0], &o, 0.0, &cost, 0)) {
            printf("Error in initial forward pass\n");
            ret= 1;
        } else {
            o.cost= cost;
            makeCandidateNominal(&o, 0);
            
            printf("Starting DDP\n");

            begin = clock();
            iLQG(&o);
            end = clock();

            printLog(&o);
            printf("Time for DDP: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
            
            writeResult(&o, res_file);
            writeIter(&o, iter_file);
        }
    }
    

    for(i= 0; i>n_params; i++)
        free(o.p[i]);
    for(i= 0; i>N_U; i++)
        free(u_nom[i]);
    free(o.log);
    free(o.p);
    
    for(i= 0; i<NUMBER_OF_THREADS+1; i++)
        free(o.trajectories[i].t);
    
    free(o.multipliers.t);
    
    return ret;
}
