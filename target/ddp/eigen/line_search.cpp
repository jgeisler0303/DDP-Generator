#include <cmath>

#include "ddp.h"
#include "line_search.h"

#ifndef PARALLEL_LINESEARCH
    #define PARALLEL_LINESEARCH 0
#endif
 
int line_search(tOptSet *o) {
    int successCandidate= -1;
    bool abort = false;
//     double start_time[MAX_ALPHA];
//     double before_ordered_time[MAX_ALPHA];
//     double in_ordered_time[MAX_ALPHA];
//     double end_time[MAX_ALPHA];
//     double total_time= omp_get_wtime();
    
    #pragma omp parallel default(none) firstprivate(o) shared(abort, successCandidate) if(PARALLEL_LINESEARCH) num_threads(NUMBER_OF_THREADS)
//     #pragma omp parallel default(none) firstprivate(o) shared(abort, successCandidate, start_time, before_ordered_time, in_ordered_time, end_time) if(1) num_threads(NUMBER_OF_THREADS)
    {
    #pragma omp for ordered schedule(dynamic, 1) 
    for(int i= 0; i < o->n_alpha; i++) {
//         start_time[i]= omp_get_wtime();
        double expected= 0.0, dcost= 0.0, cnew= 0.0, alpha= 0.0;
        int pass_ok= false;
        
        #ifdef _OPENMP
        const int useCandidate= omp_get_thread_num();
        #else
        const int useCandidate= 0;
        #endif
        
        #pragma omp flush (abort)
        if(!abort) {
            alpha= o->alpha[i];
            pass_ok= forward_pass(o->candidates[useCandidate], o, alpha, cnew, 0);
            if(pass_ok) {
                dcost= o->cost - cnew;
                expected= -alpha*(o->dV[0] + alpha*o->dV[1]);
            }
        }
//         before_ordered_time[i]= omp_get_wtime();
        #pragma omp ordered
        {
//             in_ordered_time[i]= omp_get_wtime();
            #pragma omp flush (abort)
            if(!abort) {
                double z= 0.0;
                if(pass_ok) {
                    
                    if(expected > 0) {
                        z = dcost/expected;
                        if(z > o->zMin) {
                            successCandidate= useCandidate;
                            abort= true;
                            #pragma omp flush (abort)
                        }
                    } else
                        if(o->log_line) o->log_line->neg_exp_red++;
                } else {
                    successCandidate= -2;
                    abort= true;
                    #pragma omp flush (abort)
                }
                
                if(i==o->n_alpha-1 || successCandidate>=0) {
                    o->new_cost= cnew;
                    o->dcost= dcost;
                    o->last_z= z;
                    o->n_ls= i+1;
                    if(o->log_line) {
                        o->log_line->n_line_searches= i+1;
                        o->log_line->alpha= alpha;
                        o->log_line->z= z;
                        o->log_line->cost= cnew;
                        o->log_line->dcost= dcost;
                        o->log_line->expected_red= expected;
                    }
                }
            } // !abort
        } // ordered
//         end_time[i]= omp_get_wtime();
    } // for
    } // parallel
//     total_time= omp_get_wtime() - total_time;
    
//     for(int i= 0; i < o->n_alpha; i++) {
//         double pass_time= before_ordered_time[i] - start_time[i];
//         double wait_time= in_ordered_time[i] - before_ordered_time[i];
//         double ordered_time= end_time[i] -  in_ordered_time[i];
//         double start_delay= 0.0;
//         double ordered_delay= 0.0;
//         if(i>0) {
//             start_delay= start_time[i] - start_time[i-1];
//             ordered_delay= in_ordered_time[i] - in_ordered_time[i-1];
//         }
//         printf("%2d: %12.5f %12.5f %12.5f %12.5f %12.5f\n", i, pass_time*1000.0, wait_time*1000.0, ordered_time*1000.0, start_delay*1000.0, ordered_delay*1000.0);
//     }
//     printf("total: %12.5f for %d with %d\n\n", total_time*1000.0, o->n_ls, omp_get_num_procs());
    return successCandidate;
}
