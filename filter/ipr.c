#include "fd.h"

void IPR(double *psitot, double *newpsitot, double *eval, par_st par, ipr_st *ipr, long_st *ist){
    /*** Computes the inverse participation ratio: sum_n |psi(x_n)|^4 ***/
    FILE *pipr, *piprc;
    long jms, jgrid, count;
    double sum, *rho;
    ipr_st *iprcut;
    rho = calloc(ist->ngrid,sizeof(double));

    omp_set_dynamic(0);
    omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(jms)
    for (jms = 0; jms < ist->mstot; jms++){
        // For each quasiparticle state
        sum = 0.0;
        for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
            // Sum up the IPR function value at all the grid points
            sum += sqr(psitot[jms*ist->ngrid+jgrid])*sqr(psitot[jms*ist->ngrid+jgrid]);
        }
        ipr[jms].val = sum;
        ipr[jms].index = jms;
    }
    
    
    printf("Sorting psi by IPR:\n"); fflush(0);
    

    qsort(ipr, ist->mstot, sizeof(ipr_st), compare_vals);

    pipr = fopen("IPR.dat", "w");
    for (jms = 0; jms < ist->mstot; jms++){
    fprintf(pipr, "%ld %.8f\n", jms, ipr[jms].val);
    }

    printf("\tjms    ipr    idx\n");
    for (jms = 0; jms < ist->mstot; jms++){
        if (jms == 0){
            printf("This is the delocalized index: %ld\n", ipr[jms].index);
            for (jgrid = 0; jgrid<ist->ngrid; jgrid++){
                rho[jgrid] = psitot[ipr[jms].index*ist->ngrid + jgrid];
            }
            writeCubeFile(rho, par, *ist, "deloc-wfn.cube");
        }
        if (jms == ist->iprcut-1){
            printf("This is the localized index: %ld\n", ipr[jms].index);
            for (jgrid = 0; jgrid<ist->ngrid; jgrid++){
                rho[jgrid] = psitot[ipr[jms].index*ist->ngrid + jgrid];
            }
            writeCubeFile(rho, par, *ist, "loc-wfn.cube");
        }
        if (jms % 20 == 0){
            printf("\t%ld %.6f %ld\n", jms, ipr[jms].val, ipr[jms].index );
        }
    }
    if (( iprcut = calloc(ist->iprcut, sizeof(ipr_st))) == NULL){
        printf("OUT OF MEMORY: iprcut\n");
        nerror("iprcut");
    }
 
    for (jms = 0; jms < ist->iprcut; jms++){
        iprcut[jms].val = ipr[jms].val;
        iprcut[jms].index = ipr[jms].index;
    }
    printf("Sorting the delocalized psi by index:\n");
    qsort(iprcut, ist->iprcut, sizeof(ipr_st), compare_index);

    printf("The IPR cutoff = %ld\n", ist->iprcut); fflush(0);
    printf("IPR of max delocalization = %g\n", (double) 1/ist->ngrid);fflush(0);
    printf("Larger values correspond to less delocalized states.\n"); fflush(0);

    
    piprc = fopen("IPRcut.dat", "w");
    for (jms = 0; jms < ist->iprcut; jms++ ){
        //printf("jms: %ld\n", jms); fflush(0);
        for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
            //printf("index: %ld\n", iprcut[jms].index);fflush(0);
            newpsitot[jms*ist->ngrid + jgrid] = psitot[iprcut[jms].index*ist->ngrid + jgrid];
        }
        eval[jms] = eval[iprcut[jms].index];
        fprintf(piprc, "%ld %ld %.8f\n", jms, iprcut[jms].index, iprcut[jms].val);fflush(0);
    }
    printf("\nDelocalized states obtained!\n");fflush(0);
    //printf("New mstot = %ld\t Nstates removed = %ld\n", i, ist->mstot - i);
    fclose(pipr); fclose(piprc); free(iprcut); free(rho);
}

// int compare_function(const void *a,const void *b) {
// double *x = (double *) a;
// double *y = (double *) b;
// if (*x < *y) return -1;
// else if (*x > *y) return 1; return 0;
// }

int compare_vals(const void *a, const void *b)
{
    struct ipr_st *x1 = (struct ipr_st *)a;
    struct ipr_st *x2 = (struct ipr_st *)b;
    if ((*x1).val < (*x2).val)
        return -1;
    else if ((*x1).val > (*x2).val)
        return 1;
    else
        return 0;
}

int compare_index(const void *a, const void *b)
{
    struct ipr_st *x1 = (struct ipr_st *)a;
    struct ipr_st *x2 = (struct ipr_st *)b;
    if ((*x1).index < (*x2).index)
        return -1;
    else if ((*x1).index > (*x2).index)
        return 1;
    else
        return 0;
}