#include "fd.h"

void getDeloc(double *psitot, double *newpsitot, long nstates, double *vx, double *vy, double *vz, double *eval, par_st par, long_st *ist){
    /*** Computes the average square deviation of wavefunctions: <dr^2> ***/
    FILE *pdr2, *pdr2cut, *pchi;
    long jms, jgrid, jx, jy, jz, jyz, count, iatom;
    double sumx, sumy, sumz, sum2, tmp, avg2, dis, maxdis_1;
    double cutoff;
    double *rho;
    dr_st *dr2, *dr2cut;

    count = 0;
    dr2 = calloc(nstates, sizeof(dr_st));
    rho = calloc(ist->ngrid,sizeof(double));
    cutoff = 0.15*(par.xmax-par.xmin) + 0.15*(par.ymax-par.ymin) + 0.15*(par.zmax-par.zmin);
    printf("This is the cutoff we are using for delocalization %g Bohr\n", cutoff);

    omp_set_dynamic(0);
    omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(jms)
    for (jms = 0; jms < nstates; jms++){
        // For each quasiparticle state
        //First compute the average of r
        sumx = 0.0;
        sumy = 0.0;
        sumz = 0.0;
        sum2 = 0.0;
        for (jz = 0; jz < ist->nz; jz++){
            for (jy = 0; jy < ist->ny; jy++){
                jyz = ist->nx * (ist->ny * jz + jy);
                for (jx = 0; jx < ist->nx; jx++){
                    jgrid = jyz + jx;
                    // Sum up the dr2 function value at all the grid points
                    tmp = sqr(psitot[jms*ist->ngrid+jgrid]);
                    
                    sumx += tmp*vx[jx]*par.dx;
                    sumy += tmp*vy[jy]*par.dy;
                    sumz += tmp*vz[jz]*par.dz; 
                    
                    sum2 += tmp*sqr(vx[jx])*par.dx;
                    sum2 += tmp*sqr(vy[jy])*par.dy;
                    sum2 += tmp*sqr(vz[jz])*par.dz;

                    //if (jgrid % 20000 == 0) printf("For state %ld tmp: %g\n", jms, tmp);
                    //if (jgrid % 20000 == 0) printf("For state %ld tmp: %g\n", jms, tmp);
                }
            }
        }
        dr2[jms].val = sum2 - (sqr(sumx)+sqr(sumy)+sqr(sumz)); // This is <dr^2>
        dr2[jms].index = jms;
        if (jms%20 == 0){
        printf("This is sum2 = %g\n", sum2);
        printf("This is sum1^2 = %g\n", (sqr(sumx)+sqr(sumy)+sqr(sumz)));
        printf("This is dr2[%ld] = %g\n", jms, dr2[jms].val);
        }
        if (dr2[jms].val > cutoff){
            count++;
        }
    }

    ist->delocCut = count;
    
    printf("Sorting psi by delocalization:\n"); fflush(0);
    qsort(dr2, nstates, sizeof(dr_st), comp_vals);

    pdr2 = fopen("dr2.dat", "w");
    for (jms = 0; jms < nstates; jms++){
    fprintf(pdr2, "%ld %.8f\n", jms, dr2[jms].val);
    }

    printf("\tjms    dr2    idx\n");
    for (jms = 0; jms < nstates; jms++){
        if (jms == 0){
            printf("This is the delocalized index: %ld\n", dr2[jms].index);
            printf("\t%ld %.6f %ld\n", jms, dr2[jms].val, dr2[jms].index );
            for (jgrid = 0; jgrid<ist->ngrid; jgrid++){
                rho[jgrid] = psitot[dr2[jms].index*ist->ngrid + jgrid];
            }
            writeCubeFile(rho, par, *ist, "deloc-wfn.cube");
        }
        if (jms == count-1){
            printf("This is the localized index: %ld\n", dr2[jms].index);
            printf("\t%ld %.6f %ld\n", jms, dr2[jms].val, dr2[jms].index );
            for (jgrid = 0; jgrid<ist->ngrid; jgrid++){
                rho[jgrid] = psitot[dr2[jms].index*ist->ngrid + jgrid];
            }
            writeCubeFile(rho, par, *ist, "loc-wfn.cube");
        }
        if (jms % 20 == 0){
            printf("\t%ld %.6f %ld\n", jms, dr2[jms].val, dr2[jms].index );
        }
    }

    //Allocate memory for the delocalized wavefunctions and their deloc parameters
    // if ((*newpsitot = calloc(count*ist->ngrid, sizeof(double))) == NULL){
    //     printf("OUT OF MEMORY: newpsitot\n");
    //     nerror("newpsitot");
    //}
    if (( dr2cut = calloc(count, sizeof(dr_st))) == NULL){
        printf("OUT OF MEMORY: dr2cut\n");
        nerror("dr2cut");
    }

    //Populate the dr2cut array with the Ncount parameters of highest delocalization
    for (jms = 0; jms < count; jms++){
        dr2cut[jms].val = dr2[jms].val;
        dr2cut[jms].index = dr2[jms].index;
    }

    printf("Sorting the delocalized psi by index:\n");
    qsort(dr2cut, count, sizeof(dr_st), comp_index); //Put the wavefunctions back in energy order.

    printf("Number of delocalized states = %ld\n", count); fflush(0);
    
    pdr2cut = fopen("dr2cut.dat", "w");
    for (jms = 0; jms < count; jms++ ){
        //printf("jms: %ld\n", jms); fflush(0);
        for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
            //if (jgrid % 20000 == 0) printf("Made it to populate delocalized states jms = %ld jgrid = %ld", jms, jgrid);fflush(0);
            //printf("index: %ld\n", dr2cut[jms].index);fflush(0);
            newpsitot[jms*ist->ngrid + jgrid] = psitot[dr2cut[jms].index*ist->ngrid + jgrid];
        }
        eval[jms] = eval[dr2cut[jms].index];
        fprintf(pdr2cut, "%ld %ld %.8f\n", jms, dr2cut[jms].index, dr2cut[jms].val);fflush(0);
    }
    printf("\nDelocalized states obtained!\n");fflush(0);
    //printf("New mstot = %ld\t Nstates removed = %ld\n", i, ist->mstot - i);
    fclose(pdr2); fclose(pdr2cut); free(dr2cut); free(dr2); free(rho); 

}

int comp_vals(const void *a, const void *b)
{
    struct dr_st *x1 = (struct dr_st *)a;
    struct dr_st *x2 = (struct dr_st *)b;
    if ((*x1).val > (*x2).val)
        return -1;
    else if ((*x1).val < (*x2).val)
        return 1;
    else
        return 0;
}

int comp_index(const void *a, const void *b)
{
    struct dr_st *x1 = (struct dr_st *)a;
    struct dr_st *x2 = (struct dr_st *)b;
    if ((*x1).index < (*x2).index)
        return -1;
    else if ((*x1).index > (*x2).index)
        return 1;
    else
        return 0;
}