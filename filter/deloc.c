#include "fd.h"

void getOverlap(double *psitot, double *newpsitot, double *eval, double *vx, double *vy, double *vz, double *rx, double *ry, double *rz, par_st *par, long_st *ist){
    /*** Computes the overlap of the true wavefunction with a "perfectly delocalized" wavefunction ***/
    FILE *pdeloc, *pdelocc;
    long jms, jgrid, jx, jy, jz, jyz, count, iatom;
    double sum, dis;
    double *rho, *gausspsi;
    deloc_st *deloc, *delocCut;
    omp_set_dynamic(0);
    omp_set_num_threads(ist->nthreads);

    if ((deloc = calloc(ist->mstot, sizeof(deloc_st)))==NULL){
    printf("\nOUT OF MEMORY: delocCut\n");fflush(0);
    nerror("delocCut");
    }
    if ((delocCut = calloc(ist->delocCut, sizeof(deloc_st)))==NULL){
    printf("\nOUT OF MEMORY: delocCut\n");fflush(0);
    nerror("delocCut");
    }
    if ((gausspsi = calloc(ist->ngrid, sizeof(double)))==NULL){
    printf("\nOUT OF MEMORY: Gauss Psi\n");fflush(0);
    nerror("gausspsi");
    }
    if ((rho = calloc(ist->ngrid,sizeof(double)))==NULL){
    printf("\nOUT OF MEMORY: Gauss Psi\n");fflush(0);
    nerror("gausspsi");
    }

    //Define an "ideally delocalized" wavefunction: Gaussian centered around each atom.

    for (jz = 0; jz < ist->nz; jz++){
        for (jy = 0; jy < ist->ny; jy++){
            jyz = ist->nx * (ist->ny * jz + jy);
            for (jx = 0; jx < ist->nx; jx++){
                sum = 0.0;
                jgrid = jyz + jx;
                for (iatom = 0; iatom < ist->natom; iatom++){
                    dis = sqrt(sqr(vx[jx] - rx[iatom]) + sqr(vy[jy] - ry[iatom]) + sqr(vz[jz] - rz[iatom]));
                    sum += exp(-(0.03125 * sqr(dis) ));
                    if(jgrid == 0){
                        printf("C %g %g %g\n", rx[iatom], ry[iatom], rz[iatom] );
                    } 
                }
                gausspsi[jgrid] = sum;
                if (jgrid % 2700 == 0.0) printf("\nDistance for grid point x = %g y = %g z = %g: %.8f\n",vx[jx],vy[jy],vz[jz], dis);fflush(0); 
                if (jgrid % 2700 == 0.0) printf("gausspsi[%ld] at atom %ld = %.12f\n", jgrid, iatom, gausspsi[jgrid]); fflush(0);
                
            }
        }
    }

    //Normalize the wavefunction
    long ms = 1;
    normalize_all(gausspsi, par->dv, ms, ist->ngrid, ist->nthreads);
    for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
        rho[jgrid] = gausspsi[jgrid]; // For writing the wavefunction cube file
        if (jgrid % 2500 == 0.0) printf("\nNorm gausspsi[jgrid]: %.12f\n", gausspsi[jgrid]);fflush(0);      
    }
    printf("Visualizing the gaussian wavefunction\n"); fflush(0);
    writeCubeFile(gausspsi, *par, *ist, "gausspsi.cube");


    //Compute the wavefuntion overlaps
#pragma omp parallel for private(jms)
    for (jms = 0; jms < ist->mstot; jms++){
        // For each quasiparticle state
        sum = 0.0;
        for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
            // Sum up the deloc function value at all the grid points
            sum += gausspsi[jgrid] * psitot[jms*ist->ngrid+jgrid];
        }
        deloc[jms].val = sum * par->dv;
        deloc[jms].index = jms;
    }
    
    printf("Sorting psi by wavefunction overlap:\n"); fflush(0);
    
    qsort(deloc, ist->mstot, sizeof(deloc_st), compare_vals);

    pdeloc = fopen("deloc.dat", "w");
    for (jms = 0; jms < ist->mstot; jms++){
    fprintf(pdeloc, "%ld %.8f\n", jms, deloc[jms].val); fflush(0);
    }

    printf("\tjms    deloc    idx\n"); fflush(0);
    for (jms = 0; jms < ist->mstot; jms++){
        if (jms == 0){
            printf("This is the delocalized index: %ld\n", deloc[jms].index);
            for (jgrid = 0; jgrid<ist->ngrid; jgrid++){
                rho[jgrid] = psitot[deloc[jms].index*ist->ngrid + jgrid];
            }
            writeCubeFile(rho, *par, *ist, "deloc-wfn.cube");
        }
        if (jms == ist->delocCut-1){
            printf("This is the localized index: %ld\n", deloc[jms].index); fflush(0);
            for (jgrid = 0; jgrid<ist->ngrid; jgrid++){
                rho[jgrid] = psitot[deloc[jms].index*ist->ngrid + jgrid];
            }
            writeCubeFile(rho, *par, *ist, "loc-wfn.cube");
        }
        if (jms % 20 == 0){
            printf("\t%ld %.6f %ld\n", jms, deloc[jms].val, deloc[jms].index );
        }
    }
    
    for (jms = 0; jms < ist->delocCut; jms++){
        delocCut[jms].val = deloc[jms].val;
        delocCut[jms].index = deloc[jms].index;
    }
    printf("Sorting the delocalized psi by index:\n");
    qsort(delocCut, ist->delocCut, sizeof(deloc_st), compare_index);

    printf("The deloc cutoff = %ld\n", ist->delocCut); fflush(0);
    printf("deloc of max delocalization = %g\n", (double) 1/ist->ngrid);fflush(0);
    printf("Larger values correspond to less delocalized states.\n"); fflush(0);

    
    pdelocc = fopen("delocCut.dat", "w");
    for (jms = 0; jms < ist->delocCut; jms++ ){
        //printf("jms: %ld\n", jms); fflush(0);
        for (jgrid = 0; jgrid < ist->ngrid; jgrid++){
            //printf("index: %ld\n", delocCut[jms].index);fflush(0);
            newpsitot[jms*ist->ngrid + jgrid] = psitot[delocCut[jms].index*ist->ngrid + jgrid];
        }
        eval[jms] = eval[delocCut[jms].index];
        fprintf(pdelocc, "%ld %ld %.8f\n", jms, delocCut[jms].index, delocCut[jms].val);fflush(0);
    }
    printf("\nDelocalized states obtained!\n");fflush(0);
    //printf("New mstot = %ld\t Nstates removed = %ld\n", i, ist->mstot - i);
    fclose(pdeloc); fclose(pdelocc); free(delocCut); free(rho);
}

// int compare_function(const void *a,const void *b) {
// double *x = (double *) a;
// double *y = (double *) b;
// if (*x < *y) return -1;
// else if (*x > *y) return 1; return 0;
// }

int compare_vals(const void *a, const void *b)
{
    struct deloc_st *x1 = (struct deloc_st *)a;
    struct deloc_st *x2 = (struct deloc_st *)b;
    if ((*x1).val > (*x2).val)
        return -1;
    else if ((*x1).val < (*x2).val)
        return 1;
    else
        return 0;
}

int compare_index(const void *a, const void *b)
{
    struct deloc_st *x1 = (struct deloc_st *)a;
    struct deloc_st *x2 = (struct deloc_st *)b;
    if ((*x1).index < (*x2).index)
        return -1;
    else if ((*x1).index > (*x2).index)
        return 1;
    else
        return 0;
}