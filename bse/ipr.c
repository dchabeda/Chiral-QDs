#include "fd.h"

void IPR(double *psi, double *ipr, long nstates, par_st par, long_st ist){
    /*** Computes the inverse participation ratio: sum_n |psi(x_n)|^4 ***/
    FILE *pipr, *piprc;
    long jms, jgrid, count,i,nhomo, lumoctr, *nlumo;
    double sum, iprcut, *psipr, *rho;

    iprcut = (double) 5/ist.ngrid; //This cutoff selects states that are at least 20% delocalized
    

//    omp_set_dynamic(0);
//    omp_set_num_threads(ist.nthreads);
//#pragma omp parallel for private(jms)
    for (jms = 0, count = 0; jms < nstates; jms++){
        // For each quasiparticle state
        sum = 0.0;
        for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
            // Sum up the IPR function value at all the grid points
            sum += sqr(psi[jms*ist.ngrid+jgrid])*sqr(psi[jms*ist.ngrid+jgrid]);
        }
        ipr[jms] = sum;
        if (sum < iprcut){
            count++;
        }
    }
    pipr = fopen("IPR.dat", "w");
    for (jms = 0; jms < nstates; jms++){
    fprintf(pipr, "%ld %.8f\n", jms, ipr[jms]);
    }
    
    psipr = (double*) calloc(ist.ngrid*count, sizeof(double));
    nlumo = (long*) calloc(ist.totallumo, sizeof(long));
    rho = (double*) calloc(ist.ngrid, sizeof(double));

    piprc = fopen("IPRcut.dat", "w");
    for (jms = 0, i = 0, nhomo = 0, lumoctr = 0; jms < nstates; jms++ ){
        if(ipr[jms] < iprcut){
            if (jms <= ist.nhomo){nhomo++;}
            if (jms > ist.nhomo){nlumo[lumoctr] = i; lumoctr++;}

            for (jgrid = 0; jgrid<ist.ngrid; jgrid++){ 
                psipr[i*ist.ngrid + jgrid] = psi[jms*ist.ngrid + jgrid];
                ipr[i] = ipr[jms];
                fprintf(piprc, "%ld %.8f\n", i, ipr[i]);
            }
            i++;
        }
    }
    printf("\nDelocalized states obtained!\nNew mstot = %ld\t Nstates removed = %ld\n", i, nstates - i);
    
    
    /*** Write cube files ***/
    int ncubes = 4;
    char str[50];
    printf("Delocalized nhomo = %ld nlumo = %ld\n", nhomo, nlumo[0]);

    for (i=0; i < ncubes; i++){
        sprintf(str,"deloc-homo-%ld.cube",i);
        for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
        rho[jgrid] = psipr[(nhomo-i)*ist.ngrid+jgrid];
        }
        writeCubeFile(rho,par,ist,str);
        
        sprintf(str,"deloc-lumo+%ld.cube",i);
        for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
        rho[jgrid] = psipr[(nlumo[0]+i)*ist.ngrid+jgrid];
        }
        writeCubeFile(rho,par,ist,str);
    }


    free(rho); free(nlumo); free(psipr);
    fclose(pipr); fclose(piprc);
}


