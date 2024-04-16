#include "fd.h"

void chiralpot(atm_st *atm, long ntot){
    int ndiv = 360;
    double theta;
    double PI = 3.1415926
    double z[3] = {0.0,0.0,1.0}
    double *rx, *ry, *rz;
    rx = (double*) calloc(ntot, sizeof(double));
    ry = (double*) calloc(ntot, sizeof(double));
    rz = (double*) calloc(ntot, sizeof(double));
    
    FILE *pf;
    printf("\n\nCHIRALPOT TEST\n\n")
    pf = fopen("conf.dat", "r")
    for (i = 0; i < ntot; i++){
        fscanf(pf, "%s %lf %lf %lf", atm[i].atyp, &rx[i], &ry[i], &rz[i]);
        printf("%s %lf %lf %lf\n", atm[i].atyp, &rx[i], &ry[i], &rz[i]);
    }


    for (i = 0; i < ndiv; i++){
        theta = 2*PI/ndiv * i
        printf("%g\n", theta)
    }

    printf("END CHIRALPOT TEST")

}