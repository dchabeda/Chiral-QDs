#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

void writeCubeFile(double *rho, char *fileName, long nx, long ny, long nz, double dx, double dy, double dz);

int main(int argc, char *argv[]){
    //usage: writeCube.x nx dx ny dy nz dz nhomo nlumo ncubeh ncubel
    FILE *pf;
    double *rho;
    long i;
    long nx = atol(argv[1]); double dx = atof(argv[2]);
    long ny = atol(argv[3]); double dy = atof(argv[4]);
    long nz = atol(argv[5]); double dz = atof(argv[6]);
    long ngrid = nx*ny*nz;
    long nhomo = atol(argv[7]);
    long nlumo = atol(argv[8]);
    long ncubeh = atol(argv[9]);
    long ncubel = atol(argv[10]);
    char str[25];

    rho = (double*) calloc(ngrid,sizeof(double));

    printf("nhomo = %ld nlumo = %ld\n", nhomo, nlumo);fflush(0);
    pf = fopen("psi.dat", "r");

    for (i=0; i < ncubeh; i++){
        sprintf(str,"homo-%ld.cube",i);
        fseek(pf,(nhomo-i-1)*ngrid*sizeof(double), SEEK_SET);
        fread(rho, ngrid, sizeof(double), pf);
        writeCubeFile(rho, str, nx, ny, nz, dx, dy, dz);
    }
    for (i=0; i < ncubel; i++){
        sprintf(str,"lumo+%ld.cube",i);
        fseek(pf,(nlumo+i)*ngrid*sizeof(double), SEEK_SET);
        fread(rho, ngrid, sizeof(double), pf);
        writeCubeFile(rho, str, nx, ny, nz, dx, dy, dz);
        printf("Done calcing cube %ld\n", i); fflush(0);
    }
    return 0;    
}

void writeCubeFile(double *rho, char *fileName, long nx, long ny, long nz, double dx, double dy, double dz){

    FILE *pf, *pConfFile;
    long iGrid, iX, iY, iZ, iYZ, nAtoms, atomType;
    double x, y, z;
    char line[80], atomSymbol[10];

    pConfFile = fopen("conf.dat", "r");
    fscanf(pConfFile, "%ld", &nAtoms);
    pf = fopen(fileName, "w");
    fprintf(pf, "CUBE FILE\n");
    fprintf(pf, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");

    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nAtoms, -33.0, -35.0, -34.0);
    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nx, dx, 0.0, 0.0);
    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ny, 0.0, dy, 0.0);
    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nz, 0.0, 0.0, dz);

    fgets(line, 80, pConfFile); 
    while(fgets(line, 80, pConfFile) != NULL) {
        sscanf(line, "%3s %lf %lf %lf", &atomSymbol, &x, &y, &z);
        if (! strcmp(atomSymbol, "Cd")) { 
            atomType = 48;
        }
        else if (! strcmp(atomSymbol, "S")) {  
            atomType = 16;
        }
        else if (! strcmp(atomSymbol, "Se")) { 
            atomType = 34;
        }
        else if (! strcmp(atomSymbol, "Zn")) {
            atomType = 30;
        }
        else if (! strcmp(atomSymbol, "Te")) {
            atomType = 52;
        }
        else if (! strcmp(atomSymbol, "C")) {
            atomType = 6;
        }
        else if (! strcmp(atomSymbol, "Si")) {
            atomType = 14;
        }
        else if (! strcmp(atomSymbol, "Cs")) {
            atomType = 55;
        }
        else if (! strcmp(atomSymbol, "Pb")) {
            atomType = 82;
        }
        else if (! strcmp(atomSymbol, "I")) {
            atomType = 53;
        }
        else if (! strcmp(atomSymbol, "In")) {
            atomType = 49;
        }
        else if (! strcmp(atomSymbol, "As")) {
            atomType = 33;
        }
        else if (! strcmp(atomSymbol, "Ga")) {
            atomType = 31;
        }
        else if (! strcmp(atomSymbol, "P")) {
            atomType = 15;
        }
        else if (! strcmp(atomSymbol, "P1")) {
            atomType = 84;
        }
        else if (! strcmp(atomSymbol, "P2")) {
            atomType = 85;
        }
        else if (! strcmp(atomSymbol, "P3")) {
            atomType = 86;
        }
        else if (! strcmp(atomSymbol, "PC5")) {
            atomType = 87;
        }
        else if (! strcmp(atomSymbol, "PC6")) {
            atomType = 88;
        }
        else { 
            atomType = 1; 
        }
        fprintf(pf, "%5li%12.6f%12.6f%12.6f%12.6f\n", atomType, 0.0, x, y, z);
    }
    
    for (iX = 0; iX < nx; iX++) {
        for (iY = 0; iY < ny; iY++) {
            for (iZ = 0; iZ < nz; iZ++) {
                iYZ = nx * (ny * iZ + iY);
                iGrid = iYZ + iX;
                fprintf(pf, "%12.5f ", rho[iGrid]);
                if (iZ % 6 == 5) {
                    fprintf(pf, "\n");
                }
            }
            fprintf(pf, "\n");
        }
    }
    /***for (iZ = 0; iZ < nz; iZ++) {
        for (iY = 0; iY < ny; iY++) {
            iYZ = nx * (ny * iZ + iY);
            for (iX = 0; iX < nx; iX++) {
                iGrid = iYZ + iX;
                fprintf(pf, "%g ", rho[iGrid]);
                if (iX % 6 == 5) {
                    fprintf(pf, "\n");
                }
            }
            fprintf(pf, "\n");
        }
    }***/
    fclose(pConfFile);
    fclose(pf);

    return;
}
