/****************************************************************************/
/* This file does the printing for the program */

#include "fd.h"

/****************************************************************************/
// prints out current time to stdout 
 
void writeCurrentTime(FILE *pf) {
  time_t startTime;

  startTime = time(NULL);
  fprintf(pf, ctime(&startTime));

  return;
}

/****************************************************************************/

void writeSeparation(FILE *pf) {
  fprintf(pf, "\n******************************************************************************\n\n");
  
  return;
}

/****************************************************************************/
/****************************************************************************/
void writeCubeFile(double *rho, par_st par, long_st ist, char *fileName) {

    FILE *pf, *pConfFile;
    long iGrid, iX, iY, iZ, iYZ, nAtoms, atomType;
    double x, y, z;
    char line[80], atomSymbol[10];

    pConfFile = fopen("conf.dat", "r");
    fscanf(pConfFile, "%ld", &nAtoms);
    pf = fopen(fileName, "w");
    fprintf(pf, "CUBE FILE\n");
    fprintf(pf, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");

    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", nAtoms, par.xmin, par.ymin, par.zmin);
    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.nx, par.dx, 0.0, 0.0);
    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.ny, 0.0, par.dy, 0.0);
    fprintf(pf, "%5li%12.6f%12.6f%12.6f\n", ist.nz, 0.0, 0.0, par.dz);

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
    
    for (iX = 0; iX < ist.nx; iX++) {
        for (iY = 0; iY < ist.ny; iY++) {
            for (iZ = 0; iZ < ist.nz; iZ++) {
                iYZ = ist.nx * (ist.ny * iZ + iY);
                iGrid = iYZ + iX;
                fprintf(pf, "%12.5f ", rho[iGrid]);
                if (iZ % 6 == 5) {
                    fprintf(pf, "\n");
                }
            }
            fprintf(pf, "\n");
        }
    }
    /***for (iZ = 0; iZ < ist.nz; iZ++) {
        for (iY = 0; iY < ist.ny; iY++) {
            iYZ = ist.nx * (ist.ny * iZ + iY);
            for (iX = 0; iX < ist.nx; iX++) {
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
