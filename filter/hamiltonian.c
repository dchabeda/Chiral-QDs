/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

void hamiltonian(zomplex *phi,zomplex *psi,double *potl,double *ksqr,long_st ist,par_st par,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long i;

  //memcpy(&fftwpsi[0],&phi[0],ist.ngrid*sizeof(fftwpsi[0]));
  // printf("\nhamiltonian: phi[146] (init) = %.8g + %.8gi", phi[146].re, phi[146].im);
  memcpy(&psi[0],&phi[0],ist.ngrid*sizeof(psi[0]));
  // printf("\nhamiltonian: psi[146] (init) = %.8g + %.8gi", psi[146].re, psi[146].im);
  kinetic(phi,ksqr,planfw,planbw,fftwpsi,ist);
  for (i = 0; i < ist.ngrid; i++){
    // if (i == 146){printf("\nhamiltonian: phi[146] (kinetic) = %.8g + %.8gi", phi[i].re, phi[i].im); }
    phi[i].re += (potl[i] * psi[i].re);
    phi[i].im += (potl[i] * psi[i].im);
    // if (i == 146){printf("\nhamiltonian: potl[146] * psi[146] = %.8g + %.8gi", (potl[i] * psi[i].re), (potl[i] * psi[i].im));}
    // if (i == 146){printf("\nhamiltonian: phi[146] (total) = %.8g + %.8gi\n", phi[i].re, phi[i].im); }
  }
  return;
}

/*****************************************************************************/

void kinetic(zomplex *psi,double *ksqr,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long_st ist)
{
  long j, flags=0;
  // FILE *pf;

  // pf = fopen("ksqr.dat", "w");
  memcpy(&fftwpsi[0],&psi[0],ist.ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  // //fftwnd_one(planfw,psi,fftwpsi);
  // printf("\nksqr[146] = %.6g\n", ksqr[146]);
  for (j = 0; j < ist.ngrid; j++){
    fftwpsi[j][0] *= ksqr[j];
    fftwpsi[j][1] *= ksqr[j];
    // fprintf(pf, "%ld %.8g\n", j, ksqr[j]);
  }
  fftw_execute(planbw);
  //fftwnd_one(planbw,fftwpsi,psi);
  // fclose(pf);
  memcpy(&psi[0],&fftwpsi[0],ist.ngrid*sizeof(psi[0]));
  return;
}

/*****************************************************************************/

