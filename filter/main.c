/*****************************************************************************/
//
//
//
/*****************************************************************************/

#include "fd.h"

/*****************************************************************************/

int main(int argc, char *argv[]) {
  FILE *ppsi;  zomplex *psi, *phi, *an; long idum;
  fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;
  par_st par;   long_st  ist; xyz_st dipole; atm_st *atm;
  double *eval, *ksqr, *vx, *vy, *vz, *zn, *potl, *rho;
  double *psitot, *el, *sige, *rx, *ry, *rz, tci, twi;
  long jgrid, jms, jns, mssav, mstot, flags=0, tid;
  printf("\nRUNNING PROGRAM: FILTER DIAGONALIZATION\n");
  /*** read initial setup from input.par ***/
  writeCurrentTime(stdout);
  writeSeparation(stdout);
  printf("INITIALIZING JOB SPECIFICATIONS\n\n");
  init_size(argc,argv,&par,&ist);
  
  /*************************************************************************/
  /*** allocating memory ***/

  fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid);
  //if ((fftwpsi  = (fftw_complex*)calloc(ist.ngrid,sizeof(fftw_complex)))==NULL)nerror("fftwpsi");
  if ((psi   = (zomplex *)calloc(ist.ngrid,sizeof(*psi)))==NULL)nerror("psi");
  if ((phi   = (zomplex *)calloc(ist.ngrid,sizeof(*phi)))==NULL)nerror("phi");
  /*** the pseudopotential stored on the grid ***/
  if ((potl  = (double*)calloc(ist.ngrid,sizeof(double)))==NULL)nerror("potl");
  /*** the kinetic energy stored on the grid ***/
  if ((ksqr = (double*)calloc(ist.ngrid,sizeof(double)))==NULL)nerror("ksqr");
  /*** the grid in the x, y, and z directions ***/
  if ((vx = (double*)calloc(ist.nx,sizeof(double)))==NULL)nerror("vx");
  if ((vy = (double*)calloc(ist.ny,sizeof(double)))==NULL)nerror("vy");
  if ((vz = (double*)calloc(ist.nz,sizeof(double)))==NULL)nerror("vz");
  /*** the positions of the atoms in the x, y, and z directions ***/
  if ((rx = (double*)calloc(ist.natom,sizeof(double)))==NULL)nerror("rx");
  if ((ry = (double*)calloc(ist.natom,sizeof(double)))==NULL)nerror("ry");
  if ((rz = (double*)calloc(ist.natom,sizeof(double)))==NULL)nerror("rz");
  if ((atm = (atm_st*)calloc(ist.natom,sizeof(atm_st)))==NULL)nerror("atm");

  /*** all wavefunctions in the energy range  ***/
  if ((psitot = (double*)calloc(ist.mstot*ist.ngrid,sizeof(double)))==NULL)nerror("psitot");
  /*** the filtered energies ***/
  if ((eval  = (double*)calloc(ist.mstot,sizeof(double)))==NULL)nerror("eval");
  if ((sige = (double*)calloc(ist.mstot,sizeof(double)))==NULL)nerror("sige");
  /*** the energy of the filters ***/
  if ((el   = (double*)calloc(ist.ms,sizeof(double)))==NULL)nerror("el");
  if ((rho  = (double *) calloc(ist.ngrid,sizeof(double)))==NULL) nerror("rho");
  /**************************************************************************/
  /*** initialize the pseudopotential, the kinetic energy, ***/
  /*** the fourier transorm, the grid, and the energy window ***/
  writeSeparation(stdout);
  printf("INITIALIZING GRID AND POTENTIAL\n\n");fflush(0);
  init(potl,vx,vy,vz,ksqr,rx,ry,rz,atm,&par,el,&ist,&planfw,&planbw,fftwpsi);
  
  /*** Write local potential cube file for visualization***/
  writeCubeFile(potl, par, ist, "localPot.cube");

  /*** initialization for the fast Fourier transform ***/
  planfw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_FORWARD,flags);
  planbw = fftw_plan_dft_3d(ist.nz,ist.ny,ist.nx,fftwpsi,fftwpsi,FFTW_BACKWARD,flags);
  
  /**************************************************************************/
  /*** calcualte the energy range of the hamitonian ***/
  writeSeparation(stdout);
  printf("CALCULATING ENERGY RANGE OF HAMILTONIAN\n\n");fflush(0);
  tci = (double)clock(); twi = (double)time(NULL);
  get_energy_range(psi,phi,potl,vx,vy,vz,ksqr,&par,ist,planfw,planbw,fftwpsi);
  printf("done calculate energy range, CPU time (sec) %g, wall run time (sec) %g\n",
    ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); fflush(0);
  
  writeSeparation(stdout);
  printf("FILTERING\n\n");fflush(0);
  /**************************************************************************/
  /*** set parameters for the newton interpolation ***/
  par.dt = sqr((double)(ist.nc) / (2.5 * par.dE));
  printf ("nChebyshev = %ld dt = %g dE = %g\n",ist.nc,par.dt,par.dE); fflush(0);
  an = (zomplex*)calloc(ist.nc*ist.ms,sizeof(zomplex));
  zn = (double*)calloc(ist.nc,sizeof(double));
  coefficient(an,zn,el,par,ist);

  /**************************************************************************/
  /*** start filtering loop.  we run over ns cycles and calculate ***/
  /*** ms filtered states at each cycle ***/
  Randomize();  idum = -random();
  printf("seed = %ld\n",idum);  fflush(0);

  for (jns = 0; jns < ist.ns; jns++) {
    init_psi(psi,ist,par,&idum);
    for (jms = 0; jms < ist.ms; jms++) {
      for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
        psitot[jns*ist.ms*ist.ngrid+ist.ngrid*jms+jgrid] = psi[jgrid].re;
	  }
	}
  }
  
  tci = (double)clock(); twi = (double)time(NULL);
  omp_set_dynamic(0);
  omp_set_num_threads(ist.nthreads);
#pragma omp parallel for private(jns)
  for (jns = 0; jns < ist.ns; jns++){
    tid = omp_get_thread_num();
    printf("Iteration %ld\n", jns);	
    filtering(&psitot[jns*ist.ms*ist.ngrid],potl,ksqr,an,zn,el,ist,par,tid,jns,&idum);
  } 
  printf("done calculating filter, CPU time (sec) %g, wall run time (sec) %g\n",
    ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); fflush(0);
  
  //if (par->saveChk == 1) printf("\n\nChk1\n\n");
  /*************************************************************************/
  /*** read all filtered states ***/
  /*ppsi = fopen("psi-filt.dat" , "r");
  fread (psitot,sizeof(double),ist.mstot*ist.ngrid,ppsi);
  fclose(ppsi);*/

  /*** orthogonalize and normalize the filtered states using an svd routine ***/
  writeSeparation(stdout);
  printf("ORTHGONALIZING FILTERED STATES\n\n");fflush(0);
  tci = (double)clock(); twi = (double)time(NULL);
  ist.mstot = portho(psitot,par.dv,ist);
  printf("mstot ortho = %ld\n", ist.mstot); fflush(0);
  normalize_all(psitot,par.dv,ist.mstot,ist.ngrid,ist.nthreads);
  printf("done calculating ortho, CPU time (sec) %g, wall run time (sec) %g\n",
    ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); fflush(0);
    
  /***********************************************************************/
  /*** diagonalize the hamiltonian in the subspace spanned by the ***/
  /*** orthogonal filtered states. generate the eigenstates of the ***/
  /*** hamitonian within the desired energy reange ***/
  writeSeparation(stdout);
  printf("DIAGONALIZING HAMILTONIAN\n\n");fflush(0);
  tci = (double)clock(); twi = (double)time(NULL);
  Hmatreal(psi, phi, psitot, potl, ksqr, eval,ist, par, planfw, planbw, fftwpsi);
  normalize_all(psitot, par.dv, ist.mstot, ist.ngrid, ist.nthreads);
  printf("done calculating Hmat, CPU time (sec) %g, wall run time (sec) %g\n",
    ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi);
  fflush(0);

  // ppsi = fopen("eval.dat", "w");
  // for (jms = 0; jms < ist.mstot; jms++) { 
  //   fprintf(ppsi, "%ld %.16g\n", jms, eval[jms]); 
  // }
  // fclose(ppsi);
  /*** If getDelocFlag is one, only output states with an inverse participation ratio
   * Greater than some threshold, par->iprcut ***/

  if (par.getDelocFlag == 1){
    long mstot;
    double *selectpsitot;
    long VBmindex, CBmaxdex, nstates;
    writeSeparation(stdout);
    tci = (double)clock(); twi = (double)time(NULL);
    printf("\nCalculating delocalization as <dr^2>...\n"); fflush(0);
    
    //Select only states that are within the targeted energy range. No false states
    for (jms = 0; jms < ist.mstot; jms++){
      if (eval[jms] < (double)(par.VBmin)){
        //printf("This is mindex eval[%ld] = %g\n",jms,eval[jms]);
        VBmindex = jms;
      }
    }
    for (jms = 0; jms < ist.mstot; jms++){
      if (eval[jms] > (double)(par.CBmax + 0.05)){
        //printf("This is maxdex eval[%ld] = %g\n",jms,eval[jms]);
        CBmaxdex = jms;
        break;
      }
    }
    printf("The VBmindex = %ld CBmaxdex = %ld\n", VBmindex, CBmaxdex); fflush(0);
    nstates = CBmaxdex - VBmindex;
    printf("There are %ld states within the energy range\n", nstates); fflush(0);

    if ((selectpsitot = calloc(nstates*ist.ngrid, sizeof(double)))==NULL){
      printf("\nOUT OF MEMORY: DELOC NEW PSI\n");fflush(0);
      nerror("selectpsitot");
    }

    memcpy(&selectpsitot[0],&psitot[VBmindex*ist.ngrid],ist.ngrid*nstates*sizeof(double));
    // for (jms = 0; jms < nstates; jms++){ //I use jms to index the original psi and jns to index the new psi.
    //   if (jms % 20 == 0){
    //       printf("This is jms: %ld\n", jms);
    //       printf("The selectpsitot[%ld] = %g\n", jms, selectpsitot[jms*ist.ngrid + (long)(20000)]);
          
    //     }
    // }

    free(psitot);
    if((psitot = calloc(nstates * ist.ngrid, sizeof(double))) == NULL){
      printf("OUT OF MEMORY: psitot in deloc\n");
      nerror("psitot in deloc");
    }
    
    getDeloc(selectpsitot, psitot, nstates, vx, vy, vz, eval, par, &ist);
    
    ist.mstot = ist.delocCut;
    printf("This is ist.mstot %ld\n", ist.mstot);fflush(0);
    //memcpy(&psitot[0], &selectpsitot[0], ist.ngrid*ist.mstot*sizeof(double));

    printf("done calculating delocalization, CPU time (sec) %g, wall run time (sec) %g\n",
    ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi); fflush(0);
  }

  writeSeparation(stdout);
  printf("WRITING PSI, EVAL, CUBE FILES\n\n");fflush(0);
  /*** write the eigenstates to a file ***/
  ppsi = fopen("psi.dat", "w");
  fwrite(psitot, sizeof(double), ist.mstot*ist.ngrid, ppsi);
  fclose(ppsi);
  printf("Done writing psi\n"); fflush(0);
  //if (par->saveChk == 1) printf("\n\nChk2\n\n");
  /*** calculate the standard deviation of these states ***/
  /*** this is used to check if there are ghost states ***/
  calc_sigma_E(psi, phi, psitot, potl, ksqr, sige, ist, par, planfw, planbw, fftwpsi);
  printf("Done calc sigE\n"); fflush(0);
  /*** write the eigenvalues and the standard deviation of the eigenstates ***/
  ppsi = fopen("eval.dat", "w");
  for (jms = 0; jms < ist.mstot; jms++) { 
    fprintf(ppsi, "%ld %.16g %g\n", jms, eval[jms], sige[jms]); 
  }
  fclose(ppsi);  
  printf("Done printing eval.\n"); fflush(0);
  /*** Write homo and lumo cube files ***/
  long i, a, ieof, nval;
  int ncubes = 4;
  char str[50];
  double evalloc, deloc;
  FILE *pf;

  par.deps = 0.0001;
  ist.nhomo = ist.nlumo = 0;
  pf = fopen("eval.dat" , "r");
  for (i = ieof = 0; ieof != EOF; i++){
    ieof = fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
    if (deloc < par.deps && evalloc < par.fermiEnergy) ist.nhomo = i;
  }
  fclose(pf);

  nval = i - 1;
  pf = fopen("eval.dat" , "r");
  for (i = 0; i <= ist.nhomo; i++) fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
  for (i = ist.nhomo+1; i < nval; i++) {
    fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
    if (deloc < par.deps) {
      ist.nlumo = i;
      break;
    }
  }
  fclose(pf);

  printf("nhomo = %ld nlumo = %ld\n", ist.nhomo, ist.nlumo);fflush(0);
  tci = (double)clock(); twi = (double)time(NULL);

  for (i=0; i < ncubes; i++){
    sprintf(str,"homo-%ld.cube",i);
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
    rho[jgrid] = psitot[(ist.nhomo-i)*ist.ngrid+jgrid];
    }
    writeCubeFile(rho,par,ist,str);
    printf("Done calcing cube %ld\n", i); fflush(0);
    sprintf(str,"lumo+%ld.cube",i);
    for (jgrid = 0; jgrid < ist.ngrid; jgrid++){
    rho[jgrid] = psitot[(ist.nlumo+i)*ist.ngrid+jgrid];
    }
    writeCubeFile(rho,par,ist,str);
  }
  

  free(rho);
  printf("done calculating cubes, CPU time (sec) %g, wall run time (sec) %g\n",
    ((double)clock()-tci)/(double)(CLOCKS_PER_SEC), (double)time(NULL)-twi);fflush(0);
  
  /*************************************************************************/
  /*** free memeory ***/
  free(psitot); free(phi); free(psi); free(potl); free(eval); 
  free(el); free(ksqr); free(vx); free(vy); free(an); free(zn);
  free(vz); free(rx); free(ry); free(rz); free(sige); free(atm);

  fftw_destroy_plan(planfw);
  fftw_destroy_plan(planbw);
  fftw_free(fftwpsi);

  /*************************************************************************/
  /*** Print out the time that the computation finished ***/
  writeSeparation(stdout);
  printf("EXITING PROGRAM: FILTER DIAGONALIZATION\n");fflush(0);
  writeCurrentTime(stdout);

  /*************************************************************************/
  /*** Exit the program ***/
  exit(0);
}
