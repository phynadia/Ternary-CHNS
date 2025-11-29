#include "global.h"
/*!
***********************************************************************************************
 * 
 * CHNS-2D-MPI
 *
 * Author: Nadia Bihari Padhan (<nadia_bihari.padhan@tu-dresden.de>).
 *
 * Institute: Institute of Scientific Computing, TU Dresden, Germany.
 *
 * DNSs of coupled Cahn-Hilliard-Navier-Stokes equations using pseudospectral 
 * method.
 *
 * This is the main program which calls multiple functions.
 * 
 * Time integration: etd2rk method.
 *
 * The program involves calculations for Navier-Stokes equations, Binary flows,
 * and ternary flows in 2D. Binary and ternary flows are employed using phase-field method.
 * It is developed in C programming language and MPI-platform.
 * 
***********************************************************************************************/

int main (int argc, char *argv[])
{

    system("mkdir -p spectra");
    system("mkdir -p vort");
    system("mkdir -p init_field");
    system("mkdir -p data");
#if defined BINARY || defined TERNARY
    system("mkdir -p phase_field");
#endif

    int n;

    MPI_Init (&argc, &argv);
    fftw_mpi_init ();
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);


    alloc_local = fftw_mpi_local_size_2d (Nx, Ny/2+1,
                                          MPI_COMM_WORLD,
                                          &local_nx,
                                          &local_x_start);

    
    fu = fopen("data/ken_t.dat","w");

    n = local_nx * Ny2;
    navg = maxiter / nfile;

    factor_x = (2 * pi) / lx;
    factor_y = (2 * pi) / ly;

    kxmax = (Nx/2.0)*factor_x;
    kymax = (Ny/2.0)*factor_y;

    nshell = round(((double)Nx/sqrt(2.0)) + 1.0);

    kxm = round(sqrt(kxmax*kxmax));
    kym = round(sqrt(kymax*kymax));

if(rank == master)
{
  printf("Allocating memories........\n");
}
    int local_size = local_nx * Ny;
    work = fftw_malloc(local_size * sizeof(double)); 

    local_x = Nx / size;
    local_size_x = (Nx / size) + 1;
//printf("local_x = %d\t local_nx = %ld\n", local_x, local_nx);
    u  = fftw_alloc_real (2 * alloc_local);
    ux = fftw_alloc_real (2 * alloc_local);
    uy = fftw_alloc_real (2 * alloc_local);
    wz = fftw_alloc_real (2 * alloc_local);
    tux = fftw_alloc_real (2 * alloc_local);
    tuy = fftw_alloc_real (2 * alloc_local);
    twz = fftw_alloc_real (2 * alloc_local);
    nlz = fftw_alloc_real (2 * alloc_local);    

    cu  = (fftw_complex *) u;
    cux = (fftw_complex *) ux;
    cuy = (fftw_complex *) uy;
    cwz = (fftw_complex *) wz;
    tcux = (fftw_complex *) tux;
    tcuy = (fftw_complex *) tuy;
    tcwz = (fftw_complex *) twz;
    cnlz = (fftw_complex *) nlz;

    cnl_z = fftw_alloc_complex (alloc_local);
    
    fz = (double *) tux;
    cfz = (fftw_complex *) fz;

#ifdef VPM
    chi = fftw_alloc_real (2 * alloc_local);
    sponge = fftw_alloc_real (2 * alloc_local);
#endif

#ifdef BINARY
    sph = (double *) malloc(sizeof(double) * nshell);

    phi  = fftw_alloc_real (2 * alloc_local);
    tphi = fftw_alloc_real (2 * alloc_local);
    nphi = fftw_alloc_real (2 * alloc_local);
    nlphi = fftw_alloc_real (2 * alloc_local);

    cphi  = (fftw_complex *) phi;
    tcphi = (fftw_complex *) tphi;
    ncphi = (fftw_complex *) nphi;
    cnlphi = (fftw_complex *) nlphi;
    cnl_phi = fftw_alloc_complex (alloc_local);

    ph = fftw_malloc(local_size_x * (Ny+2) * sizeof(double));
    aph = fftw_malloc(Ny * sizeof(double));
    bph = fftw_malloc(Ny * sizeof(double));

    coeff = (sigma / epsilon) * lambda;

    fph = fopen("data/volume_droplet.dat","w");
    fe = fopen("data/free_en_t.dat", "w");
    fdef = fopen("data/gamma.dat", "w");

#endif

#ifdef TERNARY
    sph1 = (double *) malloc(sizeof(double) * nshell);
    sph2 = (double *) malloc(sizeof(double) * nshell);
    sp12 = (double *) malloc(sizeof(double) * nshell);
    ph1 = fftw_malloc(local_size_x * (Ny+2) * sizeof(double));
    ph2 = fftw_malloc(local_size_x * (Ny+2) * sizeof(double));
    aph1 = fftw_malloc(Ny * sizeof(double));
    aph2 = fftw_malloc(Ny * sizeof(double));
    bph1 = fftw_malloc(Ny * sizeof(double));
    bph2 = fftw_malloc(Ny * sizeof(double));

    phi1 = fftw_alloc_real (2 * alloc_local);
    phi2 = fftw_alloc_real (2 * alloc_local);
    tphi1 = fftw_alloc_real (2 * alloc_local);
    tphi2 = fftw_alloc_real (2 * alloc_local);
    muF1 = fftw_alloc_real (2 * alloc_local);
    muF2 = fftw_alloc_real (2 * alloc_local);
    nlphi1 = fftw_alloc_real (2 * alloc_local);
    nlphi2 = fftw_alloc_real (2 * alloc_local);

    cphi1 = (fftw_complex *) phi1;
    cphi2 = (fftw_complex *) phi2;
    tcphi1 = (fftw_complex *) tphi1;
    tcphi2 = (fftw_complex *) tphi2;
    cmuF1 = (fftw_complex *) muF1;
    cmuF2 = (fftw_complex *) muF2;
    cnlphi1 = (fftw_complex *) nlphi1;
    cnlphi2 = (fftw_complex *) nlphi2;

    cnl_phi1 = fftw_alloc_complex (alloc_local);
    cnl_phi2 = fftw_alloc_complex (alloc_local);

    fph = fopen("data/volume_phase1.dat","w");
    fe = fopen("data/free_en_t.dat", "w");
    fdef1 = fopen("data/gamma1.dat", "w");
    fdef2 = fopen("data/gamma2.dat", "w");

    gamma1 = (sigma_12 - sigma_23 + sigma_13);
    gamma2 = (sigma_12 - sigma_13 + sigma_23);
    gamma3 = (sigma_23 - sigma_12 + sigma_13);
    del = (6.0 * gamma1 * gamma2 * gamma3);
    del /= (gamma1*gamma2 + gamma1*gamma3 + gamma2*gamma3);

#endif

#ifdef WRITE_ARRAY
    if(rank == 0)
    {
      u1 = (double *) malloc(Nx*(Ny+2) * sizeof(double));
    }
#endif
    kx = (double *) malloc(local_nx  * sizeof(double));
    ky = (double *) malloc(Ny2  * sizeof(double));
    ek = (double *) malloc(nshell * sizeof(double));
    ensk = (double *) malloc(nshell * sizeof(double));
    ekx = (double *) malloc(kxm * sizeof(double));
    eky = (double *) malloc(kym * sizeof(double));
    
//    printf("local_nx = %ld \t local_x_start = %ld \t rank = %d\n", local_nx, local_x_start, rank);
if(rank == master)
{
    printf("\nCreating plans for FFTW........\n");
}

    pfor = fftw_mpi_plan_dft_r2c_2d(Nx, Ny, u, cu, MPI_COMM_WORLD, FFTW_ESTIMATE);
    pinv = fftw_mpi_plan_dft_c2r_2d(Nx, Ny, cu, u, MPI_COMM_WORLD, FFTW_ESTIMATE);

if(rank == master)
{
  printf("\nSetting initial conditions.........\n");
}


//..................................................................
  if(rank == master)
  {	
    printf("gamma1 = %lf\n", gamma1);
    printf("gamma2 = %lf\n", gamma2);
    printf("gamma3 = %lf\n", gamma3);
    printf("del = %lf\n", del);
  }
    kxky();
    init();
    init_forcing();
//..................................................................
    
//    it = 0; ifile = 0; fl = 0;
 
#ifdef BINARY
    free_energy(cphi);
    kinetic_energy(cwz, cphi);
    perimeter(cphi);
    spectrum(cwz, cphi);
#elif defined TERNARY
    free_energy(cphi1, cphi2);
    kinetic_energy(cwz, cphi1, cphi2);
    perimeter(cphi1, cphi2);
    spectrum(cwz, cphi1, cphi2);
#else
    kinetic_energy(cwz);
    spectrum(cwz);
#endif

    #ifdef WRITE_ARRAY
      copy(cwz, tcwz, n); ifft(tcwz, twz);

      sprintf(filename, "vort/wz.%d", fl);
      write_array_real(twz, filename);

      #ifdef BINARY
        copy(cphi, tcphi, n); ifft(tcphi, tphi);
        sprintf(filename, "phase_field/phi.%d", fl);
        write_array_real(tphi, filename);
      #elif defined TERNARY
        copy(cphi1, tcphi1, n);  copy(cphi2, tcphi2, n);
	ifft(tcphi1, tphi1);     ifft(tcphi2, tphi2);
	sprintf(filename, "phase_field/phi1.%d", fl);
        write_array_real(tphi1, filename);
	sprintf(filename, "phase_field/phi2.%d", fl);
        write_array_real(tphi2, filename);
      #endif	
     fl++;
    #endif

if(rank == master)
{
    printf("\nTime stepping start.............\n");
}
//********************************Time Starts******************************************************
//    t = dt;
      //copy(cwz, cnlz, n);
    for(it = 1; it <= maxiter; it++)
    {

t += dt;

    #ifdef BINARY
      etd2rk(cwz, cphi);
      free_energy(cphi);
      kinetic_energy(cwz, cphi);
      perimeter(cphi);

      if(it % (navg/2) == 0)
      {
        spectrum(cwz, cphi);
      }
    #elif defined TERNARY
      etd2rk(cwz, cphi1, cphi2);
      free_energy(cphi1, cphi2);
      kinetic_energy(cwz, cphi1, cphi2);
      perimeter(cphi1, cphi2);

      if(it % (navg/2) == 0)
      {
        spectrum(cwz, cphi1, cphi2);
      }
    #else 
      etd2rk(cwz);
      kinetic_energy(cwz);

      if(it % (navg/2) == 0)
      { 
        spectrum(cwz);
      }
    #endif
  

  #ifdef WRITE_ARRAY

  tmpi = 0.0;
  t1 = MPI_Wtime();   
 
    if(it % navg == 0)
    {

      #ifdef DUMP_REAL_SERIAL

        copy(cwz, tcwz, n);
        ifft(tcwz, twz);
 
        sprintf(filename, "vort/wz.%d", fl);
        write_array_real(twz, filename);

      #ifdef BINARY
        copy(cphi, tcphi, n);
        ifft(tcphi, tphi);

        sprintf(filename, "phase_field/phi.%d", fl);
        write_array_real(tphi, filename);
      #elif defined TERNARY
        copy(cphi1, tcphi1, n); copy(cphi2, tcphi2, n);
        ifft(tcphi1, tphi1);  ifft(tcphi2, tphi2);

        sprintf(filename, "phase_field/phi1.%d", fl);
        write_array_real(tphi1, filename);
        sprintf(filename, "phase_field/phi2.%d", fl);
        write_array_real(tphi2, filename);

      #endif
      
        fl++;
      
   
   #elif defined DUMP_REAL_MPI
      copy(cwz, tcwz, n); truncation(tcwz);
      ifft(tcwz, twz);

      sprintf(filename, "vort/wz.%d", fl);
      dump_real_mpi_IO(twz, filename);
      
      #ifdef BINARY
        copy(cphi, tcphi, n); truncation(tcphi);
        ifft(tcphi, tphi);

        sprintf(filename, "phase_field/phi.%d", fl);
        dump_real_mpi_IO(tphi, filename);
      #elif defined TERNARY
        copy(cphi1, tcphi1, n); copy(cphi2, tcphi2, n);
        truncation(tcphi1); truncation(tcphi2);
        ifft(tcphi1, tphi1);  ifft(tcphi2, tphi2);

        sprintf(filename, "phase_field/phi1.%d", fl);
        dump_real_mpi_IO(tphi1, filename);

        sprintf(filename, "phase_field/phi2.%d", fl);
        dump_real_mpi_IO(tphi2, filename);

      #endif

      fl++;
   #endif //DUMP_REAL_SERIAL || DUMP_REAL_MPI
   } //navg........

  t2 = MPI_Wtime();
  tmpi += (t2 - t1);

  #endif //WRITE_ARRAY...............     

 //     t += dt;
    }//...Time loop ends here.........................
//***************************************************************************************************

  #ifdef WRITE_ARRAY
    MPI_Reduce(rank == 0? MPI_IN_PLACE:&tmpi, &tmpi, 1, MPI_DOUBLE, MPI_SUM, master, MPI_COMM_WORLD);

    if(rank == master)
    {
      printf("MPI time for file writing = %lf\n", tmpi);
    }
  #endif
//..........Free the allocated memories.................................   
   fftw_destroy_plan(pfor);
   fftw_destroy_plan(pinv);
   fclose(fu); free(work);
   fftw_free(cux); fftw_free(cuy);
   fftw_free(cwz); 
   fftw_free(tcux); fftw_free(tcuy);
   fftw_free(cnlz);  fftw_free(cnl_z);
   free(ek); free(ensk); free(ekx); free(eky);   
if(rank == master)
{
   fftw_free(u1);
}
#ifdef VPM
   fftw_free(chi); fftw_free(sponge);
#endif
#ifdef BINARY
    fclose(fph); fclose(fe); fclose(fdef);
    fftw_free(cphi); fftw_free(tcphi); fftw_free(ncphi);
    fftw_free(cnlphi); fftw_free(cnl_phi);
    free(sph); fftw_free(ph);
    fftw_free(aph); fftw_free(bph);
#elif defined TERNARY
    fclose(fph); fclose(fe); fclose(fdef1); fclose(fdef2);
    fftw_free(cphi1); fftw_free(cphi2);
    fftw_free(tcphi1); fftw_free(tcphi2);
    fftw_free(muF1); fftw_free(muF2);
    fftw_free(nlphi1); fftw_free(nlphi2);
    fftw_free(cnl_phi1); fftw_free(cnl_phi2);
    free(sph1); free(sph2); free(sp12);
    fftw_free(ph1); fftw_free(ph2);
    fftw_free(aph1); fftw_free(aph2);
    fftw_free(bph1); fftw_free(bph2);
#endif
   MPI_Finalize();
return 0;

}
//...........End of the main function..............................................
