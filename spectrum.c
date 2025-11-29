#include "global.h"

void spectrum(fftw_complex *cwz
#ifdef BINARY
, fftw_complex *cphi
#elif defined TERNARY
, fftw_complex *cphi1, fftw_complex *cphi2
#endif
)
{
    int ik, ikx, iky;
    double ksqr, kxsqr, kysqr, eu, ens;
    double eph, eph1, eph2, eph12, scale;
    double kmod, delta_k;

    delta_k = 2.0 * pi / lx;

    scale = 1.0 / (Nx * Ny);

//...................................................
    get_velocity(cwz, cux, cuy);
    truncation(cux); truncation(cuy);
#ifdef BINARY
    truncation(cphi);
#elif defined TERNARY
    truncation(cphi1); truncation(cphi2);
#endif    
//...................................................

    for(ik = 0; ik < nshell; ik++)
    {
      ek[ik] = 0.0;
      ensk[ik] = 0.0;
#ifdef BINARY
      sph[ik] = 0.0;
#elif defined TERNARY      
      sph1[ik] = 0.0;
      sph2[ik] = 0.0;
      sp12[ik] = 0.0;
#endif
    }
    for(ikx = 0; ikx < kxm; ikx++)
    {
      ekx[ikx] = 0.0;
    }
    for(iky = 0; iky < kym; iky++)
    {
      eky[iky] = 0.0;
    }

    
    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {

        Idx = IDXC(i, j);
        kxsqr = kx[i] * kx[i];
        kysqr = ky[j] * ky[j];
        ksqr = kxsqr + kysqr;

	kmod = sqrt(ksqr);
        ik  = (int) ((kmod/delta_k) + 0.5);
	ikx = (int) ((sqrt(kxsqr)/delta_k) + 0.5);
	iky = (int) ((sqrt(kysqr)/delta_k) + 0.5);

        eu = cux[Idx] * conj(cux[Idx]) + cuy[Idx] * conj(cuy[Idx]);
        eu *= scale * scale;
        ens = eu * ksqr;
#ifdef BINARY
        eph = cphi[Idx] * conj(cphi[Idx]);
        eph *= scale * scale;

	if(ky[j] == 0.0 || ky[j] == kymax)
        {
          sph[ik] += eph;
        }else
          {
            sph[ik] += 2.0 * eph;
          }

#elif defined TERNARY
        eph1 = cphi1[Idx] * conj(cphi1[Idx]);
        eph1 *= scale * scale;
	eph2 = cphi2[Idx] * conj(cphi2[Idx]);
        eph2 *= scale * scale;

	eph12 = creal(cphi1[Idx]) * creal(cphi2[Idx])
	      + cimag(cphi1[Idx]) * cimag(cphi2[Idx]);
	eph12 *= scale * scale;
        eph12 = fabs(eph12);

	if(ky[j] == 0.0 || ky[j] == kymax)
        {
          sph1[ik] += eph1;
	  sph2[ik] += eph2;
	  sp12[ik] += eph12;
        }else
          {
            sph1[ik] += 2.0 * eph1;
	    sph2[ik] += 2.0 * eph2;
	    sp12[ik] += 2.0 * eph12;
          }

#endif
        if(ky[j] == 0.0 || ky[j] == kymax)
        {
          ek[ik] += 0.5 * eu;
          ensk[ik] += 0.5 * ens;
	  ekx[ikx] += 0.5 * eu;
	  eky[iky] += 0.5 * eu;
        }else
          {
            ek[ik] += eu;
            ensk[ik] += ens;
	    ekx[ikx] += eu;
	    eky[iky] += eu;
          }

      }
    }

    MPI_Reduce(rank == 0? MPI_IN_PLACE:ek, ek, nshell, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(rank == 0? MPI_IN_PLACE:ensk, ensk, nshell, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(rank == 0? MPI_IN_PLACE:ekx, ekx, kxm, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(rank == 0? MPI_IN_PLACE:eky, eky, kym, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);

#ifdef BINARY
    MPI_Reduce(rank == 0? MPI_IN_PLACE:sph, sph, nshell, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);
#elif defined TERNARY
    MPI_Reduce(rank == 0? MPI_IN_PLACE:sph1, sph1, nshell, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(rank == 0? MPI_IN_PLACE:sph2, sph2, nshell, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(rank == 0? MPI_IN_PLACE:sp12, sp12, nshell, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);
#endif

     if(rank == 0)
    {

      //  if(it % navg == 0)
      //  {
//.......................................................................		
          sprintf(filename, "spectra/spectrum.%d", ifile);
          fek = fopen(filename, "w");

          for(ik = 1; ik < nshell; ik++)
          {
           double ik_new;
           ik_new = ik * (2 * pi / lx);	   
           #ifdef BINARY
            fprintf(fek, "%.5e \t %.12e \t %.12e \t %.12e\n", ik_new, ek[ik], ensk[ik], sph[ik]);
           #elif defined TERNARY
            fprintf(fek, "%.5e \t %.12e \t %.12e \t %.12e \t %.12e \t %.12e\n", 
			    ik_new, ek[ik], ensk[ik], sph1[ik], sph2[ik], sp12[ik]);
           #else		  
            fprintf(fek, "%.5e \t %.12e \t %.12e\n", ik_new, ek[ik], ensk[ik]);
           #endif 
          }

	  fclose(fek);
//........................................................................
	  sprintf(filename, "spectra/Xspectrum.%d", ifile);
          fek = fopen(filename, "w");

          for(ikx = 1; ikx < kxm; ikx++)
          { 
            fprintf(fek, "%d \t %.12e\n", ikx, ekx[ikx]);
          }

          fclose(fek);
//........................................................................
          sprintf(filename, "spectra/Yspectrum.%d", ifile);
          fek = fopen(filename, "w");

          for(iky = 1; iky < kym; iky++)
          {
            fprintf(fek, "%d \t %.12e\n", iky, eky[iky]);
          }

          fclose(fek);
//........................................................................	  
          ifile++;
        //}

    }//...rank.........

}//End of the function......................................................
