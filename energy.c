#include"global.h"
/*!
 ****************************************************************************
 * This function calculates the kinetic-energy spectra directly in k-space.
 *****************************************************************************/

void energy(fftw_complex *cwz)
{
    double dealias, dax, day;
    double eu, ksqr, scale;
    int ik;

    dax = (Nx/3.0) * factor_x;
    day = (Ny/3.0) * factor_y;
    dealias = dax*dax + day*day;

    eu  = 0.0; 
    scale = 1.0 / (Nx * Ny);

    for(ik = 0; ik < nshell; ik++)
    {
      Ek[ik] = 0.0; ek[ik] = 0.0;
    }
//..........................................................
    get_velocity(cwz, cux, cuy);
//..........................................................
    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
      
        Idx = IDXC(i, j);

        ksqr = kx[i] * kx[i] + ky[j] * ky[j];
	ik = round(sqrt(ksqr));
          
        eu = cux[Idx] * conj(cux[Idx]) + cuy[Idx] * conj(cuy[Idx]);
	eu *= scale * scale;

	if(ky[j] == 0.0 || ky[j] == kymax)
	{
          ek[ik] += 0.5 * eu; 
	}else
	  {
	    ek[ik] += eu;
	  }

        if(ksqr >= dealias)
          ek[ik] = 0.0;

      }
    }

    MPI_Reduce(ek, Ek, nshell, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
	Et = 0.0;
        for(ik = 0; ik < nshell; ik++)
        {
           Et += Ek[ik]; 
        }
	fprintf(fe, "%lf \t %.12e\n", t, Et);
//        printf("%lf \t %lf\n", t, Et);

        if(it % navg == 0)
	{
	  sprintf(filename, "spectra/spectrum.%d", ifile);
          fek = fopen(filename, "w");

	  for(ik = 1; ik < nshell; ik++)
	  {
	    fprintf(fek, "%d \t %.12e\n", ik, Ek[ik]);
	  }

          ifile++;
	  fclose(fek);
	}
    }

}//End of the function energy().
