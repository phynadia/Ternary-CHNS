#include "global.h"

void free_energy(
#ifdef BINARY		
  fftw_complex *cphi
#elif defined TERNARY
  fftw_complex *cphi1, fftw_complex *cphi2
#endif  
)
{
    double c, c1, c2, c3, scale;
    double c1x, c1y, c2x, c2y;    
    double ephi, egphi;
    int n;

    n = local_nx * Ny2;
    scale = 1.0 / (Nx * Ny);
    ephi = 0.0;

#ifdef BINARY
    copy(cphi, tcphi, n);
    ifft(tcphi, tphi);
#elif defined TERNARY
    copy(cphi1, tcphi1, n);
    copy(cphi2, tcphi2, n);
    ifft(tcphi1, tphi1);
    ifft(tcphi2, tphi2);
#endif

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
        
	Idx = IDXR(i, j);
       #ifdef BINARY
         c = tphi[Idx];

	 tphi[Idx] = c * c * (1.0 - c) * (1.0 - c);
	 tphi[Idx] *= (12.0 * (sigma/epsilon));

	 ephi += (tphi[Idx] * scale);

       #elif defined TERNARY
         c1 = tphi1[Idx]; c2 = tphi2[Idx];
	 c3 = 1.0 - c1 - c2;

	 tphi1[Idx] = gamma1 * c1 * c1 * (1.0 - c1) * (1.0 - c1)
		    + gamma2 * c2 * c2 * (1.0 - c2) * (1.0 - c2)
		    + gamma3 * c3 * c3 * (c1 + c2) * (c1 + c2);

	 tphi1[Idx] *= (6.0 / epsilon);
      
         ephi += (tphi1[Idx] * scale);
       #endif

      }	      
    }	    
    
      MPI_Reduce(rank == 0? MPI_IN_PLACE:&ephi, &ephi, 1, MPI_DOUBLE,
                           MPI_SUM, 0, MPI_COMM_WORLD);
     

#ifdef BINARY

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);
        
	tcux[Idx] = I * kx[i] * cphi[Idx];
	tcuy[Idx] = I * ky[j] * cphi[Idx];
      }	      
    }	    

    ifft(tcux, tux);
    ifft(tcuy, tuy);

    egphi = 0.0;

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {

        Idx = IDXR(i, j);

        c1 = tux[Idx]; c2 = tuy[Idx];

        egphi += (c1 * c1 + c2 * c2);

      }
    }

    egphi *= (0.75 * sigma * epsilon) * scale;
#endif // BINARY

#ifdef TERNARY

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

        tcux[Idx] = I * kx[i] * cphi1[Idx];
        tcuy[Idx] = I * ky[j] * cphi1[Idx];
	
	tcphi1[Idx] = I * kx[i] * cphi2[Idx];
	tcphi2[Idx] = I * ky[j] * cphi2[Idx];

      }
    }

    ifft(tcux, tux); ifft(tcuy, tuy);

    ifft(tcphi1, tphi1);
    ifft(tcphi2, tphi2);

    egphi = 0.0;

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {

        Idx = IDXR(i, j);

        c1x = tux[Idx]; c1y = tuy[Idx];
	c2x = tphi1[Idx]; c2y = tphi2[Idx];

        egphi += (gamma1 * (c1x * c1x + c1y * c1y) 
		+ gamma2 * (c2x * c2x + c2y * c2y) 
                + gamma3 * (c1x + c2x) * (c1x + c2x)
		+ gamma3 * (c1y + c2y) * (c1y + c2y));

      }
    }

    egphi *= (0.375 * epsilon * scale);
#endif    

    MPI_Reduce(rank == 0? MPI_IN_PLACE:&egphi, &egphi, 1, MPI_DOUBLE, 
		          MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
      fprintf(fe,"%lf \t %.12lf \t %.12lf\n", t, ephi, egphi);

      if(ephi >= 100)
      {
        MPI_Abort(MPI_COMM_WORLD, rank);
      }
    }

}//End of the function..................................................	
