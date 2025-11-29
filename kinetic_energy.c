#include "global.h"

void kinetic_energy(fftw_complex *cwz
#ifdef BINARY
  ,fftw_complex *cphi
#elif defined TERNARY
  ,fftw_complex *cphi1, fftw_complex *cphi2
#endif  
)
{	
    double cx, cy, px, py, ph, ph1, ph2;
    double scale, eu, eph, eph1, eph2;
    int count, count1, count2, n;

    n = local_nx * Ny2;
//................................................................
    get_velocity(cwz, cux, cuy);
    truncation(cux); truncation(cuy);    
    ifft(cux, ux);
    ifft(cuy, uy);

#ifdef BINARY
    copy(cphi, tcphi, n);
    truncation(tcphi);
    ifft(tcphi, tphi);
#elif defined TERNARY
    copy(cphi1, tcphi1, n);
    copy(cphi2, tcphi2, n);
    truncation(tcphi1); truncation(tcphi2);
    ifft(tcphi1, tphi1);
    ifft(tcphi2, tphi2);
#endif
//................................................................

    scale = 1.0 / (Nx * Ny);
    eu = 0.0; eph = 0.0;
    eph1 = 0.0; eph2 = 0.0;
    count = 0; count1 = 0; count2 = 0;

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
        Idx = IDXR(i, j);

        cx = ux[Idx]; cy = uy[Idx];
	px = cx; py = cy;

	eu += 0.5 * (cx * cx + cy * cy);

    #ifdef BINARY	
	ph = tphi[Idx];

	if(ph >= 0.5)
	{
	  ph = 1.0;
	  count += 1;
	}else{
	       ph = 0.0;
	     }	
	cx *= ph; cy *= ph;

	eph += 0.5 * (cx * cx + cy * cy);

    #elif defined TERNARY
        ph1 = tphi1[Idx];
	ph2 = tphi2[Idx];

	if(ph1 >= 0.5)
	{
	  ph1 = 1.0;
	  count1 += 1;
	}else{
	       ph1 = 0.0;
	     }	
        if(ph2 >= 0.5)
        {
          ph2 = 1.0;
	  count2 += 1;
        }else{
               ph2 = 0.0;
             }

        cx *= ph1; cy *= ph1;
	px *= ph2; py *= ph2;

	eph1 += 0.5 * (cx * cx + cy * cy);
	eph2 += 0.5 * (px * px + py * py);

    #endif	
	
      }	      
    }	   

   eu *= scale; 
   MPI_Reduce(rank == 0? MPI_IN_PLACE:&eu, &eu, 1, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD); 

#ifdef BINARY
   MPI_Reduce(rank == 0? MPI_IN_PLACE:&count, &count, 1, MPI_INT,
                          MPI_SUM, 0, MPI_COMM_WORLD);
  
   MPI_Reduce(rank == 0? MPI_IN_PLACE:&eph, &eph, 1, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);
   if(rank == 0 && count != 0)
   {	   
     eph /= count;

     fprintf(fu,"%lf \t %.12lf \t %.12lf\n", t, eu, eph);
   }
#elif defined TERNARY
   MPI_Reduce(rank == 0? MPI_IN_PLACE:&count1, &count1, 1, MPI_INT,
                          MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(rank == 0? MPI_IN_PLACE:&count2, &count2, 1, MPI_INT,
                          MPI_SUM, 0, MPI_COMM_WORLD);

   MPI_Reduce(rank == 0? MPI_IN_PLACE:&eph1, &eph1, 1, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);

   MPI_Reduce(rank == 0? MPI_IN_PLACE:&eph2, &eph2, 1, MPI_DOUBLE,
                          MPI_SUM, 0, MPI_COMM_WORLD);
   if(rank == 0)
   {
     if (count1 != 0)	
     {   
       eph1 /= count1;
     }

     
     if (count2 != 0)
     {
       eph2 /= count2;
     }  
  
     fprintf(fu,"%lf \t %.12lf \t %.12lf \t %.12lf\n", t, eu, eph1, eph2);
   }


#else
   if(rank == 0)
     fprintf(fu,"%lf \t %.12lf\n", t, eu);
#endif

}//End of the function..........................	
