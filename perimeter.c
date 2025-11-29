#include "global.h"

void perimeter(
#ifdef BINARY		
fftw_complex *cphi
#elif defined TERNARY
fftw_complex *cphi1, fftw_complex *cphi2
#endif
)
{
   int count1, count2, Idxi, Idxj, c1, c2;
   double Sp_0, Sp_t, def;
   double Sp1_0, Sp1_t, def1;
   double Sp2_0, Sp2_t, def2;
   int n = local_nx * Ny2;

   int left_proc = rank - 1;
   int right_proc = rank + 1;
   int last_rank = size - 1;

   if(left_proc < 0)
      left_proc = MPI_PROC_NULL;
   if(right_proc >= size)
      right_proc = MPI_PROC_NULL;

#ifdef BINARY
    copy(cphi, tcphi, n);
    ifft(tcphi, tphi);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
        Idx = IDXR(i, j);
      	
        ph[Idx] = tphi[Idx];

	if(ph[Idx] >= 0.5)
	{
	  ph[Idx] = 1.0;
	}else
	    {
	      ph[Idx] = -1.0;
	    }

      }
    }

    for(j = 0; j < Ny; j++)
    {
      aph[j] = ph[j];
    }

    if(rank % 2 == 0)
    {	    
      MPI_Send(aph, Ny, MPI_DOUBLE, left_proc, 23, MPI_COMM_WORLD);
      MPI_Recv(bph, Ny, MPI_DOUBLE, right_proc, 25, MPI_COMM_WORLD, &status); 
    }
    if(rank % 2 == 1)
    {
      MPI_Recv(bph, Ny, MPI_DOUBLE, right_proc, 23, MPI_COMM_WORLD, &status);
      MPI_Send(aph, Ny, MPI_DOUBLE, left_proc, 25, MPI_COMM_WORLD);
    }
    
    if(rank == 0)
    {
      MPI_Send(aph, Ny, MPI_DOUBLE, last_rank, 24,MPI_COMM_WORLD);
    }
    if(rank == last_rank)
    {
      MPI_Recv(bph, Ny, MPI_DOUBLE, 0, 24, MPI_COMM_WORLD, &status);
    }

    for(j = 0; j < Ny; j++)
    {
      Idx = IDXR(local_nx, j);
      
      ph[Idx] = bph[j];
    }

    count1 = 0; count2 = 0;

    for(i = 0; i < local_nx; i++)
    {
      int ip = i + 1;
      	
      for(j = 0; j < Ny; j++)
      { 
        int jp = j + 1;
	if(jp >= Ny)
	{
          jp = 0;
	}
        Idx = IDXR(i, j);
	Idxi = IDXR(ip, j);
	Idxj = IDXR(i, jp);

	if((ph[Idx] * ph[Idxi]) < 0.0 || (ph[Idx] * ph[Idxj]) < 0.0)
	{
	  count1 += 1;
	}
	if(ph[Idx] > 0.0)
	{
	  count2 += 1;
	}

      }
    }

   MPI_Reduce(rank == 0? MPI_IN_PLACE:&count1, &count1, 1, MPI_INT,
                          MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(rank == 0? MPI_IN_PLACE:&count2, &count2, 1, MPI_INT,
                          MPI_SUM, 0, MPI_COMM_WORLD);

   if(rank == 0)
   {

    Sp_t = count1 * dx;
    Sp_0 = 2.0 * sqrt(pi * count2 * dx * dy);

    def = (Sp_t/Sp_0) - 1.0;

    fprintf(fdef,"%g \t %g \t %g \t %g\n", t, Sp_0, Sp_t, def);
  }


#elif defined TERNARY
    copy(cphi1, tcphi1, n);
    copy(cphi2, tcphi2, n);
    ifft(tcphi1, tphi1);
    ifft(tcphi2, tphi2);
//..........................................................
    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
        Idx = IDXR(i, j);

        ph1[Idx] = tphi1[Idx];
	ph2[Idx] = tphi2[Idx];

        if(ph1[Idx] >= 0.5)
        {
          ph1[Idx] = 1.0;
        }else
            {
              ph1[Idx] = -1.0;
            }

	if(ph2[Idx] >= 0.5)
        {
          ph2[Idx] = 1.0;
        }else
            {
              ph2[Idx] = -1.0;
            }


      }
    }

    for(j = 0; j < Ny; j++)
    {
      aph1[j] = ph1[j];
      aph2[j] = ph2[j];
    }

    if(rank % 2 == 0)
    {
      MPI_Send(aph1, Ny, MPI_DOUBLE, left_proc, 23, MPI_COMM_WORLD);
      MPI_Send(aph2, Ny, MPI_DOUBLE, left_proc, 24, MPI_COMM_WORLD);
      MPI_Recv(bph1, Ny, MPI_DOUBLE, right_proc, 25, MPI_COMM_WORLD, &status);
      MPI_Recv(bph2, Ny, MPI_DOUBLE, right_proc, 26, MPI_COMM_WORLD, &status);
    }
    if(rank % 2 == 1)
    {
      MPI_Recv(bph1, Ny, MPI_DOUBLE, right_proc, 23, MPI_COMM_WORLD, &status);
      MPI_Recv(bph2, Ny, MPI_DOUBLE, right_proc, 24, MPI_COMM_WORLD, &status);
      MPI_Send(aph1, Ny, MPI_DOUBLE, left_proc, 25, MPI_COMM_WORLD);
      MPI_Send(aph2, Ny, MPI_DOUBLE, left_proc, 26, MPI_COMM_WORLD);
    }

    if(rank == 0)
    {
      MPI_Send(aph1, Ny, MPI_DOUBLE, last_rank, 27,MPI_COMM_WORLD);
      MPI_Send(aph2, Ny, MPI_DOUBLE, last_rank, 28,MPI_COMM_WORLD);
    }
    if(rank == last_rank)
    {
      MPI_Recv(bph1, Ny, MPI_DOUBLE, 0, 27, MPI_COMM_WORLD, &status);
      MPI_Recv(bph2, Ny, MPI_DOUBLE, 0, 28, MPI_COMM_WORLD, &status);
    }

    for(j = 0; j < Ny; j++)
    {
      //Idx = IDXR(local_nx, j);
      Idx = j + (Ny+2) * local_nx; 
      ph1[Idx] = bph1[j];
      ph2[Idx] = bph2[j];
    }

    count1 = 0; count2 = 0;
    c1 = 0; c2 = 0;

    for(i = 0; i < local_nx; i++)
    {
      int ip = i + 1;
      for(j = 0; j < Ny; j++)
      {
        int jp = j + 1;
        if(jp >= Ny)
        {
          jp = 0;
        }
        Idx = IDXR(i, j);
        Idxi = IDXR(ip, j);
        Idxj = IDXR(i, jp);

        if((ph1[Idx] * ph1[Idxi]) < 0.0 || (ph1[Idx] * ph1[Idxj]) < 0.0)
        {
          count1 += 1;
        }
        if(ph1[Idx] > 0.0)
        {
          count2 += 1;
        }

	if((ph2[Idx] * ph2[Idxi]) < 0.0 || (ph2[Idx] * ph2[Idxj]) < 0.0)
        {
          c1 += 1;
        }
        if(ph2[Idx] > 0.0)
        {
          c2 += 1;
        }

      }
    }

   MPI_Reduce(rank == 0? MPI_IN_PLACE:&count1, &count1, 1, MPI_INT,
                          MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(rank == 0? MPI_IN_PLACE:&count2, &count2, 1, MPI_INT,
                          MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(rank == 0? MPI_IN_PLACE:&c1, &c1, 1, MPI_INT,
                          MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(rank == 0? MPI_IN_PLACE:&c2, &c2, 1, MPI_INT,
                          MPI_SUM, 0, MPI_COMM_WORLD);


   if(rank == 0)
   {

    Sp1_t = count1 * dx;
    Sp1_0 = 2.0 * sqrt(pi * count2 * dx * dy);

    def1 = (Sp1_t/Sp1_0) - 1.0;

    Sp2_t = c1 * dx;
    Sp2_0 = 2.0 * sqrt(pi * c2 * dx * dy);

    def2 = (Sp2_t/Sp2_0) - 1.0;


    fprintf(fdef1,"%g \t %g \t %g \t %g\n", t, Sp1_0, Sp1_t, def1);
    fprintf(fdef2,"%g \t %g \t %g \t %g\n", t, Sp2_0, Sp2_t, def2);
  }
 
#endif   


}//End of the function..................................
