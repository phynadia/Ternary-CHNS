#include "global.h"
/*!
 ********************************************************************************
 
 This function introduce forcing to the flows.
 The following forcing functions can be used:
 
 1. Taylor-Green forcing.
 2. Constant forcing, which can be introduced in many selected wave numbers.
    (<Direct numerical simulation of homogeneous turbulence with hyperviscosity
     A. G. Lamorgese, D. A. Caughey, and S. B. Pope>)

 ********************************************************************************/

#define kk2_low (k_begin * k_begin)
#define kk2_high pow((k_begin+k_range), 2.0)

#ifndef DOXYGEN_SHOULD_SKIP_THIS

int Nforced, *IDXforce;
double ee, eeall, famp;
double ksqr, *in, fac;
double x, y;
int i1;

#endif

void init_forcing()
{

#ifdef CONSTANT_FORCING
    Nforced = 0;

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        ksqr = kx[i]*kx[i] + ky[j]*ky[j];

	if((kk2_low <= ksqr) && (ksqr <= kk2_high))
	{
	  Nforced++;
	}

      }
    }

    in = (double*) malloc (Nforced * sizeof(double));
    IDXforce = (int*) malloc (Nforced * sizeof(int));

    Nforced = 0;

    for (i = 0; i < local_nx; i++)
    {	    
      for (j = 0; j < Ny2; j++)
      {	            
        Idx = IDXC (i, j);
        ksqr = kx[i]*kx[i] + ky[j]*ky[j];
        
	fac = 1.0;

       
          if (ky[j] == 0.0 || ky[j] == Ny/2)
            fac = 1.0;
          else
            fac = 2.0;


          if ((kk2_low <= ksqr) && (ksqr <= kk2_high)) 
	  {
            in[Nforced] = fac;
            IDXforce[Nforced] = Idx;
            Nforced++;
          }

      }
    }

#endif

#ifdef TG_FORCING

    for(i = 0; i < local_nx; i++)
    {
      i1 = i + local_x_start;
      x = i1 * dx;      
      for(j = 0; j < Ny; j++)
      {
	y = j * dy;

	Idx = IDXR(i , j);

	fz[Idx] = forcing_amplitude * (cos(kf*y) + sin(kf*x));
	
      }
    }

    dfft(fz, cfz);

#endif

}


void forcing(fftw_complex *cux, fftw_complex *cuy, fftw_complex *cuz)
{

#ifdef CONSTANT_FORCING	
  double scale;

  scale = 1.0 / (Nx*Ny*Nz); 
  ee = 0.0; eeall = 0.0; 

    for(i = 0; i < Nforced; i++)
    {
      Idx = IDXforce[i];

      ee += in[i] * (cux[Idx] * conj(cux[Idx]) + cuy[Idx] * conj(cuy[Idx]));
    } 

    ee *= scale * scale;

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Allreduce (&ee, &eeall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(eeall > 0.0)
    {
      famp = dt * forcing_amplitude / (2.0*eeall);
    }else
      famp = 0.0;

    for(i = 0; i < Nforced; i++)
    {
      Idx = IDXforce[i];

      cux[Idx] += famp * cux[Idx];
      cuy[Idx] += famp * cuy[Idx];
      cuz[Idx] += famp * cuz[Idx];
    }
#endif

#ifdef TG_FORCING
    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
        for(k = 0; k < Nz2; k++)
	{
	  Idx = IDXC(i, j, k);

         cux[Idx] += dt * cfx[Idx];
         cuy[Idx] += dt * cfy[Idx];
         cuz[Idx] += dt * cfz[Idx];

	}
      }
    }
#endif

}

/*******************************************************************************************/
