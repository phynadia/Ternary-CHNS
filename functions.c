#include"global.h"
/*!
 ********************************************************************************
 * Program: CHNS-2D-MPI
 * 
 * 
 * 
 * This file contains many small useful functions called by other functions. 
 *
 ********************************************************************************/
void kxky()
{

  int i1;

  for(i = 0; i < local_nx; i++)
  { 
    i1 = i + local_x_start;	  
    if(i1 <= Nx/2)
    {
      kx[i] = factor_x * i1;
    }else
        {
          kx[i] = factor_x * (i1 - Nx);
        }
  }

  for(i = 0; i < Ny2; i++)
  {
    ky[i] = factor_y * i;
  }

}//End of the function kxkykz().

/************************************************************************************/
void copy(fftw_complex *in, fftw_complex *out, int size)
{
  int i;

  for(i = 0; i < size; i++)
  {
    out[i] = in[i];
  }
}

/*************************************************************************************/
void truncation(fftw_complex *cv)
{
   
   double dax, day;
   double ksqr, dealias;
   double ksqrmax, pp;

   ksqrmax = (double)(kxmax*kxmax + kymax*kymax);
   
   dax = (Nx / 4.0) * factor_x;
   day = (Ny / 4.0) * factor_y;

   dealias = dax*dax;// + day*day;

   for(i = 0; i < local_nx; i++)
   {
     for(j = 0; j < Ny2; j++)
     {
         
        Idx = IDXC(i, j);

	ksqr = kx[i]*kx[i] + ky[j]*ky[j];

        #ifdef SMOOTH_FILTER
          pp = ksqr / ksqrmax;
          cv[Idx] *= exp(-36.0 * pow(pp,18));
          cv[Idx] *= exp(-36.0 * pow(pp,18));
        #else
          if(ksqr > dealias)
          {
            cv[Idx] = 0.0;
          }
        #endif
       
     }
   }

}//End of the function truncation().

/***************************************************************************************/
void divergence(fftw_complex *cux, fftw_complex *cuy, 
                fftw_complex *div)
{
    fftw_complex c1, c2;

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

        c1 = kx[i] * cux[Idx];
        c2 = ky[j] * cuy[Idx];

        div[Idx] = I * (c1 + c2);
        
      }
    }

}
/**************************************************************************************/
