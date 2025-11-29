#include "global.h"

void get_velocity(fftw_complex *cwz, fftw_complex *cux, fftw_complex *cuy)
{

    double ksqr;
    int n = local_nx * Ny2;

    copy(cwz, tcwz, n);

  //......Poisson solver psi_k = wz_k / ksqr............
    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

        ksqr = kx[i]*kx[i] + ky[j]*ky[j];

	if(ksqr > 1e-12)
        {
	  tcwz[Idx] /= ksqr;
        }else{
                tcwz[Idx] = 0.0;
             }

      }	      
    }  

  //....Get velocity from stream-function psi. u = nabla x psi

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

	cux[Idx] =  I * ky[j] * tcwz[Idx];
        cuy[Idx] = -I * kx[i] * tcwz[Idx];

      }	      
    }	    

}//End of the function get_velocity.
