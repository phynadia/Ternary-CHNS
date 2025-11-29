#include "global.h"
/*!
 ***********************************************************************
 ***********************************************************************
 * Program: CHNS2D-MPI
 *
 * This curl function takes velocity fields and return the curl of the
 * fields and stores in the same arrays.
 *
 * The calculations happen in k-space.
 *
 * 1. \f$ \bf{\omega_k} = \iota k \times \bf{u_k} \f$
 *
 ***********************************************************************/

void curl(fftw_complex *cux, fftw_complex *cuy, fftw_complex *cwz)
{

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {

	Idx = IDXC(i, j);

	cwz[Idx] = I * (kx[i] * cuy[Idx] - ky[j] * cux[Idx]);
          
      }
    }

}//End of the function curl.
