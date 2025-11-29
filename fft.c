#include "global.h"
/*!
 ******************************************************************************************
 Program: CHNS-2D-MPI

 This functions evaluate forward and backward FFT. 
 Note that the ifft function gives normalized values.

 ******************************************************************************************/

void dfft(double *u, fftw_complex *cu)
{
    fftw_mpi_execute_dft_r2c (pfor, u, cu);
}	

void ifft(fftw_complex *cu, double *u)
{
    double scale;
    int N;
    N = local_nx * (Ny+2);

    scale = 1.0 / (Nx * Ny);

    fftw_mpi_execute_dft_c2r (pinv, cu, u);

    for(i = 0; i < N; i++)
    {
      u[i] *= scale;
    }

}

/*******************************************************************************************/
