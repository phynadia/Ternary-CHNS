#include "global.h"

void runge_kutta4(fftw_complex *cwz
#ifdef BINARY
, fftw_complex *cphi
#elif defined TERNARY
, fftw_complex *cphi1, fftw_complex *cphi2
#endif
)
{

    double dth, nuksqr, ksqr;
    dth = 0.5 * dt; 
/*********************************************************************
                        Step - 1
**********************************************************************/
    #ifdef BINARY
      nonlin(cwz, cnlz, cphi, cnlphi);
    #elif defined TERNARY
      potential(cphi1, cphi2, cmuF1, cmuF2);
      nonlin(cwz, cnlz, cphi1, cphi2, cnlphi1, cnlphi2);
    #else
      nonlin(cwz, cnlz);
    #endif

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

          ksqr = kx[i] * kx[i] + ky[j] * ky[j];          
          nuksqr = nu * ksqr;
 
          crk1w[Idx] = cwz[Idx] + dth * cnlz[Idx]
                                - dth * nuksqr * cwz[Idx];
          #ifdef BINARY

            crk1phi[Idx] = cphi[Idx] + dth * cnlphi[Idx];

          #elif defined TERNARY

            crk1phi1[Idx] = cphi1[Idx] + dth * cnlphi1[Idx];
            crk1phi2[Idx] = cphi2[Idx] + dth * cnlphi2[Idx];

          #endif

      }
    }

//if(rank == master)
//printf("Test-1...........\n");
/*********************************************************************
                        Step - 2
**********************************************************************/
    #ifdef BINARY

      nonlin(crk1w, cnlz, crk1phi, cnlphi);

    #elif defined TERNARY

      potential(crk1phi1, crk1phi2, cmuF1, cmuF2);
      nonlin(crk1w, cnlz, crk1phi1, crk1phi2, cnlphi1, cnlphi2);

    #else

      nonlin(crk1w, cnlz);

    #endif

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

          ksqr = kx[i] * kx[i] + ky[j] * ky[j]; 
          nuksqr = nu * ksqr;

          crk2w[Idx] = cwz[Idx] + dth * cnlz[Idx]
                                - dth * nuksqr * crk1w[Idx];

          #ifdef BINARY

            crk2phi[Idx] = cphi[Idx] + dth * cnlphi[Idx];

          #elif defined TERNARY

            crk2phi1[Idx] = cphi1[Idx] + dth * cnlphi1[Idx];
            crk2phi2[Idx] = cphi2[Idx] + dth * cnlphi2[Idx];

          #endif
      }
    }
/***********************************************************************
                       Step - 3
************************************************************************/
    #ifdef BINARY

      nonlin(crk2w, cnlz, crk2phi, cnlphi);

    #elif defined TERNARY

      potential(crk2phi1, crk2phi2, cmuF1, cmuF2);
      nonlin(crk2w, cnlz, crk2phi1, crk2phi2, cnlphi1, cnlphi2);

    #else

      nonlin(crk2w, cnlz);

    #endif

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

          ksqr = kx[i] * kx[i] + ky[j] * ky[j]; 
          nuksqr = nu * ksqr;

          crk3w[Idx] = cwz[Idx] + dt * cnlz[Idx]
                                - dt * nuksqr * crk2w[Idx];

          #ifdef BINARY
            
            crk3phi[Idx] = cphi[Idx] + dt * cnlphi[Idx]; 

          #elif defined TERNARY

            crk3phi1[Idx] = cphi1[Idx] + dt * cnlphi1[Idx];
            crk3phi2[Idx] = cphi2[Idx] + dt * cnlphi2[Idx];

          #endif

      }
    }
/***********************************************************************
                      Step - 4 
(cnlx, cnly, cnlz) contain 4th step.
************************************************************************/
    #ifdef BINARY

      nonlin(crk3w, cnlz, crk3phi, cnlphi);
 
    #elif defined TERNARY

      potential(crk3phi1, crk3phi2, cmuF1, cmuF2);
      nonlin(crk3w, cnlz, crk3phi1, crk3phi2, cnlphi1, cnlphi2);

    #else

      nonlin(crk3w, cnlz);

    #endif
/***********************************************************************/

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

          ksqr = kx[i]*kx[i] + ky[j]*ky[j];
          nuksqr = nu * ksqr;

          cwz[Idx] = (1.0/3.0) * (-cwz[Idx] + crk1w[Idx] + 2.0 * crk2w[Idx] + crk3w[Idx] 
                                  - dth * (nuksqr * crk3w[Idx] - cnlz[Idx]));

          #ifdef BINARY
 
            cphi[Idx] = (1.0 / 3.0) * (-cphi[Idx] + crk1phi[Idx] + 2.0 * crk2phi[Idx]
                                       + crk3phi[Idx] + dth * cnlphi[Idx]);

          #elif defined TERNARY
            cphi1[Idx] = (1.0 / 3.0) * (-cphi1[Idx] + crk1phi1[Idx] + 2.0 * crk2phi1[Idx]
                                        + crk3phi1[Idx] + dth * cnlphi1[Idx]);
        
            cphi2[Idx] = (1.0 / 3.0) * (-cphi2[Idx] + crk1phi2[Idx] + 2.0 * crk2phi2[Idx]
                                        + crk3phi2[Idx] + dth * cnlphi2[Idx]);
          #endif

      }
    }

/***********************************************************************/
}//End of the function runge_kutta1=4
