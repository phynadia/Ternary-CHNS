#include "global.h"
/*!
 **************************************************************************************
 * This function evaluates the surface tension term present in Navier-Stokes equations.
 *
 * Surface tension:
 * 1. For binary:  \f$ \mu \bf{\nabla} c\f$
 *
 * 2. For ternary: \f$ \Sigma_i {\mu}_i c_i\f$
 *
 **************************************************************************************/

void surface_tension()
{
   double ksqr;
   int n;
   n = local_nx * Ny2;

#ifdef BINARY
   copy(cphi, ncphi, n);
   truncation(ncphi);
   ifft(ncphi, nphi);

   for(i = 0; i < local_nx; i++)
   {
     for(j = 0; j < Ny; j++)
     {
       Idx = IDXR(i, j);

       nphi[Idx] = nphi[Idx] * (1.0 - nphi[Idx]) 
                             * (1.0 - 2.0 * nphi[Idx]);       
     }
   }

   dfft(nphi, ncphi);

   for(i = 0; i < local_nx; i++)
   {
     for(j = 0; j < Ny2; j++)
     {     
       Idx = IDXC(i, j);

       ksqr = kx[i]*kx[i] + ky[j]*ky[j];

       ncphi[Idx] = 1.5 * sigma * epsilon * ksqr * cphi[Idx]
                  + 24.0 * (sigma / epsilon) * ncphi[Idx];       
     }
   }


    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

        tcux[Idx] = I * kx[i] * cphi[Idx];
        tcuy[Idx] = I * ky[j] * cphi[Idx];

          /**************************************
	   * 
	   * Evaluation of grad(c), which contain
	   * three components.
	   *
	   **************************************/
       
      }
    }
   
    truncation(tcux); truncation(tcuy); 
    truncation(ncphi);

    ifft(tcux, tux);
    ifft(tcuy, tuy);
    ifft(ncphi, nphi);
  
    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {        
         Idx = IDXR(i, j);
         
         #ifdef VPM
          tux[Idx] *= nphi[Idx] * (1.0 - chi[Idx]);
          tuy[Idx] *= nphi[Idx] * (1.0 - chi[Idx]);
         #else
          tux[Idx] *= nphi[Idx];
          tuy[Idx] *= nphi[Idx];
         #endif   
	  /**************************************
	   * Evaluation of surface tension:
	   *  mu grad(c).
	   *
	   **************************************/        
      }
    }
    
    dfft(tux, tcux);
    dfft(tuy, tcuy);

    curl(tcux, tcuy, tcwz);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {     
         Idx = IDXC(i, j);

        cnlz[Idx] += tcwz[Idx];
        
      }
    } //..........BINARY.................

#elif defined TERNARY

    copy(cmuF1, tcphi1, n);
    copy(cmuF2, tcphi2, n);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {        
        Idx = IDXC(i, j);
       
        ksqr = kx[i]*kx[i] + ky[j]*ky[j];       
 
        tcphi1[Idx] += 0.75 * gamma1 * epsilon * ksqr * cphi1[Idx];
        tcphi2[Idx] += 0.75 * gamma2 * epsilon * ksqr * cphi2[Idx];

        tcux[Idx] = I * kx[i] * cphi1[Idx];
        tcuy[Idx] = I * ky[j] * cphi1[Idx];

      }
    }
   
    truncation(tcux); truncation(tcuy); 
    truncation(tcphi1); truncation(tcphi2);

    ifft(tcux, tux); ifft(tcuy, tuy);
    ifft(tcphi1, tphi1); ifft(tcphi2, tphi2);
   
    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
         Idx = IDXR(i, j);

        tux[Idx] *= ((1.0 + (gamma3/gamma1)) * tphi1[Idx] + (gamma3/gamma2) * tphi2[Idx]);
        tuy[Idx] *= ((1.0 + (gamma3/gamma1)) * tphi1[Idx] + (gamma3/gamma2) * tphi2[Idx]);

        #ifdef VPM
          tux[Idx] *= (1.0 - chi[Idx]);
          tuy[Idx] *= (1.0 - chi[Idx]);
        #endif

      }
    }

    dfft(tux, tcux);
    dfft(tuy, tcuy);

    curl(tcux, tcuy, tcwz);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

         cnlz[Idx] += tcwz[Idx];

      }
    }

/********************************************************************/
    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {      
        Idx = IDXC(i, j);

        tcux[Idx] = I * kx[i] * cphi2[Idx];
        tcuy[Idx] = I * ky[j] * cphi2[Idx];

      }
    }
   
    truncation(tcux); truncation(tcuy);

    ifft(tcux, tux);
    ifft(tcuy, tuy);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {        
        Idx = IDXR(i, j);

        tux[Idx] *= ((gamma3/gamma1) * tphi1[Idx] + (1.0 + (gamma3/gamma2)) * tphi2[Idx]);
        tuy[Idx] *= ((gamma3/gamma1) * tphi1[Idx] + (1.0 + (gamma3/gamma2)) * tphi2[Idx]);

        #ifdef VPM
          tux[Idx] *= (1.0 - chi[Idx]);
          tuy[Idx] *= (1.0 - chi[Idx]);
        #endif        
      }
    }

    dfft(tux, tcux);
    dfft(tuy, tcuy);

    curl(tcux, tcuy, tcwz);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

          cnlz[Idx] += tcwz[Idx];
      
      }
    }

#endif //TERNARY

}//..End of the function surface_tension()........................................................


