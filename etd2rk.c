#include"global.h"
/*!
 *******************************************************************************
 
  Program: CHNS2D-MPI

  
  Most of the calculations are done in this function.

  ETD2RK time integration scheme is used(Ref. Cox & Mathew).

 *******************************************************************************/

void etd2rk(fftw_complex *cwz
#ifdef BINARY
, fftw_complex *cphi
#elif defined TERNARY
, fftw_complex *cphi1, fftw_complex *cphi2
#endif
)
{
    double ksqr, eta, c, zeta;
    int n = local_nx * Ny2;

#ifdef BINARY
    double beta;
    
    nonlin(cwz, cnlz, cphi, cnlphi);
/*************************************************
 *
 *  (cwx, cwy, cwz) = u x w + surface tension.
 *  tcphi = advection + 3rd order nonlinear term.
 *
 *************************************************/
    copy(cnlz, cnl_z, n);
    copy(cnlphi, cnl_phi, n);
#elif defined TERNARY
    double beta;

    potential(cphi1, cphi2, cmuF1, cmuF2);  

/********************************************
 * 
 * muF1 and muF2 contain only the nonlinear 
 * parts of the chemical potentials. muF1 and 
 * muF2 are evaluated using tphi1 and tphi2
 * without altering them.
 *
 ********************************************/

    nonlin(cwz, cnlz, cphi1, cphi2, cnlphi1, cnlphi2);  

/***********************************************
 * 
 * At this point, tcphi1 and tcphi2 contain the 
 * entire nonlinear terms (advection + 3rd order 
 * nonlinearity) of CH equations. 
 *
 * (cwx, cwy, cwz) = u x w + surface-tension.
 ***********************************************/    
    
    copy(cnlz, cnl_z, n);
    
    copy(cnlphi1, cnl_phi1, n);
    copy(cnlphi2, cnl_phi2, n);
/************************************************
 *
 * tcphi1 and tcphi2 are copied into temporary 
 * arrays for the calculations of the second step 
 * of time integration(ETD2RK).
 *
 ************************************************/

#else
    nonlin(cwz, cnlz);
    
    copy(cnlz, cnl_z, n);
#endif
    

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {

	  Idx = IDXC(i, j);

	  ksqr = kx[i] * kx[i] + ky[j] * ky[j];
	  
	  c = - nu * ksqr;
          eta = exp(-nu * dt * ksqr);
	  if(fabs(c) > 1e-12)
	  {
	    zeta = (eta - 1.0) / c;
	  }else{
	    zeta = dt;
	  }

	  cwz[Idx] = cwz[Idx] * eta + cnlz[Idx] * zeta;

    #ifdef BINARY
      c = -1.5 * lambda * sigma * epsilon * ksqr * ksqr;
      
      beta = exp(dt * c);

      if(fabs(c) > 1e-12)
      { 
        zeta = (beta - 1.0) / c;
      }else{
        zeta = dt;
      }

      cphi[Idx] = cphi[Idx] * beta + cnlphi[Idx] * zeta;
      
    #elif defined TERNARY

      c = -0.75 * lambda * epsilon;
      c *= ksqr * ksqr;
      beta = exp(dt * c);

      if(fabs(c) > 1e-12)
      {
        zeta = (beta - 1.0) / c;
      }
      else{
            zeta = dt;
          }

      cphi1[Idx] = cphi1[Idx] * beta + cnlphi1[Idx] * zeta;
      cphi2[Idx] = cphi2[Idx] * beta + cnlphi2[Idx] * zeta;

    #endif
	
      }
    }

/***********************************************************
 *
 * The beginning of the second step of the time-integration.
 * 
 * All the steps are equivalent to the counterpart of step-1
 * of the time integration.
 *
 ***********************************************************/

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

          c = - nu * ksqr;
          eta = exp(-nu * dt * ksqr);

          if(fabs(c) > 1e-12)
          {
            zeta = (eta - 1.0 - c * dt) / (dt * c * c);
          }else{
            zeta = 0.5 * dt;
          }

          cwz[Idx] = cwz[Idx] + (cnlz[Idx] - cnl_z[Idx]) * zeta;

     #ifdef BINARY
          c = -1.5 * lambda * sigma * epsilon * ksqr * ksqr;
          
          beta = exp(dt * c);

          if(fabs(c) != 0.0)
          {
            zeta = (beta - 1.0 - c * dt) / (dt * c * c);
          }else{
            zeta = 0.5 * dt;
          }

          cphi[Idx] = cphi[Idx] + (cnlphi[Idx] - cnl_phi[Idx]) * zeta;
    
     #elif defined TERNARY
         c = -0.75 * lambda * epsilon;
         c *= ksqr * ksqr;
         beta = exp(dt * c);

         if(fabs(c) > 1e-12)
         {
           zeta = (beta - 1.0 - c * dt);
           zeta /= (dt * c * c);
         }
         else
            {
              zeta = 0.5 * dt;
            }
        
         cphi1[Idx] = cphi1[Idx] + (cnlphi1[Idx] - cnl_phi1[Idx])* zeta;
         cphi2[Idx] = cphi2[Idx] + (cnlphi2[Idx] - cnl_phi2[Idx])* zeta;

     #endif

       
      }
    }

#ifdef FORCING
    forcing(cux, cuy, cuz);
#endif

#ifdef GRAVITY
  #ifdef TERNARY
    
  #endif
#endif

}//End of the function ETD2RK.......................................................
