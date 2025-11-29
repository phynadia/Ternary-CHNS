#include"global.h"
/*!
 ***************************************************************************
 This function evaluates all the nonlinear terms present in the equations.

 1. Navier-Stokes(NS):  \f$ u \times \omega\f$
 2. CHNS: Surface tension added in NS.
 3. Third order nonlinearity present in Cahn-Hilliard equation (CH).
 4. Advection term for order parameters.

 ***************************************************************************/

void nonlin(fftw_complex *cwz, fftw_complex *cnlz
#ifdef BINARY
, fftw_complex *cphi, fftw_complex *cnlphi
#elif defined TERNARY
, fftw_complex *cphi1, fftw_complex *cphi2, fftw_complex *cnlphi1, fftw_complex *cnlphi2
#endif
)
{  
    int n;
    double c1, c2, c3;
    double  ksqr;
    double Mchi;

    n = local_nx * Ny2;

//....Evaluating the nonlinear term -(u . nabla) wz.
    
    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

        tcux[Idx] = I * kx[i] * cwz[Idx];
        tcuy[Idx] = I * ky[j] * cwz[Idx];
      }
    }

//...Get velocity from the vorticity..........................        
    get_velocity(cwz, cux, cuy);
//............................................................
    truncation(cux); truncation(cuy); 
    truncation(tcux); truncation(tcuy);

    ifft(cux, ux); ifft(tcux, tux);
    ifft(cuy, uy); ifft(tcuy, tuy);

    #ifdef SPONGE
      copy(cwz, tcwz, n);
      truncation(tcwz);
      ifft(tcwz, twz);
    #endif

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
        Idx = IDXR(i, j);

        nlz[Idx] = ux[Idx] * tux[Idx] + uy[Idx] * tuy[Idx];
        #ifdef VPM
          ux[Idx] *= (chi[Idx] * etainv);
          uy[Idx] *= (chi[Idx] * etainv);
        #endif

        #ifdef SPONGE
          twz[Idx] *= sponge[Idx];
        #endif
      }
    }
 
    dfft(nlz, cnlz);

    #ifdef VPM
      dfft(ux, cux);
      dfft(uy, cuy);
      curl(cux, cuy, tcux);
    #endif

    #ifdef SPONGE
      dfft(twz, tcwz);
    #endif

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);
        
          cnlz[Idx] *= (-1.0);
        #ifdef VPM
          cnlz[Idx] -= tcux[Idx];
        #endif

        #ifdef SPONGE
          cnlz[Idx] -= tcwz[Idx];
        #endif 
      }
    }

    #if defined BINARY || defined TERNARY
      surface_tension();      
    #endif

#ifdef GRAVITY
    #ifdef BINARY
      copy(cphi, tcphi, n);
      ifft(tcphi, tphi);

      for(i = 0; i < local_nx; i++)
      {
        for(j = 0; j < Ny; j++)
        {
          Idx = IDXR(i, j);

          tphi[Idx] = (rho2 + (rho1 - rho2) * tphi[Idx]) - 1.0;
          tphi[Idx] *= (-1.0) * g;
       
        }
      }

      dfft(tphi, tcphi);

      for(i = 0; i < local_nx; i++)
      {
        for(j = 0; j < Ny2; j++)
        {     
          Idx = IDXC(i, j);
 
          tcphi[Idx] *= (I * kx[i]);

          cnlz[Idx] += tcphi[Idx];

        }
      }
      
    #elif defined TERNARY
      copy(cphi1, tcphi1, n); copy(cphi2, tcphi2, n);
      ifft(tcphi1, tphi1); ifft(tcphi2, tphi2);

      for(i = 0; i < local_nx; i++)
      {
        for(j = 0; j < Ny; j++)
        {
          Idx = IDXR(i, j);
            
          c1 = tphi1[Idx]; c2 = tphi2[Idx];
          c3 = 1.0 - c1 - c2;

          tphi1[Idx] = (rho3 + (rho1 - rho3) * c1 + (rho2 - rho3) * c2) - 1.0;
          tphi1[Idx] *= (-1.0) * g;

        }
      }

      dfft(tphi1, tcphi1);

      for(i = 0; i < local_nx; i++)
      {
        for(j = 0; j < Ny2; j++)
        {
          Idx = IDXC(i, j);

          tcphi1[Idx] *= (I * kx[i]);
          cnlz[Idx] += tcphi1[Idx];
        }
      }
    #endif //BINARY || TERNARY...........
#endif //GRAVITY.........

//.........Now evaluating R.H.S of Cahn-Hilliard equation.
#ifdef BINARY
    fftw_complex nln;

    get_velocity(cwz, cux, cuy);
    copy(cux, tcux, n);
    copy(cuy, tcuy, n);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

        tcphi[Idx]  = I * kx[i] * cphi[Idx];
        ncphi[Idx]  = I * ky[j] * cphi[Idx];

      }
    }


    truncation(tcphi); truncation(ncphi);
    truncation(tcux); truncation(tcuy);

  
    ifft(tcphi, tphi); ifft(ncphi, nphi); 
    ifft(tcux, tux); ifft(tcuy, tuy);   

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {       
        Idx = IDXR(i, j);

        nlphi[Idx] = (tux[Idx] * tphi[Idx] + tuy[Idx] * nphi[Idx]); 

        #ifdef VPM
          nlphi[Idx] *= (1.0 - chi[Idx]);
        #endif
      }
    }
    
    dfft(nlphi, cnlphi);
  
    copy(cphi, tcphi, n); truncation(tcphi); ifft(tcphi, tphi);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
        Idx = IDXR(i, j);

        tphi[Idx] = tphi[Idx] * (1.0 - tphi[Idx]) * (1.0 - 2.0 * tphi[Idx]);

      }
    }

    dfft(tphi, tcphi);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {  
        Idx = IDXC(i, j);

        ksqr = kx[i]*kx[i] + ky[j]*ky[j];

        tcphi[Idx] = 24.0 * (sigma/epsilon) * ksqr * tcphi[Idx];

      }
    }

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

        cnlphi[Idx] = -cnlphi[Idx] - lambda * tcphi[Idx];

      }
    }
                                                           
#elif defined TERNARY

    get_velocity(cwz, cux, cuy);
    copy(cux, tcux, n);
    copy(cuy, tcuy, n);

    truncation(tcux); truncation(tcuy);
 
    ifft(tcux, tux); ifft(tcuy, tuy);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

          tcphi1[Idx] = I * kx[i] * cphi1[Idx];
          tcphi2[Idx] = I * ky[j] * cphi1[Idx];

      }
    }
/******************************************************
 *(tcphi1, tcphi2) = grad(cphi1)
 ******************************************************/


    truncation(tcphi1); truncation(tcphi2);

    ifft(tcphi1, tphi1); ifft(tcphi2, tphi2);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
        Idx = IDXR(i, j);

        nlphi1[Idx] = tphi1[Idx] * tux[Idx]
                    + tphi2[Idx] * tuy[Idx]; 
        #ifdef VPM
          nlphi1[Idx] *= (1.0 - chi[Idx]);
        #endif
      }
    }

    dfft(nlphi1, cnlphi1);
   

/********************************************************
 * advection1 = -[(1-chi)u].grad(phi1)
 *
 ********************************************************/

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

        tcphi1[Idx] = I * kx[i] * cphi2[Idx];
        tcphi2[Idx] = I * ky[j] * cphi2[Idx];
      }
    }
/*******************************************************
 *
 * (tcphi1, tcphi2) = grad(cphi2)
 *
 *******************************************************/

    truncation(tcphi1); truncation(tcphi2);

    ifft(tcphi1, tphi1); ifft(tcphi2, tphi2);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
        Idx = IDXR(i, j);

          nlphi2[Idx] = tphi1[Idx] * tux[Idx]
                      + tphi2[Idx] * tuy[Idx];

        #ifdef VPM
          nlphi2[Idx] *= (1.0 - chi[Idx]);
        #endif
      }
    }
/*********************************************************
 *  
 *  advection2 = -[(1-chi)u].grad(phi2)
 *   
 *********************************************************/

    dfft(nlphi2, cnlphi2);

/********************************************************************************
 * 
 * At this point, tcphi1 and tcphi2 contain the advection terms of Cahn-Hilliard 
 * equations.
 * Note: tcphi1 and tcphi2 are the temporary variables which is helpful for 
 * calculating different quantities at different point of times.
 *
 ********************************************************************************/
/********************************************************************
 ********************************************************************/

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {   
        Idx = IDXC(i, j);

        ksqr = kx[i]*kx[i] + ky[j]*ky[j];

        cnlphi1[Idx] = -cnlphi1[Idx] - (lambda/gamma1) * ksqr * cmuF1[Idx];
	cnlphi2[Idx] = -cnlphi2[Idx] - (lambda/gamma2) * ksqr * cmuF2[Idx];
         
	 /******************************************************************
	  * tcphi(1 & 2) = RHS of CH-equations
	  ******************************************************************/

      }
    }

#endif

}//...End of the function nonlin()..................................................
