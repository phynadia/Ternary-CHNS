#include "global.h"
/*!
 ***********************************************************************
 * Program: CHNS2D-MPI
 *
 * This function calculates the nonlinear terms of chemical potentials.
 *
 *
 ***********************************************************************/

#ifdef TERNARY
void potential(fftw_complex *cphi1, fftw_complex *cphi2, fftw_complex *cmuF1, fftw_complex *cmuF2)
{

    double G1, G2, G3, c1, c2, c3;
    double psi1, psi2, psi3, ksqr;
    double fc1, fc2, fc3;
    double gamma1_inv, gamma2_inv, gamma3_inv;
    double c1c2c3;

    gamma1_inv = 1.0 / gamma1;
    gamma2_inv = 1.0 / gamma2;
    gamma3_inv = 1.0 / gamma3;

    int n = local_nx * Ny2;

    copy(cphi1, tcphi1, n);    
    copy(cphi2, tcphi2, n);

    truncation(tcphi1); truncation(tcphi2);
    ifft(tcphi1, tphi1);
    ifft(tcphi2, tphi2);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny; j++)
      {
          Idx = IDXR(i, j);

          c1 = tphi1[Idx]; c2 = tphi2[Idx]; c3 = (1.0 - c1 - c2);
	  c1c2c3 = c1 * c2 * c3;
/*
          psi1 = c1*c1 / pow((1.0+c1*c1), alpha);
          psi2 = c2*c2 / pow((1.0+c2*c2), alpha);
          psi3 = c3*c3 / pow((1.0+c3*c3), alpha);

          G1 = (-2.0 * alpha * psi1 * c1 * c2*c2 * c3*c3) / (1.0 + c1*c1)
             + 2.0 * c1 * (c3*c3 * psi2 + c2*c2 * psi3);
          G1 *= Xlambda;

          G2 = 2.0 * c2 * (c3*c3 * psi1 + c1*c1 * psi3)
             - (2.0 * alpha * psi2 * c1*c1 * c2 * c3*c3) / (1.0 + c2*c2);
          G2 *= Xlambda;

          G3 = 2.0 * c3 * (c2*c2 * psi1 + c1*c1 * psi2)
             - (2.0 * alpha * psi3 * c1*c1 * c2*c2 * c3) / (1.0 + c3*c3);
          G3 *= Xlambda;

          muF1[Idx] = gamma1*c1*(1.0 - c1)*(1.0 - 2.0*c1) - del * c1*c2*c3;
          muF1[Idx] *= (12.0 / epsilon);

          muF2[Idx] = gamma2*c2*(1.0 - c2)*(1.0 - 2.0*c2) - del * c1*c2*c3;
          muF2[Idx] *= (12.0 / epsilon);

          G1 = G1 * (1.0/gamma2 + 1.0/gamma3) - G2 / gamma2 - G3 / gamma3;
          G1 *= (2.0 * del) / epsilon;

          G2 = G2 * (1.0/gamma1 + 1.0/gamma3) - G1 / gamma1 - G3 / gamma3;
          G2 *= (2.0 * del) / epsilon;

          muF1[Idx] += G1;
          muF2[Idx] += G2;
*/
          fc1 = c1 * (1.0 - c1) * (1.0 - 2.0 * c1);
	  fc2 = c2 * (1.0 - c2) * (1.0 - 2.0 * c2);
	  fc3 = c3 * (1.0 - c3) * (1.0 - 2.0 * c3);

	  G1 = gamma1 * fc1 + 6.0 * Xlambda * c1c2c3 * c2 * c3;
	  G2 = gamma2 * fc2 + 6.0 * Xlambda * c1c2c3 * c1 * c3;
	  G3 = gamma3 * fc3 + 6.0 * Xlambda * c1c2c3 * c1 * c2;

	  muF1[Idx] = (G1 - G2) * gamma2_inv + (G1 - G3) * gamma3_inv;
	  muF2[Idx] = (G2 - G1) * gamma1_inv + (G2 - G3) * gamma3_inv;

	  muF1[Idx] *= 2.0 * del / epsilon;
	  muF2[Idx] *= 2.0 * del / epsilon;
      }
    }

    dfft(muF1, cmuF1);
    dfft(muF2, cmuF2);


}//End of the function............................................................

#endif
