#include "global.h"

/*!
******************************************************************************
 * 
 * Function to initialize velocity and phase fields.
 * 
 * There are three types of initial conditions for velocity fields.
 * 1. Taylor-Green inital velocity fields
 *    
 *    \f$u_x =  sin(x)cos(y)cos(z)\f$   
 *    \f$u_y = -cos(x)sin(y)cos(z)\f$
 *
 *    \f$u_z = 0\f$
 *    
 * 2. Zero initial velocity fields.
 * 3. Constant velocity fields initialized in k-space
      (Follow: Deterministic forcing of homogeneous, isotropic turbulence
               Neal P. Sullivan, Shankar Mahalingam, and Robert M. Kerr)   
 * 4. Reading from files.
 *
 * For binary and ternary flows, the phase fields are initialized with 
 *  hyperbolic tangent profile.            
 *******************************************************************************/

void init()
{

  int *den_state, *den_state1;
  int ncube, kmax, ik;
  fftw_complex v1, v2, v3;
  double p11, p12, p22;
  double ee, vk, r1, r2, ksqr;
  double x, y, z;

  it = 0; ifile = 0; fl = 0;
  t = 0.0;

#ifdef TG_INIT	

    srand((unsigned)time(NULL));

    for(i = 0; i < local_nx; i++)
    {
      x = (i + local_x_start) * dx;
      for(j = 0; j < Ny; j++)
      {
	y = j * dy;

	Idx = IDXR(i, j);
    
        r1 = (rand() % 2) + 1;
        r2 = (rand() % 2) + 1;

	ux[Idx] =  sin(x) * cos(y);
	uy[Idx] = -cos(x) * sin(y);

	ux[Idx] *= uamp; uy[Idx] *= uamp;
          
      }	      
    }

    dfft(ux, cux); 
    dfft(uy, cuy);

    curl(cux, cuy, cwz);
    truncation(cwz);

#elif defined INIT_READ_DATA
    it = 0; ifile = 0; fl = 0;
    t = 0.0;

//    sprintf(filename, "init_field/wz.in");
//    read_array_real(wz, filename);
//    dfft(wz, cwz);
//    truncation(cwz);

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

        cwz[Idx] = 0.0;

      }
    }

#elif defined ZERO_INIT

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      { 
        Idx = IDXC(i, j);

        cwz[Idx] = 0.0;
      
      }
    }

#else

    ncube = Nx * Ny;
    kmax = round(sqrt(kxmax*kxmax + kymax*kymax));
    
    den_state =  malloc((kmax+1) * sizeof(long int));
    den_state1 = malloc((kmax+1) * sizeof(long int));

    for(i = 0; i <= kmax; i++)
    {
      den_state[i] = 0;
    }

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

	ksqr = kx[i]*kx[i] + ky[j]*ky[j];
	ik = (int) rint(sqrt(ksqr));

	if(ky[j] == 0 || ky[j] == kymax)
	  den_state[ik] += 1;	  
	else
          den_state[ik] += 2;
	
      }
    }

    MPI_Allreduce(den_state, den_state1, (kmax+1), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    srand((unsigned)time(NULL));

    for(i = 0; i < local_nx; i++)
    {
      for(j = 0; j < Ny2; j++)
      {
        Idx = IDXC(i, j);

        ksqr = kx[i]*kx[i] + ky[j]*ky[j];
        ik = round(sqrt(ksqr));
        ee = ksqr * ksqr * exp(-2.0*ksqr);
	  
	if(den_state1[ik] != 0)
	{
	  vk = sqrt(ee / (3.0*den_state1[ik]));
	}

	  vk *= ncube;

	  r1 = rand() / (double)(RAND_MAX);
	  r2 = rand() / (double)(RAND_MAX);

	  v1 = uamp * vk * exp(I * 2.0 * pi * r1);
	  v2 = uamp * vk * exp(I * 2.0 * pi * r2);

	  if(ksqr > 1e-5)
	  {

	    p11 = 1.0 - kx[i]*kx[i] / (double) ksqr;
	    p22 = 1.0 - ky[j]*ky[j] / (double) ksqr;
	    p12 = -kx[i]*ky[j] / (double) ksqr;

	  }else
	        {
	          p11 = 1.0 / 2.0;
		  p22 = p11;
		  p12 = -1.0 / 2.0;
	        }

	    cux[Idx] = p11 * v1 + p12 * v2;
	    cuy[Idx] = p12 * v1 + p22 * v2;

      }
    }

    curl(cux, cuy, cwz);

    truncation(cwz);
#endif  

#ifdef VPM

    double eps, eta_inv;
    double yd, yu, xd, xu;
    yd = 6.0 * dy;
    xd = 6.0 * dx;
    yu = ly - yd;
    xu = lx - xd;
    eps = 3.0 * dx;
    eta_inv = 1.0 / eta;
    
    for(i = 0; i < local_nx; i++)
    {
      x = (i + local_x_start) * dx;
      for(j = 0; j < Ny; j++)
      {
        y = j * dy;

        Idx = IDXR(i, j);
        
        chi[Idx] = 0.5 * (1.0 - tanh(2.0*(x - xd)/eps))
                 + 0.5 * (1.0 + tanh(2.0*(x - xu)/eps));
          
        #ifdef SPONGE
        sponge[Idx] = 0.5 * (1.0 - tanh(2.0*(y - yd)/eps))
                    + 0.5 * (1.0 + tanh(2.0*(y - yu)/eps));

        sponge[Idx] *= eta_inv;
        #endif

      }
    }

#endif

#ifdef BINARY

   #ifdef INIT_READ_DATA
      sprintf(filename, "init_field/phi.in");
      read_array_real(phi, filename);
      dfft(phi, cphi);
      truncation(cphi);
   #else
    for(i = 0; i < local_nx; i++)
    {
      x = (i + local_x_start) * dx;
      for(j = 0; j < Ny; j++)
      {
        y = j * dy;
        
        Idx = IDXR(i, j);

        phi[Idx] = 1.0 - tanh((sqrt((x-lx/2.0)*(x-lx/2.0)+(y-ly/2.0)
                              *(y-ly/2.0))-(pi/4.5)) / (0.5*epsilon));   
	
//	phi[Idx] = (1-tanh(((x-pi)*(x-pi)/(pi*pi/4.0)
//                 + (y-pi)*(y-pi)/(pi*pi/25.0) - 1.0)/(0.5*epsilon)));

        phi[Idx] *= 0.5;
       
        
      }
    }

      dfft(phi, cphi);
      truncation(cphi);
   #endif
#elif defined TERNARY

     #ifdef INIT_READ_DATA
      sprintf(filename, "init_lens/phi1.in");
      read_array_real(phi1, filename);
      dfft(phi1, cphi1);
      truncation(cphi1);

      sprintf(filename, "init_lens/phi2.in");
      read_array_real(phi2, filename);
      dfft(phi2, cphi2);
      truncation(cphi2);
     #else

    srand((unsigned)time(NULL));
    double rn1, rn2;

int seed = 10 + rank * 2;
srand(seed);

    for(i = 0; i < local_nx; i++)
    {
      x = (i + local_x_start) * dx;
      for(j = 0; j < Ny; j++)
      {
        y = j * dy;

        Idx = IDXR(i, j);
        
	rn1 = (double) rand() / (double)(RAND_MAX);
	rn2 = (double) rand() / (double)(RAND_MAX);

        phi1[Idx] = 0.0; phi2[Idx] = 0.0;

        phi1[Idx] = (1.0 / 4.0) + 0.001 * rn1;
        phi2[Idx] = (1.0 / 4.0) + 0.001 * rn2; 
          
      }
    }

      dfft(phi1, cphi1);
      dfft(phi2, cphi2);
      truncation(cphi1);
      truncation(cphi2);
    #endif
#endif

}// End of the function init().	
