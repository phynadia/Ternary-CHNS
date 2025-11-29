#ifndef __GLOBAL__
#define __GLOBAL__

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include <mpi.h>
#include <fftw3-mpi.h>
#include <time.h>
#include "macro.h"

#define IDXC(i, j) (j + Ny2 * i)
#define IDXR(i, j) (j + (Ny+2) * i)

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define master 0

#define pi (M_PI)
#define lx (4.0*pi)
#define ly (4.0*pi)

#define Nx (1024)
#define Ny (1024)
#define Ny2 (Ny/2+1)

#define dx (lx/Nx)
#define dy (ly/Ny)

#define nu (0.01)
#define dt (1e-4)
#define maxiter (2e5)
#define nfile (1000)

#ifdef CONSTANT_FORCING
  #define k_begin 1
  #define k_range 2
#endif

#ifdef TG_FORCING
  #define kf 3
#endif

#define uamp 1.0
#define forcing_amplitude 0.1

#ifdef VPM
  #define eta (1.5 * dt)
  #define etainv (1.0 / eta)
  extern double *chi, *sponge;
#endif

    extern fftw_plan pfor, pinv;
    extern int rank, size;
    extern MPI_Status status;
    extern ptrdiff_t alloc_local, local_nx, local_x_start;

    extern int i, j, Idx, nshell, it, ifile, navg, fl;
    extern int local_x, local_size_x;
    extern double *ux, *uy, *u, *u1;
    extern double *wz, *tux, *tuy, *twz;
    extern double *fz, *nlz;
    extern double *kx, *ky, *ek, *ekx, *eky, *ensk;
    extern double kxmax, kymax;
    extern int kxm, kym;
    extern double t, Et, factor_x, factor_y;
    extern double t1, t2, tmpi;
    extern double *work;
    extern FILE *fu, *fek;
    extern char filename[100];

    extern fftw_complex *cux, *cuy, *cu;
    extern fftw_complex *cwz, *tcwz;
    extern fftw_complex *tcux, *tcuy;
    extern fftw_complex *cfz;
    extern fftw_complex *cnlz;
    extern fftw_complex *cnl_z;

    void init();
    void dfft();
    void ifft();
    void copy();
    void curl();
    void kxky();
    void nonlin();
    void runge_kutta4();
    void energy();
    void forcing();
    void init_forcing();
    void truncation();
    void write_array_real();
    void read_array_real();
    void runge_kutta4();
    void etd2rk();
    void dump_real_mpi_IO();
    void get_velocity();
    void kinetic_energy();
    void spectrum();

#ifdef BINARY

#define sigma 2.0
#define epsilon (3.0*dx)
#define lambda 1e-6
#define radius (pi/4.0)

#ifdef GRAVITY
  #define rho1 1.044
  #define rho2 1.0
  #define g 3.0
#endif

    extern FILE *fph, *fe, *fdef;
    extern double coeff, gamma1, gamma2, gamma3, del;
    extern double *phi, *tphi, *nphi, *nlphi, *sph;
    extern double *ph, *aph, *bph;
    extern fftw_complex *cphi, *tcphi, *ncphi, *cnlphi;
    extern fftw_complex *cnl_phi;    

    void divergence();
    void surface_tension();
    void free_energy();
    void perimeter();

#elif defined TERNARY

#define sigma_12 (0.2)
#define sigma_23 (0.45)
#define sigma_13 (0.2)

#define epsilon (3.0*dx)
#define lambda (1e-4)
#define Xlambda (2.0)
#define alpha (0.0)

#ifdef GRAVITY
  #define rho1 0.957
  #define rho2 1.044
  #define rho3 1.0
  #define g 3.0
#endif

#define radius (pi/4.0)

    extern FILE *fph, *fe, *fdef1, *fdef2;
    extern double gamma1, gamma2, gamma3, del;
    extern double *phi1, *phi2, *tphi1, *tphi2;
    extern double *muF1, *muF2, *tphi3, *ph1, *ph2;
    extern double *aph1, *aph2, *bph1, *bph2;
    extern double *nlphi1, *nlphi2, *sph1, *sph2, *sp12;
    extern fftw_complex *cphi1, *cphi2, *tcphi1, *tcphi2;
    extern fftw_complex *cmuF1, *cmuF2;
    extern fftw_complex *cnl_phi1, *cnl_phi2;
    extern fftw_complex *cnlphi1, *cnlphi2;

    void divergence();
    void surface_tension();
    void potential();
    void free_energy();
    void perimeter();

#endif //TERNARY & BINARY

#endif //DOXYGEN_SKIP

#endif //__GLOBAL__
