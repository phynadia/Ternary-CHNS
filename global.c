#include "global.h"
#ifdef VPM
  double *chi, *sponge;
#endif

    fftw_plan pfor, pinv;
    int rank, size;
    MPI_Status status;
    ptrdiff_t alloc_local, local_nx, local_x_start;

    int i, j, Idx, nshell, it, ifile, navg, fl;
    int local_x, local_size_x;
    double *ux, *uy, *u, *u1;
    double *wz, *tux, *tuy, *twz;
    double *fz, *nlz;
    double *kx, *ky, *ek, *ekx, *eky, *ensk;
    double kxmax, kymax;
    int kxm, kym;
    double t, Et, factor_x, factor_y;
    double t1, t2, tmpi;
    double *work;
    FILE *fu, *fek;
    char filename[100];

    fftw_complex *cux, *cuy, *cu;
    fftw_complex *cwz, *tcwz;
    fftw_complex *tcux, *tcuy;
    fftw_complex *cfz;
    fftw_complex *cnlz;
    fftw_complex *cnl_z;
#ifdef BINARY
    FILE *fph, *fe, *fdef;
    double coeff, gamma1, gamma2, gamma3, del;
    double *phi, *tphi, *nphi, *nlphi, *sph;
    double *ph, *aph, *bph;
    fftw_complex *cphi, *tcphi, *ncphi, *cnlphi;
    fftw_complex *cnl_phi;    
#elif defined TERNARY
    FILE *fph, *fe, *fdef1, *fdef2;
    double gamma1, gamma2, gamma3, del;
    double *phi1, *phi2, *tphi1, *tphi2;
    double *muF1, *muF2, *tphi3, *ph1, *ph2;
    double *aph1, *aph2, *bph1, *bph2;
    double *nlphi1, *nlphi2, *sph1, *sph2, *sp12;
    fftw_complex *cphi1, *cphi2, *tcphi1, *tcphi2;
    fftw_complex *cmuF1, *cmuF2;
    fftw_complex *cnl_phi1, *cnl_phi2;
    fftw_complex *cnlphi1, *cnlphi2;
#endif //TERNARY & BINARY

