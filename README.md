  Source code for CHNS-2D-MPI
 
  Author: Nadia Bihari Padhan (<nadia_bihari.padhan@tu-dresden.de>).
 
  Institute: Institute of Scientific Computing, TU Dresden, Germany.
  
  cite at: 1. Hyperuniformity in ternary fluid mixtures: the role of wetting and hydrodynamics
              Nadia Bihari Padhan† and Axel Voigt, arXiv:2506.22647
           2. Unveiling the spatiotemporal evolution of liquid-Lens coalescence: Self-similarity, vortex quadrupoles, and turbulence in a three-phase fluid system
              Nadia Bihari Padhan† and Rahul Pandit
              Physics of Fluids, 35, 112105 (2023)
           3. The Cahn-Hilliard-Navier-Stokes Framework for Multiphase Fluid Flows: Laminar, Turbulent, and Active
              Nadia Bihari Padhan and Rahul Pandit
              Journal of Fluid Mechanics, vol. 1010, P1,(2025)
  
  DNSs of coupled Cahn-Hilliard-Navier-Stokes equations using pseudospectral 
  method.
 
  main program: main.c
  parameters and global variables: global.h
  all macros: macro.h
  
  Time integration: etdRK-2 method.
 
  The program involves calculations for Navier-Stokes equations, Binary flows,
  and ternary flows in 2D. Binary and ternary flows are employed using phase-field method.
  It is developed in C programming language and MPI-platform.

  Install FFTW library compatible with MPI.
  In Makefile provide proper directory name for FFTW library.
  After each changes in .h files, use make clean before compiling with make.

  Uncomment the necessary macros in macro.h.
