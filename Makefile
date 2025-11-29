CC = mpicc
LIBS = -lfftw3_mpi -lfftw3 -lm

all:chns2d.exe

chns2d.exe: global.o main.o init.o fft.o functions.o curl.o nonlin.o etd2rk.o forcing.o surface_tension.o write_array.o read_array.o chem_potential.o velocity.o free_energy.o kinetic_energy.o spectrum.o perimeter.o
	$(CC) -o $@ $+ $(LIBS)

clean:
	$(RM) *.o *.exe a.out 

scrub: clean
	$(RM) *.dat
