include ../make.inc

MPI_DIR=/usr/local/Cluster-Apps/openmpi/gnu/1.8.8/

MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../FFTXlib $(MOD_FLAG)../LAXlib $(MOD_FLAG)../Modules $(MOD_FLAG).

# FLAGS for c++ with IntelMPI when QE was compiled with Intel Fortran
MPICXX=mpicxx
MPICXXFLAGS=-DOMPI_SKIP_MPICXX=1 -O2 -g -fopenmp -fPIC -I../../src -I${MPI_DIR}/include -fpermissive
MPILIBS=-L${MPI_DIR}/lib -lmpi_usempi -lmpi_mpifh -lmpi -lgfortran

# location of required libraries
# part 1: hi-level libraries for building pw.x
PWOBJS = \
../PW/src/libpw.a \
../Modules/libqemod.a

# part 2: lo-level libraries for all of Q-E
LIBOBJS = \
../FFTXlib/libqefft.a \
../LAXlib/libqela.a   \
../clib/clib.a   \
../iotk/src/libiotk.a

# part 3: local HT library and progams
OBJ=comm.o libpwscf-ht.o pwinput.o

all : ht-qe.x

ht-qe.x : prepinp.o $(OBJ) $(PWOBJS) $(LIBOBJS) 
	$(MPICXX) $(LDFLAGS) -o $@ $^ $(PWOBJS) $(LIBOBJS) $(LIBS) $(MPILIBS)

%.o: %.cpp
	$(MPICXX) -c -D__USE_MPI $(MPICXXFLAGS) $< -o $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L

# explicit dependencies

comm.o: comm.h 
