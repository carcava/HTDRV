include ../make.sys

MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../FFTXlib $(MOD_FLAG)../LAXlib $(MOD_FLAG)../Modules $(MOD_FLAG).

# FLAGS for c++ with IntelMPI when QE was compiled with Intel Fortran
MPICXX=icpc
MPICXXFLAGS=-DOMPI_SKIP_MPICXX=1 -O2 -Wall -g -fPIC -I../../src -I${INTELMPI_HOME}/include64
MPILIBS=-openmp  -L${INTELMPI_HOME}/lib64  -lz  -lifcore -L$(I_MPI_ROOT)/lib64 -lmpi -lmpiif

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
OBJ=comm.o libpwscf-ht.o

all : ht-qe.x

ht-qe.x : prepinp.o $(OBJ) $(PWOBJS) $(LIBOBJS) 
	$(MPICXX) $(LDFLAGS) -o $@ $^ $(MPILIBS) $(LIBS)

%.o: %.cpp
	$(MPICXX) -c -D__USE_MPI $(MPICXXFLAGS) $< -o $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L

# explicit dependencies

comm.o: comm.h 
