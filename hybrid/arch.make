# arch.make for ifort compiler with mkl libraries

SIESTA_ARCH=ifc-mkl
FC=ifort
FFLAGS=  -w -O2 -mp

LIBS= -L/opt/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64 -lmkl_core -lpthread -lmkl_intel_lp64 -lmkl_sequential

SYS=bsd
RANLIB=ranlib
.F.o:
	$(FC) -c $(FFLAGS) $(DEFS) $<
.f.o:
	$(FC) -c $(FFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(FREE_F90) $<
.F90.o:
	$(FC) -c $(FFLAGS) $(DEFS) $(FREE_F90_CPP) $<
