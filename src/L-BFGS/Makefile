default: lbfgsb

FC =
#FFLAGS =
#OBJPATH := .
LIO_PATH = ../../lio
include ${LIO_PATH}/lioamber/Makefile.optchck
include ${LIO_PATH}/lioamber/Makefile.options

#FC = gfortran
FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized -c

LINPACK = linpack.f
BLAS    = blas.f
TIMER   = timer.f

lbfgsb: $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(LINPACK) $(BLAS) $(TIMER) 

clean: 
	rm -f *.o  

