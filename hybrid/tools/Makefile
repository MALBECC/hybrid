.DEFAULT_GOAL = all

SOURCES := $(wildcard *.f90)
TARGETS := $(SOURCES:%.f90=%)
CLTARGETS := $(TARGETS:%=../bin/%)

FC =
FFLAGS =
LIO_PATH = ../..
include ${LIO_PATH}/lioamber/Makefile.optchck
include ${LIO_PATH}/lioamber/Makefile.options
include ../src/Makefile.rulz
LFLAGS += -L$(LIO_PATH)/g2g


all: $(TARGETS)
%: %.f90 
	$(FC) $(LFLAGS) -o ../bin/$@ $< $(LIBS)

clean:
	rm -r $(CLTARGETS)

