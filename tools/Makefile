.DEFAULT_GOAL = all

SOURCES := $(wildcard *.f90)
TARGETS := $(SOURCES:%.f90=%)
CLTARGETS := $(TARGETS:%=../bin/%)

FC =
include ../src/Makefile.options
include ../src/fdf/Makefile.rulz

all: $(TARGETS)
%: %.f90 
	$(FC) $(FFLAGS) -o ../bin/$@ $< $(LIBS)
clean:
	rm -r $(CLTARGETS)

