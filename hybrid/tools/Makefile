.DEFAULT_GOAL = all

SOURCES := $(wildcard *.f90)
TARGETS := $(SOURCES:%.f90=%)
CLTARGETS := $(TARGETS:%=../bin/%)

FC =
HFLAGS =
LIO_PATH = ../..
include ${LIO_PATH}/lioamber/Makefile.optchck
include ${LIO_PATH}/lioamber/Makefile.options
include ../src/Makefile.rulz
HFLAGS += -L$(LIO_PATH)/g2g


all: $(TARGETS)
%: %.f90 
	$(FC) $(HFLAGS) -o ../bin/$@ $< $(LIBS)
clean:
	rm -r $(CLTARGETS)

