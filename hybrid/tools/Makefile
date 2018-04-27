.DEFAULT_GOAL = all
include ../Makefile.arch

SOURCES := $(wildcard *.f90)
TARGETS := $(SOURCES:%.f90=%)
CLTARGETS := $(TARGETS:%=../bin/%)


all: $(TARGETS)
%: %.f90 
	$(FC) -o ../bin/$@ $<

clean:
	rm -r $(CLTARGETS)

