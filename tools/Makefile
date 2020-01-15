.DEFAULT_GOAL = all

SOURCES := $(wildcard *.f90)
TARGETS := $(SOURCES:%.f90=%)
CLTARGETS := $(TARGETS:%=../bin/%)

FC =
intel     ?= 0
analytics ?= 0
################################################################################
ifeq ($(intel),0)
  FC      = gfortran
  LIBS   += -llapack
  DEBUG01 = -g -Wno-tabs -fcheck=all
  DEBUG02 = -g -Wall -Wextra -fcheck=all -fbacktrace -Wtabs
else ifeq ($(intel),1)
  FC      = ifort
  LIBS   += -llapack
  DEBUG01 = -g -check bounds -traceback
  DEBUG02 = -g -check bounds -traceback -check all -fp-stack-check
else ifeq ($(intel),2)
  FC      = ifort
  LC      = ifort
  FFLAGS += -fpp -module $(OBJPATH)
  LFLAGS += -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include
  LIBS   += -lmkl_blas95_lp64 -lmkl_lapack95_lp64
  LIBS   += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
  LIBS   += -liomp5 -lpthread
  DEBUG01 = -g -check bounds -traceback
  DEBUG02 = -g -check bounds -traceback -check all -fp-stack-check
else
  $(error Invalid value in option intel = $(intel))
endif

################################################################################
ifeq ($(analytics),0)
  FFLAGS += -O3
else ifeq ($(analytics),1)
  DEFINE += -DPROFILING
  FFLAGS += -pg
else ifeq ($(analytics),2)
  DEFINE += -DDEBUGGING -DPRINT_MATRICES
  FFLAGS += $(DEBUG01)
else ifeq ($(analytics),3)
  DEFINE += -DDEBUGGING -DPRINT_MATRICES
  FFLAGS += $(DEBUG02)
else
  $(error Invalid value in option analytics = $(analytics))
endif

all: $(TARGETS)
%: %.f90 
	$(FC) $(FFLAGS) -o ../bin/$@ $< $(LFLAGS) $(LIBS)
clean:
	rm -r $(CLTARGETS)

