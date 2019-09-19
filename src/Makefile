default: hybrid

FC =
FFLAGS =
LFLAGS =
LIBS =
SRCDIRS :=
OBJPATH = obj

include Makefile.options
include Makefile.depends

MAKEFILES =  Makefile
MAKEFILES += Makefile.depends
MAKEFILES += Makefile.options


SRCDIRS += L-BFGS

vpath %.f   $(SRCDIRS)
vpath %.F   $(SRCDIRS)
vpath %.f90 $(SRCDIRS)
vpath %.F90 $(SRCDIRS)
vpath %.o   $(OBJPATH)
vpath %.mod $(OBJPATH)
vpath %.a   $(OBJPATH)


$(OBJPATH):
	mkdir $@


options:
	@echo
	@echo "Compilation options used: "
ifeq ($(intel),2)
	@echo "Compiler : ${FC} and the MKL library for lapack"
else
	@echo "Compiler : ${FC}"
endif
ifeq (${precision},0)
	@echo "precision: ${precision} = simple"
else ifeq (${precision},1)
	@echo "precision: ${precision} = double"
else ifeq (${precision},2)
	@echo "precision: ${precision} = quad"
endif
	@echo "precision not full implemented yet"

	@echo "analytics: ${analytics}"
ifeq (${qm_level},0)
	@echo "qm_level: ${qm_level}, no QM program defined"
else ifeq (${qm_level},1)
	@echo "qm_level: ${qm_level}, using Lio as QM level"
	@echo "LIOHOME= $${LIOHOME}"
endif
	@echo
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@cp version.f90 compinfo.F90
	@($(MAKE) -s "FC=$(FC)" "FFLAGS=$(FFLAGS)" obj/version.o)
	@rm -f temp0.f compinfo.F90
	@echo

$(FDF): $(OBJPATH)
	(cd fdf ; $(MAKE) module)


#RULZ
$(OBJPATH)/%.o :  %.f  $(FDF) $(MAKEFILES) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LFLAGS) $(LIBS)

$(OBJPATH)/%.o :  %.F  $(FDF) $(MAKEFILES) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LFLAGS) $(LIBS)

$(OBJPATH)/%.o :  %.f90  $(FDF)  $(MAKEFILES) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LFLAGS) $(LIBS)

$(OBJPATH)/%.o :  %.F90  $(FDF)  $(MAKEFILES) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LFLAGS) $(LIBS)

%.mod : %.o $(FDF) | $(OBJPATH)
	@touch $@



hybrid: $(OBJECTS:%.o=$(OBJPATH)/%.o) $(OBJPATH)/$(FDF) | version options
	@mkdir -p ../bin
	$(FC) $(FFLAGS) -fopenmp hybrid.f -o ../bin/$@ $^ $(LFLAGS) $(LIBS)
	@echo
	@echo "Compilation Succesfull "

clean:
	(cd fdf ; $(MAKE) clean)
	rm -r obj ../bin *.mod

