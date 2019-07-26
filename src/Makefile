default: hybrid

FC =
FFLAGS =
LIBS =
OBJPATH := .
LIO_PATH = 
OBJPATH = obj

include Makefile.options
include Makefile.depends

OBJECTS: $(MOD_OBJS:%.mod=$(OBJPATH)/%.mod)

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
endif
	@echo

version: $(OBJPATH)
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed  "s/SIESTA_ARCH/${SIESTA_ARCH}/g" version.f90 > temp0.f
	@sed  "s/FFLAGS/${FC} ${FFLAGS}/g" temp0.f > compinfo.F90
	@($(MAKE) -s "FC=$(FC)" "FFLAGS=$(FFLAGS)" obj/version.o)
	@rm -f temp0.f compinfo.F90
	@echo

$(FDF): $(OBJPATH)
	(cd fdf ; $(MAKE) module)


#RULZ
$(OBJPATH)/%.o :  %.f  $(FDF) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)

$(OBJPATH)/%.o : L-BFGS/%.f $(FDF) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)

$(OBJPATH)/%.o :  %.F  $(FDF) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)

$(OBJPATH)/%.o :  %.f90  $(FDF) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)

$(OBJPATH)/%.o :  %.F90  $(FDF) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)

$(OBJPATH)/%.mod :  %.f90  $(FDF) | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)

$(OBJPATH)/%.mod :  %.f  | $(OBJPATH)
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS)


hybrid: $(OBJECTS:%.o=$(OBJPATH)/%.o) $(MOD_OBJS:%.mod=$(OBJPATH)/%.mod) $(FDF) version options
	@mkdir -p ../bin
	$(FC) -cpp hybrid.f -o ../bin/$@ $(OBJECTS:%.o=$(OBJPATH)/%.o) $(MOD_OBJS:%.mod=$(OBJPATH)/%.mod) $(OBJPATH)/$(FDF) $(LIBS)
	@echo
	@echo "Compilation Succesfull "

clean:
	(cd fdf ; $(MAKE) clean)
	rm -r obj ../bin *.mod

