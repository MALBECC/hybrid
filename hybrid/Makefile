#
# Makefile for HYBRID 
#
.SUFFIXES: .f .F .o .a  .f90 .F90
#
default: what hybrid 
#
include arch.make
#
#AGREGAR ESTA LINEA SI DA ERROR DE M2C
%.o: %.mod
#------------------- 
what:
	@echo
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 1

# Note that machine-specific files are now in top Src directory.
OBJS =  chkdim.o conjgr.o dynamics.o fixed.o iofa.o ofc.o functions.o \
	inver.o iorho.o ioxv.o iovext.o memory.o paste.o reclat.o \
	reinit.o timer.o timestamp.o vmb.o volcel.o cgvc.o iocg.o \
	pcpot.o centerall.o ran3.o dot.o readall.o subconstr.o \
	subfree.o writeall.o qmmm_lst_blk.o solv_assign.o \
	mmforce.o ljef.o solv_ene_fce.o linkatom.o recipes.o \
	assignation.o mulliken.o elecfield.o hybrid.o \

# FDF Module
FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) module)

# Routines using fdf calls.
iorho.o ioxv.o iovext.o iofa.o : $(FDF)
cgvc.o iocg.o ioxvconstr.o ofc.o : $(FDF)
fixed.o reinit.o readall.o writeall.o : $(FDF)
subconstr.o subfree.o : $(FDF)
linkatom.o solv_assign.o elecfield.o : $(FDF)
asign.o : $(FDF)

# Rest of Modules
# Note : precision must be the first module
MOD_OBJS=precision.o \
	 ionew.o \
         memoryinfo.o \
	 parsing.o \
	 parallel.o \
	 sys.o \
	 fsiesta.o \
	 version.o \
	 neutralatom.o \

COM_OBJS=$(OBJS) $(SYSOBJ)
ALL_OBJS=$(MOD_OBJS) $(COM_OBJS)

# Routines using modules
memory.o: memoryinfo.o
solv_assign.o timer.o vmb.o iorho.o iovext.o ioxv.o : precision.o
qmmm_lst_blk.o readall.o cgvc.o chkdim.o dynamics.o : precision.o
iofa.o fixed.o ioxvconstr.o assign.o hybrid.o : precision.o
ionew.o chkdim.o solv_assign.o hybrid.o readall.o qmmm_lst_blk.o fixed.o assign.o : sys.o
subconstr.o subfree.o :sys.o
iovext.o iorho.o : parallel.o
hybrid.o : fsiesta.o
hybrid.o mmforce.o : neutralatom.o
neutralatom.o hybrid.o solv_assign.o :ionew.o
reinit.o readall.o subfree.o writeall.o : ionew.o
qmmm_lst_blk.o mulliken.o iovext.o iorho.o : ionew.o
iofa.o iocg.o ioxv.o ioxvconstr.o ofc.o : ionew.o
assign.o : parse.o

version:
	@echo
	@echo "==> Incorporating information about present compilation (compiler and flags)"
	@sed  "s/SIESTA_ARCH/${SIESTA_ARCH}/g" version.f > temp0.f
	@sed  "s/FFLAGS/${FC} ${FFLAGS}/g" temp0.f > compinfo.F90
	@($(MAKE) -s "FC=$(FC)" "FFLAGS=$(FFLAGS)" compinfo.o)
	@rm -f temp0.f compinfo.F90
	@mv compinfo.o version.o
	@echo

hybrid: version $(FDF) $(COMP_LIBS) $(ALL_OBJS)
	$(FC) -o hybrid \
               $(LDFLAGS) $(ALL_OBJS) $(FDF) $(COMP_LIBS) $(LIBS)

siesta:
	@( cd ../../../Src/; $(MAKE) siesta) 

clean: 
	@echo "==> Cleaning object, library, and executable files"
	rm -f hybrid *.o  *.a *.pcl *.pc fdf.Standard
	rm -f *.mod
	(cd fdf ; $(MAKE) clean)
