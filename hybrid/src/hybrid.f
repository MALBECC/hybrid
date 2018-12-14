      Program HYBRID 
!****************************************************************************
! Program for QM-MM calculations. 
! Original version used the SIESTA code 
! (https://departments.icmab.es/leem/siesta/CodeAccess/index.html) to treat
! at DFT level the QM subsystem. 
! This version(2017) has been modified for use Lio code to treat at DFT level
! the QM-MM subsystem.
! It uses an own implementation of Amber99 force field parametrization to treat
! the MM subsystem. 
! Original idea by D. Scherlis, D. Estrin and P. Ordejon. 2001.
! Original QM-MM interfase with Siesta by D. Scherlis and P. Ordejon. 2001.
! Original subroutines by D. Scherlis. 2001/2002.
! Solvent implementation using the Amber99 force field parametrization
! by A. Crespo and M. Marti. 2002/2003.
! LinkAtom subrouitnes by M. Marti. 2003.
! Further modifications of QM-MM Siesta by A. Crespo and P. Ordejon. 2004.
! Developing of HYBRID: program for QM-MM calculations using the SIESTA
! as SUBROUTINE by A. Crespo and P. Ordejon. 2004.
! This code uses several subroutines of SIESTA  as well as the 
! SiestaAsSubroutine and FDF packages. Crespo, May 2005. 
! Modified vesion for use Lio as QM-MM level and debugs, N. Foglia 2017
!*****************************************************************************

! Modules
      use precision, only: dp
      use sys, only: die
      use fdf
      use ionew, only: io_setup, IOnode    

      implicit none
! General Variables
      integer :: istep, inicoor,fincoor !actual, initial and final number of move step for each restrain
      integer :: idyn !kind of movent, idyn=0 (minimization) only avalable at this moment
      integer :: istp !number of move step for each restrain starting on 1
      integer :: nmove !max number of move step for each restrain
      integer :: na_u !number of QM atoms
      integer :: nesp !number of QM species
      integer :: step !total number of steps in a MMxQM movement (not tested in this version)
      integer :: ifmax(2) !number of atom with max force on each cicle
      real(dp) :: fmax !max force on each cicle 
      integer :: icfmax(2) !number of atom with max force on each cicle after constrain
      real(dp) :: cfmax !max force on each cicle after constrain
      real(dp) :: cstress(3,3) !Stress tensor, may be included in CG minimizations 
      real(dp) :: strtol !Maximum stress tolerance
      real(dp) :: Etot ! QM + QM-MM interaction energy
      real(dp) :: fres !sqrt( Sum f_i^2 / 3N )
      real(dp) :: ftol !MAx force tol criteria (in Ry/Bohr)
      real(dp) :: ftot(3) ! Total force in system
      real(dp) :: dxmax !Maximum atomic displacement in one CG step (Bohr)
      real(dp) :: tp !Target pressure
      real(dp) :: volume !Cell volume
      real(dp), dimension(:,:), allocatable :: xa, fa !position and forces of QM atoms
      integer, dimension(:), allocatable :: isa,iza !Chemical Specie Label, and atomic charge
      character, dimension(:), allocatable :: atsym*2 !atomic symbol
      logical :: usesavexv, foundxv !control for coordinates restart
      logical :: usesavecg !control for restart CG
      logical :: varcel !true if variable cell optimization
      logical :: relaxd ! True when CG converged
      logical :: qm, mm ! True when system have a subsystem QM,MM
      character :: slabel*20 ! system label, name of outputs
      character :: sname*150 !name of system
      character :: paste*25
      external :: paste
! Solvent (MM) General variables
      integer :: nac !number of MM atoms
      character*4,  dimension(:), allocatable, save :: atname,aaname,attype,qmattype


! atom, j, atname(i), aaname(i), ch, resnum(i), rclas(1:3,na_u+i)
! ATOM     12  CE  MET H   1       7.246  15.955   0.940


! Outputs names
      character :: xyzlabel*20

! Conversion factors
      real(dp) :: Ang !r_in_ang=r_in_bohr * Ang
      real(dp) :: eV !E_in_Hartree=E_in_eV * eV
      real(dp) :: kcal ! E_in_kcal_mol-1 = E_in_eV * kcal
! Auxeliars
      integer :: i, ia, iunit, ix, j, k
!comentando variables, arriba


! Others that need check
!!!! General Variables
      real(dp) :: dt
      real(dp) :: ucell(3,3)
      real(dp) :: volcel
      external :: volcel
!!!! Solvent General variables



C General Variables
      external
     .  chkdim, cgvc,  ioxv, fixed1, assign, 
     .  iofa, ofc, reinit, read_qm, 
     .  read_md, fixed2

C Solvent General variables


      character*5, dimension(:), allocatable, save::
     .  bondtype

      character*8, dimension(:), allocatable, save::
     .  angletype

      character*11, dimension(:), allocatable, save::
     .  dihetype,imptype

      real(dp), dimension(:), allocatable, save ::
     .  kbond,bondeq,kangle,angleeq,kdihe,diheeq,
     .  perdihe,kimp,impeq,perimp,desplaz

      integer, dimension(:,:,:), allocatable, save ::
     .  atange,atangm,atdihe,atdihm,atimp

      integer, dimension(:,:), allocatable, save ::
     .  ng1

      integer, dimension(:), allocatable, save ::
     .  aanum,multidihe,multiimp, bondxat,angexat,
     .  dihexat,dihmxat,angmxat,impxat,
     .  scalexat,nonbondedxat,izs, qmstep

      integer natot,nfree,nparm,nfce, mmsteps, imm, 
     .  nroaa, nbond, nangle, ndihe, nimp , wricoord,
     .  atxres(20000), totcoor

      double precision  Etot_amber,sfc,timestep,
     .  Etots,rcorteqm,rcortemm,Elj,Eelec

      double precision, dimension(:,:), allocatable, save::
     .  fce_amber, qmpc, fdummy, cfdummy, rclas, vat

      double precision, dimension(:), allocatable, save::
     .  qmpcf, Em, Rm, pc, masst

      integer, dimension(:,:,:), allocatable, save::
     .  evaldihe,evaldihm
      integer, dimension(:,:), allocatable, save::
     .  ng1type,angetype,angmtype,dihety,dihmty,impty,
     .  nonbonded,scale
      logical, dimension(:,:), allocatable, save::
     .  evaldihelog,evaldihmlog

      logical firstpc,constropt,actualiz, free, deuter 

C Link Atom variables
        logical   linkatom
        integer numlink,linkqm(15,4),linkmm(15,4),linkat(15),
     .  linkmm2(15,4,3),parametro(15,22,4)
        double precision pclinkmm(15,15),Emlink(15,4),
     .  Elink,distl(15),kch(15),rch(15)
        character linkqmtype(15,4)*4

C ConstrOpt variables
        integer nconstr,typeconstr(20),atmsconstr(20,20)
        integer nstepconstr,istepconstr, ndists(20)
        double precision kforce(20),ro(20),rt(20),rini,rfin,dr,
     .                   coef(20,10)
 
C Cut Off QM-MM variables
        double precision rcorteqmmm, radbloqmmm
        integer, dimension(:), allocatable, save::
     .  blocklist,blockqmmm,listqmmm
c nuevas, Nick
	double precision, allocatable, dimension(:,:) :: r_cut_QMMM,
     .  F_cut_QMMM
	double precision, allocatable, dimension(:) :: Iz_cut_QMMM
	double precision :: r12
	integer :: at_MM_cut_QMMM, r_cut_pos
	integer, allocatable, dimension(:) :: r_cut_list_QMMM
	integer :: i_qm, i_mm
        logical :: done, done_freeze, done_QMMM
        logical :: do_SCF, do_QM_forces
	logical, allocatable, dimension(:) :: MM_freeze_list
c cantidad de e desapareados
	double precision :: spin 
C variables para esquema de optimizacion
        integer :: opt_scheme
        integer :: optimization_lvl
C 1- solo considera restrain
C 2- solo opimiza atomos MM
C 3- optimiza todo
	integer ::  itest, atot_block

C activa-desactiva calculo de propiedades en lio
	logical :: do_properties

C more variables
        logical water,foundcrd,foundvat,writeipl 

C Solvent external variables
        external
     .  solv_assign, solv_ene_fce, qmmm_lst_blk, wrtcrd,
     .  centermol, centerdyn, link1, link2, link3, ljef,
     .  mmForce, readcrd,  prversion, ioxvconstr,  
     .  wripdb, wriene, wrirtc, subconstr1, subconstr2, subconstr3

C agragados Nick para Lio
	integer :: charge, inick, jnick

C para centrado
	double precision, dimension(3) :: Rcm
	double precision :: Tmass
C parche para omitir bonds en extremos terminales
        double precision :: radblommbond

C--------------------------------------------------------------------

	opt_scheme=1
	do_SCF=.true.
	do_QM_forces=.true.
	rcorteqmmm=0.d0
	radbloqmmm=0.d0
	do_properties=.false.

C Initialise some variables
        relaxd=.false.
        varcel=.false.
        tp=0.0
        cstress=0.0
        strtol=0.0
        qm=.true.
        mm=.true.
C added init, Nick
	Etot_amber=0.d0
	charge=0
	spin=0.d0

C Initialize IOnode
      call io_setup   

C Print version information
      call prversion
      write(6,'(/,a,i4,a)') 'hybrid: Running in serial mode' 

C Start time counter 
      call timestamp('Start of run')      

C Factors of conversion to internal units 
      Ang    = 1._dp / 0.529177_dp
      eV     = 1._dp / 27.211396132_dp
      kcal   = 1.602177E-19_dp * 6.02214E23_dp / 4184.0_dp

C Initialise read 
      call reinit(slabel, sname)

C Read the number of QM atoms
      na_u=fdf_integer('NumberOfAtoms',0)
      if (na_u.eq.0) then
        write(6,'(/a)') 'hybrid: Running with no QM atoms'
        qm=.false.
      endif
C Read the number of MM atoms
      nac = fdf_integer('NumberOfSolventAtoms',0)
      if (nac.eq.0) then
        write(6,'(/a)') 'hybrid: Running with no MM atoms'
        mm=.false.
      endif
C Read the number of species
      nesp = fdf_integer('NumberOfSpecies',0)
      if(qm.and.nesp.eq.0) then
      call die("hybrid: You must specify the number of species")
      endif

      allocate(xa(3,na_u))
      allocate(fa(3,na_u))
      allocate(isa(na_u))
      allocate(iza(na_u))
      allocate(atsym(nesp))
 
C Read QM coordinates and labels
      write(6,*)
      write(6,"('read:',71(1h*))")
      if(qm) then
        call read_qm(na_u,nesp,isa,iza,xa,atsym,charge, spin)
      endif !qm

C Allocation of solvent variables
      natot = nac + na_u

	allocate(r_cut_list_QMMM(nac))
c referencia posiciones de atomos MM con los vectores cortados

      nparm = 500
c nparm en el numero de tipos de bonds q tiene definido el amber.parm. NO DEBERIA ESTAR fijo, Nick

      allocate(izs(natot))
      allocate(Em(natot))
      allocate(Rm(natot))
      allocate(pc(0:nac))
ccambiado, Nick
      allocate(rclas(3,natot))

      allocate(MM_freeze_list(natot))

      allocate(masst(natot))
      allocate(vat(3,natot))
      allocate(cfdummy(3,natot))
      allocate(fdummy(3,natot))
      allocate(qmattype(na_u))
      allocate(attype(nac))
      allocate(atname(nac))
      allocate(aaname(nac))
      allocate(aanum(nac))
      allocate(ng1(nac,6))
      allocate(blocklist(natot))
      allocate(blockqmmm(nac))
      allocate(listqmmm(nac))
      allocate(fce_amber(3,nac))
      allocate(ng1type(nac,6))
      allocate(angetype(nac,25))
      allocate(angmtype(nac,25))
      allocate(evaldihe(nac,100,5))
      allocate(evaldihm(nac,100,5))
      allocate(dihety(nac,100))
      allocate(dihmty(nac,100))
      allocate(impty(nac,25))
      allocate(nonbonded(nac,100))
      allocate(scale(nac,100))
      allocate(evaldihelog(nac,100))
      allocate(evaldihmlog(nac,100))
      allocate(qmpc(na_u,2))
      allocate(qmpcf(na_u))
      allocate(desplaz(nac))
      allocate(scalexat(nac))
      allocate(nonbondedxat(nac))
      allocate(kbond(nparm),bondeq(nparm),
     .           bondtype(nparm))
      allocate(kangle(nparm),angleeq(nparm),
     .           angletype(nparm))
      allocate(kdihe(nparm),diheeq(nparm),
     .           dihetype(nparm),multidihe(nparm),
     .            perdihe(nparm))
      allocate(kimp(nparm),impeq(nparm),
     .           imptype(nparm),multiimp(nparm),
     .            perimp(nparm))
      allocate(atange(nac,25,2))
      allocate(atangm(nac,25,2))
      allocate(atdihe(nac,100,3))
      allocate(atdihm(nac,100,3))
      allocate(bondxat(nac))
      allocate(angexat(nac))
      allocate(dihexat(nac))
      allocate(dihmxat(nac))
      allocate(angmxat(nac))
      allocate(impxat(nac))
      allocate(atimp(nac,25,4))

C Some definitions 
      ucell=0.d0
      fa=0.d0
      fdummy = 0.d0
      cfdummy=0.d0
      vat=0.d0
      Elj=0.d0
      Elink=0.d0
      Etot=0.d0
      Etots=0.d0
      Eelec=0.d0
      nfree=natot
      blocklist = 0
      blockqmmm = 0
      listqmmm = 0
      step = 0
      nstepconstr = 0
      numlink = fdf_integer('LinkAtoms',0)
      linkatom = .false.
      if(numlink.ne.0) linkatom = .true.
      constropt = .false.
      constropt = fdf_block('ConstrainedOpt',iunit)
      free = .false.
      free = fdf_block('FreeEnergy',iunit)
      foundvat = .false.
      writeipl = fdf_boolean('WriIniParLas',.false.)

C Read and assign Solvent variables 
       if(mm) then
       call solv_assign(na_u,natot,nac,nroaa,Em,Rm,attype,pc,
     .  ng1,bondxat,angexat,atange,angmxat,atangm,dihexat,atdihe,
     .  dihmxat,atdihm,impxat,atimp,
     .  nbond,kbond,bondeq,bondtype,
     .  nangle,kangle,angleeq,angletype,
     .  ndihe,kdihe,diheeq,dihetype,multidihe,perdihe,
     .  nimp,kimp,impeq,imptype,multiimp,perimp,
     .  nparm,aaname,atname,aanum,qmattype,rclas,
     .  rcorteqmmm,rcorteqm,rcortemm,sfc,timestep,
     .  radbloqmmm,atxres,radblommbond)
       endif !mm
c cambio cutoff a unidades atomicas, Nick
      rcorteqmmm=rcorteqmmm*Ang
      radbloqmmm=radbloqmmm*Ang


       rclas(1:3,1:na_u) = xa(1:3,1:na_u)

C Read simulation data 
       call read_md( idyn, nmove,
     .               dt, dxmax, ftol, 
     .               usesavecg, usesavexv , na_u,  
     .               natot, nfce, wricoord, mmsteps )

C assignation of masses and species 
       call assign(na_u,nac,atname,iza,izs,masst)

C Read cell shape and atomic positions from a former run
       call ioxv( 'read', natot, ucell, rclas, vat, foundxv, foundvat ) 
       if (foundxv) xa(1:3,1:na_u)=rclas(1:3,1:na_u)

C Reading LinkAtom variables
      if(qm.and.mm) then
        if (linkatom) then
        call link1(numlink,linkat,linkqm,linkmm,linkqmtype,
     .            linkmm2,ng1,nac,na_u,qmattype,rclas)
c sets to zero Em for HL and CQM
        do i=1,numlink
          Em(linkat(i))=0.d0
          Em(linkqm(i,1:1))=0.d0
        enddo
        endif !LA
      endif !qm & mm

C Read a crd file from a former run
      call readcrd(na_u,nac,masst,linkatom,linkat,numlink,
     .             rclas,vat,foundcrd,foundvat)
      if(foundcrd) xa(1:3,1:na_u)=rclas(1:3,1:na_u)



C Sets LinkAtoms' positions
      if(qm.and.mm) then
        if (linkatom) then
        call link3(numlink,linkat,linkqm,linkmm,rclas,
     .    natot,na_u,nac,distl)
        xa(1:3,1:na_u)=rclas(1:3,1:na_u)
C Sets to zero pc and Em for MMLink 
          do i=1,numlink
            pclinkmm(i,1:4)=pc(linkmm(i,1:4))
            pc(linkmm(i,1:1))=0.d0
            pc(linkmm(i,2:4))=pc(linkmm(i,2:4))+pclinkmm(i,1)/3.d0
C two options: Em(linkmm(i,1:1))=0.0 or Em(linkmm(i,1:4))=0.0
            Emlink(i,1:1)=Em(na_u+linkmm(i,1:1))
            Em(na_u+linkmm(i,1:1))=0.d0
          enddo
        endif !LA
      endif !qm & mm

C Dump initial coordinates to output 
      if (writeipl) then
        if(qm) then
          write(6,'(/a)') 'hybrid: Atomic coordinates (Ang) and species'
          write(6,"(i6,2x,3f10.5,2x,i3)")
     .         (ia, (xa(ix,ia)/Ang, ix=1,3), isa(ia), ia=1,na_u)
        endif !qm
C Dump initial Solvent coordinates to output
        if(mm) then      
          write(6,*)
          write(6,"('hybrid: Number of solvent atoms:',i6)") nac
          write(6,"('hybrid: Number of residues:',i4)") nroaa
          write(6,'(/a)')'hybrid: Solvent coordinates (Ang)'
          write(6,"(i6,2x,A4,2x,A4,2x,i5,2x,3f12.6)")
     .     (i,atname(i),aaname(i),aanum(i),
     .     (rclas(ix,i+na_u)/Ang, ix=1,3),i=1,nac)
          write(6,'(/a)') 'hybrid: Atoms parameters'
        endif !mm

        if(qm.and.mm) then
          write(6,"(i6,2x,A4,2x,f12.7,2x,f12.7)")
     .     (i,qmattype(i),Em(i),Rm(i),i=1,na_u)
        endif !qm & mm

        if(mm) then
          write(6,"(i6,2x,A4,2x,f9.6,2x,f9.6,2x,f9.6)")
     .     (i,attype(i),Em(i+na_u),Rm(i+na_u),pc(i),i=1,nac)
        endif !mm
      endif !writeipl
      call flush(6)

c Initialize .xyz
	xyzlabel = paste( slabel, '.xyz' )
	open(unit=34, file=xyzlabel)




C Initialize siestaslabel, filerho and filevqm 
      if(qm) then
! First call to siesta_forces to initialise grid
        Etot=0.d0
        fa=0.d0

! initialize lio 
        call init_lio_hybrid(na_u, nac, charge, iza, spin)

! calculate cell volume
        volume = volcel( ucell )

C writes cell
        write(6,'(/,a,3(/,a,3f12.6))')
     .       'hybrid: unit cell vectors (Ang) from siesta run:',
     .       ('    ' , (ucell(ix,ia)/Ang,ix=1,3), ia =1,3)

C Center system 
        if (.not.foundxv) call centermol(na_u,xa,rclas,ucell,natot)
      endif !qm

C Calculate Rcut & block list QM-MM 
      if(qm.and.mm) then
c esto hay q volarlo y hacer otra lista de cuts, Nick
c        call qmmm_lst_blk(na_u,nac,natot,nroaa,atxres,rclas,
c     .  rcorteqmmm,radbloqmmm,blockqmmm,listqmmm,rcorteqm,slabel)
      endif !qm & mm

C Read fixed atom constraints
      call fixed1(na_u,nac,natot,nroaa,rclas,blocklist,
     .            atname,aaname,aanum,water)




c########################################################################################
c########################################################################################
c########################################################################################


C Start loop over constrained optimization steps
        if(constropt) then
         call subconstr1(nconstr,typeconstr,kforce,nstepconstr,
     .   rini,rfin,atmsconstr,dr,ro,ndists,coef,constropt)
        endif
        allocate(qmstep(nstepconstr+1))
        qmstep=0

        do istepconstr=1,nstepconstr+1   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< RESTRAIN
c istepconstr marcA LA POSICION DEL RESTRAIN, NICK

	optimization_lvl=3

	if (opt_scheme .eq. 1 .and. .false.) then
	  optimization_lvl=1
c	  do_SCF=.false.
c	  do_QM_forces=.false.
	end if

        if(constropt) then
          write(6,*)
          write(6,'(A)')    '*******************************'
          write(6,'(A,i5)') '  Constrained Step : ', istepconstr
          write(6,'(A)')    '*******************************'
	  write(666,*) "paso ", istepconstr
        endif

C Begin of coordinate relaxation iteration ============================
      if (idyn .eq. 0) then
        inicoor = 0
        fincoor = nmove
      endif
        totcoor = fincoor - inicoor + 1





	at_MM_cut_QMMM=nac

C Aca empieza a mover para 1 valor de constrain, Nick

C Start loop over coordinate changes
      istp = 0
      do istep = inicoor,fincoor    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Paso opt
      istp = istp + 1
 
        write(6,'(/2a)') 'hybrid:                 ',
     .                    '=============================='

        if (idyn .ne. 0 ) STOP 'only minimization avalable'

        write(6,'(28(" "),a,i6)') 'Begin CG move = ',istep
        write(6,'(2a)') '                        ',
     .                    '=============================='


	write(*,*) "nivel de optimizacion: ", optimization_lvl


! Calculate Energy and Forces using Lio as Subroutine
      if(qm) then

c cutofff QM-MM
c	if (.true.) then !activa el cutoff

	  if (istep.eq.inicoor) then ! define lista de interacciones en el primer paso
!falta una condicion q se acive si pido cuts

!mod(istep,10).eq.1) then !actualiza cuts cada 10 pasos

ca.true.) then !activa el cutoff
cistep
cif (mod(i,10).eq.1) write(*,*) "i, mod ", i, mod(i,10)

	   write(*,*) "defino cuts "
	    if (allocated(r_cut_QMMM)) deallocate(r_cut_QMMM)
	    if (allocated(F_cut_QMMM)) deallocate(F_cut_QMMM)
	    if (allocated(Iz_cut_QMMM)) deallocate(Iz_cut_QMMM)
	    r_cut_list_QMMM=0
	    r_cut_pos=0
	    at_MM_cut_QMMM=0
	    if (istepconstr.eq.1) MM_freeze_list=.false.

	    do i_mm=1, nac
	      i_qm=0
	      done=.false.
              done_freeze=.false.
             done_QMMM=.false.
	      do while (i_qm .lt. na_u .and. .not. done)
	        i_qm=i_qm+1
                r12=(rclas(1,i_qm)-rclas(1,i_mm+na_u))**2.d0 +
     .              (rclas(2,i_qm)-rclas(2,i_mm+na_u))**2.d0 +
     .              (rclas(3,i_qm)-rclas(3,i_mm+na_u))**2.d0

	        r12=sqrt(r12) !distancia atomo MM a aqomo QM
c		write(*,*) "dist r12", i_mm, r12
	        if(r12 .lt. rcorteqmmm .and. .not. done_QMMM) then !luego cambiar cut, cut esta en bohrs, Cuidado Nick
	          done_QMMM=.true.
	          at_MM_cut_QMMM=at_MM_cut_QMMM+1
	          r_cut_pos=r_cut_pos+1
	          r_cut_list_QMMM(i_mm)=r_cut_pos
	        end if
c aca agregar freeze list

		if (istepconstr.eq.1) then !define lista de movimiento para la 1er fotoa
c		MM_freeze_list=.false.
	          if(r12 .gt. radbloqmmm .and. .not. done_freeze) then
	             MM_freeze_list(i_mm+na_u)=.true.
	             done_freeze=.true.
	          end if
		end if
		done=done_QMMM .and. done_freeze
	      end do
	    end do


	  allocate (r_cut_QMMM(3,at_MM_cut_QMMM+na_u),
     .    F_cut_QMMM(3,at_MM_cut_QMMM+na_u),
     .    Iz_cut_QMMM(at_MM_cut_QMMM+na_u))

	  r_cut_QMMM=0.d0
	  F_cut_QMMM=0.d0
	  Iz_cut_QMMM=0

c	end if! fin lista solo paso 1


c		write(*,*) "pc", pc
	  end if


c	if (optimization_lvl.gt.1) then !acivo para esquema de optimizacion
	  if (.true.) then !activo solo si quiero cutoff

!copy position and nuclear charges to cut-off arrays
            do i=1,natot !barre todos los atomos
              if (i.le.na_u) then !QM atoms
                r_cut_QMMM(1:3,i)= rclas(1:3,i)
c               Iz_cut_QMMM(i)= pc(i)
              else if (r_cut_list_QMMM(i-na_u) .ne. 0) then !MM atoms inside cutoff
		r_cut_QMMM(1:3,r_cut_list_QMMM(i-na_u)+na_u) = rclas(1:3,i)
                Iz_cut_QMMM(r_cut_list_QMMM(i-na_u))= pc(i-na_u)
              end if
            end do


 
	   call SCF_hyb(na_u, at_MM_cut_QMMM, r_cut_QMMM, Etot, 
     .     F_cut_QMMM,
     .           Iz_cut_QMMM, do_SCF, do_QM_forces, do_properties) !fuerzas lio, Nick


c return forces to fullatom arrays
            do i=1, natot
	      if (i.le.na_u) then !QM atoms
	         fdummy(1:3,i)=F_cut_QMMM(1:3,i)
	      else if (r_cut_list_QMMM(i-na_u).ne.0) then !MM atoms in cut-off
                 fdummy(1:3,i)=
     .           F_cut_QMMM(1:3,r_cut_list_QMMM(i-na_u)+na_u)
	      end if
            end do

c	stop "chekeo de cutoff"

	  else   !sin cut de movimiento
  	    if (nac .gt. 0) then
c luego cambiar con el cut de sdistancia
             call SCF_hyb(na_u, nac, rclas, Etot, fdummy, pc(1:nac),
     .            do_SCF, do_QM_forces, do_properties) !fuerzas lio, Nick
c atomos totales, atomos MM, posiciones, enerfia, fuerza, carga nuclear
	    else
	      call SCF_hyb(na_u, nac, rclas, Etot, fdummy, pc(nac),do_SCF,
     .             do_QM_forces, do_properties)
	    end if
	  end if
c	end if

          qmstep(istepconstr) = qmstep(istepconstr) + 1
        endif !qm

C Start MMxQM loop
      do imm=1,mmsteps    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Pasos MMxQM
      step = step +1
      if(mmsteps.ne.1) then
        write(6,*)
        write(6,'(A)')    '*******************************'
        write(6,'(A,i5)') '   MM x QM Step : ', imm 
        write(6,'(A)')    '*******************************'
      endif


C Calculation of last QM-MM interaction: LJ Energy and Forces only 
      if((qm.and.mm)) then
        call ljef(na_u,nac,natot,rclas,Em,Rm,fdummy,Elj,listqmmm)
      endif !qm & mm


C LinkAtom: set again linkmm atoms parameters
      if(qm.and.mm) then
        if(linkatom) then
          do i=1,numlink
            pc(linkmm(i,1:4))=pclinkmm(i,1:4)
            Em(na_u+linkmm(i,1:1))=Emlink(i,1:1)
          enddo
        endif !LA
      endif !qm & mm

C Calculate pure Solvent energy and forces 
      if(mm) then
      call solv_ene_fce(natot,na_u,nac,ng1,rclas,Em,Rm,pc(1:nac),
     .    Etot_amber,fce_amber,attype,
     .    nbond,nangle,ndihe,nimp,multidihe, multiimp,kbond,bondeq,
     .    kangle,angleeq,kdihe,diheeq,kimp,impeq,perdihe,perimp,
     .    bondtype,angletype,dihetype,imptype,
     .    bondxat,angexat,angmxat,dihexat,dihmxat,impxat,atange,
     .    atangm,atdihe,atdihm,atimp,
     .    ng1type,angetype,angmtype,evaldihe,evaldihm,
     .    dihety,dihmty,impty,nonbonded,scale,nonbondedxat,scalexat,
     .    evaldihelog,evaldihmlog,step,nparm,
     .    actualiz,rcortemm,
     .    atname,aaname,sfc,dt,
     .    water,masst,radblommbond)
      endif !mm

C converts fdummy to Kcal/mol/Ang  
      fdummy(1:3,1:natot)=fdummy(1:3,1:natot)*Ang/eV*kcal
 
C add famber to fdummy  
      if(mm) then
      fdummy(1:3,na_u+1:natot)=fdummy(1:3,na_u+1:natot)
     .       +fce_amber(1:3,1:nac)
      endif !mm

C Calculation of LinkAtom Energy and Forces
      if(qm.and.mm ) then
        if(linkatom) then
        call link2(numlink,linkat,linkqm,linkmm,linkmm2,rclas,
     .  natot,na_u,nac,fdummy,ng1,attype,nparm,
     .  nbond,nangle,ndihe,nimp,multidihe,multiimp,kbond,bondeq,
     .  kangle,angleeq,kdihe,diheeq,kimp,impeq,perdihe,perimp,
     .  bondtype,angletype,dihetype,imptype,linkqmtype,
     .  bondxat,Elink,parametro,step,kch,rch)
C Set again link atmos parameters to zero for next step  
        do i=1,numlink
        pclinkmm(i,1:4)=pc(linkmm(i,1:4))
        pc(linkmm(i,1:1))=0.d0
        pc(linkmm(i,2:4))=pc(linkmm(i,2:4))+pclinkmm(i,1)/3.d0
        Em(na_u+linkmm(i,1:1))=0.d0
        enddo
        endif ! LA
      endif !qm & mm

	if (.true.) then
	   do itest=1, natot
	   write(969,423) itest, rclas(1,itest),rclas(2,itest),
     .   rclas(3,itest),fdummy(1,itest)*0.5d0,fdummy(2,itest)*0.5d0,
     .   fdummy(3,itest)*0.5d0
	   end do
	end if



        if(optimization_lvl.eq.1) fdummy=0.d0 !slo mueve por restrain


C Calculation of Constrained Optimization Energy and Forces 
      if(imm.eq.1) then
        if(constropt) then
        call subconstr2(nconstr,typeconstr,kforce,rini,rfin,ro,rt,
     .  nstepconstr,atmsconstr,natot,rclas,fdummy,istp,istepconstr,
     .  ndists,coef)
        endif 
      endif !imm

C Converts fdummy to Ry/Bohr 
      fdummy(1:3,1:natot)=fdummy(1:3,1:natot)/Ang*eV/kcal

C Writes final energy decomposition
       Etots=2.d0*Etot+Elj+((Etot_amber+Elink+Eelec)/kcal*eV)
C estan duplicados los pesos de algunas energias
       Etots=0.5d0*Etots

       write(6,*)
       write(6,'(/,a)') 'hybrid: Energy Decomposition (eV):'
       if(qm) write(6,'(a,2x,F16.6)')           'Esiesta:',Etot/eV      !esta bien comparada con G

       if(qm.and.mm) write(6,'(a,2x,F16.6)')    'Elj:    ',Elj/eV       ! da la mitad si 1 atomo en QM y otro MM q si los 2 son MM
       if(mm) write(6,'(a,2x,F16.6)')      'Esolv:  ',Etot_amber/kcal   !esta esta bien comparada con Amber
       if(Elink.ne.0.0) write(6,'(a,2x,F16.6)') 'Elink:  ',Elink/kcal
       if(Eelec.ne.0.0) write(6,'(a,2x,F16.6)') 'Eefield:',Eelec/kcal
       if(qm.and.mm) write(6,'(a,2x,F16.6)')    'Etots:  ',Etots/eV
       call flush(6)

C Sets fdummy to zero inside mmxqm step
       if(qm.and.mm) then
         if(imm.ne.1) then
           fdummy(1:3,1:na_u) = 0.d0
           vat(1:3,1:na_u) = 0.d0
           if(linkatom) then
             do i=1,numlink
               fdummy(1:3,linkmm(i,1)+na_u)=0.d0
               vat(1:3,linkmm(i,1)+na_u)=0.d0
             enddo
           endif !LA
         endif !imm
       endif !qm & mm

C Impose constraints to atomic movements by changing forces
       call fixed2(na_u,nac,natot,nfree,blocklist,blockqmmm,
     .             fdummy,cfdummy,vat)

C Accumulate coordinates in PDB/CRD file for animation
      call wripdb(na_u,slabel,rclas,natot,step,wricoord,nac,atname,
     .            aaname,aanum,nesp,atsym,isa,listqmmm,blockqmmm)

c freeza QM atom
        if (optimization_lvl .eq. 2)  then !congela QM, luego ponerle una variable
        write(*,*) "cuidado congelo QM"
        do inick=1, na_u
            cfdummy(1:3,inick) = 0.d0
        end do
        end if

c freeza MM atom
c falta un if q active el freezado
        do inick=1, natot
          if(MM_freeze_list(inick)) then
            cfdummy(1:3,inick) = 0.d0
          end if
        end do


c escribo xyz
          write(34,*) natot
          write(34,*)
          do inick=1,natot
              if (inick.le.na_u) then
                  write(34,345) iza(inick), rclas(1:3,inick)/Ang
              else
                  write(34,346) pc(inick-na_u), rclas(1:3,inick)/Ang
              end if
          enddo

C Write atomic forces 
      fmax = 0.0_dp
      cfmax = 0.0_dp
      fres = 0.0_dp
      ftot = 0.0_dp 
      ifmax = 0
      icfmax = 0
      do ix=1,3
        do ia = 1,natot
          ftot(ix) = ftot(ix) + cfdummy(ix,ia) 
          fres = fres + fdummy(ix,ia)*fdummy(ix,ia) 
        enddo
      enddo
      fres = dsqrt( fres / (3.0_dp*natot ) )
      fmax=maxval(dabs(fdummy))
      ifmax=maxloc(dabs(fdummy))
      cfmax=maxval(dabs(cfdummy))
      icfmax=maxloc(dabs(cfdummy))

      write(6,'(/,a)') 'hybrid: Atomic forces (eV/Ang):'
      write(6,'(i6,3f12.6)') (ia,(cfdummy(ix,ia)*Ang/eV,ix=1,3),
     .                        ia=1,nfce)
      write(6,'(43(1h-),/,a4,3f12.6)') 'Tot',(ftot(ix)*Ang/eV,ix=1,3)
      write(6,'(43(1h-),/,a4,f12.6,a,i5)') 'Max',fmax*Ang/eV,
     .                       '  free, atom  ',ifmax(2)
      write(6,'(a4,f12.6,a)')'Res',fres*Ang/eV,
     .                       '  sqrt( Sum f_i^2 / 3N )'
      write(6,'(43(1h-),/,a4,f12.6,a,i5)') 'Max',cfmax*Ang/eV,
     .                       '  cons, atom  ',icfmax(2)
      if(nfce.ne.natot) call iofa(natot,cfdummy)


	write(123,*) "fuerzas antes de CG"
	write(123,*) "paso y nivel de optim", istep, optimization_lvl
        do inick=1, natot
            write(123,*) inick, cfdummy(1:3,inick)
        end do

C Move atoms 
      if (idyn .eq. 0 ) then 
          call cgvc( natot, rclas, cfdummy, ucell, cstress, volume,
     .             dxmax, tp, ftol, strtol, varcel, relaxd, usesavecg )
      endif


c #################################################################################
c centro traj respecto al cm, solo en el caso de sistema solo qm
      if (qm .and. .not. mm) then
	write(*,*) "centrando sistema"
        Rcm=0.d0
        Tmass=0.d0
        do i=1, natot
        do j=1,3
           Rcm(j)=Rcm(j)+rclas(j,i)*masst(i)
         end do
           Tmass=Tmass+masst(i)
        end do
        Rcm=Rcm/Tmass

        do i=1, natot
          do j=1,3
            rclas(j,i)=rclas(j,i)-Rcm(j)
          end do
        end do
      end if
c #################################################################################

C      iunit = 2
c      ntcon = 0
c      if(imm.ne.1) ntcon = na_u

c       if(idyn .ne. 0 .and. idyn .ne. 5) then
c      if(qm) call centerdyn(na_u,rclas,ucell,natot)
c       endif

C Write Energy in file
      call wriene(step,slabel,idyn,Etots,cfmax)

c sets variables for next siesta cycle
      fa = 0.d0
      fdummy = 0.d0 
      cfdummy = 0.d0
      xa(1:3,1:na_u)=rclas(1:3,1:na_u)
      call flush(6)

C Calculation Hlink's New positions 
      if(qm.and.mm) then
        if(linkatom) then
        call link3(numlink,linkat,linkqm,linkmm,rclas,
     .    natot,na_u,nac,distl)
        xa(1:3,1:na_u)=rclas(1:3,1:na_u)
        endif !LA
      endif !qm & mm

C Save last atomic positions and velocities
      call ioxv( 'write', natot, ucell, rclas, vat, foundxv, foundvat )

C write atomic constraints each step
      call wrtcrd(natot,rclas)

C Exit MMxQM loop
      enddo !imm                          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Pasos MMxQM
      if((mmsteps.ne.1).and.(imm.ne.1)) relaxd = .false.

C Exit coordinate relaxation loop
	if (relaxd) then
	  if (opt_scheme.eq.1) then
	  if (optimization_lvl.eq.1) then
	      optimization_lvl=2
c              do_SCF=.true.
c              do_QM_forces=.true.
	  elseif (optimization_lvl.eq.2) then
              optimization_lvl=3
c	      do_SCF=.true.
          else
	      optimization_lvl=1
	      goto 10
	  end if
	  end if
c	else
c	  if (optimization_lvl.eq.2) do_SCF=.false.
	end if
!      if (relaxd) goto 10 

       enddo                              !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Paso opt
 10   continue
C End of coordinate-relaxation loop ==================================

C End of constrain optimization loop

	write(*,956) rt(1), Etot/eV, Elj/eV, Etot_amber/kcal, 
     .  Elink/kcal, Eelec/kcal, Etots/eV

      if(constropt) then
       call subconstr3(ro(1),rt(1),dr,Etots)
c aca escribe el rce
       call wrirtc(slabel,Etots,rt(1),istepconstr,na_u,nac,natot,
     .             rclas,atname,aaname,aanum,nesp,atsym,isa)
       call ioxvconstr( natot, ucell, rclas, vat, istepconstr )
      endif


C calculo de propiedades en lio
      do_properties=.true.
      call SCF_hyb(na_u, at_MM_cut_QMMM, r_cut_QMMM, Etot,
     .     F_cut_QMMM,
     .     Iz_cut_QMMM, do_SCF, do_QM_forces, do_properties) !fuerzas lio, Nick
      do_properties=.false.

      enddo !istepconstr                  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< RESTRAIN



	close(34)

      if(qm) call lio_finalize()
C Dump last coordinates to output
      if (writeipl) then
      if(qm) then
      write(6,'(/a)')'hybrid: Last atomic coordinates (Ang) and species'
      write(6,"(i6,2x,3f10.5,2x,i3)")
     .         (ia, (xa(ix,ia)/Ang, ix=1,3), isa(ia), ia=1,na_u)
      endif !qm
C Dump last Solvent coordinates to output
      if(mm) then
      write(6,'(/a)')'hybrid: Last solvent coordinates (Ang)'
      write(6,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     . ('ATOM',i,atname(i),aaname(i),'A',aanum(i),
     . (rclas(ia,na_u+i)/Ang,ia=1,3), i=1,nac)
      endif !mm
      endif !writeipl
      call flush(6)

C Print final date and time
      call timestamp('End of run')
 345  format(2x, I2,    2x, 3(f10.6,2x))
 346  format(2x, f10.6, 2x, 3(f10.6,2x))
 956  format(2x, "Econtribution", 7(f18.6,2x))
 423  format(2x, I6,2x, 6(f20.10,2x))
      end program


