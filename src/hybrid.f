      Program HYBRID 
!****************************************************************************
! Program for QM-MM calculations. 
!
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
! Modified version for use Lio as QM-MM level and debugs, N. Foglia 2017
! Icluded NEB method N. Foglia 2018
!*****************************************************************************

! Modules
      use precision, only: dp
      use sys, only: die
      use fdf
      use ionew, only: io_setup, IOnode    
      use scarlett, only: istep, nmove, nesp, inicoor,fincoor, idyn, 
     . natot,na_u,qm, mm, atsym, nparm, 
     . xa, fa,
     . masst,isa, iza, pc, 
     . nac, atname, aaname, attype, qmattype, aanum, ng1, bondtype,
     . kbond,bondeq, bondxat, angletype, kangle,angleeq, angexat,
     . angmxat,dihetype, kdihe,diheeq, perdihe, multidihe, dihexat, 
     . dihmxat, imptype, kimp,impeq,perimp, multiimp, impxat, 
     . scalexat, scale, nonbondedxat, nonbonded, Em, Rm,
     . fce_amber, fdummy, cfdummy, ng1type, angetype, angmtype,
     . dihety,dihmty,impty, evaldihelog, evaldihmlog,
     . atange, atangm,
     . atdihe,atdihm,atimp,
     . rclas, vat, aat, izs, evaldihe,evaldihm, 
     . linkatom, numlink, linkat, linkqm, linkmm, linkmm2, parametro,
     . linkqmtype, Elink, distl, pclinkmm, Emlink, frstme, pi,
!cutoff
     . r_cut_list_QMMM,blocklist,blockqmmm, blockall,
     . listqmmm,MM_freeze_list, natoms_partial_freeze, 
     . natoms_partial_freeze, coord_freeze, 
!NEB
     . NEB_Nimages, 
     . NEB_firstimage, NEB_lastimage,  
     . aclas_BAND_old,
     . rclas_BAND,
     . vclas_BAND, fclas_BAND, Energy_band,
     . NEB_distl,
     . ucell,
     . ftol,
     . Ang, eV, kcal, 
!Dynamics
     . Ekinion, tempion, tempinit, tt, tauber, tempqm, kn, vn, mn,
!FIRE
     . time_steep, Ndescend, time_steep_max, alpha,
!Lio
     . charge, spin,
!outputs
     . writeRF, slabel, traj_frec

      implicit none
! General Variables
      integer :: replica_number !auxiliar
      integer :: istp !number of move step for each restrain starting on 1
      integer :: step !total number of steps in a MMxQM movement (not tested in this version)
      integer :: ifmax(2) !number of atom with max force on each cicle
      real(dp) :: fmax !max force on each cicle 
      integer :: icfmax(2) !number of atom with max force on each cicle after constrain
      real(dp) :: cfmax !max force on each cicle after constrain
      real(dp) :: cstress(3,3) !Stress tensor, may be included in CG minimizations 
      real(dp) :: strtol !Maximum stress tolerance
      real(dp) :: Etot ! QM + QM-MM interaction energy
      real(dp) :: fres !sqrt( Sum f_i^2 / 3N )
      real(dp) :: ftot(3) ! Total force in system
      real(dp) :: dxmax !Maximum atomic displacement in one CG step (Bohr)
      real(dp) :: tp !Target pressure
      real(dp) :: volume !Cell volume
      logical :: usesavexv, foundxv !control for coordinates restart
      logical :: usesavecg !control for restart CG
      logical :: varcel !true if variable cell optimization
      logical :: relaxd ! True when CG converged
      character :: paste*25
      external :: paste
      logical :: actualiz!MM interaction list control
!!!!
      integer :: nfree !number of atoms without a contrain of movement
      integer :: mmsteps !number of MM steps for each QM step, not tested
      integer :: nroaa !number of residues
      integer :: nbond, nangle, ndihe, nimp !number of bonds, angles, dihedrals and impropers defined in amber.parm
      integer, dimension(:), allocatable :: atxres !number ot atoms in residue i, no deberia estar fija la dimension
      double precision :: Etot_amber !total MM energy
      double precision :: Elj !LJ interaction (only QMMM)
      double precision :: Etots !QM+QMMM+MM energy
      integer :: ntcon !number of restrained degrees of freedom JOTA 
      integer :: cmcf !Center of Mass Coordinates Fixed JOTA

! ConstrOpt variables
      logical :: constropt !activate restrain optimizaion
      integer :: nconstr !number of constrains
      integer :: nstepconstr !numero de pasos en los que barrera la coordenada de reaccion (limitado a 1-100 hay q cambiar esto, Nick)
      integer, dimension(:), allocatable :: typeconstr !type of cosntrain (1 to 8)
      double precision, dimension(:), allocatable :: kforce !force constant of constrain i
      double precision :: rini,rfin  !initial and end value of reaction coordinate
      double precision, dimension(:), allocatable :: ro ! fixed value of constrain in case nconstr > 1 for contrains 2+
      double precision, dimension(:), allocatable :: rt ! value of reaction coordinate in constrain i
      double precision :: dr !change of reaction coordinate in each optimization dr=(rfin-rini)/nstepconstr
      double precision, dimension(:,:), allocatable :: coef !coeficients for typeconstr=8
      integer, allocatable :: atmsconstr(:,:), ndists(:) !atomos incluidos en la coordenada de reaccion
      integer :: istepconstr !auxiliar

! Cut Off QM-MM variables
      double precision :: rcorteqm ! not used
      double precision :: rcorteqmmm ! distance for QM-MM interaction
      double precision :: rcortemm ! distance for LJ & Coulomb MM interaction
      double precision :: radbloqmmm ! distance that allow to move MM atoms from QM sub-system
      double precision :: radblommbond !parche para omitir bonds en extremos terminales, no se computan bonds con distancias mayores a radblommbond
      double precision :: radinnerbloqmmm !distance that not allow to move MM atoms from QM sub-system
      integer :: res_ref ! residue that is taken as reference to fix atoms by radial criteria in full MM simulations JOTA  
      logical ::  recompute_cuts
! Lio
      logical :: do_SCF, do_QM_forces !control for make new calculation of rho, forces in actual step
      logical :: do_properties !control for lio properties calculation

! Optimization scheme
      integer :: opt_scheme ! turn on optimization scheme
      integer :: optimization_lvl ! level of movement in optimization scheme:
! 1- only QM atoms with restrain
! 2- only MM atoms
! 3- all

!variables para centrado del sistema, mejorar esto agregando el tensor de inercia, Nick
      logical :: Nick_cent !activa centrado y conservacion de los ejes de inercia
      integer :: firstcent ! marca el 1er paso
      double precision, dimension(9) :: Inivec !inertia tensor autovectors


!variables para cuts
      double precision, allocatable, dimension(:,:) :: r_cut_QMMM,
     .  F_cut_QMMM
      double precision, allocatable, dimension(:) :: Iz_cut_QMMM
      integer :: at_MM_cut_QMMM, r_cut_pos
      double precision :: r12 !auxiliar
      integer :: i_qm, i_mm ! auxiliars
      logical :: done, done_freeze, done_QMMM !control variables

! restarts
      logical :: foundcrd
      logical :: foundvat !found velocities in restart

! Outputs variables
      character :: xyzlabel*20 !*.xyz name
      integer  :: nfce !number of atoms for whom the forces will be writen, wrifces=0,1,2 => nfce = 0, na, nat
      integer ::  wricoord !number of steps for write coordinates
      logical :: writeipl


! Auxiliars
      integer :: i, ia, imm, iunit, ix, j, k, inick, jnick, itest

! Others that need check
!!!! General Variables
      real(dp) :: dt
      real(dp) :: volcel
      external :: volcel

      external
     .  chkdim, cgvc,  ioxv, fixed1, assign,
     .  iofa, ofc, reinit, read_qm,
     .  read_md, fixed2

!!!! Solvent General variables
      double precision  :: sfc
      logical :: water

! Solvent external variables
       external
     . solv_assign, solv_ene_fce, qmmm_lst_blk, wrtcrd,
     . centermol, centerdyn, link1, link2, link3, ljef,
     . mmForce, readcrd,  prversion, ioxvconstr,  
     . wripdb, wriene, wrirtc, subconstr1, subconstr2, subconstr3


!--------------------------------------------------------------------
!need to move this to an initializacion subroutine
      allocate(atxres(20000))
      allocate(typeconstr(20), kforce(20), ro(20), rt(20), coef(20,10))
      allocate(atmsconstr(20,20), ndists(20))

! Initialise some variables
      relaxd=.false.
      varcel=.false.
      tp=0.0
      cstress=0.0
      strtol=0.0
      qm=.true.
      mm=.true.
      Etot_amber=0.d0
      charge=0
      spin=0.d0
      do_SCF=.true.
      do_QM_forces=.true.
      rcorteqmmm=0.d0
      radbloqmmm=0.d0
      do_properties=.false.
      Nick_cent=.false.
      foundxv=.false.
      recompute_cuts=.true.
      cmcf = 0
! Initialize IOnode
      call io_setup   

! Print version information
      call prversion

! Start time counter 
      call timestamp('Start of run')      

! Factors of conversion to internal units 
      call init_hybrid('Constants')

! Initialise read 
      call reinit(slabel) !, sname)

! Read and initialize basics variables 
      call init_hybrid('Jolie')

! Some definitions 
      ucell=0.d0
      fa=0.d0
      fdummy = 0.d0
      cfdummy=0.d0
      vat=0.d0
      Elj=0.d0
      Elink=0.d0
      Etot=0.d0
      Etots=0.d0
      nfree=natot
      blocklist = 0
      blockqmmm = 0
      blockall = 0 !JOTA
      listqmmm = 0
      step = 0
      nstepconstr = 0
      numlink = fdf_integer('LinkAtoms',0)
      linkatom = .false.

      if(numlink.ne.0) linkatom = .true.

      opt_scheme=0
      opt_scheme = fdf_integer('OptimizationScheme',0)
      write(*,*) "opt_scheme vale", opt_scheme

      constropt = .false.
      constropt = fdf_block('ConstrainedOpt',iunit)
      foundvat = .false.
      writeipl = fdf_boolean('WriIniParLas',.false.)

! Read and assign Solvent variables 
      if(mm) then
        call solv_assign(na_u,natot,nac,nroaa,Em,Rm,attype,pc,
     .  ng1,bondxat,angexat,atange,angmxat,atangm,dihexat,atdihe,
     .  dihmxat,atdihm,impxat,atimp,
     .  nbond,kbond,bondeq,bondtype,
     .  nangle,kangle,angleeq,angletype,
     .  ndihe,kdihe,diheeq,dihetype,multidihe,perdihe,
     .  nimp,kimp,impeq,imptype,multiimp,perimp,
     .  nparm,aaname,atname,aanum,qmattype,rclas,
     .  rcorteqmmm,rcorteqm,rcortemm,sfc,
     .  radbloqmmm,atxres,radblommbond,radinnerbloqmmm,res_ref)
      endif !mm

! changing cutoff to atomic units
      rcorteqmmm=rcorteqmmm*Ang
      rcorteqmmm=rcorteqmmm**2 !we will compare square in cutoff

      rclas(1:3,1:na_u) = xa(1:3,1:na_u)

! Read simulation data 
      call read_md( idyn, nmove, dt, dxmax, ftol, 
     .              usesavecg, usesavexv , Nick_cent, na_u,  
     .              natot, nfce, wricoord, mmsteps, tempinit,
     .              tt, tauber,mn)

! Assignation of masses and species 
      call assign(na_u,nac,atname,iza,izs,masst)
      if (idyn .ne. 1) then ! Read cell shape and atomic positions from a former run
        call ioxv('read',natot,ucell,rclas,vat,foundxv,foundvat,'X',-1) 
        if (foundxv) xa(1:3,1:na_u)=rclas(1:3,1:na_u)
      else 
	call init_hybrid("NEB")
	call NEB_make_initial_band(usesavecg)
      end if


! Reading LinkAtom variables
      if(qm.and.mm) then
        if (linkatom) then
          call link1(numlink,linkat,linkqm,linkmm,linkqmtype,
     .               linkmm2,ng1,nac,na_u,qmattype,rclas)
! sets to zero Em for HL and CQM
          do i=1,numlink
            Em(linkat(i))=0.d0
            Em(linkqm(i,1:1))=0.d0
          enddo
        endif !LA
      endif !qm & mm


! Read a crd file from a former run
      call readcrd(na_u,nac,masst,linkatom,linkat,numlink,
     .             rclas,vat,foundcrd,foundvat)
      if(foundcrd) xa(1:3,1:na_u)=rclas(1:3,1:na_u)

     
! Sets LinkAtoms' positions
      if(qm.and.mm) then
        if (linkatom) then
	  if (idyn .ne. 1) then
            call link3(numlink,linkat,linkqm,linkmm,rclas,
     .               natot,na_u,nac,distl)
            xa(1:3,1:na_u)=rclas(1:3,1:na_u)
	  else !NEB case
	    NEB_distl=1.09
	    do replica_number = NEB_firstimage, NEB_lastimage 
	      rclas(1:3,1:natot)=rclas_BAND(1:3,1:natot,replica_number)
	      distl(1:15)=NEB_distl(1:15,replica_number)
	      frstme=.true.
	      call link3(numlink,linkat,linkqm,linkmm,rclas,
     .               natot,na_u,nac,distl)
	      rclas_BAND(1:3,1:na_u,replica_number)=rclas(1:3,1:na_u)
	      NEB_distl(1:15,replica_number)=distl(1:15)
	    end do
	   end if

	  ! Sets to zero pc and Em for MMLink 
            do i=1,numlink
              pclinkmm(i,1:4)=pc(linkmm(i,1:4))
              pc(linkmm(i,1:1))=0.d0
              pc(linkmm(i,2:4))=pc(linkmm(i,2:4))+pclinkmm(i,1)/3.d0
! two options: Em(linkmm(i,1:1))=0.0 or Em(linkmm(i,1:4))=0.0
              Emlink(i,1:1)=Em(na_u+linkmm(i,1:1))
              Em(na_u+linkmm(i,1:1))=0.d0
            enddo
        endif !LA
      endif !qm & mm


! Dump initial coordinates to output 
      if (writeipl) then
        if(qm) then
          write(6,'(/a)') 'hybrid: Atomic coordinates (Ang) and species'
          write(6,"(i6,2x,3f10.5,2x,i3)")
     .         (ia, (xa(ix,ia)/Ang, ix=1,3), isa(ia), ia=1,na_u)
        endif !qm

! Dump initial Solvent coordinates to output
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
     .         (i,qmattype(i),Em(i),Rm(i),i=1,na_u)
        endif !qm & mm

        if(mm) then
          write(6,"(i6,2x,A4,2x,f9.6,2x,f9.6,2x,f9.6)")
     .     (i,attype(i),Em(i+na_u),Rm(i+na_u),pc(i),i=1,nac)
        endif !mm
      endif !writeipl
      call flush(6)


! Initialize .xyz
      xyzlabel = paste( slabel, '.xyz' )
      open(unit=34, file=xyzlabel)

! Initialize Lio 
      if(qm) then
        Etot=0.d0
        fa=0.d0
        call init_lio_hybrid(na_u, nac, charge, iza, spin)

! calculate cell volume
        volume = volcel( ucell )

! Center system 
        if (.not.foundxv)call centermol(na_u,xa,rclas,ucell,natot)
      endif !qm

C Calculate Rcut & block list QM-MM 

      if(qm.and.mm) then
        call qmmm_lst_blk(na_u,nac,natot,nroaa,atxres,rclas,
     .  rcorteqmmm,radbloqmmm,blockqmmm,listqmmm,rcorteqm,slabel,
     .  radinnerbloqmmm,blockall)
      endif !qm & mm
        

! Calculate blockall for full-MM simulation JOTA
      if(mm.and. .not. qm) then

        call fixed0(res_ref,natot,nroaa,atxres,rclas,blockall,
     .              radbloqmmm,radinnerbloqmmm)

      endif
! nrjota hardcodeado = 1

! Read fixed atom constraints
      call fixed1(na_u,nac,natot,nroaa,rclas,blocklist,
     .            atname,aaname,aanum,water)

! Count fixed degrees of freedom
      call fixed3(natot,blockall,ntcon)

! Build initial velocities according to Maxwell-Bolzmann distribution

        if ((idyn .gt. 3) .and. (.not. foundvat))
     .  call vmb(natot,tempinit,masst,rclas,0,vat,cmcf,blockall,ntcon)
!tempinit




!########################################################################################
!#####################################  MAIN LOOPS  #####################################
!########################################################################################


! Start loop over constrained optimization steps
      if(constropt) then
        call subconstr1(nconstr,typeconstr,kforce,nstepconstr,
     .        rini,rfin,atmsconstr,dr,ro,ndists,coef,constropt)
      endif

      do istepconstr=1,nstepconstr+1   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< RESTRAIN Loop
! istepconstr marca la posicion del restrain
        optimization_lvl=3
        if (opt_scheme .eq. 1) then
          optimization_lvl=1
        end if

        if(constropt) then
          write(6,*)
          write(6,'(A)')    '*******************************'
          write(6,'(A,i5)') '  Constrained Step : ', istepconstr
          write(6,'(A)')    '*******************************'
        endif

! Begin of coordinate relaxation iteration ============================
        if (idyn .lt. 7 ) then ! case 0 1 2 3
          inicoor = 0
          fincoor = nmove
        endif
        at_MM_cut_QMMM=nac


! Start loop over coordinate changes
        istp = 0
        do istep = inicoor,fincoor    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CG optimization cicle
          istp = istp + 1
 
          write(6,'(/2a)') 'hybrid:                 ',
     .                    '=============================='

          if (idyn .ge. 7) 
     .    STOP 'only CG, QM, FIRE or NEB minimization available'

          write(6,'(28(" "),a,i6)') 'Begin move = ',istep
          write(6,'(2a)') '                        ',
     .                    '=============================='
          write(6,*) "Optimization level: ", optimization_lvl

!start loot over NEB images
	do replica_number = NEB_firstimage, NEB_lastimage      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Band Replicas
	  if (idyn .eq.1) then
	    rclas(1:3,1:natot)=rclas_BAND(1:3,1:natot,replica_number)
	  end if


! Calculate Energy and Forces using Lio as Subroutine
          if(qm) then
            if (allocated(r_cut_QMMM)) deallocate(r_cut_QMMM)
            if (allocated(F_cut_QMMM)) deallocate(F_cut_QMMM)
            if (allocated(Iz_cut_QMMM)) deallocate(Iz_cut_QMMM)
 
            call compute_cutsqmmm(at_MM_cut_QMMM,istepconstr,radbloqmmm,
     .      rcorteqmmm,nroaa,atxres)

            allocate (r_cut_QMMM(3,at_MM_cut_QMMM+na_u), 
     .      F_cut_QMMM(3,at_MM_cut_QMMM+na_u), 
     .      Iz_cut_QMMM(at_MM_cut_QMMM+na_u))
            r_cut_QMMM=0.d0
            F_cut_QMMM=0.d0
            Iz_cut_QMMM=0


c	  recompute_cuts=.true.
c	  if (istep.eq.inicoor) recompute_cuts=.true.
c	  if (replica_number.gt.1) recompute_cuts=.false.

c	  if (recompute_cuts) then ! define lista de interacciones en el primer paso de cada valor del restrain
c	    if (allocated(r_cut_QMMM)) deallocate(r_cut_QMMM)
c	    if (allocated(F_cut_QMMM)) deallocate(F_cut_QMMM)
c	    if (allocated(Iz_cut_QMMM)) deallocate(Iz_cut_QMMM)
c	    r_cut_list_QMMM=0
c	    r_cut_pos=0
c	    at_MM_cut_QMMM=0

c	    if (istepconstr.eq.1) then
c		MM_freeze_list=.true.
c		do i_qm=1,na_u
c		  MM_freeze_list(i_qm)=.false.
c		end do
c	    end if

c	    do i_mm=1, nac !MM atoms
c	      i_qm=0
c	      done=.false.
c              done_freeze=.false.
c              done_QMMM=.false.
c	      do while (i_qm .lt. na_u .and. .not. done) !QM atoms
c	        i_qm=i_qm+1
c                r12=(rclas(1,i_qm)-rclas(1,i_mm+na_u))**2.d0 +
c     .              (rclas(2,i_qm)-rclas(2,i_mm+na_u))**2.d0 +
c     .              (rclas(3,i_qm)-rclas(3,i_mm+na_u))**2.d0
c
c	        if(r12 .lt. rcorteqmmm .and. .not. done_QMMM) then
c	          done_QMMM=.true.
c	          at_MM_cut_QMMM=at_MM_cut_QMMM+1
c	          r_cut_pos=r_cut_pos+1
c	          r_cut_list_QMMM(i_mm)=r_cut_pos
c	        end if



c		if (istepconstr.eq.1) then !define lista de movimiento para la 1er foto
c	          if(r12 .lt. radbloqmmm .and. .not. done_freeze) then
c	             MM_freeze_list(i_mm+na_u)=.false.
c	             done_freeze=.true.
c	          end if
c		end if
c


c		done=done_QMMM .and. done_freeze
c	      end do
c	    end do


c	  allocate (r_cut_QMMM(3,at_MM_cut_QMMM+na_u),
c     .    F_cut_QMMM(3,at_MM_cut_QMMM+na_u),
c     .    Iz_cut_QMMM(at_MM_cut_QMMM+na_u))

c	  r_cut_QMMM=0.d0
c	  F_cut_QMMM=0.d0
c	  Iz_cut_QMMM=0
c	  end if


!copy position and nuclear charges to cut-off arrays
          do i=1,natot !barre todos los atomos
            if (i.le.na_u) then !QM atoms
              r_cut_QMMM(1:3,i)= rclas(1:3,i)
            else if (r_cut_list_QMMM(i-na_u) .ne. 0) then !MM atoms inside cutoff
	      r_cut_QMMM(1:3,r_cut_list_QMMM(i-na_u)+na_u) = rclas(1:3,i)
              Iz_cut_QMMM(r_cut_list_QMMM(i-na_u))= pc(i-na_u)
            end if
          end do
	  call SCF_hyb(na_u, at_MM_cut_QMMM, r_cut_QMMM, Etot, 
     .     F_cut_QMMM,
     .     Iz_cut_QMMM, do_SCF, do_QM_forces, do_properties) !fuerzas lio, Nick


c return forces to fullatom arrays
          do i=1, natot
	    if (i.le.na_u) then !QM atoms
	      fdummy(1:3,i)=F_cut_QMMM(1:3,i)
	    else if (r_cut_list_QMMM(i-na_u).ne.0) then !MM atoms in cut-off
              fdummy(1:3,i)=
     .        F_cut_QMMM(1:3,r_cut_list_QMMM(i-na_u)+na_u)
	    end if
          end do
        endif !qm    termino el if(qm)
! here Etot in Hartree, fdummy in Hartree/bohr

! Start MMxQM loop
          do imm=1,mmsteps    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MMxQM Steps
            step = step +1
            if(mmsteps.ne.1) then
              write(6,*)
              write(6,'(A)')    '*******************************'
              write(6,'(A,i5)') '   MM x QM Step : ', imm 
              write(6,'(A)')    '*******************************'
            endif
	
! Calculation of last QM-MM interaction: LJ Energy and Forces only 
            if((qm.and.mm)) then
              call ljef(na_u,nac,natot,rclas,Em,Rm,fdummy,Elj,listqmmm)
            endif !qm & mm

! LinkAtom: set again linkmm atoms parameters
            if(qm.and.mm) then
              if(linkatom) then
                do i=1,numlink
                  pc(linkmm(i,1:4))=pclinkmm(i,1:4)
                  Em(na_u+linkmm(i,1:1))=Emlink(i,1:1)
                enddo
              endif !LA
            endif !qm & mm

! Calculate pure Solvent energy and forces
! Forces in Hartree/bohr JOTA
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


! fce_amber in kcal/(molAng) 

! converts fdummy to Kcal/mol/Ang          
            fdummy(1:3,1:natot)=fdummy(1:3,1:natot)*Ang*kcal/eV

! here Etot in Hartree, fdummy in kcal/mol Ang

! add famber to fdummy  
            if(mm) then
              fdummy(1:3,na_u+1:natot)=fdummy(1:3,na_u+1:natot)
     .        +fce_amber(1:3,1:nac)
            endif !mm

! here fdummy in kcal/mol/Ang

! Calculation of LinkAtom Energy and Forces
            if(qm.and.mm ) then
              if(linkatom) then
        call link2(numlink,linkat,linkqm,linkmm,linkmm2,rclas,
     .  natot,na_u,nac,fdummy,ng1,attype,nparm,
     .  nbond,nangle,ndihe,nimp,multidihe,multiimp,kbond,bondeq,
     .  kangle,angleeq,kdihe,diheeq,kimp,impeq,perdihe,perimp,
     .  bondtype,angletype,dihetype,imptype,linkqmtype,
     .  bondxat,Elink,parametro,step)
! Set again link atmos parameters to zero for next step  
                do i=1,numlink
        pclinkmm(i,1:4)=pc(linkmm(i,1:4))
        pc(linkmm(i,1:1))=0.d0
        pc(linkmm(i,2:4))=pc(linkmm(i,2:4))+pclinkmm(i,1)/3.d0
        Em(na_u+linkmm(i,1:1))=0.d0
                enddo
              endif ! LA
            endif !qm & mm

        if(optimization_lvl.eq.1) fdummy=0.d0 !only move atoms with restrain

! Calculation of Constrained Optimization Energy and Forces 
        if(imm.eq.1) then
          if(constropt) then
        call subconstr2(nconstr,typeconstr,kforce,rini,rfin,ro,rt,
     .  nstepconstr,atmsconstr,natot,rclas,fdummy,istp,istepconstr,
     .  ndists,coef)
	  if(idyn .ge. 4) call subconstr4(istep,rt(1),slabel)
          endif 
        endif !imm


! Converts fdummy Hartree/bohr
        fdummy(1:3,1:natot)=fdummy(1:3,1:natot)*eV/(Ang*kcal)

! here Etot in Hartree, fdummy in Hartree/bohr

! Writes final energy decomposition

	Etots=Etot+1.d0*(Elj+(Etot_amber+Elink)*eV/kcal) !Elj in Hartree, Etot_amber and Elink in kcal/mol
       
! here Etot in Hartree
	write(6,*)
	write(6,'(/,a)') 'hybrid: Energy Decomposition (eV):'
	if(qm) write(6,999)           'Elio :',Etot/eV      
	if(qm.and.mm) write(6,999)    'Elj:    ',Elj/eV       
	if(mm) write(6,999)      'Esolv:  ',Etot_amber/kcal   
	if(Elink.ne.0.0) write(6,999) 'Elink:  ',Elink/kcal
	write(6,999)    'Etots:  ',Etots/eV
       call flush(6)
! saque if para Etots

! Sets fdummy to zero inside mmxqm step
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

! here Etot in Hartree, fdummy in Hartree/bohr

! Impose constraints to atomic movements by changing forces
       call fixed2(na_u,nac,natot,nfree,blocklist,blockqmmm,
     .             fdummy,cfdummy,vat,optimization_lvl,blockall)
! from here cfdummy is the reelevant forces for move system
! here Etot in Hartree, cfdummy in Hartree/bohr


! Accumulate coordinates in PDB/CRD file for animation
      call wripdb(na_u,slabel,rclas,natot,step,wricoord,nac,atname,
     .            aaname,aanum,nesp,atsym,isa,listqmmm,blockqmmm)

! freeze QM atom   Jota, meter todo esto en fixed 2
c        if (optimization_lvl .eq. 2)  then
c          do inick=1, na_u
c            cfdummy(1:3,inick) = 0.d0
c          end do
c        end if

! freeze MM atom
c	if (qm) then
c        do inick=1, natot
c          if(MM_freeze_list(inick)) then
c            cfdummy(1:3,inick) = 0.d0
c          end if
c        end do
c	endif !jota

! partial freeze
c	do inick=1,natoms_partial_freeze
c	  do jnick=1,3
c	    if (coord_freeze(inick,1+jnick) .eq. 1) then
c	      cfdummy(jnick,coord_freeze(inick,1))=0.d0
c	      vat(jnick,coord_freeze(inick,1))=0.d0 !JOTA
c	    end if
c	  end do
c	enddo   JOTA lo pusimos en fixed2


! write xyz, hay q ponerle un if para escribir solo cuando se necesita
c	call write_xyz(natot, na_u, iza, pc, rclas)

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

! here Etot in Hartree, cfdummy in Hartree/bohr

      if (mn .eq. 0 .and. idyn .eq. 6) then
        mn=dble(3*natot-ntcon-cmcf)*tt*8.617d-5*(50.d0*dt)**2
        write(6,'(/,a)') 'Calculating Nose mass as Ndf*Tt*KB*(50dt)**2'
        write(6,999) "mn =", mn
      endif

      if (idyn .ne. 1 ) then !Move atoms with a CG algorithm

        if (writeRF .eq. 1) then!save coordinates and forces for integration 
           do itest=1, natot
	      write(969,423) itest, rclas(1:3,itest)*Ang,
c     .        cfdummy(1:3,itest)*kcal/(eV *Ang)  ! Ang, kcal/ang mol
     .        cfdummy(1:3,itest)*kcal*Ang/eV  ! Ang, kcal/ang mol JOTA saco *0.5 
           end do
        end if

	Ekinion=0.d0
 
	if (idyn .eq. 0) then !Conjugated Gradient
	  call cgvc( natot, rclas, cfdummy, ucell, cstress, volume,
     .             dxmax, tp, ftol, strtol, varcel, relaxd, usesavecg )
	elseif (idyn .eq. 2) then !Quick Minimization
	  call check_convergence(relaxd, cfdummy)
	  if (.not. relaxd) call quick_min(natot, rclas, cfdummy, aat,
     .    vat, masst)
	elseif (idyn .eq. 3) then !FIRE
	  call check_convergence(relaxd, cfdummy)
	  if (.not. relaxd) call FIRE(natot, rclas,cfdummy, aat, vat, 
     .    masst, time_steep,Ndescend, time_steep_max, alpha)
	elseif (idyn .eq. 4) then
	  call verlet2(istp, 3, 0, natot, cfdummy, dt,
     .        masst, ntcon, vat, rclas, Ekinion, tempion, nfree, cmcf)
!iquench lo dejamos como 0, luego cambiar
        elseif (idyn .eq. 5) then
          call berendsen(istp,3,natot,cfdummy,dt,tauber,masst,
     .        ntcon,vat,rclas,Ekinion,tempion,tt,nfree,cmcf)
        elseif (idyn .eq. 6) then
          call nose(istp,natot,cfdummy,tt,dt,masst,mn,ntcon,vat,rclas,
     .        Ekinion,kn,vn,tempion,nfree,cmcf)

!tauber, tt
!iunit fijado en 3

	else
	  STOP "Wrong idyn value"
	end if

	write(6,999)
     .  'hybrid: Temperature Antes:', tempion, ' K' 
       if(idyn .gt. 3) then
         call calculateTemp(Ekinion,tempion,tempqm,vat,ntcon,
     . nfree,cmcf)

        write(6,999)
     .  'hybrid: Kinetic Energy (eV):',Ekinion/eV
        write(6,999)
     .  'hybrid: Total Energy + Kinetic (eV):',(Etots+Ekinion)/eV
        write(6,999)
     .  'hybrid: System Temperature:', tempion, ' K'
        if(qm) write(6,999)
     .  'hybrid: QM SubSystem Temperature:', tempqm, ' K'
        if(idyn .eq. 6) write(6,999)
     .   'Kinetic energy of Nose variable:', kn, ' eV'
        if(idyn .eq. 6) write(6,999)
     .   'Potential energyy of Nose var:', vn, 'eV'

!      if(qm) call centerdyn(na_u,rclas,ucell,natot)
	if (MOD((istp - inicoor),traj_frec) .eq. 0)
     .  call wrirtc(slabel,Etots,dble(istp),istp,na_u,nac,natot,
     .      rclas,atname,aaname,aanum,nesp,atsym,isa)
       endif



!Nick center
        if (qm .and. .not. mm .and. Nick_cent) then
          write(*,*) "Centrando Nick"
          firstcent=0
          if (istepconstr.eq.1 .and. istep.eq.inicoor ) firstcent=1
            call center_rotation(natot, masst, rclas, firstcent, Inivec)
        end if

!Write Energy in file
        call wriene(step,slabel,idyn,Etots,cfmax)

! sets variables for next cycle
        fa = 0.d0
        fdummy = 0.d0
        cfdummy = 0.d0
        xa(1:3,1:na_u)=rclas(1:3,1:na_u) !Bohr
        call flush(6)

! Calculation Hlink's New positions 
        if(qm.and.mm) then
          if(linkatom) then
            call link3(numlink,linkat,linkqm,linkmm,rclas,
     .      natot,na_u,nac,distl)
            xa(1:3,1:na_u)=rclas(1:3,1:na_u)
          endif !LA
        endif !qm & mm

! Save last atomic positions and velocities
        call ioxv( 'write',natot,ucell,rclas,vat,foundxv,foundvat,'',-1)
! write atomic constraints each step
        call wrtcrd(natot,rclas)


      elseif (idyn.eq.1) then !Save forces and energy for a NEB optimization
        fclas_BAND(1:3,1:natot,replica_number)=cfdummy(1:3,1:natot) !Hartree/bohr
	fa = 0.d0
        fdummy = 0.d0
        cfdummy = 0.d0
        Energy_band(replica_number)=Etots ! eV
      endif

! Exit MMxQM loop
      enddo !imm                          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Pasos MMxQM

      end do!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Band Replicas



      if (idyn .eq. 1 ) then !Move atoms in a NEB scheme

	if (writeRF .eq. 1) then!save coordinates and forces for integration 

	  open(unit=969, file="Pos_forces.dat")
	  do replica_number = NEB_firstimage, NEB_lastimage ! Band Replicas
	    do itest=1, natot
	      write(969,423) itest,
     .        rclas_BAND(1:3,itest,replica_number)/Ang,
     .        fclas_BAND(1:3,itest,replica_number)*Ang/eV*kcal ! Ang, kcal/ang mol
	    end do
	    write(969,423)
	  end do
	  close(969)
	end if


	call NEB_save_traj_energy()
	call NEB_steep(istep, relaxd) 

! Calculation Hlink's New positions 
        if(qm.and.mm) then
          if(linkatom) then
	    do replica_number = NEB_firstimage, NEB_lastimage
              rclas(1:3,1:natot)=rclas_BAND(1:3,1:natot,replica_number)
	      distl(1:15)=NEB_distl(1:15,replica_number)
	      call link3(numlink,linkat,linkqm,linkmm,rclas,
     .               natot,na_u,nac,distl)
	      rclas_BAND(1:3,1:na_u,replica_number)=rclas(1:3,1:na_u)
	      NEB_distl(1:15,replica_number)=distl(1:15)
	    end do
	  endif !LA
	endif !qm & mm
      end if


      if((mmsteps.ne.1).and.(imm.ne.1)) relaxd = .false.

! Exit coordinate relaxation loop
	if (relaxd) then

	  if (opt_scheme.eq.1) then
	    if (optimization_lvl.eq.1) then
	      optimization_lvl=2
	    elseif (optimization_lvl.eq.2) then
              optimization_lvl=3
            else
	      optimization_lvl=1
	      goto 10
	    end if
	  else
	    goto 10
	  end if
	end if
      enddo                              !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CG optimization cicle
 10   continue


!Write optimiced structure(s) and energy(es)
	if (idyn.eq.1) then
	   do replica_number = 1, NEB_Nimages
		Etots=Energy_band(replica_number) -Energy_band(1) ! Hartree
		rclas=rclas_BAND(1:3,1:natot,replica_number)
		vat=vclas_BAND(1:3,1:natot,replica_number)
!guardo en rce y rcg.
           call wrirtc(slabel,Etots,dble(replica_number),replica_number,
     .             na_u,nac,
     .             natot,
     .             rclas,atname,aaname,aanum,nesp,atsym,isa)
           call ioxvconstr(natot,ucell,rclas,vat,replica_number)
           end do
	  call NEB_restart(2, .false.) !write restart full precision
	else
	  write(*,956) rt(1), Etot/eV, Elj/eV, Etot_amber/kcal,
     .  Elink/kcal, Etots/eV
          if(constropt) then
           call subconstr3(ro(1),rt(1),dr,Etots)
! write .rce 
           call wrirtc(slabel,Etots,rt(1),istepconstr,na_u,nac,natot,
     .             rclas,atname,aaname,aanum,nesp,atsym,isa)
           call ioxvconstr( natot, ucell, rclas, vat, istepconstr )
          endif
	end if




! properties calculation in lio for optimized geometry
      if (idyn .ne. 1 .and. qm) then
      do_properties=.true.
      call SCF_hyb(na_u, at_MM_cut_QMMM, r_cut_QMMM, Etot,
     .     F_cut_QMMM,
     .     Iz_cut_QMMM, do_SCF, do_QM_forces, do_properties)
      do_properties=.false.
      end if
      enddo !istepconstr                 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< RESTRAIN Loop



      close(34)

      if(qm) call lio_finalize()

! Dump last coordinates to output
      if (writeipl) then
        if(qm) then
      write(6,'(/a)')'hybrid: Last atomic coordinates (Ang) and species'
      write(6,"(i6,2x,3f10.5,2x,i3)")
     .         (ia, (xa(ix,ia)/Ang, ix=1,3), isa(ia), ia=1,na_u)
        endif !qm

! Dump last Solvent coordinates to output
        if(mm) then
      write(6,'(/a)')'hybrid: Last solvent coordinates (Ang)'
      write(6,'(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)')
     . ('ATOM',i,atname(i),aaname(i),'A',aanum(i),
     . (rclas(ia,na_u+i)/Ang,ia=1,3), i=1,nac)
        endif !mm
      endif !writeipl
      call flush(6)

! Print final date and time
      call timestamp('End of run')

 345  format(2x, I2,    2x, 3(f10.6,2x))
 346  format(2x, f10.6, 2x, 3(f10.6,2x))
 956  format(2x, "Econtribution", 7(f18.6,2x))
 423  format(2x, I6,2x, 6(f20.10,2x))
 444  format(2x, I6,2x, 3(f20.10,2x))
 999  format(a,2x,F30.18)
      end program HYBRID


