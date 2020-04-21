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
      use ionew, only: io_setup!, IOnode
      use scarlett, only: istep, nmove, nesp, inicoor,fincoor, idyn,
     . natot,na_u,nroaa,qm, mm, atsym,
     . xa, fa,
     . masst,isa, iza, pc,
     . nac, atname, aaname, atxres, attype, qmattype, aanum, ng1,
     . bondxat, angexat,
     . angmxat, dihexat,
     . dihmxat, impxat,
     . Em, Rm,
     . fdummy, cfdummy,
     . rclas, vat, aat, izs,
     . linkatom, numlink, linkat, linkqm, linkmm, linkmm2,
     . linkqmtype, Elink, distl, pclinkmm, Emlink, frstme, !pi,
!cutoff
     . blocklist,blockqmmm, blockall,
     . listqmmm,
!NEB
     . NEB_Nimages,
     . NEB_firstimage, NEB_lastimage,
     . rclas_BAND,
     . vclas_BAND, fclas_BAND, Energy_band,
     . NEB_distl,
     . ucell,
     . ftol,
     . Ang, eV, kcal,
!Type 9 restraint
     . rref,rshiftm,rshiftm2,fef,rshiftsd,rclas_cut,natmsconstr,fef_cut,
     . Steep_change, rshxrshm, cov_matrix, cov_matrix_inverted, feopt,
!Dynamics
     . Ekinion, tempion, tempinit, tt, tauber, tempqm, kn, vn, mn,
!!FIRE
     . time_steep, Ndescend, time_steep_max, alpha,
!! LBFGS
     . lbfgs_verbose, lbfgs_num_corr,
!!Lio
     . charge, spin,
!!outputs
     . writeRF, slabel, traj_frec,
!! QM LEVEL
     . qm_level


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
!      logical :: actualiz!MM interaction list control
!!!!
      integer :: nfree !number of atoms without a contrain of movement
      integer :: mmsteps !number of MM steps for each QM step, not tested
      integer :: nbond, nangle, ndihe, nimp !number of bonds, angles, dihedrals and impropers defined in amber.parm
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
      logical :: foundxvr, foundvatr !control for retraint type 9 coordinates restart
      double precision, dimension(:,:), allocatable :: vatr !auxiliary variable for retraint type 9 calculation

! Cut Off QM-MM variables
      double precision :: rcorteqm ! not used
      double precision :: rcorteqmmm ! distance for QM-MM interaction
      double precision :: rcortemm ! distance for LJ & Coulomb MM interaction
      double precision :: radbloqmmm ! distance that allow to move MM atoms from QM sub-system
      double precision :: radblommbond !parche para omitir bonds en extremos terminales, no se computan bonds con distancias mayores a adblommbond
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
      integer :: at_MM_cut_QMMM!, r_cut_pos

! restarts
      logical :: foundcrd
      logical :: foundvat !found velocities in restart

! Outputs variables
      character :: xyzlabel*20 !*.xyz name
      integer  :: nfce !number of atoms for whom the forces will be writen, wrifces=0,1,2 => nfce = 0, na, nat
      integer ::  wricoord !number of steps for write coordinates
      logical :: writeipl


! Auxiliars
      integer :: i, ia, imm, iunit, ix, itest, inneri, j, at1, k !, k, j, jnick, inick
      integer ::  at2 ,k1, k2 ! cov matrix construction
      logical :: leqi

! Others that need check
!!!! General Variables
      real(dp) :: dt !time step
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


! L_BFGS variables
!	integer,  parameter    :: m = 5, iprint = 1
	real(dp), parameter    :: factr  = 1.0d+7, pgtol  = 1.0d-99
	character(len=60)      :: task, csave
	logical                :: lsave(4)
	integer                :: isave(44)
	real(dp)               :: dsave(29)
	integer,  allocatable  :: nbd(:), iwa(:)
	real(dp), allocatable  :: wa(:), limlbfgd(:)
	integer :: looplb

!	rcorteqmmm=0.d0
! Free energy gradient
      logical :: rconverged ! True when rshiftm converged
      double precision :: maxforce
      integer :: maxforceatom, auxiliarunit, auxiliaruniti, i12, j12
      integer :: INFO_inver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!--------------------------------------------------------------------
!need to move this to an initializacion subroutine
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
      do_properties=.false.
      Nick_cent=.false.
      foundxv=.false.
      recompute_cuts=.true.
      cmcf = 0
      imm=1
! Initialize IOnode
      call io_setup

! Print version information
      call prversion

! Start time counter
      call timestamp('Start of run')

! Factors of conversion to internal units
      call init_hybrid('Constants')

! Initialise read
      call reinit(slabel)

! Read and initialize basics variables
      call init_hybrid('Jolie')

! Some definitions
      ucell=0.d0
      fa=0.d0
      fdummy = 0.d0
      cfdummy=0.d0
      vat=0.d0
      aat=0.d0
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
        call MM_atoms_assign(nac, na_u, natot, atname, aaname, rclas,
     . nroaa, aanum, qmattype, rcorteqm, rcorteqmmm, rcortemm,
     . radbloqmmm, radblommbond, radinnerbloqmmm, res_ref, nbond,
     . nangle, ndihe, nimp, attype, pc, Rm, Em, ng1, bondxat,angexat,
     . angmxat, dihexat, dihmxat, impxat)
	else
	rcorteqmmm=0.d0
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
      if (idyn .ne. 1 ) then ! Read cell shape and atomic positions from a former run
      if (usesavexv) then
        call ioxv('read',natot,ucell,rclas,vat,foundxv,foundvat,'X',-1)
        if (foundxv) xa(1:3,1:na_u)=rclas(1:3,1:na_u)
      end if
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
     .               natot,nac,distl)
            xa(1:3,1:na_u)=rclas(1:3,1:na_u)
	  else !NEB case
	    NEB_distl=1.09
	    do replica_number = NEB_firstimage, NEB_lastimage
	      rclas(1:3,1:natot)=rclas_BAND(1:3,1:natot,replica_number)
	      distl(1:15)=NEB_distl(1:15,replica_number)
	      frstme=.true.
	      call link3(numlink,linkat,linkqm,linkmm,rclas,
     .               natot,nac,distl)
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
!      open(unit=34, file=xyzlabel)

! Initialize Lio
      if(qm) then
        Etot=0.d0
        fa=0.d0
#ifdef LIO
        call init_lio_hybrid(1,na_u, nac, charge, iza, spin)
#else
	if (leqi(qm_level,'orca')) then
	else
	  stop 'No QM program defined'
	end if
#endif

! calculate cell volume
        volume = volcel( ucell )

! Center system
!        if (.not.foundxv)call centermol(na_u,xa,rclas,ucell,natot)
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
     .  call vmb(natot,tempinit,masst,vat,cmcf,blockall,ntcon)
!tempinit

!L-BFGS initialization
        if (idyn .eq.-2) then
	write(*,*) "inicializo"
          allocate ( nbd(3*natot),limlbfgd(3*natot))
          allocate ( iwa(9*natot) )
          allocate ( wa(6*lbfgs_num_corr*natot + 15*natot +
     .    11*lbfgs_num_corr**2 + 8*lbfgs_num_corr) )
          nbd=0
          limlbfgd=0.d0
          task = 'START'
        endif

!########################################################################################
!#####################################  MAIN LOOPS  #####################################
!########################################################################################


! Start loop over constrained optimization steps

      if(constropt) then
        call subconstr1(nconstr,typeconstr,kforce,nstepconstr,
     .        rini,rfin,atmsconstr,dr,ro,ndists,coef,constropt)
!        if(feopt .and. typeconstr(1) .ne. 9) STOP
!     .  "feopt selected without typeconstraint 9"
        if (nconstr .eq. 1 .and. typeconstr(1) .eq. 9) then
          allocate(vatr(3,natot))
	 call ioxv('read',natot,ucell,rref,vatr,foundxvr,foundvatr,'r',-1)
!cambiar ucell cuando haya caja
        endif
        if (typeconstr(1) .eq. 9 .and. feopt) then !alocatea cosas para FE
          allocate(rshiftm(3,natot),rshiftm2(3,natot),fef(3,natot),
     .    rshiftsd(3,natot))
        endif
      endif
!      endif
        do istepconstr=1,nstepconstr+1   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< RESTRAIN Loop
! istepconstr marca la posicion del restrain
        optimization_lvl=3
        if (opt_scheme .eq. 1) optimization_lvl=1

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
     .    STOP 'only STEEP, CG, QM, FIRE or NEB minimization available'

          write(6,'(28(" "),a,i6)') 'Begin move = ',istep
          write(6,'(2a)') '                        ',
     .                    '=============================='
          write(6,*) "Optimization level: ", optimization_lvl

!start loot over NEB images
	do replica_number = NEB_firstimage, NEB_lastimage      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Band Replicas
	  if (idyn .eq.1) rclas(1:3,1:natot)=
     .                 rclas_BAND(1:3,1:natot,replica_number)

	  if (NEB_Nimages.gt.1) then
          write(6,'(/2a)') '                        ',
     .           '================================================'
          write(6,'(28(" "),a,i6)') 'Force calculation on Image = ',
     .    replica_number
          write(6,'(2a)') '                        ',
     .           '================================================'
	  end if


	  do imm=1,mmsteps !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MMxQM Steps


       if (.not. feopt) then
       call do_energy_forces(rcorteqmmm, radbloqmmm, Etot,
     . do_SCF, do_QM_forces, do_properties, istp, step,
     . nbond, nangle, ndihe, nimp, Etot_amber, Elj,
     . Etots, constropt,nconstr, nstepconstr, typeconstr, kforce, ro,
     . rt, coef, atmsconstr, ndists, istepconstr, rcortemm,
     . radblommbond, optimization_lvl, dt, sfc, water,
     . imm, rini, rfin)
      else
       call fe_opt(rcorteqmmm, radbloqmmm, Etot,
     .  do_SCF, do_QM_forces, do_properties, istp, step,
     .  nbond, nangle, ndihe, nimp, Etot_amber, Elj,
     .  Etots, constropt,nconstr, nstepconstr, typeconstr, kforce, ro,
     .  rt, coef, atmsconstr, ndists, istepconstr, rcortemm,
     .  radblommbond, optimization_lvl, dt, sfc, water,
     .  imm,rini,rfin,maxforce,maxforceatom,rconverged,ntcon,
     .  nfree,cmcf)
      endif
!       enddo

      call wripdb(na_u,slabel,rclas,natot,step,nac,atname,
     .            aaname,aanum,nesp,atsym,isa,listqmmm,blockqmmm)

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
      write(6,'(43(1h-),/,a4,f22.6,a,i5)') 'Max',fmax*Ang/eV,
     .                       '  free, atom  ',ifmax(2)
      write(6,'(a4,f22.6,a)')'Res',fres*Ang/eV,
     .                       '  sqrt( Sum f_i^2 / 3N )'
      write(6,'(43(1h-),/,a4,f22.6,a,i5)') 'Max',cfmax*Ang/eV,
     .                       '  cons, atom  ',icfmax(2)
      if(nfce.ne.natot) call iofa(natot,cfdummy)

! here Etot in Hartree, cfdummy in Hartree/bohr

      if (mn .eq. 0.d0 .and. idyn .ge. 6) then
        mn=dble(3*natot-ntcon-cmcf)*tt*8.617d-5*(50.d0*dt)**2
        write(6,'(/,a)') 'Calculating Nose mass as Ndf*Tt*KB*(50dt)**2'
        write(6,999) "mn =", mn
      endif

      if (idyn .ne. 1 ) then !Move atoms with a CG algorithm

        if (writeRF .eq. 1) then!save coordinates and forces for integration
	   open(unit=969, file="Pos_forces.dat")
           do itest=1, natot
	      write(969,423) itest, rclas(1:3,itest)*Ang,
     .        cfdummy(1:3,itest)*kcal*Ang/eV  ! Ang, kcal/ang mol
           end do
	   close(969)
        end if

	Ekinion=0.d0


!Move system
	if (idyn .eq. 0) then !Conjugated Gradient
	  call cgvc( natot, rclas, cfdummy, ucell, cstress, volume,
     .             dxmax, tp, ftol, strtol, varcel, relaxd, usesavecg )
	elseif (idyn .eq. 2) then !Quick Minimization
	  call check_convergence(relaxd, natot, cfdummy)
	  if (.not. relaxd) call quick_min(natot, rclas, cfdummy, aat,
     .    vat, masst)
	elseif (idyn .eq. 3) then !FIRE
	  call check_convergence(relaxd, natot, cfdummy)
	  if (.not. relaxd) call FIRE(natot, rclas, cfdummy, vat,
     .    time_steep, Ndescend, time_steep_max, alpha)
	elseif (idyn .eq. 4) then !verlet
	  call verlet2(istp, 3, 0, natot, cfdummy, dt,
     .        masst, ntcon, vat, rclas, Ekinion, tempion, nfree, cmcf)
!iquench lo dejamos como 0, luego cambiar
        elseif (idyn .eq. 5) then !berendsen
          call berendsen(istp,3,natot,cfdummy,dt,tauber,masst,
     .        ntcon,vat,rclas,Ekinion,tempion,tt,nfree,cmcf)
        elseif (idyn .eq. 6) then !nose
          call nose(istp,natot,cfdummy,tt,dt,masst,mn,ntcon,vat,rclas,
     .        Ekinion,kn,vn,tempion,nfree,cmcf)
	elseif (idyn .eq. -1) then !steepest descend
          call check_convergence(relaxd, natot, cfdummy)
          Steep_change = .false.
	  call steep(natot, rclas, cfdummy, Etots, istep)
	elseif (idyn .eq. -2) then !L-BFGS
	  call check_convergence(relaxd, natot, cfdummy)
	  looplb=1
	  do while(looplb.eq.1)
	    call setulb (3*natot, lbfgs_num_corr, rclas, limlbfgd,
     .                   limlbfgd, nbd,
     .                   Etots,-cfdummy, factr, pgtol, wa, iwa, task,
     .                   lbfgs_verbose, csave, lsave, isave, dsave )
	    if(task(1:2).eq.'FG') looplb=0
	  end do
	else
	  STOP "Wrong idyn value"
	end if


       if(idyn .gt. 3 .or. feopt) then

       if(.not. feopt) then
	write(6,999)
     .  'hybrid: Temperature Antes:', tempion, ' K'

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
       endif
!      if(qm) call centerdyn(na_u,rclas,ucell,natot)
	if (MOD((istp - inicoor),traj_frec) .eq. 0) then
       call wrirtc(slabel,Etots,dble(istp),istp,na_u,nac,natot,
     .      rclas,atname,aaname,aanum,nesp,atsym,isa)
       endif
       endif



!Nick center
        if (qm .and. .not. mm .and. Nick_cent) then
          write(*,*) "Centrando Nick"
          firstcent=0
          if (istepconstr.eq.1 .and. istep.eq.inicoor ) firstcent=1
            call center_rotation(natot, masst, rclas, firstcent, Inivec)
        end if

!Write Energy in file
        if (.not. feopt) then
        call wriene(step,slabel,idyn,Etots,cfmax)
        endif
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
     .      natot,nac,distl)
            xa(1:3,1:na_u)=rclas(1:3,1:na_u)
          endif !LA
        endif !qm & mm

! Save last atomic positions and velocities
	call ioxv( 'write',natot,ucell,rclas,vat,foundxv,foundvat,'X',-1)
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
	enddo !imm                          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< <<< Pasos MMxQM REVISAR JOTA

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


!	call NEB_save_traj_energy(istp,slabel)
	call NEB_steep(istp, relaxd, atmsconstr)
!       enddo JOTAJOTAJOTA

	call NEB_save_traj_energy(istp,slabel)

! Calculation Hlink's New positions
        if(qm.and.mm) then
          if(linkatom) then
	    do replica_number = NEB_firstimage, NEB_lastimage
              rclas(1:3,1:natot)=rclas_BAND(1:3,1:natot,replica_number)
	      distl(1:15)=NEB_distl(1:15,replica_number)
	      call link3(numlink,linkat,linkqm,linkmm,rclas,
     .               natot,nac,distl)
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
          if(.not. feopt) then
	  write(*,956) rt(1), Etot/eV, Elj/eV, Etot_amber/kcal,
     .  Elink/kcal, Etots/eV
          if(constropt) then
           call subconstr3(ro(1),rt(1),dr,Etots)
! write .rce
           call wrirtc(slabel,Etots,rt(1),istepconstr,na_u,nac,natot,
     .             rclas,atname,aaname,aanum,nesp,atsym,isa)
           call ioxvconstr( natot, ucell, rclas, vat, istepconstr )
          endif
          endif
	end if




! properties calculation in lio for optimized geometry
!      if (idyn .ne. 1 .and. qm) then
!      do_properties=.true.
!      call SCF_hyb(na_u, at_MM_cut_QMMM, r_cut_QMMM, Etot,
!     .     F_cut_QMMM,
!     .     Iz_cut_QMMM, do_SCF, do_QM_forces, do_properties)
!      do_properties=.false.
!      end if
      enddo !istepconstr                 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< RESTRAIN Loop



!      close(34)


#ifdef LIO
        if(qm) call lio_finalize()
#endif


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

! 345  format(2x, I2,    2x, 3(f10.6,2x))
! 346  format(2x, f10.6, 2x, 3(f10.6,2x))
 956  format(2x, "Econtribution", 7(f22.8,2x))
 423  format(2x, I6,2x, 6(f20.10,2x))
! 444  format(2x, I6,2x, 3(f20.10,2x))
 999  format(a,2x,F30.18)
      end program HYBRID
