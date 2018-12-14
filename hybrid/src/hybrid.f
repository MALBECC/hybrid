      Program HYBRID 

C****************************************************************************
C Program for QM-MM calculations. 
C It uses the SIESTA code to treat at DFT level the QM subsystem. 
C It uses an own implementation of Amber99 force fiel parametrization
C to treat the MM subsystem. 
C Original idea by D. Scherlis, D. Estrin and P. Ordejon. 2001.
C Original QM-MM interfase with Siesta by D. Scherlis and P. Ordejon. 2001.
C Original subroutines by D. Scherlis. 2001/2002.
C Solvent implementation using the Amber99 force field parametrization
C by A. Crespo and M. Marti. 2002/2003.
C LinkAtom subrouitnes by M. Marti. 2003.
C Further modifications of QM-MM Siesta by A. Crespo and P. Ordejon. 2004.
C Developing of HYBRID: program for QM-MM calculations using the SIESTA
C as SUBROUTINE by A. Crespo and P. Ordejon. 2004.
C This code uses several subroutines of SIESTA  as well as the 
C SiestaAsSubroutine and FDF packages.
C Crespo, May 2005. 
C*****************************************************************************

C Modules
      use fsiesta
      use precision
      use sys
      use fdf
      use neutralatom
      use ionew, only: io_setup, IOnode    

C General Variables
      implicit none

      integer
     .  fincoor, i, ia, ianneal, idyn, ifinal, inicoor, iquench, Nodes, 
     .  istp, istart, istep, iunit, ix, j, k, nmove, nspin, ntm(3), 
     .  ntcon, na_u, nesp, ntpl, nsm,  mullipop, step, 
     .  ia1, ia2, iadispl, ixdispl, ifmax(2), icfmax(2)

      real(dp)
     .  Ang, cfmax, cftem, cstress(3,3), dt, Ekinion, Etot, eV, 
     .  fmax, fmean, fres, ftem, ftol, ftot(3), mn, Pint, strtol, 
     .  bulkm, dxmax, dx

      real(dp)
     .  tauber, taurelax, tempinit, tempion, tp, tt, ucell(3,3), kn,  
     .  vcell(3,3), vn, volcel, volume, dvol, vna, grvna(3)

      real(dp), dimension(:,:), allocatable :: 
     .  xa, fa

      real, dimension(:,:), allocatable ::
     .  DRho

      real, dimension(:), allocatable ::
     .  Vqm, Rho

      integer, dimension(:), allocatable ::
     .  isa,iza

      character, dimension(:), allocatable ::
     .  atsym*2

      logical
     .  found, foundxv, usesavecg,  usesavexv, varcel, relaxd, qm , mm

      character
     .  slabel*20, sname*150, siestaslabel*20, 
     .  filerho*30, filevqm*30, paste*25 

      parameter ( nsm = 2 )

      external
     .  anneal, berendsen, chkdim, cgvc,  ioxv, iorho, fixed1, assign, 
     .  iofa, iovext, ofc, memory, nose, paste, reinit, read_qm, 
     .  read_md, timer, verlet2, vmb, volcel, write_lab, fixed2

C Solvent General variables

      character*4,  dimension(:), allocatable, save ::
     .   atname,aaname,attype,qmattype

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

      integer nac,natot,nfree,nparm,nfce, mmsteps, imm, 
     .  nroaa, nbond, nangle, ndihe, nimp , wricoord,
     .  atxres(20000), totcoor

      double precision  Etot_amber,sfc,timestep,
     .  Etots,kcal,rcorteqm,rcortemm,Elj,Eelec

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

C more variables
        logical water,foundcrd,foundvat,writeipl 

C Solvent external variables
        external
     .  solv_assign, solv_ene_fce, qmmm_lst_blk, wrtcrd,
     .  centermol, centerdyn, link1, link2, link3, pcpot, ljef,
     .  mmForce, readcrd,  subfree, prversion, ioxvconstr, elecfield, 
     .  wripdb, wriene, wrirtc, subconstr1, subconstr2, subconstr3

C variables para esquema de optimizacion
	integer :: opt_scheme
	integer :: optimization_lvl
C 1- solo considera restrain
C 2- solo opimiza atomos MM
C 3- optimiza todo


C--------------------------------------------------------------------

C Initialise some variables
        relaxd=.false.
        varcel=.false.
        tempion=0.0
        tp=0.0
        vcell=0.0
        cstress=0.0
        strtol=0.0
        Pint=0.0
        bulkm=0.0
        ntm=0
        nspin=1
        ntpl=1
        qm=.true.
        mm=.true.

C Initialize IOnode
      call io_setup   

C Print version information
      call prversion
      write(6,'(/,a,i4,a)') 'hybrid: Running in serial mode' 

C Start time counter 
      call timestamp('Start of run')      
      call timer( 'all', 0 )
      call timer( 'Hybrid', 1 )
      call timer( 'Setup', 1 )

C Factors of conversion to internal units 
      Ang    = 1._dp / 0.529177_dp
      eV     = 1._dp / 13.60580_dp
      kcal   = 1.602177E-19_dp * 6.02214E23_dp / 4184.0_dp

C Initialise read 
      call reinit(slabel, sname)

C set processor number for the QM siesta calculation
      Nodes=fdf_integer('NumberMPInodes',1)
C Read the number of QM atoms
      na_u=fdf_integer('NumberOfAtoms',0)
      if (na_u.eq.0) then
c      call die("You must specify number of QM atoms")
      write(6,'(/a)') 'hybrid: Running with no QM atoms'
      qm=.false.
      endif
C Read the number of MM atoms
      nac = fdf_integer('NumberOfSolventAtoms',0)
      if (nac.eq.0) then
c      call die("You must specify number of solvent atoms")
      write(6,'(/a)') 'hybrid: Running with no MM atoms'
      mm=.false.
      endif
C Read the number of species
      nesp = fdf_integer('NumberOfSpecies',0)
      if(qm.and.nesp.eq.0) then
      call die("hybrid: You must specify the number of species")
      endif

      allocate(xa(3,na_u))
      call memory('A','D',3*na_u,'hybrid')
      allocate(fa(3,na_u))
      call memory('A','D',3*na_u,'hybrid')       
      allocate(isa(na_u))
      call memory('A','D',na_u,'hybrid')
      allocate(iza(na_u))
      call memory('A','D',na_u,'hybrid')
      allocate(atsym(nesp))
      call memory('A','D',nesp,'hybrid')
 
C Read QM coordinates and labels
      write(6,*)
      write(6,"('read:',71(1h*))")
      if(qm) then
      call read_qm(na_u,nesp,isa,iza,xa,atsym)
      endif !qm

C Allocation of solvent variables
      natot = nac + na_u
      nparm = 500

      allocate(izs(natot))
      call memory('A','D',nac,'hybrid')
      allocate(Em(natot))
      call memory('A','D',natot,'hybrid')
      allocate(Rm(natot))
      call memory('A','D',natot,'hybrid')
      allocate(pc(nac))
      call memory('A','D',nac,'hybrid')
      allocate(rclas(3,natot))
      call memory('A','D',3*natot,'hybrid')
      allocate(masst(natot))
      call memory('A','D',natot,'hybrid')
      allocate(vat(3,natot))
      call memory('A','D',3*natot,'hybrid')
      allocate(cfdummy(3,natot))
      call memory('A','D',3*natot,'hybrid')
      allocate(fdummy(3,natot))
      call memory('A','D',3*natot,'hybrid')
      allocate(qmattype(na_u))
      call memory('A','D',na_u,'hybrid')
      allocate(attype(nac))
      call memory('A','D',nac,'hybrid')
      allocate(atname(nac))
      call memory('A','D',nac,'hybrid')
      allocate(aaname(nac))
      call memory('A','D',nac,'hybrid')
      allocate(aanum(nac))
      call memory('A','D',nac,'hybrid')
      allocate(ng1(nac,6))
      call memory('A','D',6*nac,'hybrid')
      allocate(blocklist(natot))
      call memory('A','D',nac,'hybrid')
      allocate(blockqmmm(nac))
      call memory('A','D',nac,'hybrid')
      allocate(listqmmm(nac))
      call memory('A','D',nac,'hybrid')
      allocate(fce_amber(3,nac))
      call memory('A','D',3*nac,'hybrid')
      allocate(ng1type(nac,6))
      call memory('A','D',6*nac,'hybrid')
      allocate(angetype(nac,25))
      call memory('A','D',25*nac,'hybrid')
      allocate(angmtype(nac,25))
      call memory('A','D',25*nac,'hybrid')
      allocate(evaldihe(nac,100,5))
      call memory('A','D',500*nac,'hybrid')
      allocate(evaldihm(nac,100,5))
      call memory('A','D',500*nac,'hybrid')
      allocate(dihety(nac,100))
      call memory('A','D',100*nac,'hybrid')
      allocate(dihmty(nac,100))
      call memory('A','D',100*nac,'hybrid')
      allocate(impty(nac,25))
      call memory('A','D',25*nac,'hybrid')
      allocate(nonbonded(nac,100))
      call memory('A','D',100*nac,'hybrid')
      allocate(scale(nac,100))
      call memory('A','D',100*nac,'hybrid')
      allocate(evaldihelog(nac,100))
      call memory('A','D',100*nac,'hybrid')
      allocate(evaldihmlog(nac,100))
      call memory('A','D',100*nac,'hybrid')
      allocate(qmpc(na_u,2))
      call memory('A','D',2*na_u,'hybrid')
      allocate(qmpcf(na_u))
      call memory('A','D',na_u,'hybrid')
      allocate(desplaz(nac))
      call memory('A','D',nac,'hybrid')
      allocate(scalexat(nac))
      call memory('A','D',nac,'hybrid')
      allocate(nonbondedxat(nac))
      call memory('A','D',nac,'hybrid')
      allocate(kbond(nparm),bondeq(nparm),
     .           bondtype(nparm))
      call memory('A','D',3*nparm,'hybrid')
      allocate(kangle(nparm),angleeq(nparm),
     .           angletype(nparm))
      call memory('A','D',3*nparm,'hybrid')
      allocate(kdihe(nparm),diheeq(nparm),
     .           dihetype(nparm),multidihe(nparm),
     .            perdihe(nparm))
      call memory('A','D',5*nparm,'hybrid')
      allocate(kimp(nparm),impeq(nparm),
     .           imptype(nparm),multiimp(nparm),
     .            perimp(nparm))
      call memory('A','D',5*nparm,'hybrid')
      allocate(atange(nac,25,2))
      call memory('A','D',50*nac,'hybrid')
      allocate(atangm(nac,25,2))
      call memory('A','D',50*nac,'hybrid')
      allocate(atdihe(nac,100,3))
      call memory('A','D',300*nac,'hybrid')
      allocate(atdihm(nac,100,3))
      call memory('A','D',300*nac,'hybrid')
      allocate(bondxat(nac))
      call memory('A','D',nac,'hybrid') 
      allocate(angexat(nac))
      call memory('A','D',nac,'hybrid')
      allocate(dihexat(nac))
      call memory('A','D',nac,'hybrid')
      allocate(dihmxat(nac))
      call memory('A','D',nac,'hybrid')
      allocate(angmxat(nac))
      call memory('A','D',nac,'hybrid')
      allocate(impxat(nac))
      call memory('A','D',nac,'hybrid')
      allocate(atimp(nac,25,4))
      call memory('A','D',100*nac,'hybrid')

C Some definitions 
      ucell=0.0
      fa=0.0
      fdummy = 0.0
      cfdummy=0.0
      vat=0.0
      Elj=0.0
      Elink=0.0
      Etot=0.0
      Etots=0.0
      Eelec=0.0
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
     .  radbloqmmm,atxres)
       endif !mm
        rclas(1:3,1:na_u) = xa(1:3,1:na_u)

C Read simulation data 
       call read_md( ianneal, idyn, ifinal, iquench, istart, nmove,
     .               dt, dxmax, ftol, mn, tauber, taurelax, tempinit,
     .               tt, usesavecg, usesavexv , na_u, dx, ia1, ia2,
     .               natot, nfce, wricoord, mullipop, mmsteps )

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
          Em(linkat(i))=0.0
          Em(linkqm(i,1:1))=0.0
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
        pc(linkmm(i,1:1))=0.0
        pc(linkmm(i,2:4))=pc(linkmm(i,2:4))+pclinkmm(i,1)/3.0
C two options: Em(linkmm(i,1:1))=0.0 or Em(linkmm(i,1:4))=0.0
        Emlink(i,1:1)=Em(na_u+linkmm(i,1:1))
        Em(na_u+linkmm(i,1:1))=0.0
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

C Initialize siestaslabel, filerho and filevqm 
      if(qm) then
        siestaslabel = paste( slabel, '.siesta' )
        filerho = paste( siestaslabel, '.DRHO' )
        filevqm = paste( siestaslabel, '.VEXT' ) 

C Write the 'slabel.siesta.fdf' file to be read by siesta as subroutine
      call write_lab( slabel, siestaslabel)

! Set physical units for communication with siesta
      call siesta_units( 'Bohr', 'Ry' )

! Launch siesta proces
      call siesta_launch( siestaslabel, Nodes )

! First call to siesta_forces to initialise grid
       call siesta_forces( siestaslabel, na_u, xa, ucell, Etot, fa )
       Etot=0.0
       fa=0.0

C Read grid information
       call iorho( 'read', filerho, ucell, ntm, nsm, ntpl, nspin,
     .              0.0d0, found )
       if(.not.found) then
       call die("iorho: Unable to read grid information")
       endif 

C Read the spin polarized variable and find the number of diagonal spin values
      if(nspin.gt.2) stop 'hybrid: nspin could be only 1 or 2'
      allocate(DRho(ntpl,nspin))
      call memory('A','D',ntpl*nspin,'hybrid')
      allocate(Vqm(ntpl))
      call memory('A','D',ntpl,'hybrid')     
      allocate(Rho(ntpl))
      call memory('A','D',ntpl,'hybrid')

C calculate cell volume
      volume = volcel( ucell )
      dvol=volume/ntpl

C writes cell
         write(6,'(/,a,3(/,a,3f12.6))')
     .       'hybrid: unit cell vectors (Ang) from siesta run:',
     .       ('    ' , (ucell(ix,ia)/Ang,ix=1,3), ia =1,3)

C Read Neutral Atom potential  
       call  read_na(nesp,atsym)

C Center system 
       if (.not.foundxv) call centermol(na_u,xa,rclas,ucell,natot)
      endif !qm

C Calculate Rcut & block list QM-MM 
      if(qm.and.mm) then
        call qmmm_lst_blk(na_u,nac,natot,nroaa,atxres,rclas,
     .  rcorteqmmm,radbloqmmm,blockqmmm,listqmmm,rcorteqm,slabel)
      endif !qm & mm

C Read fixed atom constraints
      call fixed1(na_u,nac,natot,nroaa,rclas,blocklist,
     .            atname,aaname,aanum,water)

C Build initial velocities according to Maxwell-Bolzmann distribution
        if (idyn .ne. 0 .and. idyn .ne. 5 .and. (.not. foundvat) ) 
     .  call vmb(natot,tempinit,masst,rclas,0,vat)

      call timer( 'Setup', 2 )
      call printmemory( 6, 0 )

C Start loop over constrained optimization steps
        if(constropt) then
        call subconstr1(nconstr,typeconstr,kforce,nstepconstr,
     .  rini,rfin,atmsconstr,dr,ro,ndists,coef,constropt)
        endif
        allocate(qmstep(nstepconstr+1))
        call memory('A','D',nstepconstr+1,'hybrid')
        qmstep=0

        do istepconstr=1,nstepconstr+1
        if(constropt) then
        write(6,*)
        write(6,'(A)')    '*******************************'
        write(6,'(A,i5)') '  Constrained Step : ', istepconstr
        write(6,'(A)')    '*******************************'
        endif

C Begin of coordinate relaxation iteration ============================
      if (idyn .eq. 0) then
        inicoor = 0
        fincoor = nmove
      else if (idyn .eq. 5) then
        inicoor = 0
        fincoor = (ia2-ia1+1)*3*2
      else 
        inicoor = istart
        fincoor = ifinal
      endif
        totcoor = fincoor - inicoor + 1

C Start loop over coordinate changes
      istp = 0
      do istep = inicoor,fincoor
      call timer( 'IterMD', 1 )      
      istp = istp + 1
 
        write(6,'(/2a)') 'hybrid:                 ',
     .                    '=============================='
        if (idyn .eq. 0 ) 
     .   write(6,'(28(" "),a,i6)') 'Begin CG move = ',istep
        if (idyn .ne. 0 .and. idyn .ne. 5)
     .   write(6,'(28(" "),a,i6)') 'Begin MD step = ',istep
        if (idyn .eq. 5)  then
         write(6,'(28(" "),a,i6)') 'Begin FC step = ',istep
          if (istep .eq. 0) then
            write(6,'(28(" "),a)') 'Undisplaced coordinates'
          else
            iadispl = (istep-mod(istep-1,6))/6+ia1
            write(6,'(28(" "),a,i6)') 'displace atom   ',
     .        iadispl
            ix = mod(istep-1,6)+1
            ixdispl = (ix - mod(ix-1,2) +1)/2
            write(6,'(28(" "),a,i6)') 'in direction    ',
     .        ixdispl
            dx=-dx
            write(6,'(28(" "),a,f8.4,a)') 'by       ',
     .                      dx, ' Bohr'

C Displace atom by dx
            rclas(ixdispl,iadispl)=rclas(ixdispl,iadispl)+dx
            xa(1:3,1:na_u)=rclas(1:3,1:na_u)
          endif
        endif
        write(6,'(2a)') '                        ',
     .                    '=============================='

C Calculate the external potential Vqm
      if(qm.and.mm) then
      call pcpot(na_u,nac,natot,ntpl,ntm,listqmmm,rclas,ucell,
     .           pc,rcorteqm,Vqm)

C Write Vqm to file 
      call iovext( 'write', filevqm, ucell, ntm, nsm, ntpl, Vqm, found )
      endif !qm & mm

! Calculate Energy and Forces using Siesta as Subroutine
      if(qm) then
       call timer( 'Siesta', 1 )
       call siesta_forces( siestaslabel, na_u, xa, ucell, Etot, fa )
       call timer( 'Siesta', 2 )   
       qmstep(istepconstr) = qmstep(istepconstr) + 1
      endif !qm

C Read DRho from file 
      if(qm.and.mm) then
       call iorho( 'read', filerho, ucell, ntm, nsm, ntpl, nspin,
     .              DRho, found )
       if(.not.found) then
       call die("iorho: Unable to read DRho from file")
       endif

C Calculation of total density difference
      if(nspin.eq.1) then
      Rho(1:ntpl) = DRho(1:ntpl,1)
      elseif(nspin.eq.2) then
      Rho(1:ntpl) = DRho(1:ntpl,1) + DRho(1:ntpl,2)
      endif
      endif !qm & mm

C Start MMxQM loop
      do imm=1,mmsteps
      step = step +1
      if(mmsteps.ne.1) then
        write(6,*)
        write(6,'(A)')    '*******************************'
        write(6,'(A,i5)') '   MM x QM Step : ', imm 
        write(6,'(A)')    '*******************************'
      endif

C Calculation of Forces over the classical atoms from QM-MM interaction
      if(qm.and.mm) then
      call timer( 'QM-MMcoupl', 1 )
      call mmforce(na_u,nac,natot,nesp,ntpl,nspin,ntm,isa,listqmmm,
     .             rclas,fdummy,pc,ucell,rcorteqm,dvol,Rho)
      endif !qm & mm
 
C Add fa to fdummy for QM atoms
      if(qm) fdummy(1:3,1:na_u) = fa(1:3,1:na_u)

C Calculation of last QM-MM interaction: LJ Energy and Forces only 
      if(qm.and.mm) then
      call ljef(na_u,nac,natot,rclas,Em,Rm,fdummy,Elj,listqmmm)
      call timer( 'QM-MMcoupl', 2 )
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
      call timer( 'MMenergy', 1 )
      call solv_ene_fce(natot,na_u,nac,ng1,rclas,Em,Rm,pc,
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
     .    water,masst)
      call timer( 'MMenergy', 2 )
      endif !mm

C converts fdummy to Kcal/mol/Ang  
      fdummy(1:3,1:natot)=fdummy(1:3,1:natot)*Ang/eV*kcal
 
C add famber to fdummy  
      if(mm) then
      fdummy(1:3,na_u+1:natot)=fdummy(1:3,na_u+1:natot)
     .       +fce_amber(1:3,1:nac)
      endif !mm

C Calculation of LinkAtom Energy and Forces
      if(qm.and.mm) then
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
        pc(linkmm(i,1:1))=0.0
        pc(linkmm(i,2:4))=pc(linkmm(i,2:4))+pclinkmm(i,1)/3.0
        Em(na_u+linkmm(i,1:1))=0.0
        enddo
        endif ! LA
      endif !qm & mm

C Calculation of Constrained Optimization Energy and Forces 
      if(imm.eq.1) then
        if(constropt) then
        call subconstr2(nconstr,typeconstr,kforce,rini,rfin,ro,rt,
     .  nstepconstr,atmsconstr,natot,rclas,fdummy,istp,istepconstr,
     .  ndists,coef)
        endif 
      endif !imm

C Free Energy calculation using Jarzynski formula 
      if(imm.eq.1) then
        if(free) then
        call subfree(natot,rclas,fdummy,dt,free,totcoor)
        endif
      endif !imm

C call to subroutine that calculates the Electric field energy
      if(mm) then
      call elecfield(nac,rclas(1:3,na_u+1:natot),
     .     fdummy(1:3,na_u+1:natot),pc,Eelec)
      endif !mm

C Converts fdummy to Ry/Bohr 
      fdummy(1:3,1:natot)=fdummy(1:3,1:natot)/Ang*eV/kcal

C Writes final energy decomposition
       Etots=Etot+Elj+((Etot_amber+Elink+Eelec)/kcal*eV)
       write(6,*)
       write(6,'(/,a)') 'hybrid: Energy Decomposition (eV):'
       if(qm) write(6,'(a,2x,F16.6)')           'Esiesta:',Etot/eV
       if(qm.and.mm) write(6,'(a,2x,F16.6)')    'Elj:    ',Elj/eV
       if(mm) write(6,'(a,2x,F16.6)')      'Esolv:  ',Etot_amber/kcal
       if(Elink.ne.0.0) write(6,'(a,2x,F16.6)') 'Elink:  ',Elink/kcal
       if(Eelec.ne.0.0) write(6,'(a,2x,F16.6)') 'Eefield:',Eelec/kcal
       if(qm.and.mm) write(6,'(a,2x,F16.6)')    'Etots:  ',Etots/eV
       call flush(6)

C Sets fdummy to zero inside mmxqm step
       if(qm.and.mm) then
       if(imm.ne.1) then
       fdummy(1:3,1:na_u) = 0.0
       vat(1:3,1:na_u) = 0.0
       if(linkatom) then
       do i=1,numlink
       fdummy(1:3,linkmm(i,1)+na_u)=0.0
       vat(1:3,linkmm(i,1)+na_u)=0.0
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

C Write Force Constant matrix if FC calculation
        if (idyn .eq. 5) call ofc(na_u,ia1,ia2,dx,cfdummy(1:3,1:na_u))

C Move atoms 
      if (idyn .eq. 0 ) then 
          call cgvc( natot, rclas, cfdummy, ucell, cstress, volume,
     .             dxmax, tp, ftol, strtol, varcel, relaxd, usesavecg )
      endif

      Ekinion  = 0.0_dp
      vn       = 0.0_dp
      kn       = 0.0_dp
      iunit = 2
      ntcon = 0
      if(imm.ne.1) ntcon = na_u

      if (idyn .eq. 1) then
       call verlet2(istp, iunit, iquench, natot, cfdummy, dt,
     .        masst, ntcon, vat, rclas, Ekinion, tempion, nfree)
      elseif (idyn .eq. 2) then
      call nose(istp, iunit, natot, cfdummy, tt, dt, masst, mn,
     .        ntcon, vat, rclas, Ekinion, kn, vn, tempion, nfree)
      elseif (idyn .eq. 3) then
      call anneal(istp, iunit, ianneal, taurelax, bulkm,
     .        natot, cfdummy, cstress, tp, tt, dt, masst, ntcon,
     .        vat, rclas, ucell, Ekinion, tempion, Pint, nfree)
      elseif (idyn .eq. 4) then
      call berendsen(istp,iunit,natot,cfdummy,dt,tauber,masst,
     .        ntcon,vat,rclas,Ekinion,tempion,tt,nfree)       
      endif

       if(idyn .ne. 0 .and. idyn .ne. 5) then
        write(6,'(/,a,f12.6)')
     .  'hybrid: Kinetic Energy (eV):',Ekinion/eV
        write(6,'(/,a,f12.6)')
     .  'hybrid: Total Energy + Kinetic (eV):',(Etots+Ekinion)/eV
        write(6,'(/,a,f12.3,a)')
     .  'hybrid: System Temperature:', tempion, ' K'
      if(qm) call centerdyn(na_u,rclas,ucell,natot)
       endif

C Write Energy in file
      call wriene(step,slabel,idyn,Etots,cfmax,Ekinion,tempion)

c sets variables for next siesta cycle
      fa = 0.0
      fdummy = 0.0 
      cfdummy = 0.0
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

C Restore original coordinates after FC displacements
      if (idyn .eq. 5 .and. istep .ne. 0) then
       rclas(ixdispl,iadispl)=rclas(ixdispl,iadispl)-dx
      endif

C Exit MMxQM loop
      enddo !imm
      if((mmsteps.ne.1).and.(imm.ne.1)) relaxd = .false.

C Exit coordinate relaxation loop
      call timer( 'IterMD', 2 )     
      if (relaxd) goto 10 
       enddo
 10   continue
C End of coordinate-relaxation loop ==================================

C End of constrain optimization loop
      if(constropt) then
       call subconstr3(ro(1),rt(1),dr,Etots)
       call wrirtc(slabel,Etots,rt(1),istepconstr,na_u,nac,natot,
     .             rclas,atname,aaname,aanum,nesp,atsym,isa)
       call ioxvconstr( natot, ucell, rclas, vat, istepconstr )
      endif
      enddo !istepconstr

! Stop siesta proces 
      if(qm) call siesta_quit( 'all' )

C Mulliken population analysis
      if(qm) then
      if(mullipop.eq.1) then
      call mulliken(na_u,isa,iza,nesp,atsym,siestaslabel,
     .              qmstep,nstepconstr+1)
      endif
      endif !qm

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

C Stop time counter
      call timer( 'Hybrid', 2 )
      call timer( 'all', 3 )
      call printmemory( 6, 0 )

C Print final date and time
        call timestamp('End of run')

      end


