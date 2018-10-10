	subroutine init_hybrid(init_type)
!initialize variables for hybrid calculations
!kind of initialization is controlled by  init_type
! N. foglia 03/2018
	use precision, only: dp
	use fdf, only: fdf_integer, fdf_block
	use sys, only: die
	use scarlett, only: natot, aclas_BAND_old, rclas_BAND, vclas_BAND, &
	fclas_BAND, fclas_BAND_fresh, Energy_band, NEB_firstimage, NEB_lastimage, NEB_Nimages, &
	PNEB, PNEB_ini_atom, PNEB_last_atom, NEB_distl, &
	Ang, eV, kcal, na_u, qm, mm, nesp, natoms_partial_freeze, coord_freeze, &
	nac, r_cut_list_QMMM, nparm, izs, Em, Rm, pc, rclas, MM_freeze_list, &
	masst, vat, aat, cfdummy, fdummy, qmattype, attype, atname, aaname, aanum, &
	ng1, blocklist, blockqmmm, listqmmm, fce_amber, ng1type, angetype, & 
	angmtype, evaldihe, evaldihm, dihety, dihmty, impty, nonbonded, &
	scale, evaldihelog, evaldihmlog, scalexat, nonbondedxat, kbond,bondeq, &
	bondtype, kangle,angleeq,angletype, kdihe,diheeq,dihetype, multidihe, &
	perdihe, kimp,impeq, imptype,multiimp, perimp, atange, atangm, atdihe, &
	atdihm, bondxat, angexat, dihexat, dihmxat, angmxat, impxat, atimp, &
	xa, fa, isa, iza, atsym, charge, spin, writeRF, frstme, &
	Ndescend, alpha, NEB_time_steep, NEB_alpha,NEB_Ndescend, time_steep, &
	NEB_move_method, Ndamped, tempion, Nav, pi, blockall
	
	implicit none
	character(len=*), intent(in) :: init_type
	character :: XF, YF, ZF
	character*3 :: XYZF
	integer :: i 
	integer :: iunit
	
	if ( init_type == 'Jolie') then
	  write(*,*) "Hi Angi!" !most important part of code of course
	! Read the number of QM atoms
	  na_u=fdf_integer('NumberOfAtoms',0)
	  if (na_u.eq.0) then
	    write(6,'(/a)') 'hybrid: Running with no QM atoms'
	    qm=.false.
	  endif
	! Read the number of MM atoms
	  nac = fdf_integer('NumberOfSolventAtoms',0)
	  if (nac.eq.0) then
	    write(6,'(/a)') 'hybrid: Running with no MM atoms'
	    mm=.false.
	  endif
	
	  if (nac.eq.0 .and. na_u.eq.0) then
	    call die("no atoms in system")
	  end if
	
	! Read the number of species
	  nesp = fdf_integer('NumberOfSpecies',0)
	  if(qm.and.(nesp.eq.0)) then
	    call die("hybrid: You must specify the number of species")
	  endif
	
	  allocate(xa(3,na_u), fa(3,na_u), isa(na_u), iza(na_u), atsym(nesp))
	
	! Read number of partial freeze atoms
	  natoms_partial_freeze = fdf_integer('NumberOfPartialFreeze',0)
	  if ( natoms_partial_freeze .gt. 0 ) then
	    allocate(coord_freeze(natoms_partial_freeze,4))
	! Read partial freeze atoms
	    if ( fdf_block('PartialFreeze',iunit) ) then
	      coord_freeze=0
	      do i=1,natoms_partial_freeze
	        read(iunit,*,err=2,end=2) coord_freeze(i, 1:4) !atom number, xfreeze, yfreeze, zfreeze
	      enddo
	
	      do i=1,natoms_partial_freeze
	        XF=" "
	        YF=" "
	        ZF=" "
	        XYZF=""
	        if (coord_freeze(i,2) .eq.1) XF="X"
	        if (coord_freeze(i,3) .eq.1) YF="Y"
	        if (coord_freeze(i,4) .eq.1) ZF="Z"
	        write(XYZF,"(A1,A1,A1)") XF, YF, ZF
	        write(*,*) "freezing atom: ", coord_freeze(i,1), "coord(s)", XYZF
	      enddo
	    endif
	  endif
	
	
	! Read QM coordinates and labels
	  write(6,*)
	  write(6,"('read:',71(1h*))")
	  if(qm) call read_qm(na_u,nesp,isa,iza,xa,atsym,charge, spin)
	
	! Allocation of solvent variables
	  natot = nac + na_u
	  allocate(r_cut_list_QMMM(nac)) ! referencia posiciones de atomos MM con los vectores cortados
	  nparm = 500 ! numero de tipos de bonds q tiene definido el amber.parm. NO DEBERIA ESTAR fijo, Nick
	!muchos allocate tienen valores fijos. habria que reveer esto en el futuro. Nick
	  allocate(izs(natot), Em(natot), Rm(natot), pc(0:nac))
	  allocate(rclas(3,natot), MM_freeze_list(natot), masst(natot))
	  allocate(vat(3,natot),aat(3,natot), fdummy(3,natot),cfdummy(3,natot))
	  allocate(qmattype(na_u), attype(nac), atname(nac))
	  allocate(aaname(nac), aanum(nac), ng1(nac,6), blocklist(natot))
	  allocate(blockqmmm(nac), listqmmm(nac), fce_amber(3,nac))
	  allocate(blockall(natot)) !JOTA
	  allocate(ng1type(nac,6), angetype(nac,25), angmtype(nac,25))
	  allocate(evaldihe(nac,100,5), evaldihm(nac,100,5))
	  allocate(dihety(nac,100), dihmty(nac,100), impty(nac,25))
	  allocate(nonbonded(nac,100), scale(nac,100), evaldihelog(nac,100))
	  allocate(evaldihmlog(nac,100), scalexat(nac))
	  allocate(nonbondedxat(nac))
	  allocate(kbond(nparm),bondeq(nparm),bondtype(nparm))
	  allocate(kangle(nparm),angleeq(nparm),angletype(nparm))
	  allocate(kdihe(nparm),diheeq(nparm),dihetype(nparm), multidihe(nparm), perdihe(nparm))
	  allocate(kimp(nparm),impeq(nparm), imptype(nparm),multiimp(nparm), perimp(nparm))
	  allocate(atange(nac,25,2), atangm(nac,25,2), atdihe(nac,100,3))
	  allocate(atdihm(nac,100,3), bondxat(nac), angexat(nac))
	  allocate(dihexat(nac), dihmxat(nac), angmxat(nac))
	  allocate(impxat(nac), atimp(nac,25,4))
	
	  writeRF=0
	  writeRF = fdf_integer('PFIntegrationOutput',0)
	  frstme=.true.
	  Ndescend=0
	  alpha=0.1d0
	  Ndamped=0
	  tempion=0.d0
	elseif ( init_type == 'Constants') then !define constants and convertion factors
	  Ang    = 1._dp / 0.529177_dp
	  eV     = 1._dp / 27.211396132_dp
	  kcal   = 1.602177E-19_dp * 6.022140857E23_dp / 4184.0_dp !23.06055347_dp
	  pi     = DACOS(-1.d0)
	  Nav    = 6.022140857d23
	
	elseif ( init_type == 'NEB') then !initialize Nudged elastic band variables
	  if (NEB_Nimages .lt. 3) STOP 'Runing NEB with less than 3 images'
	  NEB_firstimage=1
	  NEB_lastimage=NEB_Nimages
	  allocate(rclas_BAND(3,natot,NEB_Nimages), vclas_BAND(3,natot,NEB_Nimages), &
	         fclas_BAND(3,natot,NEB_Nimages), aclas_BAND_old(3,natot,NEB_Nimages))
	  allocate (fclas_BAND_fresh(3,natot,NEB_Nimages))
	  allocate(Energy_band(NEB_Nimages))
	  allocate(NEB_distl(15, NEB_Nimages))
	  rclas_BAND=0.d0
	  vclas_BAND=0.d0
	  fclas_BAND=0.d0
	  Energy_band=0.d0
	  PNEB=fdf_integer('PNEB',0)

	  if ( PNEB .eq.1 ) then
	    PNEB_ini_atom=fdf_integer('PNEBi',1)
	    PNEB_last_atom=fdf_integer('PNEBl',natot)
	  end if
	
	  if (NEB_move_method .eq.3) then
	    allocate(NEB_time_steep(NEB_Nimages), NEB_alpha(NEB_Nimages),NEB_Ndescend(NEB_Nimages))
	    NEB_time_steep=time_steep
	    NEB_alpha=alpha
	    NEB_Ndescend=0
	  end if

	else
	  STOP "Wrong init_type"
	end if
	
	return
 2     stop 'read: problem reading partial freeze block'
	end subroutine init_hybrid

