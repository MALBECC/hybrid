!****************************************************************************
! This routine read and assign parameters, coordinates and conectivity to
! all MM and QM/MM atoms.
! modified from oroginal solv_assign subroutine
! N. Foglia 2019
!*****************************************************************************

	subroutine MM_atoms_assign(nac, na_u, natot, atname, aaname, rclas, &
	nroaa, aanum, qmattype, rcorteqm, rcorteqmmm, rcortemm, radbloqmmm, &
	radblommbond, radinnerbloqmmm, res_ref, nbond,nangle)


!        (na_u,natot,nac,nroaa,Em,Rm,attype,pc,
!     .  ng1,bondxat,angexat,atange,angmxat,atangm,dihexat,atdihe,
!     .  dihmxat,atdihm,impxat,atimp,
!     .  nbond,kbond,bondeq,bondtype,
!     .  nangle,kangle,angleeq,angletype,
!     .  ndihe,kdihe,diheeq,dihetype,multidihe,perdihe,
!     .  nimp,kimp,impeq,imptype,multiimp,perimp,
!     .  nparm,aaname,atname,aanum,qmattype,rclas,
!     .  rcorteqmmm,rcorteqm,rcortemm,sfc,
!     .  radbloqmmm,atsinres,radblommbond,radinnerbloqmmm,res_ref)


!relacion entre variables nuevas y viejas para testeo:
! fdf_at_name=atname
! fdf_res_name=aaname
! fdf_res_num=resnum


	use precision, only : dp
	use sys,only : die
	use fdf,only : fdf_block
	use scarlett, only: Ang

	implicit none
	integer, intent(in) :: nac, na_u, natot !number of QM and MM and total atoms
	character*4, dimension(nac), intent(out) :: atname, aaname !atom and residue name in .fdf
	integer, dimension(nac) :: resnum !residue number in .fdf
	double precision, dimension(3,natot), intent(inout) :: rclas !postion of all atoms
	integer :: ivalue(natot*2) !relation between atom number and atom position in .fdf
	integer :: iunit ! unit number to read
	integer, intent(out) :: nroaa !number of residues
	integer, dimension(nac), intent(out) :: aanum !residue number for each atom
	character*4, dimension(na_u), intent(out) :: qmattype(na_u)
	integer, intent(out) :: res_ref ! residue that is taken as reference to fix atoms by radial criteria in full MM simulations
	double precision, intent(out) :: rcorteqm ! not used. not removed for keep compatibility with old .fdf files
	double precision, intent(out) :: rcorteqmmm, rcortemm !cut off criteria for QM/MM and MM interacions
	double precision, intent(out) :: radbloqmmm, radinnerbloqmmm ! freeze criteria
	double precision, intent(out) :: radblommbond !cut off criteria for MM bonds
	integer :: ncon !aditional conectivities defined in .fdf
	integer, dimension(:,:), allocatable :: con     
	logical :: foundamber

!
	integer, intent(out) :: nbond, nangle
!	double precision, dimension(:), intent(out) :: kbond,bondeq
!	character*5, dimension(:), intent(out) :: bondtype


!auxiliars
	integer :: i, j, k
	character*4 :: atom, ch4
	character*1 ::  ch, ch1, exp
	integer, dimension(2,1000) :: con2

!        integer i,j,k,l,m,n,natot,nac,na_u,nroaa,iunit,                                                              
!     .  ncon,nparm,nbond,nangle,ndihe,nimp,atsinres(20000)     
!        integer ng1(nac,6),con2(2,1000)
!        double precision, dimension(:,:), allocatable, save::
!     .  qaa,Rma,Ema 
!        double precision Rm(natot),Em(natot),
!     .  pc(0:nac),rclas(3,natot)
!cagregue 0 apc
!        character ch*1,exp
! 	character*4 atom
!	character*4, dimension(:), allocatable, save::
!     .  resname
!        character*4, dimension(:,:), allocatable, save::
!     .  atnamea,attypea,atsp
!	integer, dimension(:,:), allocatable, save::
!     .	atnu
!        integer, dimension(:,:), allocatable, save::      
!     .  nataa,con
!	integer, dimension(:), allocatable, save::
!     .  atxres
!        integer atange(nac,25,2),atangm(nac,25,2),
!     .  atdihe(nac,100,3),atdihm(nac,100,3)
!        integer bondxat(nac),angexat(nac),
!     .  dihexat(nac),dihmxat(nac),angmxat(nac)
!        integer  impxat(nac),atimp(nac,25,4)
!        character*4 aanamea(200)
!        integer atomsxaa(200)
!        integer aanum(nac),resnum(nac)
!        double precision rcorteqm,rcortemm,sfc,rcorteqmmm
!        character*4 atname(nac),aaname(nac),attype(nac),qmattype(na_u)
!        character*5 bondtype(nparm)
!        character*8 angletype(nparm)
!        character*11 dihetype(nparm),imptype(nparm)
!        double precision kbond(nparm),bondeq(nparm),kangle(nparm),
!     .  angleeq(nparm),kdihe(nparm),diheeq(nparm),perdihe(nparm),
!     .  kimp(nparm),impeq(nparm),perimp(nparm)
!        integer multidihe(nparm),multiimp(nparm)
!        double precision radbloqmmm
!	logical foundamber
!        character ch1*1,ch4*4
!ccorrecion del bug de bonds extras, Nick
!	integer :: ivalue(natot*2), inick,jnick
!c parche para q no ponga bonds entre extremos terminales
!        double precision radblommbond, radinnerbloqmmm
!        integer res_ref
!	ivalue=0


! Initialice variables
	rclas = 0.d0
	atname="NICK"
	aaname="NICK"
	resnum=-1
	aanum =-1
	nroaa =-1
	res_ref=1
	rcorteqm = 1.d-06
	rcortemm = 100.d0
	rcorteqmmm=0.d0
	radbloqmmm=0.d0
	radblommbond=999999999999.99999
	radinnerbloqmmm=0.d0
	res_ref=1
	ncon=0
	con2=-1
	foundamber=.false.
! nullify some solvent vbles
!      ng1 = 0
!      nangle = 0
!      ndihe = 0
!      nimp = 0
!      sfc=2.d0

! Read MM system from .fdf file
	if ( fdf_block('SolventInput',iunit) ) then !asing unit number to iunit for read SolventInput block
	  do i=1,nac ! barre todos los atomos clasicos del sistema
	    read(iunit,err=11,end=11,fmt='(A4,I7,2x,A4,A4,A,I4,4x,3f8.3)') &
	    atom, j, atname(i), aaname(i), ch, resnum(i), rclas(1:3,na_u+i)
! read: Char*4, atom number, atom name, residue name, char*1, residue number, x,y,z positions
! ATOM     71  H   ILE H   5      21.490  24.570  40.320
	    ivalue(j)=i !relation between atom number and atom position in .fdf
	 enddo
	else
	  call die("solvent: You must specify the solvent coordinates  &
	  if NumberOfSolventAtoms not equal 0")      
	endif

! Change coordinates to atomic units
	rclas(1:3,1:natot) = rclas(1:3,1:natot) * Ang

! Assigns number of residues 
	k=1
	aanum(1)=k
	do i=2,nac !all MM atoms
	  if (resnum(i).eq.resnum(i-1)) then !check is residue number if same as atom i-1
	    aanum(i)=aanum(i-1)
	  elseif (resnum(i).ne.resnum(i-1)) then
	    k = k+1
	    aanum(i)= k
	  endif
	enddo

	nroaa = aanum(nac)
	if(nroaa.eq.0) call die("solvent: Number of residues can not be zero")


! Read QM atom type for QM/MM
	if(na_u.ne.0) then
	  if ( fdf_block('SoluteAtomTypes',iunit) ) then

	    do i=1,na_u+1
	      read(iunit,*,end=21,err=21) ch4
	        ch1=ch4(1:1)

	        if(i.le.na_u) then
		  if(ch1.eq.'%') &
		    call die('Atom types in SoluteAtomTypes are lower than na_u')
		else
		  if(ch1.eq.'%') &
		    call die('Atom types in SoluteAtomTypes are greater than na_u')
		end if

	        qmattype(i)=ch4 !type of QM atom in amber.parm
  	    enddo
	  else
	    call die('SoluteAtomTypes Block if not defined')
	  endif
	endif


! Cutoff & Freezing Radious
	if ( fdf_block('CutOffRadius',iunit) ) then
	  read(iunit,*,err=31,end=31) exp, rcorteqm ! not used
	  read(iunit,*,err=31,end=31) exp, rcorteqmmm ! distance for QM-MM interaction
	  read(iunit,*,err=31,end=31) exp, rcortemm ! distance for LJ & Coulomb MM interaction
	  read(iunit,*,err=31,end=31) exp, radbloqmmm ! distance that allow to move MM atoms from QM sub-system
	  read(iunit,*,err=31,end=31) exp, radblommbond !parche para omitir bonds en extremos terminales, no se computan bonds con distancias mayores a radblommbond
	  ! remover al incluir un bloque que ajuste conectividades
	  read(iunit,*,err=50,end=50) exp, radinnerbloqmmm !distance that not allow to move MM atoms from QM sub-system
	  read(iunit,*,err=50,end=50) exp, res_ref ! residue that is taken as reference to fix atoms by radial criteria in full MM simulations
	else
	  write(6,'(/a)') 'WARNING: CutOffRadius not defined. Using defaults'
	endif
  50    continue

	if(rcorteqmmm.le.1.e-8) then
	  call die('WARNING: QM-MM cut-off radius to close to zero')
	endif


! Aditional MM Connectivity block
	if ( fdf_block('SolventConnectivity',iunit) ) then
	i=1
	exp= 'N'
	do while (exp.ne.'%' .or. i.eq.1000)
	  read(iunit,'(a)',advance='no',err=41,end=41) exp
	  if(exp.ne.'%') then
	    read(iunit,*,err=41,end=41) exp, con2(1:2,i)
	    i=i+1
	  end if
	end do
	if(i.eq.1000) &
	  call die('SolventConnectivity block must not exeed 1000 lines')
	  if(ncon.ne.0) write(6,'(/,a)') 'Reading new connectivities block'
	  ncon=i-1
	  allocate(con(2,ncon)) !arreglo esto para conectividades correctas segun el numero de atomo en el fdf, Nick

	  do j=1, ncon
	    con(1:2,1:ncon)=con2(1:2,1:ncon)
	    write(*,*) "extra interaction added ",con(1,j)," & ",con(2,j)
          end do
	else
	  allocate(con(2,1)) !agregado nick, sino trae problemas cuando ncon=0
	endif


! Checkin MM parameters file "amber.parm"
	inquire( file="amber.parm", exist=foundamber )
	if(.not.foundamber) call die("solvent: 'amber.parm' file not found")

	write(*,*) "FF1"
!Reading MM force Field (bond, angle, dihe e imp)
	call Force_field_read(nbond, nangle)
!	call amber_union_parms(nbond,kbond,bondeq,bondtype,nangle,kangle,  &
!	angleeq,angletype,ndihe,kdihe,diheeq,dihetype,multidihe,perdihe,   &
!	nimp,kimp,impeq,imptype,multiimp,perimp,nparm)
	write(*,*) "FF2"


	Return

! Errors
 11    continue
	  write(*,*) "Error reading MM coordinates in .fdf"
	  if (i.ne.1) then
	    write(*,*) "last correct read:"
	    write(*,*) atom, j, atname(i-1), aaname(i-1), ch, resnum(i-1), rclas(1:3,na_u+i-1)
	  else
	    write(*,*) "First atoms is incorrect"
	  end if
	  stop 'Check .fdf file'
 21    continue
	  write(*,*) "Error reading SoluteAtomTypes block"
	  write(*,*) "Last correct read: ", ch4
	  stop 'Check .fdf file' 
 31    continue
	  stop 'Error reading CutOffRadius block'
 41    continue
	  write(*,*) "Error reading solvent connectivities in .fdf"
	  stop 

	end subroutine MM_atoms_assign





!**********************************************************************
! This subroutine reads bonds, angles, dihedrals and improper torsions 
! parameters of force field

	subroutine Force_field_read(nbond, nangle)
	use scarlett, only: kbond,bondeq,bondtype,angletype,kangle,angleeq

!     . nac, atname, aaname, attype, qmattype, aanum, ng1, bondtype,
!     . kbond,bondeq, bondxat, angletype, kangle,angleeq, angexat,
!     . angmxat,dihetype, kdihe,diheeq, perdihe, multidihe, dihexat,
!     . dihmxat, imptype, kimp,impeq,perimp, multiimp, impxat,
!     . scalexat, scale, nonbondedxat, nonbonded, Em, Rm,
!    . fce_amber, fdummy, cfdummy, ng1type, angetype, angmtype,
!     . dihety,dihmty,impty, evaldihelog, evaldihmlog,
!     . atange, atangm,
!     . atdihe,atdihm,atimp,
!     . rclas, vat, aat, izs, evaldihe,evaldihm,
!     . linkatom, numlink, linkat, linkqm, linkmm, linkmm2, parametro,
!     . linkqmtype, Elink, distl, pclinkmm, Emlink, frstme, pi,

!nbond,kbond,bondeq,bondtype,nangle, &
!	kangle,angleeq,angletype,ndihe,kdihe,diheeq,dihetype,multidihe,  &
!	perdihe,nimp,kimp,impeq,imptype,multiimp,perimp,nparm)
	
	use ionew 
	implicit none
!	double precision, dimension(:), allocatable, intent(out) :: kbond,bondeq
!	character*5, dimension(:), allocatable, intent(out) :: bondtype
	integer :: ui
	logical :: search
	character*12 :: option
	integer :: i
	integer, intent(out) :: nbond,nangle
!	integer nparm  
!        character*5 bondtype(nparm)
!        character*8 angletype(nparm)
!        character*11 dihetype(nparm),imptype(nparm)
!        double precision kbond(nparm),bondeq(nparm),kangle(nparm),
!     .  angleeq(nparm),kdihe(nparm),diheeq(nparm),perdihe(nparm),
!     .  kimp(nparm),impeq(nparm),perimp(nparm)
!        integer nbond,nangle,ndihe,nimp,i,j,k,multidihe(nparm),
!     .          multiimp(nparm)

! Reading Bonds
	nbond=0
	call io_assign(ui)
	open(unit=ui,file="amber.parm")
	search=.true.
	do while (search)
	  read (ui,*,err=100,end=100) option
	  if (option.eq.'bonds') then
	   read(ui,*,err=100,end=100) nbond
	   allocate(bondtype(nbond),kbond(nbond),bondeq(nbond))
	   do i=1,nbond
	     read(ui,10,err=100,end=100) bondtype(i),kbond(i),bondeq(i)
	   enddo
	   search=.false.
	 endif
	enddo
	call io_close(ui)

	nangle=0
! Reading Angles
	call io_assign(ui)
	open(unit=ui,file="amber.parm")
	search=.true.
	do while (search)
	read (ui,*,err=200,end=200) option
	if (option.eq.'angles') then
	  read(ui,*,err=200,end=200) nangle
	  allocate(angletype(nangle),kangle(nangle),angleeq(nangle))
	  do i=1,nangle
	    read(ui,20,err=200,end=200) angletype(i),kangle(i),angleeq(i)
	  enddo
	  search=.false.
	endif
	enddo
	call io_close(ui)
!
!c diherdrals
!        call io_assign(ui)
!        open(unit=ui,file="amber.parm")
!        search=.true.
!        do while (search)
!        read (ui,*,err=300,end=300) option
!        if (option.eq.'dihes') then
!        	read(ui,*,err=300,end=300) ndihe
!                if(ndihe.gt.nparm) stop 'solvent: Increase nparm'
!        	do i=1,ndihe
!        	read(ui,30,err=300,end=300) dihetype(i),multidihe(i),kdihe(i),
!     .                      diheeq(i),perdihe(i)
!        	enddo
!        search=.false.
!        endif
!        enddo
!        call io_close(ui)
!
!c impropers
!        call io_assign(ui)
!        open(unit=ui,file="amber.parm")
!        search=.true.
!        do while (search)
!        read (ui,*,err=400,end=400) option
!        if (option.eq.'imps') then
!        	read(ui,*,err=400,end=400) nimp
!                if(nimp.gt.nparm) stop 'solvent: Increase nparm'
!        	do i=1,nimp
!        	read(ui,40,err=400,end=400) imptype(i),multiimp(i),kimp(i),
!     .                      impeq(i),perimp(i)
!	        enddo
!        search=.false.
!        endif
!        enddo
!        call io_close(ui)
!
 10     format(A5,2x,F5.1,4x,F6.4)
 20     format(A8,3x,F5.1,6x,F6.2)
! 30     format(A11,3x,I1,3x,F6.3,7x,F7.3,10x,F6.3)
! 40     format(A11,3x,I1,3x,F6.3,7x,F7.3,10x,F6.3)
	return
 100    write(*,*) &
	'solvent: Problem reading bonds block in amber.parm file',i
	stop
 200    write(*,*) &
	'solvent: Problem reading angles block in amber.parm file',i
        stop
! 300    write(*,*) 
!     .  'solvent: Problem reading dihes block in amber.parm file',i
!        stop
! 400    write(*,*) 
!     .  'solvent: Problem reading imps block in amber.parm file',i
!        stop
	end subroutine Force_field_read
!**********************************************************************


