!****************************************************************************
! This routine read and assign parameters, coordinates and conectivity to
! all MM and QM/MM atoms.
! modified from oroginal solv_assign subroutine
! N. Foglia 2019
!*****************************************************************************

	subroutine MM_atoms_assign(nac, na_u, natot, atname, aaname, rclas, &
	nroaa, aanum, qmattype, rcorteqm, rcorteqmmm, rcortemm, radbloqmmm, &
	radblommbond, radinnerbloqmmm, res_ref, nbond, nangle, ndihe, nimp, &
	attype, pc, Rm, Em)


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
	integer, intent(out) :: nbond, nangle, ndihe, nimp !number of bonds, angles, sihesral and impropers in amber.parm file
	integer :: FF_residues, FF_max_at_by_res !number of residues and max number of atoms by residue in force field


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
	double precision, dimension(:,:), allocatable, save:: Rma,Ema
!     .  qaa
	double precision, dimension(0:nac), intent(out) :: pc
	double precision, dimension(natot), intent(out) :: Rm, Em
!     .  pc(0:nac),rclas(3,natot)
!cagregue 0 apc
!        character ch*1,exp
! 	character*4 atom
	character*4, dimension(:), allocatable :: resname
	character*4, dimension(:,:), allocatable :: atnamea, attypea
!,attypea,atsp
	integer, dimension(:,:), allocatable :: atnu
	integer, dimension(:,:), allocatable :: nataa
!     .  nataa,con
	integer, dimension(:), allocatable :: atxres
!        integer atange(nac,25,2),atangm(nac,25,2),
!     .  atdihe(nac,100,3),atdihm(nac,100,3)
!        integer bondxat(nac),angexat(nac),
!     .  dihexat(nac),dihmxat(nac),angmxat(nac)
!        integer  impxat(nac),atimp(nac,25,4)
	character*4, dimension(:), allocatable ::  aanamea
	integer, dimension(:), allocatable :: atomsxaa
!        integer aanum(nac),resnum(nac)
!        double precision rcorteqm,rcortemm,sfc,rcorteqmmm
	character*4, dimension(nac), intent(out) ::  attype
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
	nbond = 0
	nangle = 0
	ndihe = 0
	nimp = 0

! nullify some solvent vbles
!      ng1 = 0
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

! Reading MM force Field (bond, angle, dihe e imp)
	call Force_field_read(nbond, nangle, ndihe, nimp)
	write(*,*) "FF2"
! Reading number of residues and atoms in force field
	call FF_atoms_by_residue(FF_residues,FF_max_at_by_res)
	write(*,*) "FF3", FF_residues

	allocate(atnamea(nroaa,FF_max_at_by_res), resname(nroaa), &
	atxres(nroaa), atnu(nroaa,FF_max_at_by_res), &
	attypea(nroaa,FF_max_at_by_res), nataa(nroaa,FF_max_at_by_res))

! Read charge and atom type, and assign to .fdf. atoms in  attype & pc
	call FF_atoms_types(FF_residues,FF_max_at_by_res, nroaa, atnamea,&
	aaname,atxres,resname, nac, atname, atnu, attypea, nataa, attype, pc)
	write(*,*) "FF4"

! Read Lenard-Jones parameters and asiign to .fdf atoms using attypea
	allocate(Ema(nroaa,FF_max_at_by_res),Rma(nroaa,FF_max_at_by_res))
	call FF_lj(nroaa,FF_max_at_by_res,attypea,Ema,Rma,na_u,natot, &
	Em,Rm,qmattype,atxres)

	write(*,*) "FF4"

!aca empieza la verdadera asignacion: segun atomo xa c/aa

!	deallocate()
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

!****************************************************
! This subroutine read Lenard-Jones parameters from amber.parm and asign it
! to .fdf atoms using attypea

	subroutine FF_lj(nroaa,FF_max_at_by_res,attypea,Ema,Rma,na_u,natot, &
	Em,Rm,qmattype,atxres)
	use ionew, only : io_assign, io_close
	use scarlett, only: Ang, kcal, eV
	implicit none

	integer, intent(in) :: nroaa, FF_max_at_by_res,na_u,natot
	character*4, dimension(nroaa,FF_max_at_by_res), intent(in) :: attypea 
	double precision, dimension(nroaa,FF_max_at_by_res), intent(out) :: Ema, Rma
	double precision, dimension(natot), intent(out) :: Rm, Em
	character*4, dimension(na_u), intent(in) :: qmattype(na_u)
	integer, dimension(nroaa), intent(in) :: atxres
	integer :: ui
	logical :: search
	character*12 :: option
	double precision, dimension(:), allocatable :: pRm, pEm
	character*4, dimension(:), allocatable ::  ljtype
	logical :: asigned
	integer :: nlj
!auxiliars
	integer :: i,j,k
	Rm=0.d0
	Em=0.d0

	call io_assign(ui)
	open(unit=ui,file="amber.parm")
 
! Read amber.parm
	search=.true.
	do while (search)
	  read (ui,*,err=1,end=1) option
	  if (option.eq.'ljs') then
	    read(ui,*,err=1,end=1) nlj !cantidad de lennard jones en el amber.parm
	    allocate(ljtype(nlj),pRm(nlj),pEm(nlj))
	    do  i=1,nlj
	      read (ui,*,err=1,end=1) ljtype(i),pRm(i),pEm(i)
	    enddo
	    search=.false.
	  endif
	enddo
	call io_close(ui)

! transform unit to hybrid 
	do i=1,nlj
	  pRm(i) = (2.0d0*pRm(i)*Ang)/(2.d0**(1.d0/6.d0))
	  pEm(i) = (pEm(i)*eV/kcal) !(pEm(i)/627.5108d0)
	enddo


! asignation of LJ parameters to MM atoms using attypea 
	do i=1,nroaa
	  do j=1,atxres(i)
	    asigned=.false.

 	    do k=1,nlj
	      if (attypea(i,j).eq.ljtype(k)) then
	        Rma(i,j) = pRm(k)
	        Ema(i,j) = pEm(k)
	        asigned=.true.
	      endif
	    enddo

	    if (.not. asigned) then
	      write(*,*) "WARNING, LJ parameter not found for", attypea(i,j)
	      write(*,*) "checkin for similar atom"

	      if (attypea(i,j).eq.'C'.or. &
	      attypea(i,j).eq.'CA'.or.attypea(i,j).eq.'CM'.or. &
	      attypea(i,j).eq.'CC'.or.attypea(i,j).eq.'CV'.or. &
	      attypea(i,j).eq.'CW'.or.attypea(i,j).eq.'CR'.or. &
	      attypea(i,j).eq.'CB'.or.attypea(i,j).eq.'C*'.or. &
	      attypea(i,j).eq.'CN'.or.attypea(i,j).eq.'CK'.or. &
	      attypea(i,j).eq.'CQ'.or.attypea(i,j).eq.'CX'.or. &
	      attypea(i,j).eq.'CY'.or.attypea(i,j).eq.'CD')  then
 
	        do k=1,nlj
	          if (ljtype(k).eq.'C') then
	            Rma(i,j) = pRm(k)
	            Ema(i,j) = pEm(k)
		    write(*,*) "Using :", ljtype(k)
	          endif
	        enddo
 
	      elseif (attypea(i,j).eq.'N'.or.attypea(i,j).eq.'NA'.or. &
	      attypea(i,j).eq.'NB'.or.attypea(i,j).eq.'NC'.or. &
	      attypea(i,j).eq.'N*'.or.attypea(i,j).eq.'N2'.or. &
	      attypea(i,j).eq.'NO'.or.attypea(i,j).eq.'NP') then
 
	        do k=1,nlj
	          if (ljtype(k).eq.'N') then
	            Rma(i,j) = pRm(k)
	            Ema(i,j) = pEm(k)
	            write(*,*) "Using :", ljtype(k)
	          endif
	        enddo
	
	      else
	        STOP "didn found any LJ parameter"
	      end if
	    endif
	  enddo
	enddo

! asignation of LJ parameters to QM atoms using qmttype 
	do i=1,na_u
	  asigned=.false.

	  do k=1,nlj
	    if (qmattype(i).eq.ljtype(k)) then
	      Rm(i) = pRm(k)
	      Em(i) = pEm(k)
	      asigned=.true.
	    endif
	   end do
	    
	    if(.not. asigned) then
	      write(*,*) "WARNING, LJ parameter not found for ", qmattype(i)
	      write(*,*) "checkin for similar atom"
	      
	      if (qmattype(i).eq.'C'.or. &
	      qmattype(i).eq.'CA'.or.qmattype(i).eq.'CM'.or. &
	      qmattype(i).eq.'CC'.or.qmattype(i).eq.'CV'.or. &
	      qmattype(i).eq.'CW'.or.qmattype(i).eq.'CR'.or. &
	      qmattype(i).eq.'CB'.or.qmattype(i).eq.'C*'.or. &
	      qmattype(i).eq.'CN'.or.qmattype(i).eq.'CK'.or. &
	      qmattype(i).eq.'CQ'.or.qmattype(i).eq.'CX'.or. &
	      qmattype(i).eq.'CY'.or.qmattype(i).eq.'CD') then
 	
	        do k=1,nlj
	          if (ljtype(k).eq.'C') then
	            Rm(i) = pRm(k)
	            Em(i) = pEm(k)
	            write(*,*) "Using :", ljtype(k)
	          endif
	        enddo

	      elseif (qmattype(i).eq.'N' .or.qmattype(i).eq.'NA'.or. &
	      qmattype(i).eq.'NB'.or.qmattype(i).eq.'NC'.or. &
	      qmattype(i).eq.'N*'.or.qmattype(i).eq.'N2'.or. &
	      qmattype(i).eq.'NO'.or.qmattype(i).eq.'NP') then
 
	        do k=1,nlj
	          if (ljtype(k).eq.'N') then
	            Rm(i) = pRm(k)
	            Em(i) = pEm(k)
	            write(*,*) "Using :", ljtype(k)
	          endif
	        enddo
 	      else
	        STOP "didn found any LJ parameter"
	      endif
	    end if
	enddo
 

!Move LJ parameters to Em and Rm
	k=na_u+1
	do i=1,nroaa
	  do j=1,atxres(i)
	    Em(k)=Ema(i,j)
	    Rm(k)=Rma(i,j)     
	    k=k+1
	  enddo
	enddo

	return
 1      write(*,*) 'Problem reading LJ block in amber.parm file'
	if (i.gt.1) then
	  write(*,*) "last correct read",  ljtype(i-1),pRm(i-1),pEm(i-1)
	else
	  write(*,*) "1st LJ is incorrect"
	end if
	stop
	end subroutine FF_lj



!*************************************************************
! This subroutine Reads atom type and charges form amber.parm and
! Assing thos parameters to each atom on .fdf in  attype & pc
	subroutine FF_atoms_types(FF_residues,FF_max_at_by_res, nroaa, atnamea,&
	aaname,atxres,resname, nac, atname, atnu, attypea, nataa, attype, pc)
	use ionew, only:io_assign, io_close
	implicit none
	integer :: trash
	integer, intent(in) :: FF_residues,FF_max_at_by_res, nroaa
	integer, intent(in) :: nac
	character*4, dimension(nac), intent(in) :: atname
	integer, dimension(nroaa,FF_max_at_by_res), intent(out) :: nataa
	character*4, dimension(nac), intent(out) ::  attype
	character*4, dimension(nroaa,FF_max_at_by_res), intent(inout) :: atnamea
	character*4, dimension(nroaa), intent(out) :: resname
	character*4, dimension(nac), intent(in) :: aaname 
	integer, dimension(nroaa,FF_max_at_by_res), intent(out) :: atnu
	character*4, dimension(nroaa,FF_max_at_by_res), intent(out) ::  attypea
	double precision, dimension(0:nac), intent(out) :: pc
	double precision, dimension(:,:), allocatable :: qaa
	double precision, dimension(:,:), allocatable :: pqaa
	character*4, dimension(:), allocatable :: paanamea
	character*4, dimension(:,:), allocatable :: patnamea,pattype
	integer, dimension(:,:), allocatable :: pnataa,patmas !atnu
	integer, dimension(:), allocatable :: patxres
	logical :: search
	character*12 :: option 
	integer :: ui
	integer, dimension(nroaa), intent(out) :: atxres
	logical, dimension(nroaa, FF_max_at_by_res) :: assigned
	logical :: kill
!auxiliars
	integer :: n1, n2, n3
	integer :: i,j,k, l,m
	integer :: i2,j2,k2

	kill=.false.

	call io_assign(ui) 
	open(unit=ui,file="amber.parm")

	allocate(paanamea(FF_residues), patxres(FF_residues), &
	patnamea(FF_residues,FF_max_at_by_res), &
	pattype(FF_residues,FF_max_at_by_res), &
	pnataa(FF_residues,FF_max_at_by_res), &
	patmas(FF_residues,FF_max_at_by_res), &
	pqaa(FF_residues,FF_max_at_by_res))

	allocate(qaa(nroaa, FF_max_at_by_res))

	search=.true.
	do while (search)
	  read (ui,*,err=1,end=1) option
	  if (option.eq.'residues') then
	    read(ui,*,err=1,end=1)  trash !number of residues in amber.parm

	    do i=1,FF_residues
	      read(ui,*,err=1,end=1) paanamea(i), patxres(i) !residue name and number of atoms in it
	        do j=1,patxres(i)
	          read(ui,*,err=1,end=1) patnamea(i,j),pattype(i,j), & !atom name, atom type
	          n1,n2,n3,pnataa(i,j),patmas(i,j),pqaa(i,j) ! x, x, x, number in residue, atomic mass, charge 
	        enddo
	    enddo
	    search=.false.
	  endif
	enddo
	call io_close(ui)


! Asign atom number for each residue in .fdf
	atxres=0
	atnu=-1
	k = 1
	do i=1,nroaa !Number of residue for clculation
	  resname(i)=aaname(k) !resname =  name of residues, aaname name of residue of k-th line in .fdf
	  do l=1,FF_residues !all residues in amber.parm
	    if (resname(i).eq.paanamea(l)) then !compare name with amber.parm
	      atxres(i)=patxres(l) !asing number of atoms to residue number in calcultion
	    endif
	  enddo

	  if(atxres(i).eq.0) then
	    write(6,*) 'Wrong residue name in.fdf.:', aaname(k)
	    STOP
	  endif

	  do j=1,atxres(i)
	    atnu(i,j)=k !atom number atom j in residue i
	    atnamea(i,j)=atname(k) !atom name
	    if (aaname(k) .ne. resname(i)) write(*,*) & !check residue name for each atom
	     "WARNING, Missing atoms in ", resname(i), i 
	    k = k+1
	  enddo
	enddo



! Assign charge and atom type to fdf atoms
	qaa=-999999999.d0
	assigned=.false.
	do i=1,nroaa !residues in .fdf
	  do k=1,FF_residues !residues in amber.parm
	    if(resname(i).eq.paanamea(k)) then !compare name of residues
	      do j=1,patxres(k) !atoms in .fdf
	      do m=1,patxres(k) !atoms in amber.parm
	        if (atnamea(i,j).eq.patnamea(k,m)) then !compare atom names
	          qaa(i,j) = pqaa(k,m) !asign charge
	          attypea(i,j) = pattype(k,m) !assing type
	          nataa(i,j) = pnataa(k,m) !assign atom number
	          assigned(i,j) = .true.
	        endif
	      enddo
	      enddo
	    endif
	  enddo
	enddo

	k=1
	do i=1,nroaa
	  do j=1,atxres(i)
	    attype(k)=attypea(i,j)
	    k=k+1
	   enddo
	enddo
 
	pc=0
	k=1
	do i=1,nroaa
	  do j=1,atxres(i)
	    pc(k)=qaa(i,j) !asing charge to each atom
	    if(.not.  assigned(i,j)) then !check asignation
	      write(6,'(a,i5)') 'Wrong atom name in .fdf :', k           
	      k2=k
	      if (i.gt.1) then
  		i2=i-1
		do j2=1,atxres(i2)
		  k=k-1
		end do
		i2=i
		do j2=1,j-1
		  k=k-1
		end do

		do i2=i-1, i
		  do j2=1,atxres(i2)
		    if (i .eq. i2 .and. j2.eq.j) then
		      write(6,*) k, resname(i2),atnamea(i2,j2),"<--- this one"
! attypea(i2,j2),"<--- this one"
		    else
		      write(6,*) k, resname(i2),atnamea(i2,j2)
!attypea(i2,j2)
		    end if
		    k=k+1
		  end do
		end do

	      else
		do j2=1,j-1
		  k=k-1
		end do

		i2=i
		do j2=1,atxres(i2)
		  if (j2.eq.j) then
		    write(6,*) k, resname(i2), attypea(i2,j2),"<--- this one"
		  else
		    write(6,*) k, resname(i2), attypea(i2,j2)
		  end if
		  k=k+1
		end do
	      end if
	      kill=.true.
	      k=k2
	    endif
	    k=k+1
	  enddo
	enddo


	deallocate(paanamea, patxres, patnamea, pattype, pnataa, patmas, pqaa,&
	qaa)
	if (kill) STOP

        return
 1      write(*,*) "Problem reading residues block in amber.parm file"
	if (i.gt.1) then
	  write(*,*) "Last correct residue read"
	  write(*,*)  paanamea(i-1), patxres(i) 
	else
	  write(*,*) "First residue with problems"
	end if
	stop
	end subroutine FF_atoms_types


!************************************************************************
!  This subroutine reads number of atoms by residue in Force Field 
	subroutine FF_atoms_by_residue(FF_residues,max_at_by_res)
	use ionew, only:io_assign, io_close
	implicit none

	integer, intent(out) :: FF_residues !number of residues in force field
	integer, intent(out) :: max_at_by_res ! max number of atoms in a residue in force field
	integer :: ui
	integer :: at_by_res
	character*4 :: nada
	character*10 :: option
	logical :: search
	integer :: i,j

	FF_residues=0
	max_at_by_res=0
	call io_assign(ui)
	open(unit=ui,file="amber.parm")
	search=.true.
	do while (search)
	  read (ui,*,err=1,end=1) option
	  if (option.eq.'residues') then
	    read(ui,*,err=1,end=1) FF_residues
	    do i=1,FF_residues
	      read(ui,*,err=1,end=1) nada,at_by_res
	      if (at_by_res .gt. max_at_by_res) max_at_by_res=at_by_res
	      do j=1,at_by_res
	        read(ui,*,err=1,end=1)
	      end do
	    end do
	    search=.false.
	  end if
	end do
	call io_close(ui)

	return
 1      write(*,*) 'Problem reading residues block in amber.parm file'
        stop
	end subroutine FF_atoms_by_residue

!**********************************************************************
! This subroutine reads bonds, angles, dihedrals and improper torsions 
! parameters of force field

	subroutine Force_field_read(nbond, nangle, ndihe, nimp)
	use scarlett, only: kbond,bondeq,bondtype,angletype,kangle,angleeq, &
	dihetype,multidihe,kdihe,diheeq,perdihe,imptype,multiimp,kimp,impeq, &
	perimp
	use ionew, only : io_assign, io_close
	implicit none
	integer :: ui
	logical :: search
	character*12 :: option
	integer :: i
	integer, intent(out) :: nbond,nangle, ndihe, nimp

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
	   bondtype=""
	   kbond=0.d0
	   bondeq=0.d0
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
	    angletype=""
	    kangle=0.d0
	    angleeq=0.d0
	    do i=1,nangle
	      read(ui,20,err=200,end=200) angletype(i),kangle(i),angleeq(i)
	    enddo
	    search=.false.
	  endif
	enddo
	call io_close(ui)
!
! Reading Diherdrals
	call io_assign(ui)
	open(unit=ui,file="amber.parm")
	search=.true.
	do while (search)
	  read (ui,*,err=300,end=300) option
	  if (option.eq.'dihes') then
	    read(ui,*,err=300,end=300) ndihe
	    allocate(dihetype(ndihe),multidihe(ndihe),kdihe(ndihe),diheeq(ndihe)&
	    ,perdihe(ndihe))
	    dihetype=""
	    multidihe=0
	    kdihe=0.d0
	    diheeq=0.d0
	    perdihe=0.d0
	    do i=1,ndihe
	      read(ui,30,err=300,end=300) dihetype(i),multidihe(i),kdihe(i), &
	      diheeq(i),perdihe(i)
	    enddo
	    search=.false.
	  endif
	enddo
	call io_close(ui)
!
! Reading Impropers
	call io_assign(ui)
	open(unit=ui,file="amber.parm")
	search=.true.
	do while (search)
	  read (ui,*,err=400,end=400) option
	  if (option.eq.'imps') then
	    read(ui,*,err=400,end=400) nimp
	    allocate(imptype(nimp),multiimp(nimp),kimp(nimp),impeq(nimp), &
	    perimp(nimp))
	    imptype=""
	    multiimp=0
	    kimp=0.d0
	    impeq=0.d0
	    perimp=0.d0
	    do i=1,nimp
	      read(ui,40,err=400,end=400) imptype(i),multiimp(i),kimp(i), &
	      impeq(i),perimp(i)
	    enddo
	    search=.false.
	  endif
	enddo
	call io_close(ui)
!
 10     format(A5,2x,F5.1,4x,F6.4)
 20     format(A8,3x,F5.1,6x,F6.2)
 30     format(A11,3x,I1,3x,F6.3,7x,F7.3,10x,F6.3)
 40     format(A11,3x,I1,3x,F6.3,7x,F7.3,10x,F6.3)
	return

 100    write(*,*) 'Problem reading bonds block in amber.parm file'
	if (i.ge.2) then
	  write(*,*) 'Last correct read'
	  write(*,10) bondtype(i-1),kbond(i-1),bondeq(i-1)
	else
	  write(*,*) 'Problems on 1st bond '
	end if
	stop

 200    write(*,*) 'Problem reading angles block in amber.parm file'
	if (i.ge.2) then
	  write(*,*) 'Last correct read'
	  write(*,20) angletype(i-1),kangle(i-1),angleeq(i-1)
	else
	  write(*,*) 'Problems on 1st angle '
	end if
	stop

 300    write(*,*) 'Problem reading dihes block in amber.parm file'
	if (i.ge.2) then
	  write(*,*) 'Last correct read'
	  write(*,30) dihetype(i-1),multidihe(i-1),kdihe(i-1), &
	  diheeq(i-1),perdihe(i-1)
	else
	  write(*,*) 'Problems on 1st dihedral '
	end if
	stop

 400    write(*,*) 'Problem reading imps block in amber.parm file'
	if (i.ge.2) then
	  write(*,*) 'Last correct read'
	  write(*,40) imptype(i-1),multiimp(i-1),kimp(i-1), &
	  impeq(i-1),perimp(i-1)
	else
	  write(*,*) 'Problems on 1st dihedral '
	end if

	stop
	end subroutine Force_field_read
!**********************************************************************


