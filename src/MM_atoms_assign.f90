!****************************************************************************
! This routine read and assign parameters, coordinates and conectivity to
! all MM and QM/MM atoms.
! modified from oroginal solv_assign subroutine
! N. Foglia 2019
!*****************************************************************************

	subroutine MM_atoms_assign(nac, na_u, natot, atname, aaname, rclas, &
	nroaa, aanum, qmattype, rcorteqm, rcorteqmmm, rcortemm, radbloqmmm, &
	radblommbond, radinnerbloqmmm, res_ref, nbond, nangle, ndihe, nimp, &
	attype, pc, Rm, Em, ng1, bondxat,angexat, angmxat, dihexat, dihmxat,&
	impxat)


!        (attype,pc,
!     .  atange,atangm,atdihe,
!     .  atdihm,atimp,
!     .  kbond,bondeq,bondtype,
!     .  kangle,angleeq,angletype,
!     .  kdihe,diheeq,dihetype,multidihe,perdihe,
!     .  kimp,impeq,imptype,multiimp,perimp,
!     .  nparm,aaname,atname,aanum,qmattype,rclas,
!     .  rcorteqmmm,rcorteqm,rcortemm,sfc,
!     .  radbloqmmm,atsinres,radblommbond,radinnerbloqmmm,res_ref)


	use precision, only : dp
	use sys,only : die
	use fdf,only : fdf_block
	use scarlett, only: Ang, atxres

	implicit none
	integer, intent(in) :: nac, na_u, natot !number of QM and MM and total atoms
	character*4, dimension(nac), intent(out) :: atname, aaname !atom and residue name in .fdf
	integer, intent(out), dimension(nac) :: bondxat, angexat,angmxat, &
	dihexat,dihmxat, impxat !number of bonds, angle(extreme), angle(middle), dihedral(extreme), dihedral(middle) and impropers for each atom
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
	integer :: number_conect_ommit
	integer, dimension(:,:), allocatable :: conect_ommit
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
	integer, intent(out) :: ng1(nac,6)
!,con2(2,1000)
	double precision, dimension(:,:), allocatable :: Rma,Ema
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
!	integer, dimension(:), allocatable :: atxres
!        integer atange(nac,25,2),atangm(nac,25,2),
!     .  atdihe(nac,100,3),atdihm(nac,100,3)
!        integer  atimp(nac,25,4)
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
	ng1 = 0
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

! Read Lenard-Jones parameters and asign to .fdf atoms using attypea
	allocate(Ema(nroaa,FF_max_at_by_res),Rma(nroaa,FF_max_at_by_res))
	call FF_lj(nroaa,FF_max_at_by_res,attypea,Ema,Rma,na_u,natot, &
	Em,Rm,qmattype,atxres)

	write(*,*) "FF5"

	number_conect_ommit=0
	if ( fdf_block('SolventOmmit',iunit) ) then
	  read(iunit,*) number_conect_ommit
	  allocate(conect_ommit(number_conect_ommit,2))
	  conect_ommit=-2
	  do i=1, number_conect_ommit
	    read(iunit,*,err=51,end=51) conect_ommit(i,1), conect_ommit(i,2)
	  end do
	else
	 allocate(conect_ommit(1,2))
	 conect_ommit=-2
	end if

	write(*,*) "FF6"

!Asign conectivity
	call FF_conectivity(nac,nataa,nroaa,atxres,atnu,resname, ng1, &
	number_conect_ommit, conect_ommit, atnamea,ncon,con, FF_max_at_by_res)

	write(*,*) "FF7"

! calculate bonds, angles, dihed & improp
	call FF_bon_ang_dih_imp(nac,nroaa,FF_max_at_by_res, ng1, &
	bondxat, angexat,angmxat,dihexat,dihmxat, impxat, atnamea, atxres, &
	resname, atnu)


!change WAT name to HOH, need to remove this in future. Nick
	do i=1,nac
	  if(aaname(i).eq.'WAT') aaname(i)='HOH'
	enddo


!check correct order for atoms in water residues
	do i=1,nroaa
	  if(resname(i).eq.'HOH') then
	    do j=1,atxres(i)
	      if(atnamea(i,j).eq.'O') then
	        if(j.ne.1) then
	          write(6,*)  'solvent: Wrong order in water residue :',i
	          stop
	        endif
	      endif
	    enddo
	  endif
	enddo

	deallocate(atnamea,atnu,resname) !,atxres)
	deallocate(attypea,nataa)
	deallocate(Ema,Rma)


!c checking ST and SV parameters
	do i=1,na_u
!	  if(qmattype(i).ne.'HO'.and.qmattype(i).ne.'HW') then
	  if(Rm(i).lt.0.or.Em(i).lt.0) then
	    write(6,'(a,i6)') 'Wrong QM LJ parameter, atom:', i
	    STOP
	  endif
!	  endif
	enddo

	do i=1,nac
!a          if(attype(i).ne.'HO'.and.attype(i).ne.'HW') then
	  if(Rm(i+na_u).lt.0.or.Em(i+na_u).lt.0) then
	    write(6,'(a,i6)') 'Wrong MM LJ parameter, atom:', i
	    write(*,*) "Rm(i+na_u)", Rm(i+na_u)
	    write(*,*) "Em(i+na_u)", Em(i+na_u)
	    STOP
	  end if
	  if (pc(i).gt.900000000.d0) then
	    write(6,'(a,i6)') 'Wrong MM charge, atom:', i
	    write(*,*) "pc(i)", pc(i)
	    STOP
          endif
!          endif
	enddo   

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
 51    continue
	  write(*,*) "Error reading SolventOmmit block in .fdf"
          stop
	end subroutine MM_atoms_assign



!--------------------------------------------------

	subroutine FF_bon_ang_dih_imp(nac,nroaa,FF_max_at_by_res, ng1, &
	bondxat, angexat,angmxat,dihexat,dihmxat, impxat, atnamea, atxres, &
	resname, atnu)

!aFF_bon_ang_dih_imp(nac,nroaa,FF_max_at_by_res, ng1, &
!	bondxat, angexat,angmxat,dihexat,dihmxat, atnamea, atxres,     &
!	resname, impxat)

	use ionew, only: io_assign, io_close
	use scarlett, only: atange, atangm, atdihe,atdihm, atimp
	implicit none
	integer, intent(in) :: nac
	integer, intent(in) :: nroaa !number of residues in .fdf
	integer, intent(in) :: FF_max_at_by_res
	integer, intent(in) :: ng1(nac,6)
	integer, intent(out), dimension(nac) :: bondxat, angexat,angmxat, &
	dihexat,dihmxat, impxat
	character*4, intent(in), dimension(nroaa,FF_max_at_by_res) :: atnamea
	integer, intent(in), dimension(nroaa) ::  atxres
	integer, intent(in), dimension(nroaa,FF_max_at_by_res) :: atnu
	character*4, intent(in), dimension(nroaa) :: resname

!c       parametros asoc a la asignacion de angulos y dihedros
!        integer k,l,m,n,t,t2,
!c       parametros asoc a la asignacion de impropers
!     .  atimp(nac,25,4),nataa(nroaa,100), size
	integer :: imptot !number of impropers in .fdf
	integer :: nresid !number of improper angles in amber.parm
	character*4, dimension(:,:,:), allocatable :: impatnamea
	character*4, dimension(:), allocatable :: presname !Residue name
	integer, dimension(:), allocatable :: pimpxres !number of impropers for each residue
	integer :: max_angle_ex, max_angle_mid
	integer :: max_dihe_ex, max_dihe_mid
	integer :: max_improp, max_improp_at
	integer, dimension(:,:), allocatable :: impnum !impnum(i,j) = atom number for j-th atom in i-th improper
	character*10 :: option
	logical :: search, assignimp
	integer :: ui
	integer :: i,j, k, m, n, t, t2
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%% improper angles
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	max_improp=0
! reed improper angles & allocate variables
	call io_assign(ui)
	open(unit=ui,file="amber.parm")
	search=.true.
	do while (search)
	  read (ui,*,err=1,end=1) option
	  if (option.eq.'impropers') then
	    read(ui,*,err=2,end=2) nresid
	    allocate(presname(nresid),pimpxres(nresid))
	    do i=1,nresid
	      read(ui,*,err=3,end=3) presname(i),pimpxres(i)
	      do j=1,pimpxres(i)
	        read(ui,*,err=4,end=4)
	      end do
	      if (pimpxres(i) .gt. max_improp) max_improp=pimpxres(i)
	write(*,*) "flag 1.12"
	    enddo
	write(*,*) "flag 1.13"
	    search=.false.
	write(*,*) "flag 1.14"
	  endif
	write(*,*) "flag 1.15"
	enddo
	write(*,*) "flag 1.16"
	call io_close(ui)
	write(*,*) "flag 2"
	allocate(impatnamea(nresid,max_improp,4))

	call io_assign(ui)
	open(unit=ui,file="amber.parm")
	search=.true.
	do while (search)
	  read (ui,*,err=1,end=1) option
	  if (option.eq.'impropers') then
	    read(ui,*,err=2,end=2) nresid
	    do i=1,nresid
	      read(ui,*,err=3,end=3) presname(i),pimpxres(i)
	      do j=1,pimpxres(i)
	        read(ui,*,err=4,end=4) impatnamea(i,j,1),impatnamea(i,j,2), &
	                               impatnamea(i,j,3),impatnamea(i,j,4)
	      enddo
	    enddo
	    search=.false.
	  endif
	enddo
	call io_close(ui)
	write(*,*) "flag 3"


! asignacion segun el atomo a partir del aa(2)
! asignacion del numero de impropios imxpat
	imptot=1
	do i=1,nresid !number of improper angles in amber.parm
	do m=1,nroaa !number of residues
	  if(presname(i).eq.resname(m)) then !check residue name = amber.parm name
	    do j=1,pimpxres(i) !impropers resid i in amber.parm
	      imptot=imptot+1 !count number of impropers
	    enddo
	  endif
	enddo
	enddo

	allocate(impnum(imptot,4))
	impnum=0

	imptot=1
	do i=1,nresid !number of improper angles in amber.parm
	do m=1,nroaa !number of residues
	  if(presname(i).eq.resname(m)) then !check residue name = amber.parm name
	    do j=1,pimpxres(i) !impropers resid i in amber.parm
	    do k=1,4 !atom name involved in impropers
	      if(impatnamea(i,j,k).eq.'+M'.and. m.ne.nroaa) then !atom name = +M
	        do n=1,atxres(m+1) !atoms in next residue
	          if(atnamea(m+1,n).eq.'N') impnum(imptot,k)=atnu(m+1,n) !atom name = N on next residue => assign atom number
	        enddo
	      elseif(impatnamea(i,j,k).eq.'-M'.and. m.ne.1) then !atom name = -M
	        do n=1,atxres(m-1) !atoms in previous residue
	          if(atnamea(m-1,n).eq.'C') impnum(imptot,k)=atnu(m-1,n) !atom name = C on previous residue => assign atom number
	        enddo
	      else !impropers inside one residue
	        do n=1,atxres(m) !atoms in actual residue
	          if(atnamea(m,n).eq.impatnamea(i,j,k)) impnum(imptot,k)=atnu(m,n) ! => assign atom number
	        enddo
	      endif
	    enddo
	    imptot=imptot+1
	    enddo
	  endif
	enddo
	enddo


	write(*,*) "flag 4"
	do i=1,imptot
	do j=1,4
	 if (impnum(i,j).eq.0) then
	   do k=1,4
	     impnum(i,k)=0
	   enddo
	 endif
	enddo
	enddo
 
	max_improp_at=0
	do i=1,nac !all MM atoms
	  k=0
	  do j=1,imptot !all impropers in system
	    assignimp=(impnum(j,1).eq.i).or.(impnum(j,2).eq.i).or.(impnum(j,3).eq.i).or.(impnum(j,4).eq.i) !i atom have an improper 
	    if (assignimp) k=k+1 !impropers for i-th atom
	  enddo
	  if (k.gt.max_improp_at) max_improp_at=k !max number of impropers for any atom
	enddo

	allocate(atimp(nac,max_improp_at,4))

	do i=1,nac !all MM atoms
	  k=0
	  do j=1,imptot !all impropers in system
	    assignimp=(impnum(j,1).eq.i).or.(impnum(j,2).eq.i).or.(impnum(j,3).eq.i).or.(impnum(j,4).eq.i) !i atom have an improper 
	    if (assignimp) then
	      k=k+1
	      atimp(i,k,1)=impnum(j,1) !1st atom number for k-th improper of i-th atom
	      atimp(i,k,2)=impnum(j,2) !2nd atom number for k-th improper of i-th atom
	      atimp(i,k,3)=impnum(j,3) !3rd atom number for k-th improper of i-th atom
	      atimp(i,k,4)=impnum(j,4) !4th atom number for k-th improper of i-th atom
	    endif
	  enddo
	  impxat(i)=k !total number of impropers for i-th atom
	enddo

	deallocate(impnum,presname,pimpxres,impatnamea)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%% END improper angles
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



! count bonds for each atom
	do i=1,nac !all atoms
	  bondxat(i)=0 !number of bonds for i-th atom
	  do j=1,6 !all posible conectivity
	    if (ng1(i,j).ne.0) bondxat(i)=bondxat(i)+1 !increase number of bonds (+1) for each atom conected to i
	  enddo
	enddo

! count angles for each atom in extreme(e)
	max_angle_ex=0
	do i=1,nac !all atoms
	  k=1
	  do j=1,bondxat(i) !all bonds for atom i-th
	    t=ng1(i,j) !atom number of bondxat(i) bond
	    do m=1,bondxat(t) !all bonds for t-th atom
	      if(ng1(t,m).ne.i) k=k+1  !if i in not t there are an angle
	    enddo
	  enddo
	  if (k.gt.max_angle_ex) max_angle_ex=k
	enddo

	allocate(atange(nac,max_angle_ex,2))
	do i=1,nac !all atoms
	  k=1
	  do j=1,bondxat(i) !all bonds for atom i-th
	    t=ng1(i,j) !atom number of bondxat(i) bond
	    do m=1,bondxat(t) !all bonds for t-th atom
	      if(ng1(t,m).ne.i) then !if i in not t assign an angle 
	        atange(i,k,1)=ng1(i,j) !assing atom number for middle atom of the k-th angle for i-th atom
	        atange(i,k,2)=ng1(t,m) !assing atom number for extreme atom of the k-th angle for i-th atom 
	        k=k+1
	      endif
	    enddo
	  enddo
	  angexat(i)=k-1 !number of angles for i-th atom in extreme
	enddo
 
! count angles for each atom in middle(m)
	max_angle_mid=0
	do i=1,nac !all atoms
	  k=1
	  do j=1,bondxat(i) !all bonds for atom i-th
	    do m=1,bondxat(i) !all bonds for atom i-th again 
	      if(ng1(i,m).gt.ng1(i,j)) k=k+1 !if both bonds are conected to different atoms there are an angle
	    enddo
	  enddo
	  if (k.gt.max_angle_mid) max_angle_mid=k
	enddo

	allocate(atangm(nac,max_angle_mid,2))
	do i=1,nac !all atoms
	  k=1
	  do j=1,bondxat(i) !all bonds for atom i-th
	    do m=1,bondxat(i) !all bonds for atom i-th again 
	      if(ng1(i,m).gt.ng1(i,j)) then !if both bonds are conected to different atoms assign an angle
	        atangm(i,k,1)=ng1(i,j) !assing atom number for extreme atom 1 of the k-th angle for i-th atom
	        atangm(i,k,2)=ng1(i,m) !assing atom number for extreme atom 2 of the k-th angle for i-th atom
	        k=k+1
	      endif
	    enddo
	  enddo
	  angmxat(i)=k-1 !number of angles for i-th atom in middle
	enddo

! count dihedral angles for each atom in extreme(e)
	max_dihe_ex=0
	do i=1,nac !all atoms
	  k=1
	  do j=1,angexat(i) !all angles with i in extreme
	    t=atange(i,j,2) !other atom in extreme for j angle
	    do m=1,bondxat(t) !all bonds of j
	      t2=ng1(t,m) !atom number bount to t 
	      if(t2.ne.atange(i,j,1)) k=k+1 !if t2 is not middle atom of angle between i and t there are a dihedral
	    enddo
	  enddo
	  if (k.gt.max_dihe_ex) max_dihe_ex=k
	enddo

	allocate(atdihe(nac,max_dihe_ex,3))
	do i=1,nac !all atoms
	  k=1
	  do j=1,angexat(i) !all angles with i in extreme
	    t=atange(i,j,2) !other atom in extreme for j angle
	    do m=1,bondxat(t) !all bonds of j
	      t2=ng1(t,m) !atom number bount to t 
	      if(t2.ne.atange(i,j,1)) then !if t2 is not middle atom of angle between i and t assign a dihedral
	        atdihe(i,k,1)=atange(i,j,1)
	        atdihe(i,k,2)=t
	        atdihe(i,k,3)=t2
	        k=k+1
	      endif
	    enddo
	  enddo
	  dihexat(i)=k-1 !number of dihedral angles for i-th atom in extreme
	enddo
 
! count dihedral angles for each atom in middle(m)
	max_dihe_mid=0
	do i=1,nac!all atoms
	  k=1
	  do j=1,angexat(i) !all angles with i in extreme
	    t=atange(i,j,1) !middle atom of j-th angle (then extreme atom is guaranteed)
	    do m=1,bondxat(i) !all bonds to i
	      t2=ng1(i,m) !atom bound to i
	      if(t2.ne.t) k=k+1 ! if t2 is not t there are an angle
	    enddo
	  enddo
	  if (k.gt.max_dihe_mid) max_dihe_mid=k
	enddo

	allocate(atdihm(nac,max_dihe_mid,3))
	do i=1,nac!all atoms
	  k=1
	  do j=1,angexat(i) !all angles with i in extreme
	    t=atange(i,j,1) !middle atom of j-th angle (then extreme atom is guaranteed)
	    do m=1,bondxat(i) !all bonds to i
	      t2=ng1(i,m) !atom bound to i
	      if(t2.ne.t) then ! if t2 is not t assign and dihedral
	        atdihm(i,k,1)=t2
	        atdihm(i,k,2)=t
	        atdihm(i,k,3)=atange(i,j,2)
	        k=k+1
	      endif
	    enddo
	  enddo
	  dihmxat(i)=k-1
	enddo

	return
 1      write(*,*) 'Problem reading amber.parm in impropers'
	stop
 2      write(*,*) 'Problem reading nresid in amber.parm from impropers'
	stop
 3      write(*,*) 'Problem reading impropers block in amber.parm file residue ',i
	stop
 4      write(*,*) 'Problem reading impropers block in amber.parm file  ',i,j
        stop

	end subroutine FF_bon_ang_dih_imp


!****************************************************
! This subroutine assign conectivity to each residue using amber.parm 


	subroutine FF_conectivity(nac,nataa,nroaa,atxres,atnu,resname,ng1, &
	number_conect_ommit, conect_ommit, atnamea,ncon,con, FF_max_at_by_res)
!nac,nataa,nroaa,atxres,atnu,resname,ng1,atnamea,ncon,con)
	use ionew, only : io_assign, io_close
	use precision, only:dp
	implicit none
	integer, intent(in) :: nac,nroaa, FF_max_at_by_res
	integer, intent(in), dimension(nroaa) :: atxres
	integer, intent(in), dimension(nroaa,100) :: atnu
	integer, intent(in), dimension(nroaa,100) :: nataa
	character*4, intent(in), dimension(nroaa)  :: resname
	integer, intent(out), dimension(nac,6) :: ng1
	integer, intent(in) :: ncon
	integer, intent(in), dimension(2,ncon) :: con
	character*4, intent(in), dimension(nroaa,FF_max_at_by_res) :: atnamea
!        integer na_u,nresid,                                                                   
	integer, intent(in) :: number_conect_ommit
	integer, dimension(number_conect_ommit,2), intent(in) :: conect_ommit

	character*4, dimension(:), allocatable :: presname
	integer, dimension(:,:,:), allocatable :: png1
	integer, dimension(:), allocatable :: bondxres
	integer :: maxconect !max number of conectivities in amber.parm
	integer :: nresid !number of residuer in amber.parm
	logical :: search
	logical :: include_conect
	character*12 :: option
	character*1 :: c1
	character*2 :: c2
	character*4 :: c4*4
	integer :: i, ii, j, k, l, m, n
	integer :: ui

	maxconect=-1
	call io_assign(ui)
	open(unit=ui,file="amber.parm")
	search=.true.
	do while (search)
	  read (ui,*,err=20,end=20) option
	  if (option.eq.'connectivity') then
	    read(ui,*,err=20,end=20) nresid
	    allocate(presname(nresid),bondxres(nresid))
	    presname=""
	    bondxres=-1
	    do i=1,nresid
	      read(ui,*,err=20,end=20) presname(i),bondxres(i)
	      do j=1,bondxres(i)
	        read(ui,*,err=20,end=20)
	        if (bondxres(i).gt.maxconect) maxconect=bondxres(i)
	      enddo
	    enddo 
	    search=.false.
	  endif
	enddo
	call io_close(ui)
	allocate(png1(nresid,maxconect,2))


	call io_assign(ui)
	open(unit=ui,file="amber.parm")
	search=.true.
	do while (search)
	  read (ui,*,err=20,end=20) option
	  if (option.eq.'connectivity') then
	    read(ui,*,err=20,end=20) nresid
	    png1=0
	    presname=""
	    bondxres=-1
	    do i=1,nresid
	      read(ui,*,err=20,end=20) presname(i),bondxres(i) !residue, number of bonds
	      do j=1,bondxres(i)
	        read(ui,*,err=20,end=20) png1(i,j,1),png1(i,j,2) !number of atoms conected in amber.parm definition
	      enddo
	    enddo    
	    search=.false.
	  endif
	enddo
	call io_close(ui)


!assign conectivity na los 1eros vecinos de cada atomo      
	do i=1,nroaa !residues in .fdf
	  do k=1,nresid !residues in amber.parm
	    if (resname(i).eq.presname(k)) then !check residue name
	      do j=1,atxres(i) !atoms in i-th residue
	        n=1  
	
	        do l=1,bondxres(k) !number of bond defined in amber.parm for l-th residue
	          if (nataa(i,j).eq.png1(k,l,1)) then !nataa(i,j) atom number for j-th atom in i-th residue check 1st atom for conectivity
	
	            do m=1,atxres(i) ! all atoms in residue i-th
	              if (nataa(i,m).eq.png1(k,l,2)) then !check 2nd atom for conectivity
	                ng1(atnu(i,j),n) = atnu(i,m) !asing conectivity
	                n=n+1
	              endif
	            enddo
	
	          elseif(nataa(i,j).eq.png1(k,l,2)) then !check 1st atom for conectivity
	
	            do m=1,atxres(i)
	              if (nataa(i,m).eq.png1(k,l,1)) then !check 2nd atom for conectivity
	                ng1(atnu(i,j),n) = atnu(i,m)
	                n=n+1
	              endif
	            enddo 
	
	          endif 
	
	        enddo
	      enddo
	    endif
	  enddo
	enddo






! generate conectivity between different aminoacids
	do i=2,nroaa !all residues in .fdf
	  include_conect=.false.
	  c4=resname(i)
	  c1=c4(1:1)
	  if(c1.eq.'N') then
	    if(c4.eq.'NME') then
	      include_conect=.true.
	    endif
	  else 
	    include_conect=.true.
	  endif

	  do ii=1, number_conect_ommit
	    if (conect_ommit(ii,1).eq.i) then
	      if (conect_ommit(ii,2).eq.i-1) then !skip conectivities for TER case in pdb, not tested yed
	        include_conect=.false.
		write(*,*) "omiting 1 conectivity!, this is not tested YED"
	      end if
	    end if
	  end do

	  if (include_conect) then
	    do j=1,atxres(i)
	      if(atnamea(i,j).eq.'N') then ! se fija si hay exclusion por connectivity extra 1er num -1 segundo el N
	        do k=1,ncon
	          if(con(1,k).eq.-1.and.con(2,k).eq.atnu(i,j)) then
	            write(*,*) 'skiping N atom number',atnu(i,j)
	            goto 5    
	          endif    
	        enddo      

	        write(*,*) 'Connecting N atom number',atnu(i,j)

	        do m=1,atxres(i-1)
	          if(atnamea(i-1,m).eq.'C') then
	            ng1(atnu(i,j),3) = atnu(i-1,m)
	            ng1(atnu(i-1,m),3) = atnu(i,j)
	            write(*,*) 'with ', atnu(i-1,m)
	          endif
	        enddo

	      endif
	    enddo
	  endif
 5	enddo


! generate conectivity between different nucleotids
	do i=2,nroaa
	  c4=resname(i)
	  c1=c4(3:3)
	  include_conect=.true.
	  if(c1.eq.'5') include_conect=.false.

	  if (include_conect) then
	    do j=1,atxres(i)
	      if(atnamea(i,j).eq.'P') then
	        do m=1,atxres(i-1)
	          c4=atnamea(i-1,m)
	          c2=c4(1:2)
	          if(c2.eq.'O3') then
	            ng1(atnu(i,j),4) = atnu(i-1,m)
	            ng1(atnu(i-1,m),2) = atnu(i,j)
	          endif
	        enddo
	      endif
	    enddo
	  endif
	enddo


! include extra conectivities defined in .fdf
	do i=1,ncon !extra conectivities 
	  if(con(1,i).eq.-1) goto 10

	  do k=1,6
	    if(ng1(con(1,i),k).eq.0) then
	      ng1(con(1,i),k)=con(2,i)
	      write(*,*) "included conectivity", con(1,i), con(2,i), k
	      do j=1,6
	        if(ng1(con(2,i),j).eq.0) then        
	          ng1(con(2,i),j)=con(1,i)   
	          write(*,*) "included conectivity", con(1,i), con(2,i), j
	          goto 10
	        endif
	      enddo
	    endif
	  enddo
 10       continue                 
	enddo

	deallocate(presname,bondxres,png1)      
	return
 20     stop &
	'solvent: Problem reading connectivity block in amber.parm file'
	end subroutine FF_conectivity


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
	Rm=-1.d0
	Em=-1.d0

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
	integer, dimension(:,:), allocatable :: pnataa,patmas 
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
 
	pc=987654321.d0
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


