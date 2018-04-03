	program Restart_maker
! this program create a restart for lio for a system of nco_nuevo occupied orbitals
! using a restart generated with nco_viejo  occupied orbitals
! N. Foglia  03/2018
	implicit none
	character*100 :: inputname,outputname
	integer :: nco_input,nco_output
	logical :: change_to_open, readOK
	logical :: open
	integer :: nco_viejo, nco_nuevo, Nbases
	integer :: new_restart
	double precision, dimension (:,:), allocatable :: matrix
	integer :: i,j, jmax
	double precision :: occup_value
	
!defaults
	inputname="just_a_big_star"
	outputname="just_a_big_star"
	nco_input=-1
	nco_output=-1
	new_restart=-1
	change_to_open=.false.
	
	call read_options(inputname,outputname,nco_input,nco_output,change_to_open, new_restart)
!check read options
	readOK=.true.
	readOK=readOK .and. (nco_input .ne. -1)
	readOK=readOK .and. (nco_output .ne. -1)
	readOK=readOK .and. .not.(inputname == "just_a_big_star")
	readOK=readOK .and. .not. (outputname == "just_a_big_star" )
	
	if (.not. readOK) STOP "error in Restart_maker.in"
	
	nco_viejo=nco_input
	nco_nuevo=nco_output
	open=change_to_open
	
	
!count number of basis
	Nbases=0
	open(unit=973, file=inputname)
	DO 
	  READ (973,*, END=10) 
	  Nbases = Nbases + 1 
	END DO 
	10 CLOSE (973) 
	
	allocate(matrix(nco_nuevo,Nbases))
	
	open(unit=973, file=inputname)
	matrix=0.d0
	do i=1, Nbases
	  read(973,*) matrix(1:nco_viejo,i)
	end do
	close(973)
	
	jmax=1
	if (change_to_open) jmax=2
	
	
	if (new_restart .gt. 0) then
	  if (new_restart .eq. 1) occup_value=2.0 !close shell calculation
	  if (new_restart .eq. 2) occup_value=1.0 !open shell calculation
	  if (new_restart .gt. 2) STOP "wrong new_restart value"
	  matrix=0.d0
	  do i=1, nco_nuevo
	    matrix(i,i)=occup_value
	  end do
	end if
	
	
	open(unit=974, file=outputname)
	do j=1, jmax
	do i=1, Nbases
	  write(974,*) matrix(1:nco_nuevo,i)
	end do
	end do
	close(974)
	end program Restart_maker
	
	
	subroutine read_options(inputname,outputname,nco_input,nco_output,change_to_open, new_restart)
	implicit none
	character*100, intent(inout) :: inputname,outputname
	integer, intent(inout) :: nco_input,nco_output
	logical, intent(inout) :: change_to_open
	integer, intent(inout) :: new_restart
	integer :: ierr, ios
	logical :: fileExists
	namelist / Milla / inputname,outputname,nco_input,nco_output,change_to_open, new_restart
	inquire(file ="Restart_maker.in", exist = fileExists)
	
	if(fileExists) then
	  open(unit = 100, file="Restart_maker.in", iostat = ios)
	  read(100, nml = Milla, iostat = iErr)
	  if(ierr.gt.0) stop 'Input error in Milla namelist.'
	  close(100)
	else
	  write(*,*) 'Restart_maker.in not found:'
	  write(*,*) 'create file with correct value for this variables:'
	  write(*,*) 'inputname'
	  write(*,*) 'outputname'
	  write(*,*) 'nco_input'
	  write(*,*) 'nco_output'
	  write(*,*) 'change_to_open'
	  STOP
	end if
	RETURN
	end subroutine read_options


