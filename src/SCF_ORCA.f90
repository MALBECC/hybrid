	subroutine SCF_orca(na_u, at_MM_cut_QMMM, r_cut_QMMM, Etot, F_cut_QMMM, Iz_cut_QMMM )
	use scarlett, only:charge, spin, iza, Ang, qm_command, qm_header_lines, qm_header
	implicit none
	integer, intent(in) :: na_u, at_MM_cut_QMMM
	double precision, dimension(3,na_u+at_MM_cut_QMMM), intent(in) :: r_cut_QMMM
	double precision, dimension(3,na_u+at_MM_cut_QMMM), intent(out) :: F_cut_QMMM
	double precision,intent(out) :: Etot
	double precision, intent(in) ::  Iz_cut_QMMM(na_u+at_MM_cut_QMMM)

        character*550 :: command
	character*550 :: paste
	external :: paste

        character*33 :: temp
        integer :: i,j




!write orca.in
	open(unit=720, file='orca.in')
	do i=1, qm_header_lines
	  write(720,'(A)') trim(qm_header(i))
	end do
	  write(720,'(A)')
	if (at_MM_cut_QMMM.gt.0) write(720,'(A)') '% pointcharges "orca_pcharges.pc"'
!'
	write(720,'(A)') ''
	write(720,'(A5,x,i2,x,i2)') '* xyz',  charge, int(2*spin+1)
	write(720,'(A)') ''
	do i=1,na_u
	  write(720,*) iza(i), r_cut_QMMM(1:3,i)/Ang
	end do
	write(720,'(A)') '*'

	close(720)
!write poit charges file
	
	if (at_MM_cut_QMMM.gt.0) then
	  open(unit=721, file='orca_pcharges.pc')
	  write(721,*) at_MM_cut_QMMM
	  do i=1, at_MM_cut_QMMM
	    write(721,*) Iz_cut_QMMM(i), r_cut_QMMM(1:3,i)
	  end do
	  close(721)
	end if

!call orca SCF & forces calculation
	command=paste(qm_command,' orca.in')
	call execute_command_line(command)
!command)

!read energy & QM forces
	F_cut_QMMM=0.d0
	open(unit=722, file='orca.engrad')
	do i=1,7
	  read(722,*) temp
	end do
	read(722,*) Etot
	do i=1,3
	  read(722,*) temp
	end do

	do i=1,na_u
	  do j=1,3
	    read(722,*) F_cut_QMMM(j,i)
	  end do
	end do
	close(722)

!read MM forces
	if (at_MM_cut_QMMM.gt.0) then
	open(unit=722, file='orca.pcgrad')
	  read(722,*) temp
	do i=1,at_MM_cut_QMMM
	  read(722,*) F_cut_QMMM(1:3,i+na_u)
	end do
	close(722)
	end if

	F_cut_QMMM=-F_cut_QMMM

!	do i=1,at_MM_cut_QMMM+na_u
!	  write(*,*) F_cut_QMMM(1:3,i)
!	end do

  700 FORMAT('a200')
	end subroutine SCF_orca
