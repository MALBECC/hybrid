	program armarXV
!create a number fotos_salida XVs using a number of fotos_entrada XVs
	implicit none
	integer :: fotos_entrada, fotos_salida
	character*35 :: nombre_base, nombre_lectura, nombre_salida
	character*3 :: stri
	double precision, dimension(:,:,:), allocatable :: r_V
	double precision, dimension(:,:), allocatable :: r_V_salida
	double precision, dimension(:), allocatable :: relpos_i, relpos_f
	double precision, dimension(3,3) :: box
	integer :: atomos
	integer :: i,j, positioni, trash
	double precision :: fact_ant, fact_post
	logical :: namelist_read, fileExists
	namelist /milla/ fotos_entrada, fotos_salida, atomos, nombre_base

	fotos_entrada=-1
	fotos_salida=-1
	atomos=-1
	nombre_base=""
	
	inquire(file = "armarXV.in", exist = fileExists)
	if(fileExists) then
	  open(unit = 100, file = "armarXV.in")
	  read(100, nml=milla)
	  close(unit = 100)
	else
	  stop 'File armarXV.in not found'
	endif
	
	namelist_read=.true.
	namelist_read=namelist_read .and. (fotos_entrada.ne.-1)
	namelist_read=namelist_read .and. (fotos_salida.ne.-1)
	namelist_read=namelist_read .and. (atomos.ne.-1)
	if (.not. namelist_read) then
	  write(*,*) "armarXV.in should have:"
	  write(*,*) "&milla"
	  write(*,*) "fotos_entrada="
	  write(*,*) "fotos_salida="
	  write(*,*) "atomos="
	  write(*,*) "nombre_base="
	  write(*,*) "&end"
	  stop
	end if

	write(*,*) "interponado posiciones entre nombre_base 1 y ", fotos_entrada
	allocate (r_V(fotos_entrada, atomos, 6))
	allocate (r_V_salida(atomos, 6))

	do i=1, fotos_entrada
	  if(i.lt.10) write(stri,'(I1)') i
	  if(i.lt.100 .and. i.ge.10) write(stri,'(I2)') i
	  nombre_lectura=trim(nombre_base)//trim(stri)
	  write(*,*) "leo ", nombre_lectura
	  open(unit=954, file=nombre_lectura)
	  read(954,*) box(1,1:3)
	  read(954,*) box(2,1:3)
	  read(954,*) box(3,1:3)
	  read(954,*) trash
	  do j=1, atomos
	    read(954,*) r_V(i,j,1:6)
	  end do
	  close(954)
	end do

	allocate(relpos_i(fotos_entrada),relpos_f(fotos_salida))

	do i=1, fotos_entrada
	  relpos_i(i)=1.d0/dble(fotos_entrada-1) * (i-1)
	end do

	do i=1, fotos_salida
	  relpos_f(i)=1.d0/dble(fotos_salida-1) * (i-1)
	end do

	do i=2, fotos_salida-1
	  positioni=1
	  do while (relpos_f(i) .gt. relpos_i(positioni))
	    positioni=positioni+1
	  end do
	  fact_ant=relpos_i(positioni)-relpos_f(i)
	  fact_post=relpos_f(i)-relpos_i(positioni-1)
	  
	 if (fact_post .lt. 0.d0 .or. fact_ant .lt. 0.d0 ) STOP "error in fact"
	  fact_post=fact_post/(relpos_i(positioni)-relpos_i(positioni-1))
	  fact_ant=fact_ant/(relpos_i(positioni)-relpos_i(positioni-1))
	  if (positioni .gt. fotos_entrada ) then
	    write(*,*) "i, positioni",i, positioni
	  end if
	  r_V_salida(1:atomos,1:6)=r_V(positioni-1,1:atomos,1:6)*fact_ant + r_V(positioni,1:atomos,1:6)*fact_post
	
	  if(i.lt.10) write(stri,'(I1)') i
	  if(i.lt.100 .and. i.ge.10) write(stri,'(I2)') i
	  nombre_salida=trim(nombre_base)//"out."//trim(stri)
	  write(*,*) "escribo ", nombre_salida
	  open(unit=955, file=nombre_salida)
	  write(955,456) box(1,1:3)
	  write(955,456) box(2,1:3)
	  write(955,456) box(3,1:3)
	  write(955,*) atomos
	  do j=1, atomos
	    write(955,457) r_V_salida(j,1:6)
	  end do
	  close(955)
	end do

  456 FORMAT(3x,3(2x,f16.9))
  457 FORMAT(3x,6(2x,f16.9))

	end program armarXV
