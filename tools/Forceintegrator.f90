	program Forceintegrator
!integrate forces in Pos_forces.dat file for obtain energy contributions
!2018 N. foglia
	implicit none
	double precision, dimension(:,:,:), allocatable :: r, F
	double precision, dimension(:), allocatable :: Work
	double precision, dimension(:), allocatable :: Work_cut
	double precision :: temp, distance
	integer :: i,j,k, natoms, npict, trash, icut
	integer :: nlines, io, atnum, atnum_old
	integer, dimension(:), allocatable :: Iz
	double precision, dimension(9) :: In0vec
	logical :: file_exists

!determine npict and natoms
	open(unit=25, file="Pos_forces.dat")
	nlines=0
	atnum=0
	atnum_old=0
	natoms=0
	npict=0
	DO
	  READ(25,*,iostat=io) atnum
	  if (atnum .eq. atnum_old+1) then
	    natoms=natoms+1
	    atnum_old=atnum
	  end if 
	  IF (io/=0) EXIT
	  nlines = nlines + 1
	END DO
	npict=nlines/natoms
	close(25)
	write(*,*) "sistem with ", natoms, "atoms and ", npict, "frames"


	temp=0.d0
	allocate (r(0:npict,natoms,3), f(0:npict,natoms,3))
	allocate(Work(natoms+1), Work_cut(30))
	allocate(Iz(natoms))

	r=0.d0
	f=0.d0
	Work=0.d0
	Work_cut=0.d0
	npict=npict-1
!read
	open(unit=25, file="Pos_forces.dat")
	k=0
	do i=0,npict
	do j=1, natoms
	  k=k+1
	  read(25,*) trash, r(i,j,1:3),f(i,j,1:3)
	end do
	  k=k+1
	  read(25,*) 
	end do
	close(25)

!center traj
	inquire(file = 'IZ_input', exist = file_exists)

	if (file_exists) then
	  write(*,*) "reading IZ from IZ_input"
	  open(unit=30, file="IZ_input")

	  do j=1, natoms
	    read(30,*) Iz(j)
	  end do

	  write(*,*) "centering trajectory"
	  In0vec=0.d0
	  do i=0,npict
	    call center_rotation(natoms, Iz, r(i,1:natoms,1:3), i, In0vec,F(i,1:natoms,1:3))
	  end do
	else
	  write(*,*) "didnt find IZ_input system will be not centered"
	end if



	open(unit=26, file="atomic_work.dat")
	open(unit=27, file="cutoff_work.dat")
!integral
	do i=0,npict-1
	  do j=1, natoms
	    distance = (r(i,j,1)-r(1,j,1))**2
	    distance = distance + (r(i,j,2)-r(1,j,2))**2
	    distance = distance + (r(i,j,3)-r(1,j,3))**2

	    do k=1,3
	      Work(j)=Work(j) - (f(i,j,k)+f(i+1,j,k))*(r(i+1,j,k)-r(i,j,k))/2.d0 

	      do icut=1, 30
	        if (distance .lt. dble(icut)**2) then
		  Work_cut(icut)=Work_cut(icut)- (f(i,j,k)+f(i+1,j,k))*(r(i+1,j,k)-r(i,j,k))/2.d0  
		end if
	      end do

	    end do !k
	  end do !j

!total work
	  temp=0.d0
	  do j=1, natoms
	    temp=temp+Work(j)
	  end do
	  write(26,556) dble(i)+1.5d0, Work(1:natoms), temp
	  write(27,556) dble(i)+1.5d0, Work_cut

	end do
	close(26)
	close(27)
	write(*,*) "integration done"
	write(*,*) "results in atomic_work.dat and cutoff_work.dat" 
	
 556 format(f5.2,2x,300000(f20.10,2x))
	end program Forceintegrator


	subroutine center_rotation(atoms, Iz, r, firstcet, In0vec,F)
!center system and keeps direction of autovectors of inercia tensor
!2017 N. foglia
	implicit none
        integer, intent(in):: atoms
!	double precision, intent(in), dimension(atoms) :: mass
	integer, dimension(atoms) :: Iz
	double precision, dimension(36) :: massZ
	double precision, intent(inout), dimension(atoms,3) :: r,F
	double precision, dimension(3) :: RCM, rtemp, Ftemp
	double precision ::  Tmass
	double precision, dimension(3,3) :: Ini
	double precision, dimension(9) :: Inivec, Intempvec
	double precision, dimension(3,3) :: Inivec2
        double precision, dimension(9), intent(inout) :: In0vec
	double precision :: rfactor
	integer :: firstcet
        double precision, dimension(1000) :: work
        double precision, dimension(3) ::w
        integer, dimension(1000) :: iwork
        integer :: info
	integer :: LWORK, LIWORK, LWMAX

	double precision :: proj1, proj0
	integer :: i, j, j2
	integer :: case12, signcase, caseIn, signcaseIn
	integer :: wi

	LWMAX=1000
	massZ= (/1.d0, 4.d0, 7.d0, 9.d0, 11.d0, 12.d0, 14.d0, 16.d0, 19.d0, 20.d0, 23.d0, 24.d0, 27.d0, 28.d0, 31.d0, 32.d0, &
	35.d0, 40.d0, 39.d0, 40.d0, 45.d0, 48.d0, 51.d0, 52.d0, 55.d0, 56.d0, 59.d0, 59.d0, 64.d0, 65.d0, 70.d0, 73.d0, 75.d0,&
	79.d0, 80.d0, 84.d0/) 
!ampliar luegpo massZ hasta 118


!calculo el centro de masa y centro al sistema
        RCM=0.d0
        Tmass=0.d0
        do i=1, atoms
          do j=1, 3
            RCM(j)=RCM(j)+massZ(Iz(i))*r(i,j)
          end do
          Tmass=Tmass+massZ(Iz(i))
        end do
        RCM=RCM/Tmass

        do i=1, atoms
          do j=1, 3
            r(i,j)=r(i,j)-RCM(j)
          end do
        end do

!calcula el tensor de inercia inicial
        Ini=0.d0
        do j=1,3
          do j2=1,3
            do i=1, atoms
               rfactor=-r(i,j)*r(i,j2)
               if (j.eq.j2) rfactor=rfactor + r(i,1)*r(i,1) + r(i,2)*r(i,2) + r(i,3)*r(i,3)
               Ini(j,j2)=Ini(j,j2)+massZ(Iz(i))*rfactor
            end do
          end do
        end do

        wi=0
	Inivec2=0.d0
        do j=1,3
          do j2=j,3
            wi=wi+1
	    Inivec2(j,j2)=Ini(j,j2)
          end do
        end do

!calculo autovectores del tensor de inercia
        LWORK = -1
        LIWORK = -1
        call DSYEVD( 'V', 'U', 3, inivec2, 3, W, WORK, LWORK, IWORK, LIWORK, INFO )
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        LIWORK = MIN( LWMAX, IWORK( 1 ) )
        call DSYEVD( 'V', 'U', 3, inivec2, 3, W, WORK, LWORK, IWORK, LIWORK, INFO )
        if (info .ne. 0) stop "info ne 0, problem in dsyevd"

        wi=0
        do j=1,3
          do j2=1,3
            wi=wi+1
            inivec(wi)=Inivec2(j2,j)
          end do
        end do


        if (firstcet.eq.0) then
           in0vec=inivec
        else !ordeno ejes para q conincidan con el tensor de inercia del paso anterior
           proj0=0.d0
           do case12=1, 6
             do signcase=1, 8
               call switch(inivec, intempvec, case12, signcase)
               proj1=multi3vec(intempvec,in0vec)
               if (proj1.gt.proj0)  then
                 caseIn=case12
                 signcaseIn=signcase
                 proj0=proj1
               end if
             end do
           end do

          call switch(inivec, in0vec, caseIn, signcaseIn)
          inivec=in0vec
	end if

!paso el sistema a la base que diagonaliza el tensor de inercia
        do i=1, atoms
          rtemp=0.d0
          rtemp(1)=r(i,1)*inivec(1)+r(i,2)*inivec(2)+r(i,3)*inivec(3)
          rtemp(2)=r(i,1)*inivec(4)+r(i,2)*inivec(5)+r(i,3)*inivec(6)
          rtemp(3)=r(i,1)*inivec(7)+r(i,2)*inivec(8)+r(i,3)*inivec(9)
          Ftemp(1)=F(i,1)*inivec(1)+F(i,2)*inivec(2)+F(i,3)*inivec(3)
          Ftemp(2)=F(i,1)*inivec(4)+F(i,2)*inivec(5)+F(i,3)*inivec(6)
          Ftemp(3)=F(i,1)*inivec(7)+F(i,2)*inivec(8)+F(i,3)*inivec(9)
          do j=1,3
            r(i,j)=rtemp(j)
	    F(i,j)=Ftemp(j)
          end do
        end do

	open(unit=966, file="trj_cent.xyz")
	write(966,*) atoms
        write(966,*)
        do i=1, atoms
          write(966,345) Iz(i), r(i,1:3)
        end do

  345  format(2x, I2,    2x, 3(f10.6,2x))

	contains

	subroutine switch(VV1in, VV2, case12, signcase)
        implicit none
        double precision, dimension(9), intent(in) :: VV1in
	double precision, dimension(9) :: VV1
	double precision, dimension(9), intent(out) :: VV2
	integer, intent(in) :: case12, signcase

	if (signcase .eq. 1) then
	   VV1=VV1in
	elseif (signcase .eq. 2) then
	   VV1(1:3)=VV1in(1:3)
           VV1(4:6)=VV1in(4:6)
           VV1(7:9)=-VV1in(7:9)
        elseif (signcase .eq. 3) then
           VV1(1:3)=VV1in(1:3)
           VV1(4:6)=-VV1in(4:6)
           VV1(7:9)=VV1in(7:9)
        elseif (signcase .eq. 4) then
           VV1(1:3)=VV1in(1:3)
           VV1(4:6)=-VV1in(4:6)
           VV1(7:9)=-VV1in(7:9)
        elseif (signcase .eq. 5) then
           VV1(1:3)=-VV1in(1:3)
           VV1(4:6)=VV1in(4:6)
           VV1(7:9)=VV1in(7:9)
        elseif (signcase .eq. 6) then
           VV1(1:3)=-VV1in(1:3)
           VV1(4:6)=VV1in(4:6)
           VV1(7:9)=-VV1in(7:9)
        elseif (signcase .eq. 7) then
           VV1(1:3)=-VV1in(1:3)
           VV1(4:6)=-VV1in(4:6)
           VV1(7:9)=VV1in(7:9)
        elseif (signcase .eq. 8) then
           VV1(1:3)=-VV1in(1:3)
           VV1(4:6)=-VV1in(4:6)
           VV1(7:9)=-VV1in(7:9)
	else
           stop "bad case selected on sign"
        end if


	if (case12 .eq. 1) then
	   VV2=VV1
	elseif (case12 .eq. 2) then
	   VV2(1:3)=VV1(1:3)
	   VV2(4:6)=VV1(7:9)
	   VV2(7:9)=VV1(4:6)
        elseif (case12 .eq. 3) then
           VV2(1:3)=VV1(4:6)
           VV2(4:6)=VV1(1:3)
           VV2(7:9)=VV1(7:9)
        elseif (case12 .eq. 4) then
           VV2(1:3)=VV1(4:6)
           VV2(4:6)=VV1(7:9)
           VV2(7:9)=VV1(1:3)
        elseif (case12 .eq. 5) then
           VV2(1:3)=VV1(7:9)
           VV2(4:6)=VV1(1:3)
           VV2(7:9)=VV1(4:6)
        elseif (case12 .eq. 6) then
           VV2(1:3)=VV1(7:9)
           VV2(4:6)=VV1(4:6)
           VV2(7:9)=VV1(1:3)
	else
	   write(*,*) "llamo con", case12
	   stop "bad case selected on switch"
	end if
	return
	end subroutine switch


	double precision function multi3vec(VV1, VV2)
	implicit none
	double precision, dimension(9), intent(in) :: VV1, VV2
	integer :: k
	multi3vec=0.d0
	do k=1, 9
	  multi3vec=multi3vec+VV1(k)*VV2(k)
	end do
	return
	end function multi3vec

	end subroutine center_rotation

