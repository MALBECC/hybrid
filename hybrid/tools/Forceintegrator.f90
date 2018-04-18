	program Forceintegrator
!integrate forces in Pos_forces.dat file for obtain energy contributions
	implicit none
	double precision, dimension(:,:,:), allocatable :: r, F
	double precision, dimension(:), allocatable :: Work
	double precision, dimension(:), allocatable :: Work_cut
	double precision :: temp, distance
	integer :: i,j,k, natoms, npict, trash, icut
	integer :: nlines, io, atnum, atnum_old

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
