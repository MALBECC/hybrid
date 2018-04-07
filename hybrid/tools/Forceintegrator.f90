	program Forceintegrator
!betha version of Forceintegrator
	implicit none
	double precision, dimension(:,:,:), allocatable :: r, F
	double precision, dimension(:), allocatable :: Work
	double precision :: temp
	integer :: i,j,k, natoms, npict, trash
	npict=24
	natoms=6666

	temp=0.d0
	allocate (r(0:npict,natoms,3), f(0:npict,natoms,3))
	allocate(Work(natoms+1))
	r=0.d0
	f=0.d0
	Work=0.d0


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

	f=f*627.509391d0 !change to kcal/bohr mol
!	r=r falta checkear q este parametro este bien
!a*0.5d0
!29177d0

!integral
	do i=0,npict-1
	  do j=1, natoms
	    do k=1,3
	      Work(j)=Work(j) - (f(i,j,k)+f(i+1,j,k))*(r(i+1,j,k)-r(i,j,k))/2.d0 
	    end do !k
	  end do !j

!total work
	  temp=0.d0
	  do j=1, natoms
	    temp=temp+Work(j)
	  end do
	  write(*,556) dble(i)+1.5d0, Work(1:natoms), temp
	end do

 556 format(f5.2,2x,300000(f20.10,2x))
	end program Forceintegrator
