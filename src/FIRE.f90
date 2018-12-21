	subroutine FIRE(natot, pos, Force, acel, vel, masst, time_steep, Ndescend, time_steep_max, alpha)
!FIRE, Bitzek, et. al., Phys. Rev. Lett. 97, 170201 (2006).
!FIRE optimization algorithm, Nfoglia 05/18 
	use scarlett, only : Ndamped
	implicit none
	integer, intent(in) :: natot
	double precision, dimension(3,natot), intent(in) :: Force
	double precision, dimension(3,natot), intent(inout) :: vel, acel, pos
	double precision, dimension(natot), intent(in) :: masst
	double precision, dimension(3,natot) :: acel_new
	double precision :: Fmod, velocity_proyected, velocity_proyected_at, velocity_mod
	double precision, intent(inout) :: time_steep, alpha, time_steep_max
	integer, intent(inout) :: Ndescend
	integer :: i,j
	logical :: damp
   	damp=.false.
	Fmod=0.d0
	do i=1, natot
	  acel_new(1:3, i) = Force(1:3, i)/masst(i)
	  do j=1,3
	    Fmod=Fmod + Force(j, i)**2
	  end do
	end do
	Fmod=sqrt(Fmod)
	vel=vel+0.5d0*(acel+acel_new)*time_steep
	velocity_proyected=0.d0 !P in paper
	velocity_mod=0.d0
    
	do i=1, natot
	  velocity_proyected_at=0.d0
	  do j=1,3
	    velocity_proyected_at=velocity_proyected_at+vel(j,i)*Force(j,i)
	    velocity_mod=velocity_mod+vel(j,i)*vel(j,i)
	  end do
!	  write(*,*) "atom", i, "VP ", velocity_proyected_at
	  velocity_proyected=velocity_proyected+velocity_proyected_at
	  if (velocity_proyected_at .lt. 0.d0) damp=.true.
	end do
!		damp=.false.
!		if (velocity_proyected .lt. 0.d0) damp=.true.
	velocity_mod=sqrt(velocity_mod)

	if (.not. damp) then
	  Ndescend=Ndescend+1
	  if (Ndescend .gt. 5) then
	    time_steep=min(time_steep*1.1d0, time_steep_max)
	    alpha=alpha*0.99d0
	  end if
	  do  i=1, natot
	    Fmod=0.d0
	    velocity_mod=0.d0
	
	    do j=1,3
	      Fmod=Fmod + Force(j, i)**2
	      velocity_mod=velocity_mod+vel(j,i)*vel(j,i)
	    end do
	     velocity_mod=sqrt(velocity_mod)
	     Fmod=sqrt(Fmod)
	    vel(1:3,i)=(1.d0-alpha)*vel(1:3,i)+alpha*velocity_mod*Force(1:3,i)/Fmod
	  end do
	else
	  Ndamped=Ndamped+1
	  write(*,*) "damping system", Ndamped
	  Ndescend=0
	  alpha=0.1d0
	  time_steep=time_steep*0.5d0
	  vel=0.d0
	end if

!move images
	pos=pos+vel*time_steep+0.5d0*acel*time_steep**2
	acel=acel_new
	return
	end subroutine FIRE

