	subroutine FIRE(natot, pos, Force, vel,  time_steep, Ndescend, time_steep_max, alpha)
!Herbol, H. C.; Stevenson J.; Clancy, P., JCTC 2017 13 (7), 3250-3259
!FIRE optimization algorithm, Nfoglia 05/19 
	use scarlett, only : Ndamped
	implicit none
	double precision, dimension(3,natot) :: vel_now, step !FIRE velocity, and movement step
	double precision :: stepsize, stepmax !square distance and ma square distance in step
	integer, intent(in) :: natot !number of atoms
	double precision, dimension(3,natot), intent(in) :: Force
	double precision, dimension(3,natot), intent(inout) :: vel, pos
	double precision :: Fmod, velocity_proyected, velocity_proyected_at, velocity_mod 
	double precision, intent(inout) :: time_steep, alpha, time_steep_max
	integer, intent(inout) :: Ndescend
	logical :: damp !damping system
	integer :: i,j !auxiliar

!Assign FIRE velocities
	velocity_proyected=0.d0
	do  i=1, natot
	  Fmod=0.d0
	  velocity_mod=0.d0
	  velocity_proyected_at=0.d0
	  do j=1,3
	    Fmod=Fmod + Force(j, i)**2
	    velocity_mod=velocity_mod+vel(j,i)*vel(j,i)
	    velocity_proyected_at=velocity_proyected_at+vel(j,i)*Force(j,i)
	    if (vel(j,i).ne.vel(j,i)) STOP "NAN in VEL FIRE"
	  end do
	  velocity_mod=sqrt(velocity_mod)
	  Fmod=sqrt(Fmod)
	  velocity_proyected=velocity_proyected+velocity_proyected_at !P in paper
	  if (Fmod.ne. 0.d0) then
	    vel_now(1:3,i)=(1.d0-alpha)*vel(1:3,i)+alpha*velocity_mod*Force(1:3,i)/Fmod
	  else
	    vel_now(1:3,i)=vel(1:3,i)
	  end if
	end do

!damp system
	damp=.false.
	if (velocity_proyected .lt. 0.d0) damp=.true.
	if (.not. damp) then
	  if (Ndescend .gt. Ndamped) then
	    time_steep=min(time_steep*1.1d0, time_steep_max)
	    alpha=alpha*0.99d0
	  end if
	  Ndescend=Ndescend+1
	else
	  write(*,*) "damping system", alpha, time_steep
	  vel_now=0.d0
	  alpha=0.1d0
	  time_steep=time_steep*0.5d0
	  Ndescend=0
	end if


!Euler step
	vel=vel_now + Force*time_steep
	step=vel*time_steep
	stepmax=0.d0
	do i=1,natot
	  stepsize=0.d0
	  do j=1,3
	    stepsize=stepmax+step(j,i)**2
	  enddo
	  if (stepsize.gt.stepmax) stepmax=stepsize
	end do

	if(stepmax .gt. 0.1d0) step=0.1d0*step/sqrt(stepmax)
	pos=pos+step

	end subroutine FIRE


