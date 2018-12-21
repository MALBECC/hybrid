	subroutine quick_min(natot, pos, Force, acel, vel, masst)
!quick min optimization algorithm
!Nfoglia 05/18 
	use scarlett, only : time_steep, Ndamped
	implicit none
	integer, intent(in) :: natot !number of atoms
	double precision, dimension(3,natot), intent(in) :: Force !forces
	double precision, dimension(3,natot), intent(inout) :: vel, acel, pos !velocity, aceleration and positions
	double precision, dimension(natot), intent(in) :: masst !atomic mass
	double precision, dimension(3,natot) :: acel_new !acelerarion in this step
	integer :: i, j !auxiliar
	double precision :: Fmod, velocity_proyected, velocity_proyectedR 
	logical :: freeze
	Fmod=0.d0
        do i=1, natot
          acel_new(1:3, i) = Force(1:3, i)/masst(i)
          do j=1,3
              Fmod=Fmod + Force(j, i)**2
          end do
        end do
	Fmod=sqrt(Fmod)

	velocity_proyected=0.d0
        velocity_proyectedR=0.d0 
          
        vel=vel+0.5d0*(acel+acel_new)*time_steep

	freeze=.false.
	i=0
        do while (.not. freeze)
	  i=i+1 !atom number
          velocity_proyectedR=0.d0
          do j=1,3
            velocity_proyectedR=velocity_proyectedR+vel(j,i)*Force(j,i)
          end do

          if (velocity_proyectedR .lt. 0.d0 ) then
	    freeze=.true.
	    velocity_proyected=0.d0
	    vel=0.d0
	    Ndamped=Ndamped+1
	    write(*,*) "damping system", Ndamped
	  else
	    velocity_proyected=velocity_proyected+velocity_proyectedR
	  end if
	  if (i .eq. natot) freeze=.true.
	end do

	pos=pos+vel*time_steep+0.5d0*acel*time_steep**2
	acel=acel_new
	return
 556 format(2x,A10, 2x, i3, 13(2x,f20.8))
	end subroutine quick_min
