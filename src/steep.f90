!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Steepest Decend algorithm for geometry optimization
!
! Nicolas Foglia, 2019
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

	SUBROUTINE steep(natot, rclas, cfdummy, Energy, istep)
	use scarlett, only: Eprev, lambda, lambda_i, Steep_change

	implicit none
	integer, intent(in) :: natot
	real*8, intent(in) :: Energy !energy

	real*8 :: Fmax !max |force| on an atoms
	integer :: istep !steep number
	double precision, dimension(3, natot) :: rclas, cfdummy


	Fmax=0.d0
	if (istep.eq.1) Eprev=Energy
	if (istep.eq.1) lambda=lambda_i
	call maxforce(cfdummy,natot,Fmax) !find max force
	call move_movesteep(lambda, Fmax, natot, rclas, cfdummy)

	if (Steep_change) then
	  if (Energy.le.Eprev) then
	    lambda=lambda*1.2d0
	  else
	    lambda=lambda*0.5d0
	  end if
	  if (lambda.gt.0.1d0) lambda=0.5d0
	  Eprev=Energy
	end if
	END SUBROUTINE steep


	subroutine maxforce(cfdummy,natot,Fmax)
	implicit none
	integer, intent(in) :: natot
	double precision, intent(out) :: Fmax
	double precision, intent(in) :: cfdummy(3,natot)
	double precision :: F_i
	integer :: i
	Fmax=0.d0
	do i=1, natot
	  F_i=cfdummy(1,i)**2 + cfdummy(2,i)**2 + cfdummy(3,i)**2
	  F_i=sqrt(F_i)
	  if (F_i .gt. Fmax) Fmax=F_i
	end do
	return
	end subroutine maxforce


	subroutine move_movesteep(lambda, Fmax, natot, rclas, cfdummy) !calculate new positions
	implicit none
	integer, intent(in) :: natot
	double precision, dimension(3,natot), intent(inout) :: rclas
	double precision, dimension(3,natot), intent(in) :: cfdummy
	double precision, intent(in) :: lambda, Fmax
	INTEGER :: i
	double precision :: a
	a=lambda/Fmax
	do i=1,natot
	  rclas(1:3,i)= rclas(1:3,i)+a*cfdummy(1:3,i)
	end do
	end subroutine move_movesteep

