        SUBROUTINE check_convergence(relaxd, force)
        use scarlett, only: natot, ftol
        implicit none
        integer :: i,j !auxiliar
        logical, intent(inout) :: relaxd
	double precision, dimension(3,natot), intent(in) :: force
        double precision :: Fmod2
        relaxd=.true.
        do i=1, natot
          do j=1,3
            Fmod2=force(j,i)**2
            relaxd=relaxd .and. (Fmod2 .lt. ftol**2)
          end do
        end do
        return
        END SUBROUTINE check_convergence
