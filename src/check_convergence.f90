        SUBROUTINE check_convergence(relaxd, natot, force)
        use scarlett, only: ftol
        implicit none
        integer :: i,j !auxiliar
        logical, intent(inout) :: relaxd
	double precision, dimension(3,natot), intent(in) :: force
        double precision :: Fmod2
        relaxd=.true.
        Fmod2Max=-1d80
        do i=1, natot
          do j=1,3
            Fmod2=force(j,i)**2
            relaxd=relaxd .and. (Fmod2 .lt. ftol**2)
            if (Fmod2 .gt. Fmod2Max) Fmod2Max=Fmod2
          end do
        end do
        return
        END SUBROUTINE check_convergence
