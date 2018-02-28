	SUBROUTINE bandmove(istep, na_u,replicas,rclas_BAND,fclas_BAND, Energy_band, relaxd)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: na_u,replicas, istep
	DOUBLE PRECISION, DIMENSION(3,na_u,replicas), INTENT(INOUT) :: fclas_BAND
	DOUBLE PRECISION, DIMENSION(3,na_u,replicas), INTENT(INOUT) :: rclas_BAND
	DOUBLE PRECISION, DIMENSION(3,na_u) :: tang_vec
	DOUBLE PRECISION, DIMENSION(replicas), intent(In) :: Energy_band
	DOUBLE PRECISION, DIMENSION(replicas) :: Energy_band_old
	DOUBLE PRECISION :: NORMVEC, proyForce, MAXFmod, MAXFmod_temp
	DOUBLE PRECISION :: SZstep
	logical, intent(INOUT) :: relaxd
	INTEGER :: replica_number
	integer :: i 


	if (istep.eq.0) then
	  relaxd=.false.
	  Energy_band_old=Energy_band
	else
	  relaxd=.true.
	  do replica_number=2, replicas-1
	    relaxd=relaxd .and. (abs(Energy_band(replica_number)-Energy_band_old(replica_number)).lt. 0.1)
	  end do
	  Energy_band_old=Energy_band
	end if

	do replica_number=2, replicas-1

!calculo direccion tg
	  tang_vec(1:3,1:na_u) = rclas_BAND(1:3,1:na_u,replica_number) - rclas_BAND(1:3,1:na_u,replica_number-1)
	  do i=1, na_u
	    NORMVEC= tang_vec(1,i)**2 + tang_vec(2,i)**2 + tang_vec(3,i)**2
	    NORMVEC=sqrt(NORMVEC)
	    tang_vec(1:3,i)=tang_vec(1:3,i)/NORMVEC
	  end do
!        write(*,*) "test tangente"
!        write(*,*) "normA"
!
!          do i=1, na_u
!            WRITE(*,*) replica_number,i,tang_vec(1,i)**2 + tang_vec(2,i)**2 + tang_vec(3,i)**2
!          end do
	  
!dejo solo fperpenda
	  MAXFmod=0.d0
	  do i=1, na_u
	    proyForce=fclas_BAND(1,i,replica_number)*tang_vec(1,i)
	    proyForce=proyForce+fclas_BAND(2,i,replica_number)*tang_vec(2,i)
	    proyForce=proyForce+fclas_BAND(3,i,replica_number)*tang_vec(3,i)
	    fclas_BAND(1:3,i,replica_number)=fclas_BAND(1:3,i,replica_number)-proyForce*tang_vec(1:3,i) 
	    MAXFmod_temp=fclas_BAND(1,i,replica_number)**2+fclas_BAND(2,i,replica_number)**2+fclas_BAND(3,i,replica_number)**2
	    if (MAXFmod_temp .gt. MAXFmod) MAXFmod=MAXFmod_temp
	  end do


	  write(*,*) "muevo replica:", replica_number
	  MAXFmod=sqrt(MAXFmod)
	  SZstep=0.01/MAXFmod
          write(*,*) "maxforce", MAXFmod, "stepsize", SZstep

!crit de convergencia
	  rclas_BAND(1:3,1:na_u,replica_number)=rclas_BAND(1:3,1:na_u,replica_number)+SZstep*fclas_BAND(1:3,1:na_u,replica_number)

	end do



	END SUBROUTINE bandmove
