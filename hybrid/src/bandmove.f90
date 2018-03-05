	SUBROUTINE bandmove(istep, na_u,replicas,rclas_BAND,fclas_BAND, Energy_band, relaxd, ftol)
!Steepest descend method for nudged elastic band.
!N. Foglia 03/2018
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: na_u,replicas, istep
	DOUBLE PRECISION, DIMENSION(3,na_u,replicas), INTENT(INOUT) :: fclas_BAND
	DOUBLE PRECISION, DIMENSION(3,na_u) :: F_spring
	DOUBLE PRECISION, DIMENSION(3,na_u,replicas), INTENT(INOUT) :: rclas_BAND
	DOUBLE PRECISION, DIMENSION(3,na_u) :: tang_vec
	DOUBLE PRECISION, DIMENSION(replicas), intent(In) :: Energy_band
	DOUBLE PRECISION, DIMENSION(replicas) :: Energy_band_old
	DOUBLE PRECISION :: NORMVEC, proyForce, MAXFmod, MAXFmod_temp
	DOUBLE PRECISION :: SZstep
	DOUBLE PRECISION :: Fmod2
	DOUBLE PRECISION, INTENT(IN) :: ftol
	logical, intent(INOUT) :: relaxd
	INTEGER :: replica_number
	INTEGER :: MAX_FORCE_REPLICA, MAX_FORCE_ATOM
	DOUBLE PRECISION ::  MAXFmod_total
	DOUBLE PRECISION, save :: MAXFmod_total_old
	DOUBLE PRECISION, save :: SZstep_base
	integer :: i 

!	write(*,*) "mievo paso ", istep
	if (istep.eq.0) then
	  SZstep_base=0.01d0
	end if

!	if (istep.eq.0) then
!	  relaxd=.false.
!	  Energy_band_old=Energy_band
!	else
!convergence criteria
!	  relaxd=.true.
!	  do replica_number=2, replicas-1
!	    relaxd=relaxd .and. (abs(Energy_band(replica_number)-Energy_band_old(replica_number)).lt. 0.04)
!	  end do
!	  Energy_band_old=Energy_band
!	end if
!check forces, sacar luego Nick
	do replica_number=2, replicas-1
	  do i=1, na_u
	    Fmod2=fclas_BAND(1,i,replica_number)**2+fclas_BAND(2,i,replica_number)**2+fclas_BAND(3,i,replica_number)**2
!	    if (Fmod2.lt.1d-9) write(*,*) "fuerza 0 ", i,replica_number
	  end do
	end do

	relaxd=.true.
	MAXFmod_total=0.d0
	MAX_FORCE_REPLICA=-1
	MAX_FORCE_ATOM=-1
!	write(*,*) "flag 111"
!	if (.not. relaxd) then
	  do replica_number=2, replicas-1
!calculo direccion tg
!	write(*,*) "flag 112"

!		if (replica_number.eq. 2) write(*,*) "F221 init", fclas_BAND


	  call calculate_tg(1,na_u,replicas,replica_number,rclas_BAND,tang_vec)
!	write(*,*) "flag 113"

!remove parallel force
	 call remove_parallel(na_u,replicas, replica_number,tang_vec,fclas_BAND)
!	write(*,*) "flag 114"

!		if (replica_number.eq. 2) write(*,*) "F221 2nd", fclas_BAND



!convergence criteria
!	        ftol
	  do i=1, na_u
	    Fmod2=fclas_BAND(1,i,replica_number)**2+fclas_BAND(2,i,replica_number)**2+fclas_BAND(3,i,replica_number)**2
!	    if (Fmod2.lt.1d-9) write(*,*) "fuerza perp 0 ", i,replica_number
	    relaxd=relaxd .and. (Fmod2 .lt. ftol**2)
	    if (Fmod2 .gt. MAXFmod_total) then
	      MAXFmod_total=Fmod2
	      MAX_FORCE_REPLICA=replica_number
	      MAX_FORCE_ATOM=i
	    end if
	  end do


!		if (replica_number.eq. 2) write(*,*) "F221 3er", fclas_BAND


!Spring force
	  call calculate_spring_force(na_u, replicas, replica_number,rclas_BAND, tang_vec, F_spring)

	  fclas_BAND(1:3, 1:na_u,replica_number)=fclas_BAND(1:3, 1:na_u,replica_number)  &
          + 0.2*MAXFmod_total*F_spring(1:3,1:na_u)


!		if (replica_number.eq. 2) write(*,*) "F221 4th", fclas_BAND

	  MAXFmod=0.d0
          do i=1, na_u
            Fmod2=fclas_BAND(1,i,replica_number)**2+fclas_BAND(2,i,replica_number)**2+fclas_BAND(3,i,replica_number)**2
            if (Fmod2 .gt. MAXFmod)  MAXFmod=Fmod2
          end do


!		if (replica_number.eq. 2) write(*,*) "F221 5th", fclas_BAND

!	  MAXFmod=0.d0
!	  do i=1, na_u
!	    MAXFmod_temp=fclas_BAND(1,i,replica_number)**2+fclas_BAND(2,i,replica_number)**2+fclas_BAND(3,i,replica_number)**2
!	    if (MAXFmod_temp .gt. MAXFmod) MAXFmod=MAXFmod_temp
!	  end do
!	  write(*,*) "max force: ", MAXFmod, "comb crit", ftol
!	  write(*,*) "flag 115"


!	  write(*,*) "muevo replica:", replica_number
	  MAXFmod=sqrt(MAXFmod)
	  SZstep=SZstep_base/MAXFmod
          write(*,*) "maxforce", MAXFmod, "stepsize", SZstep, &
          "stepsize_base", SZstep_base

	  rclas_BAND(1:3,1:na_u,replica_number)=rclas_BAND(1:3,1:na_u,replica_number)+SZstep*fclas_BAND(1:3,1:na_u,replica_number)

	end do
	
        if (istep.eq.0) MAXFmod_total_old=MAXFmod_total
	if (MAXFmod_total .gt. MAXFmod_total_old) SZstep_base=SZstep_base*0.85d0
	MAXFmod_total_old=MAXFmod_total


	if (relaxd) then
	  write(*,*) "system converged"
	else
	  write(*,*) "system NOT converged"  
	end if
	write(*,*) "max force ", sqrt(MAXFmod_total), "on atom ", MAX_FORCE_ATOM,  &
        "in replica ", MAX_FORCE_REPLICA, "conv criteria", ftol
!	end if

	if (SZstep_base.lt.1d-8) then
	  relaxd=.true.
	  write(*,*) "max precision reached on atomic displacement"
	end if

!	write(*,*) "flag 11end"

	END SUBROUTINE bandmove


	SUBROUTINE calculate_tg(method,na_u,replicas,replica_number,rclas_BAND,tang_vec)
	IMPLICIT NONE
        INTEGER, INTENT(IN) :: method,na_u,replicas, replica_number
        DOUBLE PRECISION, DIMENSION(3,na_u,replicas), INTENT(IN) :: rclas_BAND
        DOUBLE PRECISION, DIMENSION(3,na_u), INTENT(INOUT) :: tang_vec
	DOUBLE PRECISION, DIMENSION(3,na_u) :: tang_vecA, tang_vecB
	DOUBLE PRECISION :: NORMVEC, NORMVECA, NORMVECB
	integer :: i

	IF (method.eq.0) then
          tang_vec(1:3,1:na_u) = rclas_BAND(1:3,1:na_u,replica_number) - rclas_BAND(1:3,1:na_u,replica_number-1)
	ELSEIF (method.eq.1) then
	  tang_vecA(1:3,1:na_u) = rclas_BAND(1:3,1:na_u,replica_number) - rclas_BAND(1:3,1:na_u,replica_number-1)
	  tang_vecB(1:3,1:na_u) = rclas_BAND(1:3,1:na_u,replica_number+1) - rclas_BAND(1:3,1:na_u,replica_number)

          do i=1, na_u
            NORMVECA= tang_vecA(1,i)**2 + tang_vecA(2,i)**2 + tang_vecA(3,i)**2
            NORMVECA=sqrt(NORMVECA)
            tang_vecA(1:3,i)=tang_vecA(1:3,i)/NORMVECA

            NORMVECB= tang_vecB(1,i)**2 + tang_vecB(2,i)**2 + tang_vecB(3,i)**2
            NORMVECB=sqrt(NORMVECB)
            tang_vecB(1:3,i)=tang_vecB(1:3,i)/NORMVECB
          end do

	   tang_vec(1:3,1:na_u)=tang_vecA(1:3,1:na_u)+tang_vecB(1:3,1:na_u)

	ELSE
	  STOP "Wrong method in calculate_tg"
	END IF


        do i=1, na_u
          NORMVEC= tang_vec(1,i)**2 + tang_vec(2,i)**2 + tang_vec(3,i)**2
          NORMVEC=sqrt(NORMVEC)
	  if (NORMVEC .lt. 1d-9) then
	write(*,*) "tangent vector null, replica: ", replica_number, "atom ", i
	    stop
	  else if (NORMVEC .ne. NORMVEC) then
	    stop "NAN in tangent vector"
	  end if
          tang_vec(1:3,i)=tang_vec(1:3,i)/NORMVEC
        end do

	RETURN
	END SUBROUTINE calculate_tg

	SUBROUTINE remove_parallel(na_u,replicas, replica_number,tang_vec,fclas_BAND) 
	IMPLICIT NONE
        INTEGER, INTENT(IN) :: na_u,replicas
        DOUBLE PRECISION, DIMENSION(3,na_u,replicas), INTENT(INOUT) :: fclas_BAND
        DOUBLE PRECISION, DIMENSION(3,na_u), INTENT(IN) :: tang_vec
        DOUBLE PRECISION ::  proyForce 
        INTEGER :: replica_number
        integer :: i

        do i=1, na_u
          proyForce=fclas_BAND(1,i,replica_number)*tang_vec(1,i)
!	  if (proyForce .ne. proyForce) then
!	    write(*,*) "proyForce x", proyForce, "i, replica_number", i, replica_number
!	    write(*,*) "fclas_BAND(1,i,replica_number), tang_vec(1,i)", fclas_BAND(1,i,replica_number), tang_vec(1,i)
!	  end if

	  proyForce=proyForce+fclas_BAND(2,i,replica_number)*tang_vec(2,i)
!	   if (proyForce .ne. proyForce) then
!	  write(*,*) "proyForce y", proyForce, "i, replica_number", i, replica_number
!	     write(*,*) "fclas_BAND(2,i,replica_number), tang_vec(2,i)", fclas_BAND(2,i,replica_number), tang_vec(2,i)
!	   end if
	  proyForce=proyForce+fclas_BAND(3,i,replica_number)*tang_vec(3,i)
!	   if (proyForce .ne. proyForce) then
!	  write(*,*) "proyForce ", proyForce, "i, replica_number", i, replica_number
!	    write(*,*) "fclas_BAND(3,i,replica_number), tang_vec(3,i)", fclas_BAND(3,i,replica_number), tang_vec(3,i)
!	   end if
	  fclas_BAND(1:3,i,replica_number)=fclas_BAND(1:3,i,replica_number)-proyForce*tang_vec(1:3,i) 
        end do

	RETURN
	END SUBROUTINE remove_parallel


	SUBROUTINE calculate_spring_force(na_u, replicas, replica_number,rclas_BAND, tang_vec, F_spring)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: na_u,replicas
        DOUBLE PRECISION, DIMENSION(3,na_u), INTENT(INOUT) :: F_spring
        DOUBLE PRECISION, DIMENSION(3,na_u,replicas), INTENT(IN) :: rclas_BAND
        DOUBLE PRECISION, DIMENSION(3,na_u), INTENT(IN)  :: tang_vec
	DOUBLE PRECISION, DIMENSION(3) :: auxvector
	DOUBLE PRECISION :: auxescalar
        INTEGER :: replica_number
        integer :: i

	F_spring=0.d0
	do i=1, na_u
	  auxvector(1:3)=rclas_BAND(1:3, i, replica_number+1)+rclas_BAND(1:3, i, replica_number-1) - 2*rclas_BAND(1:3, i, replica_number) 
	  auxescalar=auxvector(1)*tang_vec(1,i)+auxvector(2)*tang_vec(2,i)+auxvector(3)*tang_vec(3,i)
	  F_spring(1:3,i)=auxescalar*tang_vec(1:3,i)
	end do

	RETURN
	END SUBROUTINE calculate_spring_force
