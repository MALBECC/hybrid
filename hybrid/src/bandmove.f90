	SUBROUTINE bandmove(istep, na_u,replicas,rclas_BAND,fclas_BAND, Energy_band, relaxd, ftol, NEB_firstimage, NEB_lastimage)
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
	DOUBLE PRECISION ::  MAXFmod_total, NEB_maxFimage
	DOUBLE PRECISION, save :: MAXFmod_total_old
	DOUBLE PRECISION, save :: SZstep_base
	logical :: NEB_freeze
	logical, dimension(replicas) :: NEB_converged_image
	INTEGER, INTENT(INOUT) :: NEB_firstimage, NEB_lastimage
	integer :: i 


	if (istep.eq.0) then!initial max steep size
	  SZstep_base=0.01d0
	  NEB_firstimage=1
	  NEB_lastimage=replicas
	elseif (mod(istep,10).eq.0) then !recheck convergence on full band every X steps
	  NEB_firstimage=1
	  NEB_lastimage=replicas
	end if
	NEB_converged_image=.true.


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
	do replica_number=NEB_firstimage+1, NEB_lastimage-1
	  do i=1, na_u
	    Fmod2=fclas_BAND(1,i,replica_number)**2+fclas_BAND(2,i,replica_number)**2+fclas_BAND(3,i,replica_number)**2
!	    if (Fmod2.lt.1d-9) write(*,*) "fuerza 0 ", i,replica_number
	  end do
	end do

	relaxd=.true.
	MAXFmod_total=0.d0
	MAX_FORCE_REPLICA=-1
	MAX_FORCE_ATOM=-1
!	if (.not. relaxd) then
	  do replica_number=NEB_firstimage+1, NEB_lastimage-1
!calculo direccion tg
	 call calculate_tg(2,na_u,replicas,replica_number,rclas_BAND,tang_vec, Energy_band)
!remove parallel force
	 call remove_parallel(na_u,replicas, replica_number,tang_vec,fclas_BAND)


!convergence criteria
!	        ftol

	  NEB_maxFimage=0.d0
	  do i=1, na_u
	    Fmod2=fclas_BAND(1,i,replica_number)**2+fclas_BAND(2,i,replica_number)**2+fclas_BAND(3,i,replica_number)**2
	    if (Fmod2 .gt. NEB_maxFimage) NEB_maxFimage=Fmod2
!	    if (Fmod2.lt.1d-9) write(*,*) "fuerza perp 0 ", i,replica_number
	    relaxd=relaxd .and. (Fmod2 .lt. ftol**2)
	    if (Fmod2 .gt. MAXFmod_total) then
	      MAXFmod_total=Fmod2
	      MAX_FORCE_REPLICA=replica_number
	      MAX_FORCE_ATOM=i
	    end if
	  end do

	  if (NEB_maxFimage.gt.ftol**2) NEB_converged_image(replica_number)=.false.

!Spring force
	  call calculate_spring_force(2,na_u, replicas, replica_number,rclas_BAND, tang_vec, F_spring)

	  fclas_BAND(1:3, 1:na_u,replica_number)=fclas_BAND(1:3, 1:na_u,replica_number)  &
          + 0.2*MAXFmod_total*F_spring(1:3,1:na_u)


	  MAXFmod=0.d0
          do i=1, na_u
            Fmod2=fclas_BAND(1,i,replica_number)**2+fclas_BAND(2,i,replica_number)**2+fclas_BAND(3,i,replica_number)**2
            if (Fmod2 .gt. MAXFmod)  MAXFmod=Fmod2
          end do


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
	  if (NEB_firstimage.eq.1 .and. NEB_lastimage.eq.replicas) then
	    write(*,*) "system converged"
	  else
	    NEB_firstimage=1
	    NEB_lastimage=replicas
	    relaxd=.false.
	    write(*,*) "system converged on frozen band checking full band"
	  end if
	else
	  write(*,*) "system NOT converged"
 
!freeze converged band
	  NEB_freeze=.true.
	  replica_number=0
	  do while (replica_number .le. replicas .and. NEB_freeze)
	    replica_number=replica_number+1
	    NEB_freeze=NEB_freeze.and.NEB_converged_image(replica_number)
	  end do
!	  NEB_firstimage=replica_number

          NEB_freeze=.true.
          replica_number=replicas+1
          do while (replica_number .ge. 1 .and. NEB_freeze)
            replica_number=replica_number-1
            NEB_freeze=NEB_freeze.and.NEB_converged_image(replica_number)
          end do
!	  NEB_lastimage=replica_number
	  write(*,*) "freezing fists ", NEB_firstimage, "images and lasts ", NEB_lastimage, " images"

	end if
	write(*,*) "max force ", sqrt(MAXFmod_total), "on atom ", MAX_FORCE_ATOM,  &
        "in replica ", MAX_FORCE_REPLICA, "conv criteria", ftol

	if (SZstep_base.lt.1d-8) then
	  relaxd=.true.
	  write(*,*) "max precision reached on atomic displacement"
	end if
	END SUBROUTINE bandmove


	SUBROUTINE calculate_tg(method,na_u,replicas,replica_number,rclas_BAND,tang_vec, Energy_band)
	IMPLICIT NONE
        INTEGER, INTENT(IN) :: method,na_u,replicas, replica_number
        DOUBLE PRECISION, DIMENSION(3,na_u,replicas), INTENT(IN) :: rclas_BAND
        DOUBLE PRECISION, DIMENSION(3,na_u), INTENT(INOUT) :: tang_vec
	DOUBLE PRECISION, DIMENSION(3,na_u) :: tang_vecA, tang_vecB
	DOUBLE PRECISION :: NORMVEC, NORMVECA, NORMVECB
	DOUBLE PRECISION, DIMENSION(replicas), intent(In) :: Energy_band
	DOUBLE PRECISION :: E0, E1, E2, Vmax, Vmin

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

	ELSEIF (method.eq.2) then
!The Journal of Chemical Physics 113, 9978 (2000); https://doi.org/10.1063/1.1323224
	  tang_vecA(1:3,1:na_u) = rclas_BAND(1:3,1:na_u,replica_number) - rclas_BAND(1:3,1:na_u,replica_number-1)
	  tang_vecB(1:3,1:na_u) = rclas_BAND(1:3,1:na_u,replica_number+1) - rclas_BAND(1:3,1:na_u,replica_number)
	  E0=Energy_band(replica_number-1)
	  E1=Energy_band(replica_number)
	  E2=Energy_band(replica_number+1)
	  do i=1, na_u
	    if ((E2.gt.E1) .and. (E1.gt.E0)) then
	      tang_vec(1:3,i)=tang_vecB(1:3,i)
	    else if ((E0.gt.E1) .and. (E1.gt.E2)) then
	      tang_vec(1:3,i)=tang_vecA(1:3,i)
	    else
	      Vmax=max(abs(E2-E1), abs(E1-E0))
	      Vmin=min(abs(E2-E1), abs(E1-E0))
	      if (E2.gt.E0) tang_vec(1:3,i)=Vmax*tang_vecB(1:3,i)+Vmin*tang_vecA(1:3,i)
	      if (E0.gt.E2) tang_vec(1:3,i)=Vmin*tang_vecB(1:3,i)+Vmax*tang_vecA(1:3,i)
	    end if
	  end do
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
	  proyForce=proyForce+fclas_BAND(2,i,replica_number)*tang_vec(2,i)
	  proyForce=proyForce+fclas_BAND(3,i,replica_number)*tang_vec(3,i)
	  fclas_BAND(1:3,i,replica_number)=fclas_BAND(1:3,i,replica_number)-proyForce*tang_vec(1:3,i) 
        end do

	RETURN
	END SUBROUTINE remove_parallel


	SUBROUTINE calculate_spring_force(methodSF, na_u, replicas, replica_number,rclas_BAND, tang_vec, F_spring)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: na_u,replicas
        DOUBLE PRECISION, DIMENSION(3,na_u), INTENT(INOUT) :: F_spring
        DOUBLE PRECISION, DIMENSION(3,na_u,replicas), INTENT(IN) :: rclas_BAND
        DOUBLE PRECISION, DIMENSION(3,na_u), INTENT(IN)  :: tang_vec
	DOUBLE PRECISION, DIMENSION(3) :: auxvector
	DOUBLE PRECISION :: auxescalar, auxescalarA, auxescalarB
        INTEGER :: replica_number
	INTEGER, INTENT(IN) :: methodSF
        integer :: i

	F_spring=0.d0
	do i=1, na_u
	  IF (methodSF.eq.1) then
	    auxvector(1:3)=rclas_BAND(1:3, i, replica_number+1)+rclas_BAND(1:3, i, replica_number-1) - 2*rclas_BAND(1:3, i, replica_number) 
	    auxescalar=auxvector(1)*tang_vec(1,i)+auxvector(2)*tang_vec(2,i)+auxvector(3)*tang_vec(3,i)
	  ELSEIF (methodSF.eq.2) then
	    auxvector(1:3)=rclas_BAND(1:3, i, replica_number+1)-rclas_BAND(1:3, i, replica_number)
	    auxescalarA=auxvector(1)**2+auxvector(2)**2+auxvector(3)**2
	    auxvector(1:3)=rclas_BAND(1:3, i, replica_number)-rclas_BAND(1:3, i, replica_number-1)
	    auxescalarB=auxvector(1)**2+auxvector(2)**2+auxvector(3)**2
	    auxescalar=sqrt(auxescalarA)-sqrt(auxescalarB)
	  ELSE
	    STOP "wrong methodSF"
	  END IF
	    F_spring(1:3,i)=auxescalar*tang_vec(1:3,i)
	end do

	RETURN
	END SUBROUTINE calculate_spring_force
