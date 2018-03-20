	SUBROUTINE NEB_steep(istep, relaxd)
!Nudged elastic band method
!movement methods avalables are: Steepest descend and quick-min and FIRE using velocity verlet 
!N. Foglia 03/2018
	use scarlett, only: natot, NEB_Nimages, masst, NEB_firstimage, NEB_lastimage, rclas_BAND, vclas_BAND, fclas_BAND, aclas_BAND_old, Energy_band,NEB_move_method, NEB_spring_constant, ftol, NEB_steep_size, NEB_MAXFmod
	implicit none
	integer, intent(in) :: istep
	double precision, dimension(3,natot) :: F_spring
	double precision, dimension(3,natot) :: tang_vec
	double precision, dimension(NEB_Nimages) :: Energy_band_old
	double precision :: NORMVEC, proyForce, MAXFmod, MAXFmod_temp
!	double precision :: SZstep
	double precision :: Fmod2
	logical, intent(inout) :: relaxd
	integer :: replica_number
	integer :: MAX_FORCE_REPLICA, MAX_FORCE_ATOM
	double precision ::  MAXFmod_total, NEB_Ekin
!	double precision, save :: MAXFmod_total_old
	logical :: NEB_freeze
	logical, dimension(NEB_Nimages) :: NEB_converged_image
	integer :: i 
	
	if (istep.eq.0) then!initial max steep size
	  NEB_steep_size=0.1d0
	  NEB_firstimage=1
	  NEB_lastimage=NEB_Nimages
	elseif (mod(istep,10).eq.0) then !recheck convergence on full band every X steps
	  NEB_firstimage=1
	  NEB_lastimage=NEB_Nimages
	end if
	NEB_converged_image=.true.
	
	relaxd=.true.
	MAXFmod_total=0.d0
	MAX_FORCE_REPLICA=-1
	MAX_FORCE_ATOM=-1
	
	do replica_number=NEB_firstimage+1, NEB_lastimage-1
	  call NEB_calculate_tg(2,replica_number,tang_vec) !calculate tangent direccion 
	  call NEB_remove_parallel(replica_number,tang_vec) !remove force in tangent direction
	  call NEB_check_convergence(relaxd, replica_number, MAXFmod_total, MAX_FORCE_REPLICA, MAX_FORCE_ATOM, NEB_converged_image)
	  call NEB_calculate_spring_force(2, replica_number, tang_vec, F_spring) !Spring force
	  fclas_BAND(1:3, 1:natot,replica_number)=fclas_BAND(1:3, 1:natot,replica_number)  &
	  + NEB_spring_constant*F_spring(1:3,1:natot) !New force
	end do
	
	if (.not. relaxd) then
	  call NEB_movement_algorithm(NEB_move_method,istep, MAXFmod_total) !move systems
!	  if (istep.eq.0) MAXFmod_total_old=MAXFmod_total
!	  if (MAXFmod_total .gt. MAXFmod_total_old) NEB_steep_size=NEB_steep_size*0.85d0
!	  MAXFmod_total_old=MAXFmod_total
	  write(*,*) "system NOT converged"
	else
	  write(*,*) "system converged"
	end if
	
	if (NEB_move_method .eq. 1) THEN
	  if (istep.eq.0) NEB_MAXFmod=MAXFmod_total
	  if (MAXFmod_total .gt. NEB_MAXFmod) NEB_steep_size=NEB_steep_size*0.85d0
	  NEB_MAXFmod=MAXFmod_total
	
	  if (NEB_steep_size.lt.1d-5) then
	    relaxd=.true.
	    write(*,*) "max precision reached on atomic displacement"
	  end if
	
	else if (NEB_move_method .eq. 2) THEN
	  CALL NEB_calculate_T(NEB_Ekin)
	  WRITE(*,*) "TOTAL K ", NEB_Ekin, "AVERAGE K", NEB_Ekin/(dble(NEB_Nimages-2))
	END IF
	write(*,*) "max force ", sqrt(MAXFmod_total), "on atom ", MAX_FORCE_ATOM,  &
        "in replica ", MAX_FORCE_REPLICA, "conv criteria", ftol
	END SUBROUTINE NEB_steep


	SUBROUTINE NEB_calculate_tg(method,replica_number,tang_vec)
	use scarlett, only: NEB_Nimages, natot, rclas_BAND, Energy_band
	implicit none
	integer, intent(in) :: method, replica_number
	double precision, dimension(3,natot), intent(inout) :: tang_vec
	double precision, dimension(3,natot) :: tang_vecA, tang_vecB
	double precision :: NORMVEC, NORMVECA, NORMVECB
	double precision :: E0, E1, E2, Vmax, Vmin
	integer :: i
	
	IF (method.eq.0) then
	  tang_vec(1:3,1:natot) = rclas_BAND(1:3,1:natot,replica_number) - rclas_BAND(1:3,1:natot,replica_number-1)
	ELSEIF (method.eq.1) then
	  tang_vecA(1:3,1:natot) = rclas_BAND(1:3,1:natot,replica_number) - rclas_BAND(1:3,1:natot,replica_number-1)
	  tang_vecB(1:3,1:natot) = rclas_BAND(1:3,1:natot,replica_number+1) - rclas_BAND(1:3,1:natot,replica_number)
	
	  do i=1, natot
	    NORMVECA= tang_vecA(1,i)**2 + tang_vecA(2,i)**2 + tang_vecA(3,i)**2
	    NORMVECA=sqrt(NORMVECA)
	    tang_vecA(1:3,i)=tang_vecA(1:3,i)/NORMVECA
	
	    NORMVECB= tang_vecB(1,i)**2 + tang_vecB(2,i)**2 + tang_vecB(3,i)**2
	    NORMVECB=sqrt(NORMVECB)
	    tang_vecB(1:3,i)=tang_vecB(1:3,i)/NORMVECB
	  end do
	
	  tang_vec(1:3,1:natot)=tang_vecA(1:3,1:natot)+tang_vecB(1:3,1:natot)
	
	ELSEIF (method.eq.2) then !The Journal of Chemical Physics 113, 9978 (2000); https://doi.org/10.1063/1.1323224
	  tang_vecA(1:3,1:natot) = rclas_BAND(1:3,1:natot,replica_number) - rclas_BAND(1:3,1:natot,replica_number-1)
	  tang_vecB(1:3,1:natot) = rclas_BAND(1:3,1:natot,replica_number+1) - rclas_BAND(1:3,1:natot,replica_number)
	  E0=Energy_band(replica_number-1)
	  E1=Energy_band(replica_number)
	  E2=Energy_band(replica_number+1)
	  do i=1, natot
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
	  STOP "Wrong method in NEB_calculate_tg"
	END IF
	
	do i=1, natot
	  NORMVEC= tang_vec(1,i)**2 + tang_vec(2,i)**2 + tang_vec(3,i)**2
	  NORMVEC=sqrt(NORMVEC)
	  if (NORMVEC .lt. 1d-9) then
	    write(*,*) "tangent vector null, replica: ",replica_number, "atom ", i
	    stop
	  else if (NORMVEC .ne. NORMVEC) then
	    stop "NAN in tangent vector"
	  end if
	  tang_vec(1:3,i)=tang_vec(1:3,i)/NORMVEC
	end do
	
	RETURN
	END SUBROUTINE NEB_calculate_tg

	SUBROUTINE NEB_remove_parallel(replica_number,tang_vec) 
	use scarlett, only: NEB_Nimages, natot, fclas_BAND
	implicit none
	double precision, dimension(3,natot), intent(in) :: tang_vec
	double precision ::  proyForce 
	integer :: replica_number
	integer :: i
	
	do i=1, natot
	  proyForce=fclas_BAND(1,i,replica_number)*tang_vec(1,i)
	  proyForce=proyForce+fclas_BAND(2,i,replica_number)*tang_vec(2,i)
	  proyForce=proyForce+fclas_BAND(3,i,replica_number)*tang_vec(3,i)
	  fclas_BAND(1:3,i,replica_number)=fclas_BAND(1:3,i,replica_number)-proyForce*tang_vec(1:3,i) 
	end do
	return
	END SUBROUTINE NEB_remove_parallel


	SUBROUTINE NEB_calculate_spring_force(methodSF, replica_number, tang_vec, F_spring)
	use scarlett, only: NEB_Nimages, NEB_Nimages, natot, rclas_BAND
	implicit none
	double precision, dimension(3,natot), intent(inout) :: F_spring
	double precision, dimension(3,natot), intent(in)  :: tang_vec
	double precision, dimension(3) :: auxvector
	double precision :: auxescalar, auxescalarA, auxescalarB
	integer :: replica_number
	integer, intent(in) :: methodSF
	integer :: i
	
	F_spring=0.d0
	do i=1, natot
	  if (methodSF.eq.1) then
	    auxvector(1:3)=rclas_BAND(1:3, i, replica_number+1)+rclas_BAND(1:3, i, replica_number-1) - 2*rclas_BAND(1:3, i, replica_number) 
	    auxescalar=auxvector(1)*tang_vec(1,i)+auxvector(2)*tang_vec(2,i)+auxvector(3)*tang_vec(3,i)
	  elseif (methodSF.eq.2) then
	    auxvector(1:3)=rclas_BAND(1:3, i, replica_number+1)-rclas_BAND(1:3, i, replica_number)
	    auxescalarA=auxvector(1)**2+auxvector(2)**2+auxvector(3)**2
	    auxvector(1:3)=rclas_BAND(1:3, i, replica_number)-rclas_BAND(1:3, i, replica_number-1)
	    auxescalarB=auxvector(1)**2+auxvector(2)**2+auxvector(3)**2
	    auxescalar=sqrt(auxescalarA)-sqrt(auxescalarB)
	  else
	    stop "wrong methodSF"
	  endif
	    F_spring(1:3,i)=auxescalar*tang_vec(1:3,i)
	end do
	return
	END SUBROUTINE NEB_calculate_spring_force


	SUBROUTINE NEB_movement_algorithm(method,istep, MAXFmod_total)
	use scarlett, only: aclas_BAND_old, NEB_Nimages, natot,rclas_BAND,vclas_BAND,fclas_BAND, masst, NEB_Ndescend, time_steep, time_steep_max, NEB_alpha, NEB_steep_size
	implicit none
	integer :: method
	double precision, dimension(3,natot,NEB_Nimages) :: v12clas_BAND
	double precision, dimension(3,natot,NEB_Nimages) :: aclas_BAND
	double precision, intent(in) :: MAXFmod_total
	double precision :: MAXFmod
	double precision :: SZstep
	integer, intent(in) :: istep
	integer :: i, j, replica_number
	double precision :: Fmod, velocity_proyected, velocity_mod
	
	if (istep .eq. 1 ) time_steep=1.d-1 !cambiar esto luego, por ahora valor arbitrario
	Fmod=0.d0
	
	if (method.eq.1) then !steepest descend
	
	  MAXFmod=sqrt(MAXFmod_total)
	  SZstep=NEB_steep_size/MAXFmod
	  write(*,*) "maxforce", MAXFmod, "stepsize", SZstep, &
	  "stepsize_base", NEB_steep_size
	  rclas_BAND(1:3,1:natot,1:NEB_Nimages)=rclas_BAND(1:3,1:natot,1:NEB_Nimages)+SZstep*fclas_BAND(1:3,1:natot,1:NEB_Nimages)
	
	
	elseif (method.eq.2) then !quick-min using velocity verlet
	
	  do i=1, natot
	    aclas_BAND(1:3, i, 1:NEB_Nimages) = fclas_BAND(1:3, i, 1:NEB_Nimages)/masst(i)
	    do j=1,3
	      do replica_number=1,NEB_Nimages
	        Fmod=Fmod + fclas_BAND(j, i, replica_number)**2
	      end do
	    end do
	  end do
	  Fmod=sqrt(Fmod)
	  if (istep .ne. 1) then
	    vclas_BAND=vclas_BAND+0.5d0*(aclas_BAND+aclas_BAND_old)*time_steep
	    velocity_proyected=0.d0      
	    do i=1, natot
	      do j=1,3
	        do replica_number=1,NEB_Nimages
	          velocity_proyected=velocity_proyected+vclas_BAND(j,i,replica_number) * fclas_BAND(j,i,replica_number)
	        end do
	      end do
	    end do
	    velocity_proyected=velocity_proyected/Fmod
	
	    if (velocity_proyected .gt. 0.d0) then
	      vclas_BAND=velocity_proyected*fclas_BAND/Fmod
	    else
	      vclas_BAND=0.3d0*velocity_proyected*fclas_BAND/Fmod !0.3 fue arbitrario. hay q liberarlo
	    end if
	
	  end if
	!freezing 1st and last images 
	  vclas_BAND(1:3, 1:natot,1)=0.d0
	  vclas_BAND(1:3, 1:natot,NEB_Nimages)=0.d0
	  aclas_BAND(1:3, 1:natot,1)=0.d0
	  aclas_BAND(1:3, 1:natot,NEB_Nimages)=0.d0
	!move images
    	  rclas_BAND=rclas_BAND+vclas_BAND*time_steep+0.5d0*aclas_BAND*time_steep**2
	  aclas_BAND_old=aclas_BAND
	
	
	elseif (method.eq.3) then !FIRE, Bitzek, et. al., Phys. Rev. Lett. 97, 170201 (2006).
	
	  do i=1, natot
	    aclas_BAND(1:3, i, 1:NEB_Nimages) = fclas_BAND(1:3, i, 1:NEB_Nimages)/masst(i)
	    do j=1,3
	      do replica_number=1,NEB_Nimages
	        Fmod=Fmod + fclas_BAND(j, i, replica_number)**2
	      end do
	    end do
	  end do
	  Fmod=sqrt(Fmod)
	
	  if (istep .ne. 1) then
	    vclas_BAND=vclas_BAND+0.5d0*(aclas_BAND+aclas_BAND_old)*time_steep
	    velocity_proyected=0.d0 !P in paper
	    velocity_mod=0.d0
	
	    do i=1, natot
	      do j=1,3
		do replica_number=1,NEB_Nimages
		  velocity_proyected=velocity_proyected+vclas_BAND(j,i,replica_number) * fclas_BAND(j,i,replica_number)
		  velocity_mod=velocity_mod+vclas_BAND(j,i,replica_number)*vclas_BAND(j,i,replica_number)
		end do
	      end do
	    end do
	    velocity_mod=sqrt(velocity_mod)
	
	    if (velocity_proyected .gt. 0.d0) then
	      NEB_Ndescend=NEB_Ndescend+1
	      if (NEB_Ndescend .gt. 5) then
	        time_steep=min(time_steep*1.1d0, time_steep_max) 
	        NEB_alpha=NEB_alpha*0.99d0
	      end if
	      vclas_BAND=(1.d0-NEB_alpha)*vclas_BAND+NEB_alpha*velocity_mod*1.d0/Fmod*fclas_BAND
	    else
	      NEB_Ndescend=0
	      NEB_alpha=0.1d0
	      time_steep=time_steep*0.5
	      vclas_BAND=0.d0
	    end if
	  else !initialice variables in 1st steep
	    NEB_Ndescend=0
	    time_steep_max=10.d0*time_steep
	    NEB_alpha=0.1d0
	  end if
	!freezing 1st and last images 
	  vclas_BAND(1:3, 1:natot,1)=0.d0
	  vclas_BAND(1:3, 1:natot,NEB_Nimages)=0.d0
	  aclas_BAND(1:3, 1:natot,1)=0.d0
	  aclas_BAND(1:3, 1:natot,NEB_Nimages)=0.d0
	!move images
	  rclas_BAND=rclas_BAND+vclas_BAND*time_steep+0.5d0*aclas_BAND*time_steep**2
	  aclas_BAND_old=aclas_BAND
	
	
	else
	
	  STOP "Wrong method in NEB_movement_algorithm"
	end if
	END SUBROUTINE NEB_movement_algorithm


	SUBROUTINE NEB_check_convergence(relaxd, replica_number, MAXFmod_total, MAX_FORCE_REPLICA, MAX_FORCE_ATOM, NEB_converged_image)
	use scarlett, only: natot, fclas_BAND, ftol, NEB_Nimages
	implicit none
	integer :: i,j !auxiliar
	integer, intent(in) :: replica_number
	integer, intent(inout) :: MAX_FORCE_REPLICA, MAX_FORCE_ATOM
	logical, intent(inout) :: relaxd
	double precision :: NEB_maxFimage, Fmod2
	double precision, intent(inout) :: MAXFmod_total
	logical, dimension(NEB_Nimages), intent(inout) :: NEB_converged_image
	
	NEB_maxFimage=0.d0
	
	do i=1, natot
	  do j=1,3
	    Fmod2=fclas_BAND(j,i,replica_number)**2
!a+fclas_BAND(2,i,replica_number)**2+fclas_BAND(3,i,replica_number)**2
	
	    if (Fmod2 .gt. NEB_maxFimage) NEB_maxFimage=Fmod2
	    relaxd=relaxd .and. (Fmod2 .lt. ftol**2)
	
	    if (Fmod2 .gt. MAXFmod_total) then
	      MAXFmod_total=Fmod2
	      MAX_FORCE_REPLICA=replica_number
	      MAX_FORCE_ATOM=i
	    end if
	  end do
	end do
	
	if (NEB_maxFimage.gt.ftol**2) NEB_converged_image(replica_number)=.false.
	return
	END SUBROUTINE NEB_check_convergence


	SUBROUTINE NEB_calculate_T(NEB_Ekin)
!calculate kinectic energy of band
	use scarlett, only: natot, NEB_Nimages, vclas_BAND, masst
	implicit none
	double precision, intent(inout) :: NEB_Ekin
	double precision :: vel2
	integer :: replica_number, i, j
	
	NEB_Ekin=0.d0
	do replica_number=1, NEB_Nimages
	  do i=1, natot
	    vel2=0.d0
	    do j=1, 3
	      vel2=vel2+vclas_BAND(j,i, replica_number)**2
	    end do
	    NEB_Ekin=NEB_Ekin+0.5d0*masst(i)*vel2
	  end do
	end do
	END SUBROUTINE NEB_calculate_T



	subroutine NEB_make_initial_band()
!Generate initial configurations of images in NEB using, .XV.i restarts or
!using only reactives and products restarts (and TS if posible)
!N. Foglia 03/2018
	use scarlett, only: natot, ucell, rclas, vat, rclas_BAND, NEB_Nimages
	implicit none
	logical :: band_xv_found
	logical :: foundxv, foundvat !control for coordinates restart
	integer :: replica_number !auxiliar
	double precision :: BAND_slope(3), BAND_const(3)
	integer :: i, k !auxiliars
	integer :: middle_point
	
	band_xv_found=.true.
	do replica_number = 1, NEB_Nimages
	  call ioxv( 'read', natot, ucell, rclas, vat, foundxv, foundvat,'X',replica_number)
	  band_xv_found=band_xv_found .and. foundxv
	  if (foundxv) then
	    rclas_BAND(1:3,1:natot,replica_number)=rclas(1:3,1:natot)
	  end if
	end do
	
	if (band_xv_found) then !restart case
	  write(*,*) "used .XV. restarts"
	else !not restart case
	  write(*,*) "didnt found all necesary restarts .XV.i"
	  write(*,*) "with i between 1 and ", NEB_Nimages
	  write(*,*) "using .XVR and .XVP for generate initial states"
	 
!read reactive coordinates
	  call ioxv('read',natot,ucell,rclas,vat, foundxv, foundvat,'R',-1)
	  if (foundxv) then
	    rclas_BAND(1:3,1:natot,1)=rclas(1:3,1:natot)
	  else
	    stop ".XVR not found"
	  end if
	
!transition state
	  call ioxv('read',natot,ucell,rclas,vat,foundxv,foundvat,'T',-1)
	  if (foundxv) then
	    write(*,*) "using  TS"
	    middle_point=1+NEB_Nimages
	    middle_point=middle_point/2
	    rclas_BAND(1:3,1:natot,middle_point)=rclas(1:3,1:natot)
	  else
	    write(*,*) "not using  TS, initial band will be created", &
	   "interpolating reactives and products"
	    middle_point=NEB_Nimages
	  end if
	
!read products coordinates
	  call ioxv('read',natot,ucell,rclas,vat,foundxv,foundvat,'P',-1)
	  if (foundxv) then
	    rclas_BAND(1:3,1:natot,NEB_Nimages)=rclas(1:3,1:natot)
	  else
	    stop ".XVP not found"
	  end if
	
!generate initial middleimages
	  do i=1,natot
	    BAND_slope(1:3)= rclas_BAND(1:3,i,middle_point)-rclas_BAND(1:3,i,1)
	    BAND_slope=BAND_slope/(dble(middle_point) - 1.d0)
	    BAND_const=rclas_BAND(1:3,i,1)-BAND_slope(1:3)
	    do k=1, middle_point
	      rclas_BAND(1:3,i,k)=BAND_slope(1:3)*dble(k) + BAND_const(1:3)
	    end do
	
	    !usign TS state case
	    if (middle_point .ne. NEB_Nimages) then
	      BAND_slope(1:3)= rclas_BAND(1:3,i,NEB_Nimages) - rclas_BAND(1:3,i,middle_point)
	      BAND_slope=BAND_slope/(dble(NEB_Nimages) - dble(middle_point))
	      BAND_const=rclas_BAND(1:3,i,middle_point) - dble(middle_point)*BAND_slope(1:3)
	      do k=middle_point, NEB_Nimages
	        rclas_BAND(1:3,i,k)=BAND_slope(1:3)*dble(k) + BAND_const(1:3)
	      end do
	    end if
	  end do
	end if
	return
	end subroutine NEB_make_initial_band



	subroutine NEB_save_traj_energy()
	use scarlett, only: natot, na_u, iza, pc, NEB_Nimages, rclas_BAND, Energy_band, Ang
	implicit none
	integer :: replica_number, i
	integer :: unitnumber
	character*13 :: fname
	
	!save energy
	do replica_number = 1, NEB_Nimages
	  write(*,*)"Energy-band", replica_number," ",Energy_band(replica_number)
	end do
	  write(*,*)"Energy-band"
	
	!open files
	do replica_number = 1, NEB_Nimages
	  unitnumber=replica_number+500
	  if (replica_number.lt. 10) then
	    write(fname,"(A7,I1,A4)") "Replica",replica_number,".xyz"
	  elseif (replica_number.ge. 10) then
	    write(fname,"(A7,I2,A4)") "Replica",replica_number,".xyz"
	  end if
	  open(unit=unitnumber,file=fname, access='APPEND')
	end do
	

	!write .xyz
	do replica_number = 1, NEB_Nimages
	  unitnumber=replica_number+500
	  write(unitnumber,*) natot
	  write(unitnumber,*)
	
	  do i=1, natot
	    if (i.le.na_u) then
	      write(unitnumber,345) iza(i), rclas_BAND(1:3,i,replica_number)/Ang
	    else
	      write(unitnumber,346) pc(i-na_u), rclas_BAND(1:3,i,replica_number)/Ang
	    end if
	  end do
	end do
	
	!close files
	do replica_number = 1, NEB_Nimages
	  unitnumber=replica_number+500
	  close(unitnumber)
	end do
 345  format(2x, I2,    2x, 3(f10.6,2x))
 346  format(2x, f10.6, 2x, 3(f10.6,2x))
	end subroutine NEB_save_traj_energy









