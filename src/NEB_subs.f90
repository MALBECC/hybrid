	SUBROUTINE NEB_steep(istep,relaxd,atmsconstr)
!Nudged elastic band method
!movement methods avalables are: Steepest descend, quick-min and FIRE using velocity verlet
!N. Foglia 03/2018
	use scarlett, only: NEB_firstimage, NEB_lastimage,  &
        NEB_move_method, NEB_spring_constant, &
        ftol, NEB_steep_size, NEB_MAXFmod, NEB_Nimages, &
	verbose_level, rclas_BAND
	implicit none
	integer, intent(in) :: istep
	logical, intent(inout) :: relaxd
	integer :: replica_number
	integer :: MAX_FORCE_REPLICA, MAX_FORCE_ATOM
	integer, dimension(20,20), intent(in) :: atmsconstr
	double precision ::  MAXFmod_total, NEB_Ekin
	logical, dimension(NEB_Nimages) :: NEB_converged_image

	if (istep.eq.1) then!initial max steep size
	  NEB_steep_size=0.1d0
!          fclas_BAND_fresh=fclas_BAND
        else
	  MAX_FORCE_REPLICA=NEB_firstimage
	  if ( (MAX_FORCE_REPLICA .ne. 1) .and. (MAX_FORCE_REPLICA.ne. NEB_lastimage)) then
!	    fclas_BAND_fresh(1:3, 1:natot,MAX_FORCE_REPLICA)=fclas_BAND(1:3, 1:natot,MAX_FORCE_REPLICA)
!	    fclas_BAND(1:3, 1:natot,MAX_FORCE_REPLICA-1)=fclas_BAND_fresh(1:3, 1:natot,MAX_FORCE_REPLICA-1)
!	    fclas_BAND(1:3, 1:natot,MAX_FORCE_REPLICA+1)=fclas_BAND_fresh(1:3, 1:natot,MAX_FORCE_REPLICA+1)
	  end if
!	  NEB_firstimage=MAX_FORCE_REPLICA-1
!	  NEB_lastimage=MAX_FORCE_REPLICA+1
	end if

	NEB_converged_image=.true.

	relaxd=.true.
	MAXFmod_total=0.d0
	MAX_FORCE_REPLICA=-1
	MAX_FORCE_ATOM=-1

!	NEB_firstimage=2
!	NEB_lastimage=NEB_Nimages-1

	do replica_number=NEB_firstimage, NEB_lastimage
	  if ( replica_number .gt. 1 .and. replica_number .lt. NEB_Nimages) then
	  call NEB_Force(replica_number,NEB_spring_constant,atmsconstr) !Force in H/bohr
	  end if
	end do

	do replica_number=2, NEB_Nimages-1
	  call NEB_check_convergence(relaxd, replica_number, MAXFmod_total, MAX_FORCE_REPLICA, MAX_FORCE_ATOM, NEB_converged_image)
	end do
	if (verbose_level.gt.5) call NEB_status()


! select image to move in this step
!	NEB_firstimage=MAX_FORCE_REPLICA-1
!	NEB_lastimage=MAX_FORCE_REPLICA+1

	if (NEB_firstimage .lt. 1) NEB_firstimage=1
	if (NEB_lastimage .gt. NEB_Nimages) NEB_lastimage=NEB_Nimages

! full band move until more test
!	NEB_firstimage=1
!	NEB_lastimage=NEB_Nimages


	if (.not. relaxd) then
	  call NEB_movement_algorithm(NEB_move_method, MAXFmod_total, NEB_firstimage, NEB_lastimage) !move systems
	  write(*,*) "system NOT converged"
	else
	  write(*,*) "system converged"
	end if

! recompute SCF only in image that was move
!	NEB_firstimage=MAX_FORCE_REPLICA
!	NEB_lastimage=MAX_FORCE_REPLICA

! full band move until more test
!	NEB_firstimage=2
!	NEB_lastimage=NEB_Nimages-1
! JOTA
	if (NEB_move_method .eq. 1) THEN
	  if (istep.eq.1) NEB_MAXFmod=MAXFmod_total
	  if (MAXFmod_total .gt. NEB_MAXFmod) NEB_steep_size=NEB_steep_size*0.85d0
	  NEB_MAXFmod=MAXFmod_total

! Comento estoy porque con feopt los gradientes son chiquitos
!   if (NEB_steep_size.lt.1d-5) then
!	    relaxd=.true.
!	    write(*,*) "max precision reached on atomic displacement"
!	  end if

	else if (NEB_move_method .ge. 2) THEN
	  CALL NEB_calculate_T(NEB_Ekin)
	  WRITE(*,*) "TOTAL K ", NEB_Ekin, "AVERAGE K", NEB_Ekin/(dble(NEB_Nimages-2))
	END IF
	write(*,*) "max force ", sqrt(MAXFmod_total), "on atom ", MAX_FORCE_ATOM,  &
        "in replica ", MAX_FORCE_REPLICA, "conv criteria", ftol
	END SUBROUTINE NEB_steep

	SUBROUTINE NEB_Force(replica_number,NEB_spring_constant,atmsconstr)
!this subroutine obtains nundged elastic band force for each atom in replica_number
	use scarlett, only: natot, rclas_BAND, Energy_band, fclas_BAND, NEB_firstimage, NEB_lastimage, NEB_CI, NEB_Nimages
	implicit none
	integer, intent(in) :: replica_number
	double precision, intent(in) :: NEB_spring_constant
	double precision, dimension(3,natot) :: tang_vec, Fspring
	integer, dimension(20,20), intent(in) :: atmsconstr
	double precision :: dist01, dist12, Forceproj
	double precision :: Emax
	integer :: Imax
	integer :: i, j

	if (replica_number.eq.1 .or. replica_number.eq. NEB_Nimages) return

	if (NEB_CI.eq.1) then !find image with max energy
	  Emax=-9.d99
	  Imax=1
	  do i=1, NEB_lastimage
	    if (Energy_band(i).gt.Emax) then
	      Imax=i
	      Emax=Energy_band(i)
	    end if
	  end do
	end if

	call NEB_calculate_tg(0,replica_number,tang_vec,atmsconstr)

!Spring Force
	Fspring=0.d0
	do i=1, natot
	  dist01=0.d0
	  dist12=0.d0
	  do j=1,3
	    dist01=dist01+(rclas_BAND(j,i, replica_number+1)-rclas_BAND(j,i, replica_number))**2
	    dist12=dist12+(rclas_BAND(j,i, replica_number)-rclas_BAND(j,i, replica_number-1))**2
	  end do
	  dist01=sqrt(dist01)
	  dist12=sqrt(dist12)
	  Fspring(1:3,i)=NEB_spring_constant*(dist01-dist12)*tang_vec(1:3,i) !here distXX in Bohr, force in Hartree/bohr
	end do

!NEB Force
	do i=1, natot
	  Forceproj=0.d0
	  do j=1,3
	    Forceproj=Forceproj+fclas_BAND(j,i,replica_number)*tang_vec(j,i)
	  end do
	  if ((NEB_CI.eq.1) .and. (replica_number.eq.Imax)) then !Climb image case
	    write(*,*) "inv gradient in image:", Imax
	    fclas_BAND(1:3,i,replica_number)=fclas_BAND(1:3,i,replica_number)-2*Forceproj*tang_vec(1:3,i)
	  else !Normal NEB case
	    fclas_BAND(1:3,i,replica_number)=fclas_BAND(1:3,i,replica_number)-Forceproj*tang_vec(1:3,i)
	    fclas_BAND(1:3,i,replica_number)=fclas_BAND(1:3,i,replica_number)+Fspring(1:3,i)
	  end if
	end do
	END SUBROUTINE NEB_Force

	SUBROUTINE NEB_status()
!this subroutine calculate distance and angles betwee NEB replicas
	use scarlett, only: natot, rclas_BAND, NEB_lastimage
	implicit none
	double precision :: cosphi, mod_prev, mod_next
	integer :: i,j,replica_number

	do replica_number=2,NEB_lastimage-1
	  cosphi=0.d0
	  mod_prev=0.d0
	  mod_next=0.d0
	  do i=1,natot
	    do j=1,3
	      cosphi=cosphi+(rclas_BAND(j,i, replica_number+1)-rclas_BAND(j,i, replica_number))* &
	           (rclas_BAND(j,i, replica_number)-rclas_BAND(j,i, replica_number-1))
	      mod_next=mod_next+(rclas_BAND(j,i, replica_number+1)-rclas_BAND(j,i, replica_number))* &
	          (rclas_BAND(j,i, replica_number+1)-rclas_BAND(j,i, replica_number))
	      mod_prev=mod_prev+(rclas_BAND(j,i, replica_number)-rclas_BAND(j,i, replica_number-1))*&
	          (rclas_BAND(j,i, replica_number)-rclas_BAND(j,i, replica_number-1))
	    end do
	  end do
	  mod_prev=sqrt(mod_prev)
	  mod_next=sqrt(mod_next)
	  cosphi=cosphi/(mod_prev*mod_next)
	  if (cosphi .gt. 1.d0) cosphi=1.d0
	  write(*,*) "NEB status, image: ", replica_number, &
	  "dist prev, distnext, cos(angle)", mod_prev,mod_next,cosphi
	end do
	end SUBROUTINE NEB_status



	SUBROUTINE NEB_calculate_tg(method,replica_number,tang_vec,atmsconstr)
!this subroutine need checks in method=1 and method=2
	use scarlett, only: natot, rclas_BAND, Energy_band, PNEB, PNEB_ini_atom, PNEB_last_atom, & !, NEB_Nimages,
	natmsconstr, feopt
	implicit none
	integer, intent(in) :: method, replica_number
	double precision, dimension(3,natot), intent(inout) :: tang_vec
  integer, dimension(20,20), intent(in) :: atmsconstr
	double precision, dimension(3,natot) :: tang_vecA, tang_vecB
	double precision :: NORMVEC, NORMVECA, NORMVECB
	double precision :: E0, E1, E2, Vmax, Vmin
	integer :: initial_atom, last_atom, at1
	integer :: i

	initial_atom=1
	last_atom=natot

	tang_vec=0.d0

	if ( PNEB .eq.1 ) then
	  initial_atom=PNEB_ini_atom
	  last_atom=PNEB_last_atom
	  write(*,*) "Runing Partial nudged elastic band between atoms", PNEB_ini_atom, "and ", PNEB_last_atom
	  tang_vec=0.d0
	end if

	IF (method.eq.0) then
	  tang_vec(1:3,initial_atom:last_atom) =  &
	  rclas_BAND(1:3,initial_atom:last_atom,replica_number+1) &
	  - rclas_BAND(1:3,initial_atom:last_atom,replica_number-1)
	ELSEIF (method.eq.1 .or. method.eq.2) then !The Journal of Chemical Physics 113, 9978 (2000); https://doi.org/10.1063/1.1323224
	  tang_vecA(1:3,initial_atom:last_atom) = &
	  rclas_BAND(1:3,initial_atom:last_atom,replica_number) &
	  - rclas_BAND(1:3,initial_atom:last_atom,replica_number-1)

	  tang_vecB(1:3,initial_atom:last_atom) =  &
	  rclas_BAND(1:3,initial_atom:last_atom,replica_number+1) &
	  - rclas_BAND(1:3,initial_atom:last_atom,replica_number)

!luego pasar la normalizacion a una subrutina, Nick
	  do i=initial_atom,last_atom
	    NORMVECA= tang_vecA(1,i)**2 + tang_vecA(2,i)**2 + tang_vecA(3,i)**2
	    NORMVECB= tang_vecB(1,i)**2 + tang_vecB(2,i)**2 + tang_vecB(3,i)**2
	    if (NORMVECA .lt. 1d-300) then
	      write(*,*) "WARNING : tangA = 0 in ", replica_number,i
	      tang_vecA(1:3,i)=0.d0
	    else if (NORMVECA .ne. NORMVECA) then
	      stop "NAN in tangent vector A"
	    else
!	      NORMVECA=sqrt(NORMVECA)
!	      tang_vecA(1:3,i)=tang_vecA(1:3,i)/NORMVECA
	    end if
	    if (NORMVECB .lt. 1d-300) then
	      write(*,*) "WARNING : tangB = 0 in ", replica_number,i
	      tang_vecB(1:3,i)=0.d0
	    else if (NORMVECB .ne. NORMVECB) then
	      stop "NAN in tangent vector A"
	    else
!	      NORMVECB=sqrt(NORMVECB)
!	      tang_vecB(1:3,i)=tang_vecB(1:3,i)/NORMVECB
	    end if
	  end do

	  IF (method.eq.1) then
	    tang_vec(1:3,initial_atom:last_atom)= &
	    tang_vecA(1:3,initial_atom:last_atom)+tang_vecB(1:3,initial_atom:last_atom)
	  ELSE
	    E0=Energy_band(replica_number-1)
	    E1=Energy_band(replica_number)
	    E2=Energy_band(replica_number+1)
	    do i=initial_atom,last_atom
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
	  END IF
	ELSE
	  STOP "Wrong method in NEB_calculate_tg"
	END IF

	do i=initial_atom,last_atom
	  NORMVEC= tang_vec(1,i)**2 + tang_vec(2,i)**2 + tang_vec(3,i)**2
	  NORMVEC=sqrt(NORMVEC)
	  if (NORMVEC .lt. 1d-300) then
	     tang_vec(1:3,i)=0.d0
	  else if (NORMVEC .ne. NORMVEC) then
	    stop "NAN in tangent vector"
	  else
	    tang_vec(1:3,i)=tang_vec(1:3,i)/NORMVEC
	  end if
	end do

! set tg=0 for not restrained atoms on free energy calculations
        if (feopt) then
						do i=1,natmsconstr
		          at1=atmsconstr(1,i)
							tang_vec(1:3,at1)=0
		        enddo
        end if

	RETURN
	END SUBROUTINE NEB_calculate_tg


	SUBROUTINE NEB_movement_algorithm(method,MAXFmod_total,NEB_firstimage,NEB_lastimage)
	use scarlett, only: aclas_BAND_old, NEB_Nimages, natot,rclas_BAND,vclas_BAND,fclas_BAND, masst,  &
	time_steep_max, NEB_steep_size, NEB_time_steep, NEB_Ndescend, NEB_alpha
	implicit none
	integer :: method
	double precision, intent(in) :: MAXFmod_total
	double precision :: MAXFmod
	double precision :: SZstep
	integer, intent(in) :: NEB_firstimage, NEB_lastimage
	integer :: replica_number
	double precision :: Fmod

	Fmod=0.d0
	if (method.eq.1) then !steepest descend

	  MAXFmod=sqrt(MAXFmod_total)
	  SZstep=NEB_steep_size/MAXFmod
	  write(*,*) "maxforce", MAXFmod, "stepsize", SZstep, &
	  "stepsize_base", NEB_steep_size
          do replica_number=NEB_firstimage+1, NEB_lastimage-1
	  rclas_BAND(1:3,1:natot,replica_number)=rclas_BAND(1:3,1:natot,replica_number)+SZstep*fclas_BAND(1:3,1:natot,replica_number)
          end do

	elseif (method.eq.2) then !quick-min using velocity verlet

	  do replica_number=NEB_firstimage+1, NEB_lastimage-1
	    write(*,*) "moving image ", replica_number
	    call quick_min(natot, rclas_BAND(:,:,replica_number), fclas_BAND(:,:,replica_number), &
	    aclas_BAND_old(:,:,replica_number), vclas_BAND(:,:,replica_number), masst)
	  end do

	elseif (method.eq.3) then !FIRE, Bitzek, et. al., Phys. Rev. Lett. 97, 170201 (2006).
	  do replica_number=NEB_firstimage+1, NEB_lastimage-1
	    write(*,*) "moving image ", replica_number
	    call FIRE(natot, rclas_BAND(:,:,replica_number), fclas_BAND(:,:,replica_number), &
	    vclas_BAND(:,:,replica_number), NEB_time_steep(replica_number), &
	    NEB_Ndescend(replica_number), time_steep_max, NEB_alpha(replica_number))
	    write(*,*) "timesteep", NEB_time_steep(replica_number)
	    if (NEB_time_steep(replica_number) .lt. 1D-10) NEB_time_steep(replica_number)=1D-2
	  end do
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



	subroutine NEB_make_initial_band(use_restart)
!Generate initial configurations of images in NEB using, .XV.i restarts or
!using only reactives and products restarts (and TS if posible)
!N. Foglia 03/2018
	use scarlett, only: natot, ucell, rclas, vat, rclas_BAND, vclas_BAND, NEB_Nimages
	implicit none
	logical, intent(in) :: use_restart
	logical :: band_xv_found, ok_restart
	logical :: foundxv, foundvat !control for coordinates restart
	integer :: replica_number !auxiliar
	double precision :: BAND_slope(3), BAND_const(3)
	integer :: i, k !auxiliars
	integer :: middle_point
	ok_restart=.false.
	if (use_restart) call NEB_restart(1, ok_restart)

	if (.not. ok_restart) then
	band_xv_found=.true.
	do replica_number = 1, NEB_Nimages
	  call ioxv( 'read', natot, ucell, rclas, vat, foundxv, foundvat,'X',replica_number)
	  band_xv_found=band_xv_found .and. foundxv
	  if (foundxv) then
	    rclas_BAND(1:3,1:natot,replica_number)=rclas(1:3,1:natot)
	    vclas_BAND(1:3,1:natot,replica_number)=vat(1:3,1:natot)
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
	end if
	return
	end subroutine NEB_make_initial_band



	subroutine NEB_save_traj_energy(step,slabel)
	use scarlett, only: natot, na_u, iza, pc, NEB_Nimages, rclas_BAND, Energy_band, Ang, fclas_BAND, eV, idyn
	implicit none
	character*20, intent(in) :: slabel
	integer,  intent(in) :: step
	integer :: replica_number, i
	integer :: unitnumber
	character*13 :: fname
	double precision :: Emax, Fmax


	Emax=maxval(dabs(Energy_band))
	Fmax=maxval(dabs(fclas_BAND))

	call wriene(step,slabel,idyn,Emax,Fmax) !Need to verify units here, Nick

	!save energy
	  write(*,*)"Energy-band in eV"
	do replica_number = 1, NEB_Nimages
	  write(*,*)"Energy-band", replica_number," ",Energy_band(replica_number)
	end do
	  write(*,*)"Energy-band"

	  write(*,*)"conv status", Emax/eV, Fmax*Ang/eV

	if (.false.) then !needs to include a high level of verbose here
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
	end if
 345  format(2x, I2,    2x, 3(f10.6,2x))
 346  format(2x, f10.6, 2x, 3(f10.6,2x))
	end subroutine NEB_save_traj_energy

	subroutine NEB_restart(read_write, ok_restart)
	use scarlett, only: NEB_move_method, rclas_BAND, vclas_BAND,  &
	aclas_BAND_old, NEB_time_steep, NEB_Ndescend, NEB_alpha, slabel
	implicit none
	integer, intent(in) :: read_write
	logical, intent(inout) :: ok_restart
	character*24 :: fname, paste
	integer :: NEB_move_method_read
	logical :: hay_restart
	external :: paste
	fname = paste(slabel,'.NEB')
	if(read_write .eq. 1) then !read
	  INQUIRE(FILE=fname, EXIST=hay_restart)
	  if (hay_restart) then
	    ok_restart=.true.
	    write(*,*) "using NEB restart ", fname
	    open(unit=935, file=fname, STATUS='UNKNOWN', ACCESS='STREAM')
	    write(935) NEB_move_method_read
	    if (NEB_move_method_read .eq. NEB_move_method)  &
	    STOP ("NEB_move_method differs from restart")

	    read(935) rclas_BAND
	    if (NEB_move_method .eq.2 .or. NEB_move_method .eq.3) then
	      read(935) vclas_BAND
	      read(935) aclas_BAND_old
	      if (NEB_move_method .eq.3) then
	        read(935) NEB_time_steep
	        read(935) NEB_Ndescend
	        read(935) NEB_alpha
	      end if
	    end if
	    close(935)
	  else
	    write(*,*) "WARNING ", fname, " not found"
	    ok_restart=.false.
	  end if
	elseif(read_write .eq. 2) then !write
	  open(unit=935, file=fname, STATUS='UNKNOWN', ACCESS='STREAM')
	    write(935) NEB_move_method
	    write(935) rclas_BAND
	    if (NEB_move_method .eq.2 .or. NEB_move_method .eq.3) then
	      write(935) vclas_BAND
	      write(935) aclas_BAND_old
	      if (NEB_move_method .eq.3) then
	        write(935) NEB_time_steep
	        write(935) NEB_Ndescend
	        write(935) NEB_alpha
	      end if
	    end if
	  close(935)
	else
	  stop "wrong number in read_write"
	end if
	end subroutine NEB_restart
