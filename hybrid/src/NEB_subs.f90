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
	    STOP ".XVR not found"
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
	   "interpolating reactivs and products"
	    middle_point=NEB_Nimages
	  end if

!read products coordinates
	  call ioxv('read',natot,ucell,rclas,vat,foundxv,foundvat,'P',-1)
	  if (foundxv) then
	    rclas_BAND(1:3,1:natot,NEB_Nimages)=rclas(1:3,1:natot)
	  else
	    STOP ".XVP not fund"
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









