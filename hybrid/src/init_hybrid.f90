	subroutine init_hybrid(init_type)
!initialize variables for hybrid calculations
!kind of initialization is controlled by  init_type
! N. foglia 03/2018
	use precision, only: dp
	use fdf, only: fdf_integer
	use sys, only: die
	use scarlett, only: natot, aclas_BAND_old, rclas_BAND, vclas_BAND, &
	fclas_BAND, Energy_band, NEB_firstimage, NEB_lastimage, NEB_Nimages, &
	Ang, eV, kcal, na_u, qm, nesp
        implicit none
        character(len=*), intent(in) :: init_type
!ire moviendo inicializaciones a este subrutina en el futuro

	if ( init_type == 'Jolie') then
	  write(*,*) "Hi Angi!" !most important part of code of course
	! Read the number of QM atoms
	  na_u=fdf_integer('NumberOfAtoms',0)
	  if (na_u.eq.0) then
	    write(6,'(/a)') 'hybrid: Running with no QM atoms'
	    qm=.false.
	  endif
	! Read the number of MM atoms
	  nac = fdf_integer('NumberOfSolventAtoms',0)
	  if (nac.eq.0) then
	    write(6,'(/a)') 'hybrid: Running with no MM atoms'
	    mm=.false.
	  endif
	
	  if (nac.eq.0 .and. na_u.eq.0) then
	    call die("no atoms in system")
	  end if
	
	! Read the number of species
	  nesp = fdf_integer('NumberOfSpecies',0)
	  if(qm.and.(nesp.eq.0)) then
	    call die("hybrid: You must specify the number of species")
	  endif
	
	  allocate(xa(3,na_u), fa(3,na_u), isa(na_u), iza(na_u), atsym(nesp))
	
	elseif ( init_type == 'Constants') then !define constants and convertion factors
	  Ang    = 1._dp / 0.529177_dp
	  eV     = 1._dp / 27.211396132_dp
	  kcal   = 1.602177E-19_dp * 6.02214E23_dp / 4184.0_dp
	elseif ( init_type == 'NEB') then !initialize Nudged elastic band variables
	  if (NEB_Nimages .lt. 3) STOP 'Runing NEB with less than 3 images'
	  NEB_firstimage=1
	  NEB_lastimage=NEB_Nimages
	  allocate(rclas_BAND(3,natot,NEB_Nimages), vclas_BAND(3,natot,NEB_Nimages), &
	         fclas_BAND(3,natot,NEB_Nimages), aclas_BAND_old(3,natot,NEB_Nimages))
	  allocate(Energy_band(NEB_Nimages))
	  rclas_BAND=0.d0
	  vclas_BAND=0.d0
	  fclas_BAND=0.d0
	  Energy_band=0.d0
	else
	  STOP "Wrong init_type"
	end if
	
	end subroutine init_hybrid

