        subroutine init_hybrid(init_type)
!initialize variables for hybrid calculations
!kind of initialization is controlled by  init_type
! N. foglia 03/2018
      use scarlett, only: natot, aclas_BAND_old, rclas_BAND, vclas_BAND, &
	fclas_BAND, Energy_band, NEB_firstimage, NEB_lastimage, NEB_Nimages
        implicit none
        character(len=*), intent(in) :: init_type
!ire moviendo inicializaciones a este subrutina
        if ( init_type == 'Jolie') then
          write(*,*) "Hi Angi!"
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

