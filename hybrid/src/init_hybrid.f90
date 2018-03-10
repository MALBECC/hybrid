        subroutine init_hybrid(init_type)
!initialize variables for hybrid calculations
!kind of initiañlization is controlled by  init_type
! N. foglia 03/2018
      use scarlett, only: natot,replicas, aclas_BAND_old, rclas_BAND, vclas_BAND, &
	fclas_BAND, Energy_band
        implicit none
        character(len=*), intent(in) :: init_type
!ire moviendo inicializaciones a este subrutina
        if ( init_type == 'hola mundo') then
          write(*,*) "pase if 1"
          write(*,*) init_type

        elseif ( init_type == 'NEB') then

	allocate(rclas_BAND(3,natot,replicas), vclas_BAND(3,natot,replicas), &
	         fclas_BAND(3,natot,replicas), aclas_BAND_old(3,natot,replicas))
        allocate(Energy_band(replicas))
        rclas_BAND=0.d0
        vclas_BAND=0.d0
        fclas_BAND=0.d0
	Energy_band=0.d0


        else
          STOP "Wrong init_type"
        end if
        end subroutine init_hybrid

