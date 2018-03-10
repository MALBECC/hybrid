	module scarlett
!general module for hybrid
!In future i will move most variables to this module 03/2018
	implicit none

! Solvent (MM) General variables
	integer :: natot !total number of atoms

!NEB variables
	integer :: replicas !number of replicas in band methods
	double precision, allocatable, dimension(:,:,:) :: aclas_BAND_old !aceleration of previus step
	double precision, dimension(:,:,:), allocatable, save:: rclas_BAND!Position of all atoms in BAND method
	double precision, dimension(:,:,:), allocatable, save:: vclas_BAND!velocities of all atoms in BAND method
	double precision, dimension(:,:,:), allocatable, save:: fclas_BAND!Force all atoms in BAND method
	double precision, dimension(:), allocatable :: Energy_band !Energy of each image



	end module scarlett
