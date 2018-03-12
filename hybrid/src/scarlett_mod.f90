	module scarlett
!general module for hybrid
!In future i will move most variables to this module 03/2018
	use precision, only: dp
	implicit none


	real(dp) :: ftol !MAx force tol criteria (in Ry/Bohr)
	double precision, dimension(:,:), allocatable, save:: rclas !Position of all atoms
	double precision, dimension(:,:), allocatable, save:: vat !velocities of all atoms, not used for CG


! Solvent (MM) General variables
	integer :: natot !total number of atoms
	double precision, dimension(:), allocatable, save:: masst !atomic mass


!NEB variables
	integer :: NEB_Nimages !number of images in band methods
	integer :: NEB_firstimage, NEB_lastimage !first and last imagas that can be move
	integer :: NEB_move_method !1= steepest descend, 2=velocity verlet
	double precision :: NEB_spring_constant !cosntant for space images
	double precision, allocatable, dimension(:,:,:) :: aclas_BAND_old !aceleration of previus step
	double precision, dimension(:,:,:), allocatable, save:: rclas_BAND!Position of all atoms in BAND method
	double precision, dimension(:,:,:), allocatable, save:: vclas_BAND!velocities of all atoms in BAND method
	double precision, dimension(:,:,:), allocatable, save:: fclas_BAND!Force all atoms in BAND method
	double precision, dimension(:), allocatable :: Energy_BAND !Energy of each image


! Others that need check
	real(dp) :: ucell(3,3)



	end module scarlett
