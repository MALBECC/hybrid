	module scarlett
!general module for hybrid
!In future i will move most variables to this module 03/2018
	use precision, only: dp
	implicit none

! General Variables
	character*30 :: sname
	integer :: istep, inicoor,fincoor !actual, initial and final number of move step for each restrain
	integer :: idyn ! kind of calculation, idyn=0 (CG minimization), idyn=1 (NEB minimization)
	integer :: nmove !max number of move step for each restrain
	integer :: nesp !number of QM species
	integer, dimension(:), allocatable :: isa,iza !Chemical Specie Label, and atomic charge
	real(dp) :: ftol !MAx force tol criteria (in Ry/Bohr)
	real(dp), dimension(:,:), allocatable :: xa, fa !position and forces of QM atoms
	double precision, dimension(:,:), allocatable :: rclas !Position of all atoms
	double precision, dimension(:,:), allocatable :: vat !velocities of all atoms 
	double precision, dimension(:,:), allocatable :: aat !aceleration of all atoms
	double precision :: time_steep
	double precision :: time_steep_max !max value of timesteep in FIRE algorithm
	double precision :: alpha !alpha value in FIRE algorithm
	integer :: Ndamped !number of times that system was damped in QM of FIRE algorithm
	integer :: Ndescend !number of consecutive steps in which F·v >= 0
	logical :: qm, mm ! True when system have a subsystem QM,MM
	integer :: nparm !number of bond types in amber.parm. esta fijado en 500 por algun motivo, hay q arreglar esto, Nick
	character, dimension(:), allocatable :: atsym*2 !atomic symbol
	character :: slabel*20 ! system label, name of outputs

! Cut Off & freeeze variables 
	integer, allocatable, dimension(:) :: r_cut_list_QMMM
	logical, allocatable, dimension(:) :: MM_freeze_list
	integer :: natoms_partial_freeze !number of atoms with 0 force in any direction
	integer,  allocatable, dimension(:,:) :: coord_freeze !cartesian direction with force = 0


! Cut Off QM-MM variables
	integer, dimension(:), allocatable, save:: blocklist,blockqmmm,listqmmm
	integer, dimension(:), allocatable, save:: blockall !JOTA
!listas para congelar atomos, hay q reveer estas subrutinas, por ahora estoy usando mis subrutinas, nick

!Dynamics
	double precision :: Ekinion ! Kinectic energy
	double precision :: tempion ! Actual Temperature of system
        double precision :: tempqm ! Actual Temperature of QM subsystem
	double precision :: tempinit ! Starting Temperature
	double precision :: tt ! Target Temperature
	double precision :: tauber ! Bath Coupling Constant Berendsen
        double precision :: kn ! Kinetic energy of Nose variable
        double precision :: vn ! Potential energyy of Nose var
        double precision :: mn ! Mass of Nose thermostat
! Solvent (MM) General variables
	integer :: nac !number of MM atoms
	integer :: na_u !number of QM atoms
	integer :: natot !total number of atoms
	double precision, dimension(:), allocatable :: masst !atomic mass


	character*4,  dimension(:), allocatable :: atname,aaname !atom and residue name
!in *.fdf
!%block SolventInput
!ATOM      1  O   WAT H   2      53.559  34.043  71.033
!             ^    ^
!             |    |
!         atname  aaname

	character*4,  dimension(:), allocatable :: attype !atom type
!in amber.parm
!residues
!101
!PRE 24
! "O1" "oh" 0 1 131072 1 8 -0.627601
!  ^    ^
!  |    |
!atname attype -> qmattype=O


	character*4,  dimension(:), allocatable :: qmattype !QM atom type
!in *.fdf
!%block SoluteAtomTypes
!  o   O5    PRE
!  ^    
!  |    
! qmattype

	integer, dimension(:), allocatable :: aanum !aanum(i) = numero de residuo al cual pertenece el atomo i
	integer, dimension(:,:), allocatable :: ng1 !ng1(i,j) j-esimo atomo al que esta unido el atomo i-esimo


!!!! bond parameters
	character*5, dimension(:), allocatable :: bondtype !bond type name
	real(dp), dimension(:), allocatable :: kbond,bondeq ! force constant and equil distance
	integer, dimension(:), allocatable :: bondxat !bondxat(i) number of bonds of atom i
!in amber.parm
!bonds
!164
!OW-HW  553.0    0.9572    ! TIP3P water
!  ^        ^     ^
!  |        |     |
!bondtype kbond bondeq


!!!! angle parameters
	character*8, dimension(:), allocatable :: angletype !angle type name
	real(dp), dimension(:), allocatable ::  kangle,angleeq ! force constant and equil angle
	integer, dimension(:), allocatable :: angexat !angexat(i) number of angles of atom i with i in an extreme
	integer, dimension(:), allocatable :: angmxat !angmxat(i) number of angles of atom i with i in middle

!in amber.parm
!angles
!375
!HW-OW-HW    100.      104.52    TIP3P water
!    ^        ^           ^
!    |        |           |
!angletype  kangle     angleeq


!!!! dihedral parameters
	character*11, dimension(:), allocatable :: dihetype !diedral type
	real(dp), dimension(:), allocatable :: kdihe,diheeq, perdihe !
	integer, dimension(:), allocatable :: multidihe
	integer, dimension(:), allocatable :: dihexat !dihexat(i) number of dihedral of atom i with i in extreme
	integer, dimension(:), allocatable :: dihmxat !dihmxat(i) number of dihedral of atom i with i in middle

! in this case E=kdihe/multidihe * (1+ cos (perdihe*dihedral-diheeq))

!in amber.parm
!dihes
!215
!X -CA-CA-X    4   14.50        180.0             2.         
!    ^         ^      ^           ^               ^
!    |         |      |           |               |
! dihetype multidihe kdihe     diheeq          perdihe


!!!! impropers 
	character*11, dimension(:), allocatable :: imptype
	real(dp), dimension(:), allocatable ::  kimp,impeq,perimp
	integer, dimension(:), allocatable :: multiimp
	integer, dimension(:), allocatable :: impxat !impxat(i) number of impropers of atom i

! E= kimp/multiimp*(1+COS(perimp*dihedral-impeq))
!in amber.parm
!imps
!53
!X -X -CB-X    1    1.0         180.              2.
!    ^         ^      ^           ^               ^
!    |         |      |           |               |
!imptype  multiimp  kimp        impeq          perimp

!!!! LJ & coulomb
	integer, dimension(:), allocatable :: scalexat !scalexat(i) number of scaled-LJ&coul atoms for atom i
	integer, dimension(:,:), allocatable :: scale !scale(i,j) j-esimo atomo de interaccion LJ&coul escalada con el atomo i-esimo

	integer, dimension(:), allocatable :: nonbondedxat !cantidad de elementos j en nonbonded(i,j) para i
	integer, dimension(:,:), allocatable :: nonbonded !nonbonded(i,j) atomo jesimo unido de alguna forma al atomo i-esimo, por lo que NO se debe calcular LJ & coulomb

	double precision, dimension(:), allocatable :: Em, Rm !LJ parameters, in hybrid units
	double precision, dimension(:), allocatable :: pc !charge of MM atoms
!in amber.parm
!ljs
!65
!  H           0.6000  0.0157            !Ferguson base pair geom.
!  ^              ^       ^
!  |              |       |
!ljtype         ∝Rm     ∝Em

	integer, dimension(:), allocatable :: izs !atomic number of a MM atom



	double precision, dimension(:,:), allocatable :: fce_amber !Total MM force
	double precision, dimension(:,:), allocatable :: fdummy !QM+MM+QMMM force 
	double precision, dimension(:,:), allocatable :: cfdummy !QM+MM+QMMM force after constrain


	integer, dimension(:,:), allocatable :: ng1type !ng1type(i,j) number of bond type defined in amber.parm for atom i, bond j
	integer, dimension(:,:), allocatable :: angetype !angetype(i,j) number of angle type defined in amber.parm for atom i, bond j, with i in extreme
	integer, dimension(:,:), allocatable ::  angmtype !angmtype(i,j) number of angle type defined in amber.parm for atom i, bond j, with i in middle
	integer, dimension(:,:), allocatable :: dihety,dihmty,impty !same with dihedrals and impropers
	logical, dimension(:,:), allocatable :: evaldihelog !control for evaluate diedral i, j on energy/force
	logical, dimension(:,:), allocatable :: evaldihmlog !control for evaluate diedral i, j on energy/force
	integer, dimension(:,:,:), allocatable :: evaldihe,evaldihm
	integer, dimension(:,:,:), allocatable ::  atange, atangm, atdihe,atdihm,atimp


! Link Atom variables
	logical :: linkatom !control for link atoms subroutines
	integer :: numlink !number of link atoms

	integer :: linkat(15) !linkat(i) numero  de atomo QM del i-esimo link atom
	integer :: linkqm(15,4),linkmm(15,4) !linkXm(i,*) lista de atomos QM/MM vecinos al link atom i
	integer :: linkmm2(15,4,3) !lista de segundos vecinos
	integer :: parametro(15,22,4) !parametro(i,*,*) posiciones en los arrays de constantes de enlace, angulos y diedros correspondientes a las interacciones del link atom i
	character :: linkqmtype(15,4)*4 ! tipo de atomo en linkqm(i,j)
	double precision :: Elink !Energy of link atoms
	double precision :: distl(15) !distancia del link atom i al atomo QM mas cercano
	double precision :: pclinkmm(15,15),Emlink(15,4)
	logical :: frstme

! Lio
      double precision :: spin !number of unpaired electrons
      integer :: charge !charge of QM sub-system



!NEB variables
	integer :: NEB_Nimages !number of images in band methods
	integer :: NEB_firstimage, NEB_lastimage !first and last imagas that can be move
	integer :: NEB_move_method !1= steepest descend, 2=velocity verlet
	double precision :: NEB_spring_constant !cosntant for space images
	double precision, allocatable, dimension(:,:,:) :: aclas_BAND_old !aceleration of previus step
	double precision, dimension(:,:,:), allocatable :: rclas_BAND!Position of all atoms in BAND method
	double precision, dimension(:,:,:), allocatable :: vclas_BAND!velocities of all atoms in BAND method
	double precision, dimension(:,:,:), allocatable :: fclas_BAND!Force all atoms in BAND method
	double precision, dimension(:,:,:), allocatable :: fclas_BAND_fresh !Force all atoms in BAND method just for energy gradient
	double precision, dimension(:), allocatable :: Energy_BAND !Energy of each image
	double precision :: NEB_steep_size !steep size in steepest descend algorithm in NEB
	double precision :: NEB_MAXFmod !Max force in NEB optimizarion
	integer :: PNEB !enable partial nudged elastic band
	integer :: PNEB_ini_atom, PNEB_last_atom ! initial and last atom in PNEB
	double precision, dimension(:,:), allocatable :: NEB_distl !distancia del link atom i al atomo QM mas cercano
	double precision, dimension(:), allocatable :: NEB_time_steep, NEB_alpha ! time steep and alpha value for image i in FIRE NEB
	integer, dimension(:), allocatable :: NEB_Ndescend !number of consecutive steps in which F·v >= 0 for image i in FIRE NEB
!outputs
	integer :: writeRF ! force integration
	integer :: traj_frec ! Frecuency to write trayectory and Energy in .rcg and .rce files
! Conversion factors
	real(dp) :: Ang !r_in_bohr=r_in_ang * Ang
	real(dp) :: eV !E_in_Hartree=E_in_eV * eV
	real(dp) :: kcal ! E_in_kcal_mol-1 = E_in_eV * kcal

! Others that need check
	real(dp) :: ucell(3,3)

! Nuevos jota
	double precision :: Nav !Avogadro Number
	double precision :: pi  
	end module scarlett
