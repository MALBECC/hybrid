!****************************************************************************
! Subroutine do_energy_forces 
!
! Calculates energy and forces for QM & MM subsystems
! Position of all atoms in rclas in Bohrs
! Final forces will be in cfdummy in Hartree/Bohr
! Final energy will be in Etots in Hartree
! Extracted from original version on hybrid.f
! J. Semelak & N. Foglia 2019
!*****************************************************************************

	subroutine do_energy_forces(rcorteqmmm, radbloqmmm, Etot, &
	do_SCF, do_QM_forces, do_properties, istp, step, &
	nbond, nangle, ndihe, nimp, Etot_amber, Elj, &
	Etots, constropt,nconstr, nstepconstr, typeconstr, kforce, ro, &
	rt, coef, atmsconstr, ndists, istepconstr, rcortemm, &
	radblommbond, optimization_lvl, dt, sfc, water, imm, rini, rfin)

! Modules
	use precision, only: dp
	use sys, only: die
	use ionew, only: io_setup
	use scarlett, only: qm, mm, na_u, natot, nroaa, masst, pc, &
	fdummy, cfdummy, nac, rclas, Em, Rm, linkatom, numlink, pclinkmm, & 
	Emlink, istep, idyn, nparm, atname, aaname, attype, ng1, &
	bondtype, kbond, bondeq, bondxat, angletype, kangle,angleeq, &
	angexat, angmxat, dihetype, kdihe,diheeq, perdihe, multidihe, dihexat, &
	dihmxat, imptype, kimp,impeq,perimp, multiimp, impxat, &
	scalexat, scale, nonbondedxat, nonbonded, fce_amber, ng1type, &
	angetype, angmtype, dihety,dihmty,impty, evaldihelog, evaldihmlog, &
	atange, atangm, atdihe, atdihm, atimp, vat, evaldihe, evaldihm, &
	linkat, linkqm, linkmm, linkmm2, parametro, linkqmtype, Elink, &
!cutoff
	r_cut_list_QMMM,blocklist,blockqmmm, blockall, listqmmm, &
!external potential
	external_potential, &
!units
	Ang, eV, kcal, &
!outputs
	slabel

	implicit none

! General Variables
	integer, intent(inout) :: step !total number of steps in a MMxQM movement (not tested in this version)
	integer, intent(in) :: istp !number of move step for each restrain starting on 1
	integer, intent(in) :: imm !MM step by QM each step
	integer, intent(in) :: nbond, nangle, ndihe, nimp !number of bonds, angles, dihedrals and impropers defined in amber.parm
	real(dp), intent(inout) :: Etot ! QM + QM-MM interaction energy
	double precision, intent(out) :: Etot_amber !total MM energy
	double precision, intent(out) :: Elj !LJ interaction (only QMMM)
	double precision, intent(out) :: Etots !QM+QMMM+MM energy
	real(dp), intent(in) :: dt !time step

! Cut Off QM-MM variables
	double precision, intent(in) :: rcortemm ! distance for LJ & Coulomb MM interaction
	double precision, intent(in) :: radblommbond !parche para omitir bonds en extremos terminales, no se computan bonds con distancias mayores a radblommbond
	double precision :: rcorteqmmm !Distance for QM-MM interaction
	double precision, intent(in) :: radbloqmmm !Distance that allow to move MM atoms from QM sub-system
	logical :: actualiz!MM interaction list control
	integer :: at_MM_cut_QMMM !MM atoms inside QM-MM cutoff
	double precision, dimension(:,:), allocatable :: r_cut_QMMM !Positions of QM & MM atoms inside QM-MM cutoff
	double precision, dimension(:,:), allocatable :: F_cut_QMMM !Forces of QM & MM atoms inside QM-MM cutoff
	double precision, dimension(:), allocatable :: Iz_cut_QMMM !Charge of MM atoms inside QM-MM cutoff

! Lio
	logical, intent(in) :: do_SCF, do_QM_forces !control for make new calculation of rho, forces in actual step
	logical, intent(in) :: do_properties !control for lio properties calculation

! Optimization scheme
	integer, intent(in) :: optimization_lvl ! level of movement in optimization scheme:
! 1- only QM atoms with restrain
! 2- only MM atoms
! 3- all

! ConstrOpt variables
	logical, intent(in) :: constropt !activate restrain optimizaion
	integer, intent(in) :: nconstr !number of constrains
	integer, intent(in) :: nstepconstr !numero de pasos en los que barrera la coordenada de reaccion (limitado a 1-100 hay q cambiar esto, Nick)
	integer, dimension(20), intent(in) :: typeconstr !type of cosntrain (1 to 8)
	double precision, dimension(20), intent(in) :: kforce !force constant of constrain i
	double precision :: rini,rfin  !initial and end value of reaction coordinate
	double precision, dimension(20), intent(in) :: ro ! fixed value of constrain in case nconstr > 1 for contrains 2+
	double precision, dimension(20), intent(inout) :: rt ! value of reaction coordinate in constrain i
	double precision, dimension(20,10), intent(in) :: coef ! coeficients for typeconstr=8
	integer, dimension(20,20), intent(in) :: atmsconstr 
	integer, dimension(20), intent(in) :: ndists !atomos incluidos en la coordenada de reaccion
	integer, intent(in) :: istepconstr !step of restraint 

! Others that need check
! Solvent General variables
	double precision, intent(in) :: sfc
	logical, intent(in) :: water
! Auxiliars
	integer :: i

! ---------------------------------------------------------------------------------------
	at_MM_cut_QMMM = nac
! Calculate Energy and Forces using Lio as Subroutine
	if(qm .and. (imm.eq.1)) then
	  if (allocated(r_cut_QMMM)) deallocate(r_cut_QMMM)
	  if (allocated(F_cut_QMMM)) deallocate(F_cut_QMMM)
	  if (allocated(Iz_cut_QMMM)) deallocate(Iz_cut_QMMM)
 
	  call compute_cutsqmmm(at_MM_cut_QMMM,istepconstr,radbloqmmm, &
	  rcorteqmmm,nroaa)

	  allocate (r_cut_QMMM(3,at_MM_cut_QMMM+na_u), F_cut_QMMM(3,at_MM_cut_QMMM+na_u), &
	  Iz_cut_QMMM(at_MM_cut_QMMM+na_u)) 

	  r_cut_QMMM=0.d0
	  F_cut_QMMM=0.d0
	  Iz_cut_QMMM=0

!copy position and nuclear charges to cut-off arrays
	  do i=1,natot !barre todos los atomos
	    if (i.le.na_u) then !QM atoms
	      r_cut_QMMM(1:3,i)= rclas(1:3,i)
	    else if (r_cut_list_QMMM(i-na_u) .ne. 0) then !MM atoms inside cutoff
	      r_cut_QMMM(1:3,r_cut_list_QMMM(i-na_u)+na_u) = rclas(1:3,i)
	      Iz_cut_QMMM(r_cut_list_QMMM(i-na_u))= pc(i-na_u)
	    end if
	  end do


#ifdef LIO
	  call SCF_hyb(na_u, at_MM_cut_QMMM, r_cut_QMMM, Etot, &
	  F_cut_QMMM, Iz_cut_QMMM, do_SCF, do_QM_forces, do_properties) !fuerzas lio, Nick
# else
	  STOP 'NO QM program defined in do_energy_forces'
#endif


! return forces to fullatom arrays
	  do i=1, natot
	    if (i.le.na_u) then !QM atoms
	      fdummy(1:3,i)=F_cut_QMMM(1:3,i)
	    else if (r_cut_list_QMMM(i-na_u).ne.0) then !MM atoms in cut-off
	      fdummy(1:3,i)= F_cut_QMMM(1:3,r_cut_list_QMMM(i-na_u)+na_u)
	    end if
	  end do
	endif !qm

! here Etot in Hartree, fdummy in Hartree/bohr
	step = step +1
	if(imm.gt.1) then
	  write(6,*)
	  write(6,'(A)')    '*******************************'
	  write(6,'(A,i5)') '   MM x QM Step : ', imm 
	  write(6,'(A)')    '*******************************'
	endif
	
! Calculation of last QM-MM interaction: LJ Energy and Forces only 
	if((qm.and.mm)) call ljef(na_u,nac,natot,rclas,Em,Rm,fdummy,Elj,listqmmm)

! LinkAtom: set again linkmm atoms parameters
	if(qm.and.mm) then
	  if(linkatom) then
	    do i=1,numlink
	      pc(linkmm(i,1:4))=pclinkmm(i,1:4)
	      Em(na_u+linkmm(i,1:1))=Emlink(i,1:1)
	    enddo
	  endif !LA
	endif !qm & mm

! Calculate pure Solvent energy and forces
! Forces in Hartree/bohr JOTA
	if(mm) then
	  call solv_ene_fce(natot,na_u,nac,ng1,rclas,Em,Rm,pc(1:nac), &
	  Etot_amber,fce_amber,attype, &
	  nbond,nangle,ndihe,nimp,multidihe, multiimp,kbond,bondeq, &
	  kangle,angleeq,kdihe,diheeq,kimp,impeq,perdihe,perimp, &
	  bondtype,angletype,dihetype,imptype, &
	  bondxat,angexat,angmxat,dihexat,dihmxat,impxat,atange, &
	  atangm,atdihe,atdihm,atimp, &
	  ng1type,angetype,angmtype,evaldihe,evaldihm, &
	  dihety,dihmty,impty,nonbonded,scale,nonbondedxat,scalexat, &
	  evaldihelog,evaldihmlog,step,nparm, &
	  actualiz,rcortemm, &
	  atname,aaname,sfc,dt, &
	  water,masst,radblommbond)
	endif !mm


! fce_amber in kcal/(molAng) 
! converts fdummy to Kcal/mol/Ang          
	fdummy(1:3,1:natot)=fdummy(1:3,1:natot)*Ang*kcal/eV

! here Etot in Hartree, fdummy in kcal/mol Ang
! add famber to fdummy  
	if(mm) fdummy(1:3,na_u+1:natot)=fdummy(1:3,na_u+1:natot) +fce_amber(1:3,1:nac)

! here fdummy in kcal/mol/Ang

! Calculation of LinkAtom Energy and Forces
	if(qm.and.mm ) then
	  if(linkatom) then
	    call link2(numlink,linkat,linkqm,linkmm,linkmm2,rclas, &
	    natot,nac,fdummy,attype,nparm, &
	    nbond,nangle,ndihe,multidihe,kbond,bondeq, &
	    kangle,angleeq,kdihe,diheeq,perdihe, &
	    bondtype,angletype,dihetype,linkqmtype, &
	    Elink,parametro,step)

! Set again link atoms parameters to zero for next step  
	    do i=1,numlink
	      pclinkmm(i,1:4)=pc(linkmm(i,1:4))
	      pc(linkmm(i,1:1))=0.d0
	      pc(linkmm(i,2:4))=pc(linkmm(i,2:4))+pclinkmm(i,1)/3.d0
	      Em(na_u+linkmm(i,1:1))=0.d0
	    enddo
	  endif ! LA
	endif !qm & mm

	if(optimization_lvl.eq.1) fdummy=0.d0 !only move atoms with restrain

! Calculation of Constrained Optimization Energy and Forces 
	if(imm.eq.1) then
	  if(constropt) then
	    call subconstr2(nconstr,typeconstr,kforce,rini,rfin,ro,rt, &
	    nstepconstr,atmsconstr,natot,rclas,fdummy,istp,istepconstr, &
	    ndists,coef)
	    if(idyn .ge. 4) call subconstr4(istep,rt(1),slabel)
	  endif 
	endif !imm

! Converts fdummy Hartree/bohr
	fdummy(1:3,1:natot)=fdummy(1:3,1:natot)*eV/(Ang*kcal) ! here Etot in Hartree, fdummy in Hartree/bohr
! Writes final energy decomposition
	Etots=Etot+1.d0*(Elj+(Etot_amber+Elink)*eV/kcal) !Elj in Hartree, Etot_amber and Elink in kcal/mol


	if (external_potential .gt. 0) call external_bias(external_potential,natot,rclas,fdummy,Etots)
      
! here Etot in Hartree
	write(6,*)
	write(6,'(/,a)') 'hybrid: Potential Energy Decomposition (eV):'
	if(qm) write(6,999)           'Elio :',Etot/eV      
	if(qm.and.mm) write(6,999)    'Elj:    ',Elj/eV       
	if(mm) write(6,999)      'Esolv:  ',Etot_amber/kcal   
	if(Elink.ne.0.0) write(6,999) 'Elink:  ',Elink/kcal
	write(6,999)    'Etots:  ',Etots/eV
	call flush(6)

! Sets fdummy to zero inside mmxqm step
	if(qm.and.mm) then
	  if(imm.ne.1) then
	    fdummy(1:3,1:na_u) = 0.d0
	    vat(1:3,1:na_u) = 0.d0
	    if(linkatom) then
	      do i=1,numlink
	        fdummy(1:3,linkmm(i,1)+na_u)=0.d0
	        vat(1:3,linkmm(i,1)+na_u)=0.d0
	      enddo
	    endif !LA
	  endif !imm
	endif !qm & mm

! here Etot in Hartree, fdummy in Hartree/bohr
! Impose constraints to atomic movements by changing forces
	call fixed2(na_u,nac,natot,blocklist,blockqmmm, &
	     fdummy,cfdummy,vat,optimization_lvl,blockall)
! from here cfdummy is the reelevant forces for move system
! here Etot in Hartree, cfdummy in Hartree/bohr

 999  format(a,2x,F30.18)
	end subroutine do_energy_forces

