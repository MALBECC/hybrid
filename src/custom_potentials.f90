!****************************************************************************
! This subroutine reads fdf specifications regarding custom potentials.
! Additional custom potential options should be added here.
! Options available
! custompot_type = 1 Special potential for a certain torsion
!   dihe_type = 1 Cosine series: sum from i=1 to ncos V(i)*cos(dihe)**(i-1)                 
! custompot_type = 2 LEPS potential (comming soon)
!****************************************************************************

  subroutine custom_potentials_assign(custompot_type)

  use fdf
  use scarlett, only: custom_dihe, cos_weights, ncos, dihe_type
  implicit none
  integer :: i, iunit, custompot_type
  character :: exp

  if (fdf_block('CustomPotential',iunit)) then
    read(iunit,*) exp, custompot_type
    write(*,*) "custompot_type", custompot_type
    if (custompot_type .eq. 1) then ! Special potential for a certain torsion
      read(iunit,*) exp,(custom_dihe(i),i=1,4)
      read(iunit,*) exp, dihe_type
      write(6,'(/,a)')  'Custom dihedral potential is turned on'
      write(6,*)  'Atoms integrating the custom dihedral:', custom_dihe(1:4)
      write(*,*) "dihe_type", dihe_type
      if (dihe_type .eq. 1) then ! Cosine series
        write(6,'(/,a)')  'Potential type: Cosine series'
        read(iunit,*) exp, ncos
        allocate(cos_weights(ncos))
        write(*,*) "ncos", ncos
        read(iunit,*) exp, (cos_weights(i),i=1,ncos)
         write(6,*)  'Cosine weights (kcal/mol): ', cos_weights(1:ncos)
      endif
    else
      write(6,'(/,a)')  'Wrong custompot_type in CustomPotential block'
    endif
  endif
  return
  100 stop 'Custom potentials: problem reading CustomPotential block'
  end subroutine custom_potentials_assign

  subroutine custom_dihe_energy_forces(Etots)
  use scarlett, only: custom_dihe, cos_weights, ncos, rclas, fdummy, dihe_type, &
  natot, eV, Ang, kcal

  implicit none
  double precision, intent(inout) :: Etots
  integer ::  i, at1, at2, at3, at4
  double precision :: dihe, Edihe, factor, dihedro, pi
  double precision, dimension(ncos) :: Vcos
  double precision, dimension(3,10) :: fnew
  double precision, dimension(12) :: fdihe

  pi=DACOS(-1.d0)
  at1=custom_dihe(1)
  at2=custom_dihe(2)
  at3=custom_dihe(3)
  at4=custom_dihe(4)

  dihe=dihedro(rclas(1,at1),rclas(2,at1),rclas(3,at1), &
               rclas(1,at2),rclas(2,at2),rclas(3,at2), &
               rclas(1,at3),rclas(2,at3),rclas(3,at3), &
               rclas(1,at4),rclas(2,at4),rclas(3,at4))

  if(dihe.gt.180) dihe=dihe-360.d0
  if(dihe.le.180) dihe=dihe+360.d0

  dihe=dihe*pi/180.d0

  if(dihe_type .eq. 1) then
    !compute energy
    Edihe=0.d0
    do i=1,ncos
      Edihe=Edihe+cos_weights(i)*(DCOS(dihe)**(i-1))
    enddo

    !summs energy to Etots
    Etots=Etots+Edihe*eV/kcal

    fnew=0.d0
    fdihe=0.d0
    !compute forces (kcal/mol)


    call diheforce2(natot,rclas,at1,at2,at3,at4,1,1.d0,fdihe)
    fnew(1:3,1)=fdihe(1:3)
    call diheforce2(natot,rclas,at1,at2,at3,at4,2,1.d0,fdihe)
    fnew(1:3,2)=fdihe(4:6)
    call diheforce2(natot,rclas,at1,at2,at3,at4,3,1.d0,fdihe)
    fnew(1:3,3)=fdihe(7:9)
    call diheforce2(natot,rclas,at1,at2,at3,at4,4,1.d0,fdihe)
    fnew(1:3,4)=fdihe(10:12)

    factor=0.d0

    do i=1,ncos
      factor=factor+cos_weights(i)*(i-1)*DCOS(dihe)**(i-2)
    enddo
    factor=factor*(-DSIN(dihe))

    fnew=factor*fnew
    fnew=fnew*eV/(Ang*kcal) ! Converts fnew to Hartree/bohr

    !summs forces to fdummy

    if((dihe.ge.0..and.dihe.le.180.).or.(dihe.gt.360)) then
      fnew=(-1.d0)*fnew
    elseif((dihe.gt.180..and.dihe.lt.360).or.(dihe.lt.0)) then
      fnew=fnew
    else
      stop 'custom dihe: Wrong dihedral angle value'
    endif

    fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
    fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
    fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)
    fdummy(1:3,at4)= fdummy(1:3,at4)+fnew(1:3,4)

  endif


  end subroutine custom_dihe_energy_forces
