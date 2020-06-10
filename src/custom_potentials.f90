!****************************************************************************************
! This subroutine reads fdf specifications regarding custom potentials.
! Additional custom potential options should be added here.
! Options available
! custompot_type = 1 Special potential for a certain torsion
!   dihe_type = 1 Ryckaert-Belleman dihedral: sum from i=1 to ncos V(i)*cos(dihe)**(i-1)
! custompot_type = 2 London-Eyring-Polanyi-Sato (LEPS) Potential for a three-site complex
! J. Semelak 2020
!****************************************************************************************

  subroutine custom_potentials_assign(custompot_type)

  use fdf
  use scarlett, only: custom_dihe, cos_weights, ncos, dihe_type, &
  lepsmask,leps1D,leps1B,leps1R0,leps3D,leps3B,leps3R0
  implicit none
  integer :: i, iunit, custompot_type
  character :: exp

  if (fdf_block('CustomPotential',iunit)) then
    read(iunit,*) exp, custompot_type
    if (custompot_type .eq. 1) then ! Special potential for a certain torsion
      read(iunit,*) exp,(custom_dihe(i),i=1,4)
      read(iunit,*) exp, dihe_type
      write(6,'(/,a)')  'Custom dihedral potential is turned on'
      write(6,*)  'Atoms integrating the custom dihedral:', custom_dihe(1:4)
      if (dihe_type .eq. 1) then ! Ryckaert-Belleman dihedral
        write(6,'(/,a)')  'Potential type: Cosine series'
        read(iunit,*) exp, ncos
        allocate(cos_weights(ncos))
        read(iunit,*) exp, (cos_weights(i),i=1,ncos)
         write(6,*)  'Cosine weights (kcal/mol): ', cos_weights(1:ncos)
      endif
    elseif (custompot_type .eq. 2) then ! LEPS
      read(iunit,*) exp,(lepsmask(i),i=1,3)
      read(iunit,*) exp,(leps1D(i),i=1,3)
      read(iunit,*) exp,(leps1B(i),i=1,3)
      read(iunit,*) exp,(leps1R0(i),i=1,3)
      read(iunit,*) exp,(leps3D(i),i=1,3)
      read(iunit,*) exp,(leps3B(i),i=1,3)
      read(iunit,*) exp,(leps3R0(i),i=1,3)
      write(6,'(/,a)')  'LEPS potential is turned on'
      write(6,*)  'Atoms integrating thethree-site complex: ', lepsmask(1:3)
      write(6,*)  'Parameters:'
      write(6,*)  '1D ' ,leps1D(1:3), ' 3D ' ,leps3D(1:3)
      write(6,*)  '1B ' ,leps1B(1:3), ' 3B ' ,leps3B(1:3)
      write(6,*)  '1R0 ' ,leps3R0(1:3), ' 3R0 ' ,leps3R0(1:3)
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
  double precision :: dihe, Edihe, factor, dihedro2, pi
  double precision, dimension(ncos) :: Vcos
  double precision, dimension(3,10) :: fnew
  double precision, dimension(12) :: fdihe

  pi=DACOS(-1.d0)
  at1=custom_dihe(1)
  at2=custom_dihe(2)
  at3=custom_dihe(3)
  at4=custom_dihe(4)

  dihe=dihedro2(rclas(1,at1),rclas(2,at1),rclas(3,at1), &
               rclas(1,at2),rclas(2,at2),rclas(3,at2), &
               rclas(1,at3),rclas(2,at3),rclas(3,at3), &
               rclas(1,at4),rclas(2,at4),rclas(3,at4))

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

    fnew=-fnew !diheforce2 subroutine does not include the minus sign regarding F=-dV
    factor=0.d0

    do i=1,ncos
      factor=factor+cos_weights(i)*(dble(i)-1.d0)*DCOS(dihe)**(i-2)
    enddo
    factor=factor*(-DSIN(dihe))

    fnew=factor*fnew
    fnew=fnew*eV/(Ang*kcal) ! Converts fnew to Hartree/bohr

    !summs forces to fdummy
    fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
    fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
    fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)
    fdummy(1:3,at4)= fdummy(1:3,at4)+fnew(1:3,4)
  endif

  end subroutine custom_dihe_energy_forces

  subroutine leps_energy_forces(Etots)
    use scarlett, only: lepsmask,leps1D,leps1B,leps1R0,leps3D,leps3B,leps3R0, &
    rclas, fdummy, natot, eV, Ang, kcal

    implicit none
    double precision, intent(inout) :: Etots
    double precision :: dist, Jprod, Eleps, JA, JB, JC
    integer ::  i, j, at1, at2, at3
    double precision, dimension(3) :: lepsR, leps1E, leps3E, lepsQ, lepsJ
    double precision, dimension(3) :: d1Edr, d3Edr, dQdr, dJdr
    double precision, dimension(3,3) :: dr1dq, dr2dq, dr3dq, dJproddq
    double precision, dimension(3,3) :: fnew

    at1=lepsmask(1)
    at2=lepsmask(2)
    at3=lepsmask(3)

    lepsR(1)=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1), &
             rclas(1,at2),rclas(2,at2),rclas(3,at2))
    lepsR(2)=dist(rclas(1,at2),rclas(2,at2),rclas(3,at2), &
             rclas(1,at3),rclas(2,at3),rclas(3,at3))
    lepsR(3)=dist(rclas(1,at1),rclas(2,at1),rclas(3,at1), &
             rclas(1,at3),rclas(2,at3),rclas(3,at3))



    lepsR=lepsR/Ang



    !compute energy
    do i=1,3
      leps1E(i)=leps1D(i)*((1-exp(-leps1B(i)*(lepsR(i)-leps1R0(i))))**2.d0)-leps1D(i)
      leps3E(i)=leps3D(i)*((1+exp(-leps3B(i)*(lepsR(i)-leps3R0(i))))**2.d0)-leps3D(i)
    end do

    do i=1,3
      lepsQ(i)=(leps1E(i)+leps3E(i))/2.d0
      lepsJ(i)=(leps1E(i)-leps3E(i))/2.d0
    end do

    Jprod=lepsJ(1)*lepsJ(1)+lepsJ(2)*lepsJ(2)+lepsJ(3)*lepsJ(3)- &
          lepsJ(1)*lepsJ(2)-lepsJ(2)*lepsJ(3)-lepsJ(3)*lepsJ(1)
    Jprod=dsqrt(Jprod)

    Eleps=lepsQ(1)+lepsQ(2)+lepsQ(3)-Jprod

    !summs energy to Etots
    Etots=Etots+Eleps*eV/kcal



    !compute forces

    !drdq(coordinate,at)
    do i=1,3
      dr1dq(i,1)=(rclas(i,at1)-rclas(i,at2))/lepsR(1)
      dr1dq(i,2)=-dr1dq(i,1)
      dr1dq(i,3)=0.d0
      dr2dq(i,1)=0.d0
      dr2dq(i,2)=(rclas(i,at2)-rclas(i,at3))/lepsR(2)
      dr2dq(i,3)=-dr2dq(i,2)
      dr3dq(i,1)=(rclas(i,at1)-rclas(i,at3))/lepsR(3)
      dr3dq(i,2)=0
      dr3dq(i,3)=-dr3dq(i,1)
    end do

    dr1dq=dr1dq/Ang
    dr2dq=dr2dq/Ang
    dr3dq=dr3dq/Ang



     do i=1,3
       d1Edr(i)=leps1D(i)*2.d0*(1-exp(-leps1B(i)*(lepsR(i)-leps1R0(i))))*leps1B(i)*exp(-leps1B(i)*(lepsR(i)-leps1R0(i)))
       d3Edr(i)=-leps3D(i)*2.d0*(1+exp(-leps3B(i)*(lepsR(i)-leps3R0(i))))*leps3B(i)*exp(-leps3B(i)*(lepsR(i)-leps3R0(i)))
     end do

     do i=1,3
       dQdr(i)=((d1Edr(i)+d3Edr(i))/2.d0)
       dJdr(i)=((d1Edr(i)-d3Edr(i))/2.d0)
     end do

   !dJproddq(coordinate,at)

     JA=2.d0*lepsJ(1)-lepsJ(2)-lepsJ(3)
     JB=2.d0*lepsJ(2)-lepsJ(1)-lepsJ(3)
     JC=2.d0*lepsJ(3)-lepsJ(2)-lepsJ(1)

    do i=1,3
      do j=1,3
        dJproddq(i,j)= JA*dJdr(1)*dr1dq(i,j)+JB*dJdr(2)*dr2dq(i,j)+JC*dJdr(3)*dr3dq(i,j)
      end do
    end do


    do i=1,3
      do j=1,3
        fnew(i,j)=dQdr(1)*dr1dq(i,j)+dQdr(2)*dr2dq(i,j)+dQdr(3)*dr3dq(i,j)- &
                  (1.d0/(2.d0*Jprod))*dJproddq(i,j)
      end do
     end do

     fnew=-fnew
     fnew=fnew*eV/(Ang*kcal) ! Converts fnew to Hartree/bohr

     !summs forces to fdummy
     fdummy(1:3,at1)= fdummy(1:3,at1)+fnew(1:3,1)
     fdummy(1:3,at2)= fdummy(1:3,at2)+fnew(1:3,2)
     fdummy(1:3,at3)= fdummy(1:3,at3)+fnew(1:3,3)

  end subroutine leps_energy_forces
