      subroutine vmb(nat,ttemp,mass,xa,isa,va,cmcf,blockall,ntcon)
C *************************************************************************
C This routine creates a velocity distribution according to the
C Maxwell-Boltzman distribution with a given temperature.
C It also imposes the constraint of zero total velocity 
C (to avoid center of mass drifts).
C
C Temperatures in Kelvin
C Mass in atomic mass units
C Velocities in bohr/fs
C
C Writen by J. Junquera  Nov. 96,  P. Ordejon  Nov 97.
C Modified by P. Ordejon to take into account constraints, June 2003
C ************************** INPUT ****************************************
C integer nat              : Number of atoms
C real*8  ttemp            : Temperature desired for the distribution (K)
C real*8  mass(nat)        : Atomic masses (in amu)
C real*8  xa(3,nat)        : Atomic positions
C integer isa(nat)         : Species indexes
C ************************* OUTPUT ****************************************
C real*8  va(3,nat)        : Atomic Velocities
C *************************************************************************
C
C  Modules
C
      use precision
      use scarlett, only: Ang, eV, coord_freeze, natoms_partial_freeze

      implicit none

      integer 
     .  nat, isa(nat), cmcf, blockall(nat), ntcon

      double precision
     .  mass(nat), ttemp, va(3,nat), xa(3,nat)

C Internal variables ..................

      integer 
     .  iseed, i, ix, iy, Node

      double precision
     .  cell(3,3), stress(3,3), fa(3,nat), 
     .  massi, tempe, velo, vtot(3)

      external velo
C .....................
      Node = 0

C Only generate the velocities on one Node in case they are truely random
      if (Node.eq.0) then

        do ix=1,3
          vtot(ix)=0.0d0
        enddo


        iseed=-17

C Loop over atoms to assing velocities .................
        do i = 1,nat
          massi = mass(i)
           
          va(1,i) = velo(iseed,ttemp,massi)
          va(2,i) = velo(iseed,ttemp,massi)
          va(3,i) = velo(iseed,ttemp,massi)

       	  vtot(1) = vtot(1) + va(1,i)
       	  vtot(2) = vtot(2) + va(2,i)
       	  vtot(3) = vtot(3) + va(3,i)
        enddo
C ...............
C Impose constraint on fixed atoms JOTA
	do i = 1,nat
	  if(blockall(i) .eq. 1) then
	    va(1,i) = 0.d0
	    va(2,i) = 0.d0
	    va(3,i) = 0.d0
	  endif
	enddo

        do i=1,natoms_partial_freeze
          do ix=1,3
            if (coord_freeze(i,1+ix) .eq. 1) then
              va(ix,coord_freeze(i,1))=0.d0
            end if
          end do
        enddo




C Impose constraint on zero center of mass velocity ....................
	if(cmcf .ne. 0) then  !Agrego if, para que no lo haga siempre
        do i = 1,nat
          do ix=1,3
            va(ix,i) = va(ix,i) - vtot(ix)/nat
          enddo
        enddo
        endif
      endif



C Impose other constraints for MD, as specified in fixed.F ..................
C Note that fixed.F imposes constraints on forces; here we
C need those on velocities, so must transform multiplying by atomic mass
C
C Zero dummy matrices for stress and cell (needed in
C call to constrains routine)
c      do ix=1,3
c      do iy=1,3
c        cell(ix,iy)=0.0d0
c        stress(ix,iy)=0.0d0
c      enddo
c      enddo
c
c      do i=1,nat
c      do ix=1,3
c        fa(ix,i)=va(ix,i)*mass(i)
c      enddo
c      enddo
c      call fixed( cell, stress, nat, isa, mass, xa, fa,
c     .                  stress, fa, ntcon )
c      do i=1,nat
c      do ix=1,3
c      va(ix,i)=fa(ix,i)/mass(i)
c      enddo
c      enddo
C ...............

      if (nat .le. 1) return

C Correct velocity to exact temperature .....................
      call temp(2,Nat,mass,va,ntcon,tempe,cmcf)    !Esta hardcodeado con iunit=2; si iunit=3 la Ecin queda mal calculada pero la temp da ok por compensacion de errores; lo dejo así porque no escribe energía cinética en ningún lado...

      if (abs(tempe-ttemp) .ge. 1.e-4 .and. tempe .ge. 1.e-4) then
        do i = 1,Nat
          do ix=1,3
            va(ix,i) = va(ix,i) * dsqrt(ttemp/tempe)
          enddo
        enddo
      endif

      if (Node.eq.0) 
     .write(6,'(/,a)') 'vmb: creating a velocity distribution from MB' 

      return
      end

      real*8 function velo(iseed,temp,mass)
C *************************************************************************
C This function assigns velocities to atoms according to the 
C Maxwell-Boltzman distribution.
C It generates random numbers drawn from the normal distribution,
C using as input random numbers distributed uniformly from 0 to 1,  
C which are provided by routine ran3.
C
C Temperatures in Kelvin
C Mass in atomic mass units
C Velocities in bohr/fs
C
C Writen by J. Junquera  Nov. 96,  P. Ordejon  Nov 97.
C ************************** INPUT ****************************************
C integer iseed            : Seed for random number generator
C real*8  temp             : Temperature desired for the distribution (K)
C real*8  mass             : Atomic mass (in amu)
C ************************* OUTPUT ****************************************
C real*8  velo             : Velocity component 
C *************************************************************************

      implicit none

      integer 
     . iseed

      real*8
     . mass,temp

C Internal variables .................

      real*8
     .  arg1, arg2, gauss, med, pi, ran3, var
      
      external
     .  ran3

C ...........

C  For other distributions med may be different from cero.
      med = 0.0
      pi = acos(-1.0d0)
      var = sqrt(temp/mass)
C  conversion factor to bohr/fs
      var = var * 0.00172309d0
 
      arg1 = sqrt((-2.) * log(ran3(iseed)))

      arg2 = 2.0d0 * pi * ran3(iseed)
      gauss = arg1 * cos(arg2)

      velo = med + var * gauss

      return
      end

      subroutine temp(iunit,natoms,ma,va,ntcon,tempe,cmcf)
C *************************************************************************
C Subroutine to calculate instantaneous temperature
C
C Written by P.Ordejon, November'97
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Atomic masses in atomic mass units
C
C Space units depend on input option:
C
C   Option iunit = 1:
C     Distances are in Angstrom
C   Option iunit = 2:
C     Distances are in Bohr
C ************************* INPUT *********************************************
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 ma(natoms)     : Atomic masses 
C real*8 va(3,natoms)   : Atomic velocities 
C integer               : Total number of position constraints imposed
C ************************* OUTOPUT *******************************************
C real*8 tempe          : Instantaneous system temperature 
C *****************************************************************************
      use scarlett, only: Ang,eV

      implicit none

      integer 
     .   natoms,iunit,ntcon,cmcf

      double precision
     .  ma(natoms),tempe,va(3,natoms)

C Internal variables ..........................................................
 
      integer
     .  i,ia,ct

      double precision
     .  fovermp,kin

c saque Ang y eV
C ........................

      if (iunit .ne. 1 .and. iunit .ne. 2 .and. iunit .ne. 3) then
        write(6,*) 'temp: Wrong iunit option;  must be 1, 2 or 3'
        stop
      endif

C Define degrees of freedom fixed
      ct = cmcf + ntcon

C Define constants and conversion factors .....................................

C      Ang = 1.d0 / 0.529177d0
C      eV  = 1.d0 / 13.60580d0

        fovermp = 0.009648534d0
      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
c        fovermp = 0.009579038d0                      !Revise la cuenta y a mi me dio 0.009648534 Jota
      elseif (iunit .eq. 2) then
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = fovermp * Ang**2 / (eV*2.d0)
      else
C convert F/m in (Hartree/Bohr)/amu  to  Bohr/fs**2
        fovermp = 2.0d0 * fovermp * Ang**2 / (eV*2.d0)
      endif
C ........................

C Calculate kinetic energy and temperature ...................
C Kinetic energy of atoms
      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        tempe = 2.0d0*kin/dble(3*natoms-ct)/8.617d-5 
      elseif (iunit .eq. 2) then
        tempe = 2.0d0*kin/dble(3*natoms-ct)/8.617d-5/(eV*2.d0)
      else
       tempe = 2.0d0*2.0d0*kin/dble(3*natoms-ct)/8.617d-5/(eV*2.d0)
!JOTA PRUEBO SACAR EL 3
      endif




C .....................

      return
      end
    
