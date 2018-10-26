      subroutine verlet2(istep,iunit,iquench,natoms,fa,dt,ma,ntcon,va,
     .                   xa,kin,temp,nfree,cmcf)
C *************************************************************************
C Subroutine for MD simulations using the velocity-Verlet Algrithm.
C (See Allen-Tildesley, Computer Simulations of Liquids, pg. 81)
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Energies are in eV
C     Distances are in Angstrom
C   Option iunit = 2:
C     Energies are in Ry
C     Distances are in Bohr
C   Option iunit = 3:
C     Energies are in Hartree
C     Distances are in Bohr
C ************************* INPUT *********************************************
C integer istep         : Number of time step during the simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer iquench       : Option for quenching:
C                              0 = no quenching (standard dynamics)
C                              1 = power quenching (set to cero velocity
C                                 components opposite to force)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces 
C real*8 dt             : Length of the time step
C real*8 ma(natoms)     : Atomic masses 
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms)   : Atomic velocities
C                         (used only if istep = 1)
C integer cmcf          : Center of mass coordinates fixed JOTA
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms)   : Atomic coordinates 
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy at current time step 
C real*8 temp           : Instantaneous system temperature
C *****************************************************************************
C
C  Modules
C
      use precision
      use scarlett, only: Ang,eV
      implicit none

      integer 
     .   natoms,ntcon,istep,iquench,iunit,nfree,cmcf

      double precision
     .  dt,fa(3,natoms),kin,ma(natoms),
     .  va(3,natoms),xa(3,natoms)


C Internal variables ..........................................................
 
      integer
     .  ct,i,ia,Node

      double precision
     .  dot,dt2,dtby2,fovermp,temp
      
      double precision, dimension(:,:), allocatable, save ::
     .  accold,vold
C ........................


      if (iunit .ne. 1 .and. iunit .ne. 2 .and. iunit .ne. 3) then
          write(6,*) 'verlet2: Wrong iunit option;  must be 1, 2 or 3'
        stop
      endif
      
      ct = cmcf + ntcon 
      if (natoms .eq. 1) ct = 0

C Allocate local memory - only done once as data must be saved. As a
C result the memory is not deallocated at the end of the routine.
      if (.not.allocated(accold)) then
        allocate(accold(3,natoms))
C        call memory('A','D',3*natoms,'verlet2')
      endif
      if (.not.allocated(vold)) then
        allocate(vold(3,natoms))
C        call memory('A','D',3*natoms,'verlet2')
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      dtby2 = dt/2.0d0

      fovermp = 0.009648534d0
      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
      elseif (iunit .eq. 2) then
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = fovermp * Ang**2 / (eV*2.d0)
      else
C convert F/m in (Hartree/Bohr)/amu  to  Bohr/fs**2
        fovermp = 2.0d0 * fovermp * Ang**2 / (eV*2.d0)    
      endif
C .............................................................................

C Compute old accelerations and velocities
C  if the time step is the first of the simulation ...........................
      if (istep .eq. 1) then
        do ia = 1,natoms
          do i = 1,3
            accold(i,ia) = fovermp * fa(i,ia) / ma(ia)
            vold(i,ia) = va(i,ia) - dt * accold(i,ia)
          enddo
        enddo
      endif
C ....................

C Compute velocities at current time step ...................................
	do ia = 1,natoms
        do i = 1,3
          va(i,ia) = vold(i,ia) + dtby2 
     .               * (accold(i,ia) + fovermp * fa(i,ia) / ma(ia))
        enddo
      enddo

C Quench option if iquench = 0 ..............................................
      if (iquench .eq. 1) then

C Quench velocity components going uphill
         do ia = 1,natoms
           do i = 1,3
             dot = va(i,ia) * fa(i,ia)
             if (dot .lt. 0.0) va(i,ia) = 0.0
           enddo
         enddo

      endif
C ................

C Compute positions at next time step.....................................
      do ia = 1,natoms
        do i = 1,3
          xa(i,ia) = xa(i,ia) + dt * va(i,ia) 
     .                  + dt2 / 2.0d0 * fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ...................

C Save current velocities and accelerations as old ones .....................
      do i = 1,3
        do ia = 1,natoms
          vold(i,ia) = va(i,ia)
          accold(i,ia) = fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ....................

C Calculate kinetic energy and temperature at current time ...................
C Kinetic energy of atoms 
      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo
C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 2.0d0*kin/dble(3*nfree-ct)/8.617d-5
      elseif (iunit .eq. 2) then
        temp = 2.0d0*kin/dble(3*nfree-ct)/8.617d-5/(eV*2.d0)
      else
        temp = 2.0d0*2.0d0*kin/dble(3*nfree-ct)/8.617d-5/(eV*2.d0)
      endif
C .....................

      return
      end



      subroutine berendsen(istep,iunit,natoms,fa,dt,tau,ma,ntcon,va,xa,
     .                   kin,temp,ctemp,nfree,cmcf)

C *************************************************************************
C Subroutine for MD simulations using the Berendsen Algrithm.
C (See Allen-Tildesley, Computer Simulations of Liquids, pg. 231)
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Energies are in eV
C     Distances are in Angstrom
C   Option iunit = 2:
C     Energies are in Ry
C     Distances are in Bohr
C   Option iunit = 3:
C     Energies are in Hartree
C     Distances are in Bohr
C
C ************************* INPUT *********************************************
C integer istep         : Number of time step during the simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces 
C real*8 dt             : Length of the time step
C real*8 ma(natoms)     : Atomic masses 
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms)   : Atomic velocities
C                         (used only if istep = 1)
C real*8 tau                   : berensen coupling cosntant
C integer cmcf          : Center of mass coordinates fixed JOTA
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms)   : Atomic coordinates 
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy at current time step 
C real*8 temp           : Instantaneous system temperature
C *****************************************************************************
C
C  Modules
C
      use precision
      use scarlett, only: Ang,eV
      implicit none

      integer
     .   natoms,ntcon,istep,iunit,nfree,cmcf

      double precision
     .  dt,fa(3,natoms),kin,ma(natoms),
     .  va(3,natoms),xa(3,natoms)

c      external
c     .  memory

C Internal variables ..........................................................

      integer
     .  ct,i,ia,Node,regstep

      double precision
     .  dot,dt2,dtby2,fovermp,temp,ctemp,tau,chi

      double precision, dimension(:,:), allocatable, save ::
     .  accold, vold

      logical bervble
      save    bervble
      data    bervble /.false./


C ........................
C******************************
C TEST JOTA - SACAR AL TERMINAR
C      tau = tau*10000.d0
C******************************


      if (iunit .ne. 1 .and. iunit .ne. 2 .and. iunit .ne. 3) then
          write(6,*) 'verlet2: Wrong iunit option;  must be 1, 2 or 3'
        stop
      endif
      ct = cmcf + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory - only done once as data must be saved. As a
C result the memory is not deallocated at the end of the routine.
      if (.not.allocated(accold)) then
        allocate(accold(3,natoms))
c        call memory('A','D',3*natoms,'berendsen')
      endif
      if (.not.allocated(vold)) then
        allocate(vold(3,natoms))
c        call memory('A','D',3*natoms,'berendsen')
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      dtby2 = dt/2.0d0

      fovermp = 0.009648534d0
      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
      elseif (iunit .eq. 2) then
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = fovermp * Ang**2 / (eV*2.d0)
      else
C convert F/m in (Hartree/Bohr)/amu  to  Bohr/fs**2
        fovermp = 2.0d0 * fovermp * Ang**2 / (eV*2.d0)
C revisar jota <--------------------------------------------------------------------------------    
      endif
C ........................


C Compute old accelerations and velocities 
C  if the time step is the first of the simulation ...........................
      if (istep .eq. 1) then
        do ia = 1,natoms
          do i = 1,3
            accold(i,ia) = fovermp * fa(i,ia) / ma(ia)
            vold(i,ia) = va(i,ia) - dt * accold(i,ia)
          enddo
        enddo
      endif
C ....................

C Compute velocities at current time step ...................................
      do ia = 1,natoms
        do i = 1,3
          va(i,ia) = vold(i,ia) + dtby2
     .               * (accold(i,ia) + fovermp * fa(i,ia) / ma(ia))
        enddo
      enddo

C Calculate kinetic energy and temperature at current time ...................
C Kinetic energy of atoms 
      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo
C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 2.0d0*kin/dble(3*nfree-ct)/8.617d-5
      elseif (iunit .eq. 2) then
        temp = 2.0d0*kin/dble(3*nfree-ct)/8.617d-5/(eV*2.d0)
      else
        temp = 2.0d0*2.0d0*kin/dble(3*nfree-ct)/8.617d-5/(eV*2.d0)
      endif
C .....................

C Variable Berendsen Coupling
        if (istep .eq. 1) then
        if (tau.eq.0.0) then
        bervble=.true.
        endif
        endif

        if(bervble) then
        regstep=200
        if(istep.le.regstep) then
        tau=real((istep-1)*(10.-dt)/(regstep-1)+dt)
        else
        tau=10.
        endif
        endif

C Calculate the Berendsen bath coupling and modify va
      chi= sqrt(1.+((dt/tau)*((ctemp/temp)-1.)))
      do i = 1,3
        do ia = 1,natoms
          va(i,ia) = va(i,ia) * chi
        enddo
      enddo

C Compute positions at next time step.....................................
      do ia = 1,natoms
        do i = 1,3
          xa(i,ia) = xa(i,ia) + dt * va(i,ia)
     .                  + dt2 / 2.0d0 * fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ...................

C Save current velocities and accelerations as old ones .....................
      do i = 1,3
        do ia = 1,natoms
          vold(i,ia) = va(i,ia)
          accold(i,ia) = fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ....................

C Calculate kinetic energy and temperature at current time ...................
C Kinetic energy of atoms 

      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo


C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 2.0d0*kin/dble(3*nfree-ct)/8.617d-5
      elseif (iunit .eq. 2) then
        temp = 2.0d0*kin/dble(3*nfree-ct)/8.617d-5/(eV*2.d0)
      else
        temp = 2.0d0*2.0d0*kin/dble(3*nfree-ct)/8.617d-5/(eV*2.d0)
      endif
C .....................

      return
      end

C***********************************************************************
C***********************************************************************
C** Subroutine that calculates intantaneous temperature of the system **
C***********************************************************************
C***********************************************************************
C iunit = 3

      subroutine calculateTemp(kin,temp,tempqm,va,ntcon,
     . nfree,cmcf)

      use precision
      use scarlett, only: Ang,eV,na_u,natot,nac,masst,qm,mm
      implicit none

      double precision, intent(inout) :: kin
      double precision, intent(out) :: temp, tempqm
      double precision, intent(in) :: va(3,natot)
      integer, intent(in) :: ntcon, cmcf, nfree
     
C internal variables
      integer :: ia, i, ct
      double precision :: fovermp, kinqm

      ct = ntcon + cmcf
        fovermp = 2.0d0 * 0.009648534d0 * Ang**2 / (eV*2.d0)

      kinqm = 0.0d0
      do ia = 1,na_u
        do i = 1,3
          kinqm = kinqm + 0.5d0 * masst(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

      kin = kinqm
      do ia = na_u+1,natot
        do i = 1,3
          kin = kin + 0.5d0 * masst(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

        tempqm = 2.0d0*2.0d0*kinqm/dble(3*na_u)/8.617d-5/(eV*2.d0)
        temp = 2.0d0*2.0d0*kin/dble(3*nfree-ct)/8.617d-5/(eV*2.d0)
      return
      end subroutine calculateTemp

C************************************************************************

      subroutine nose(istep,natoms,fa,tt,dt,
     .                ma,mn,ntcon,va,xa,kin,kn,vn,
     .                temp,nfree,cmcf)
C *************************************************************************
C Subroutine for MD simulations with CONTROLLED TEMPERATURE.
C The temperature is controlled with a NOSE thermostat.
C The use of Nose dynamics provides trajectories which sample the 
C isothermal ensamble.
C The equations of motion are integrated with a modified Verlet 
C algorithm, with a selfconsistent set of equations for the 
C Nose variables.
C
C Written by P.Ordejon, November'96
C ************************* UNITS ******************************************
C (Como si fuera iunit = 3 JOTA)
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C Energies are in Hartree
C Distances are in Bohr
C ************************* INPUT *********************************************
C integer istep         : Number of time step during simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 fa(3,natoms)   : Atomic forces
C real*8 tt             : Target temperature
C real*8 dt             : Length of the time step 
C real*8 ma(natoms)     : Atomic masses
C real*8 mn             : Mass of Nose thermostat
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms)   : Atomic velocities
C                         (used only if istep = 1)
C integer cmcf          : Center of mass coordinates fixed JOTA
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms)   : Atomic coordinates
C                        (input: current time step; output: next time step)
C ************************* OUTOPUT *******************************************
C real*8 kin            : Kinetic energy of the atomic system
C real*8 kn             : Kinetic energy of Nose variable 
C real*8 vn             : Potential energyy of Nose var
C real*8 temp           : Instantaneous system temperature 
C *****************************************************************************
C
C  Modules
C
      use precision
      use scarlett, only: Ang,eV
      implicit none

      integer
     .   natoms,ntcon,istep,nfree,cmcf

      double precision
     .  dt,fa(3,natoms),kin,kn,
     .  ma(natoms),mn,tt,
     .  va(3,natoms),vn,xa(3,natoms)

c      external
c     .  memory
C Internal variables .........................................................

      integer
     .  ct,i,ia,Node

      save x,xold

      double precision
     .  diff,dt2,dtby2,fact,fovermp,
     .  tekin,temp,tol,twodt,
     .  x,xdot,xlast,xnew,xold

      double precision, dimension(:,:), allocatable, save ::
     .  xanew,xaold
C .............................................................................

      ct = cmcf + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory and initialise
      if (.not.allocated(xanew)) then
        allocate(xanew(3,natoms))
c        call memory('A','D',3*natoms,'nose')
      endif
      if (.not.allocated(xaold)) then
        allocate(xaold(3,natoms))
c        call memory('A','D',3*natoms,'nose')
        do ia = 1,natoms
          do i = 1,3
            xaold(i,ia)=0.0d0
          enddo
        enddo
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      dtby2 = dt/2.0d0
      twodt = dt*2.0d0
      tol   = 1.0d-12

      fovermp = 0.009648534d0*Ang**2/eV
C  fovermp converts F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2

C  convert target ionic temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in eV if Temp is in Kelvin)

c	write(7171,*) tt, nfree, ct 

      tekin = 0.5d0*dble(3*nfree-ct)*8.617d-5*eV*tt

C Initialize variables if current time step is the first of the simulation

      if (istep .eq. 1) then
        x = 0.0d0
        xold = 0.0d0
        do ia = 1,natoms
          do i = 1,3
            xaold(i,ia) = xa(i,ia) - dt * va(i,ia)
     .                + (dt2/2.0d0) * fovermp * fa(i,ia) / ma(ia)
          enddo
        enddo
      endif

C Compute uncorrected next positions .....................................
      do ia = 1,natoms
        do i = 1,3
          xanew(i,ia) =  2.0d0 * xa(i,ia) - xaold(i,ia) +
     .                   dt2 * fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo

C Compute uncorrected velocities and kinetic energy ......................
      kin = 0.d0
      do ia = 1,natoms
        do i = 1,3
          va(i,ia) = (xanew(i,ia) - xaold(i,ia)) / twodt
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

C Compute initial guess for Nose variables at next time step .............
      xnew = 2.0d0 * x - xold

C Start selfconsistency loop to calculate Nose variable ..................
10    continue

      xlast = xnew

C xdot and hdot (time derivatives at current time), and related stuff
      xdot = (xnew - xold) / twodt
      fact = (1.0/(1.0+xdot*dtby2))

C  Compute Nose variable for next iteration
      xnew = 2.0d0 * x - xold
     .       + (dt2/mn) * 2.0 * (fact**2 * kin - tekin)

C Check if selfconsistency has been reached
      diff = abs(xnew - xlast)
      if (xlast .eq. 0.0d0) then
        if (diff .gt. tol)  goto 10
      else
        if (diff/abs(xlast) .gt. tol)  goto 10
      endif
C ...................

C Calculate corrected atomic coordinates at next time step, 
C and corrected velocities and kinetic energy at current time step .........
      do ia = 1,natoms
        do i = 1,3
          xanew(i,ia) = fact * ( xanew (i,ia) +
     .                   dtby2 * xdot * xaold(i,ia))
          va(i,ia) = fact * va(i,ia)
        enddo
      enddo
      kin = kin * fact**2
C ...................

C Save current atomic positions as old ones, 
C   and next positions as current ones
      do i = 1,3
        do ia = 1,natoms
          xaold(i,ia) = xa(i,ia)
          xa(i,ia) = xanew(i,ia)
        enddo
      enddo

      xold = x
      x = xnew

C Calculate Kinetic and potential energies ................................
C Kinetic energy of Nose variable
      kn = (1.0d0 / 2.0d0) * mn * xdot**2

C Potential energy of Nose variable (in eV)
      vn = 2.0d0 * tekin * xold

C Instantaneous temperature (Kelvin)

      temp = 2.0d0*kin/dble(3*nfree-ct)/8.617d-5/eV

C .....................

      return
      end subroutine nose

   
