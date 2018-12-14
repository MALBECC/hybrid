      subroutine nose(istep,iunit,natoms,fa,tt,dt,
     .                ma,mn,ntcon,va,xa,kin,kn,vn,
     .                temp,nfree)
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
C Temperature in Kelvin
C Time in femtoseconds
C Atomic masses in atomic mass units
C
C Other units depend on input option:
C
C   Option iunit = 1:
C     Energies are in eV
C     Distances are in Angstrom
C     Nose mass in eV * fs**2
C   Option iunit = 2:
C     Energies are in Ry
C     Distances are in Bohr
C     Nose mass in Ry * fs**2
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
      implicit none

      integer 
     .   natoms,ntcon,istep,iunit,nfree

      double precision
     .  dt,fa(3,natoms),kin,kn,
     .  ma(natoms),mn,tt,
     .  va(3,natoms),vn,xa(3,natoms)

      external
     .  memory
C Internal variables .........................................................

      integer
     .  ct,i,ia,Node

      save x,xold

      double precision
     .  Ang,diff,dt2,dtby2,eV,fact,fovermp,
     .  tekin,temp,tol,twodt,
     .  x,xdot,xlast,xnew,xold

      double precision, dimension(:,:), allocatable, save ::
     .  xanew,xaold
C .............................................................................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
          write(6,*) 'nose: Wrong iunit option;  must be 1 or 2'
        stop
      endif
      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory and initialise
      if (.not.allocated(xanew)) then
        allocate(xanew(3,natoms))
        call memory('A','D',3*natoms,'nose')
      endif
      if (.not.allocated(xaold)) then
        allocate(xaold(3,natoms))
        call memory('A','D',3*natoms,'nose')
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

      Ang = 1.d0 / 0.529177d0
      eV  = 1.d0 / 13.60580d0

      if (iunit .eq. 1) then
C  convert target ionic temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in eV if Temp is in Kelvin)
         tekin = 0.5d0 * (3.d0 * nfree - ct) * 8.617d-5 * tt
C  convert F/m in (eV/Angstrom)/amu  to  Angstrom/fs**2
        fovermp = 0.009579038
      else
C  convert target temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in Ry if Temp is in Kelvin)
        tekin = eV * 0.5d0 * (3.d0 * nfree - ct) * 8.617d-5 * tt
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif
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
C ..................

C Compute uncorrected next positions .....................................
      do ia = 1,natoms
        do i = 1,3
          xanew(i,ia) =  2.0d0 * xa(i,ia) - xaold(i,ia) +
     .                   dt2 * fovermp * fa(i,ia) / ma(ia)
        enddo
      enddo
C ...................

C Compute uncorrected velocities and kinetic energy ......................
      kin = 0.d0
      do ia = 1,natoms
        do i = 1,3
          va(i,ia) = (xanew(i,ia) - xaold(i,ia)) / twodt
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo
C ..................

C Compute initial guess for Nose variables at next time step .............
      xnew = 2.0d0 * x - xold
C ...................

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
      if (iunit .eq. 1) then
        temp = kin / (0.5d0 * (3.d0 * nfree - ct) * 8.617d-5)
      else
        temp = kin / (0.5d0 * (3.d0 * nfree - ct) * 8.617d-5 * eV)
      endif
 
C .....................

      return
      end

      subroutine anneal(istep,iunit,ianneal,taurelax,bulkm,
     .               natoms,fa,stress,tp,tt,dt,
     .               ma,ntcon,va,xa,h,kin,
     .               temp,pressin,nfree)
C *************************************************************************
C Subroutine for MD simulations with a TARGET PRESSURE AND TEMPERATURE.
C The system is driven to a desired temperature and pressure in
C a given time, by rescaling the velocities and the cell shape and size.
C It needs an estimate of the bulk modulus of the system, to determine
C the rate of change of cell shape to accomodate to the target
C pressure in the required time. A wrong estimate will simply
C drive the system to the desired pressure, but in a different time
C than the especified in the input. (See Kittel for representative
C values of bulk moduli for materials).
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
C     Pressures are in eV/A**3 ( = 160.20506 GPa)
C     Energies are in eV
C     Distances are in Angstrom
C     Bulk modulus in eV/A**3
C   Option iunit = 2:
C     Pressures are in Ry/Bohr**3 
C     Energies are in Ry
C     Distances are in Bohr
C     Bulk modulus in Ry/Bohr**3
C ************************* INPUT *********************************************
C integer istep       : Number of time step during simulation
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer ianneal     : Work mode option:
C                       1 = reach target temperature only
C                       2 = reach target pressure only
C                       3 = reach target temperature and pressure
C real*8 taurelax     : Relaxation time to reach desired T and P
C real*8 bulkm        : Estimate of the Bulk Modulus of the system
C integer natoms      : Number of atoms in the simulation cell
C real*8 fa(3,natoms) : Atomic forces 
C real*8 stress(3,3)  : Stress tensor components 
C real*8 tp           : Target pressure 
C real*8 tt           : Target temperature 
C real*8 dt           : Length of the time step 
C real*8 ma(natoms)   : Atomic masses 
C integer ntcon         : Total number of position constraints imposed
C real*8 va(3,natoms) : Atomic velocities
C                       (used only if istep = 1)
C ******************* INPUT AND OUTPUT ****************************************
C real*8 xa(3,natoms) : Atomic coordinates 
C                      (input: current time step; output: next time step)
C real*8 h(3,3)       : Matrix of the vectors defining the unit
C                       cell at the current time step 
C                       h(i,j) is the ith component of jth basis vector
C                      (input: current time step; output: next time step)
C ************************* OUTPUT *******************************************
C real*8 kin            : Kinetic energy of the atomic system 
C real*8 temp           : Instantaneous system temperature 
C real*8 pressin        : Instantaneous system pressure 
C *****************************************************************************
C
C  Modules
C
      use precision
      implicit none

      integer 
     .   natoms,ntcon,istep,ianneal,iunit,nfree

      double precision
     .  bulkm,dt,fa(3,natoms),h(3,3),kin,
     .  ma(natoms),stress(3,3),taurelax,tp,tt,
     .  va(3,natoms),xa(3,natoms)

C Internal variables .............................................................

      integer
     .  ct,i,ia,info,j,k,Node

      double precision
     .  Ang,dt2,eV,
     .  fovermp,
     .  hi(3,3),hs(3),
     .  pgas,press(3,3),pressin,rfac,rfac2,
     .  tekin,temp,twodt,
     .  vol,volcel

      double precision, dimension(:,:), allocatable, save ::
     .  s,sdot,snew,sold,sunc

       external
     .  volcel, memory
C ....................................................................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
          write(6,*) 'anneal: Wrong iunit option;  must be 1 or 2'
        stop
      endif
      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory
      allocate(s(3,natoms))
      call memory('A','D',3*natoms,'anneal')
      allocate(sdot(3,natoms))
      call memory('A','D',3*natoms,'anneal')
      allocate(snew(3,natoms))
      call memory('A','D',3*natoms,'anneal')
      allocate(sunc(3,natoms))
      call memory('A','D',3*natoms,'anneal')
      if (.not.allocated(sold)) then
        allocate(sold(3,natoms))
        call memory('A','D',3*natoms,'anneal')
      endif

C Define constants and conversion factors .......................................
      dt2   = dt**2
      twodt = dt*2.0d0

      Ang = 1.d0 / 0.529177d0
      eV  = 1.d0 / 13.60580d0

      if (iunit .eq. 1) then
C  convert target ionic temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in eV if Temp is in Kelvin)
        tekin = 0.5d0 * (3.d0 * nfree - ct) * 8.617d-5 * tt
C  convert F/m in (eV/Angstrom)/amu  to  Angstrom/fs**2
        fovermp = 0.009579038
      else
C  convert target temperature into target kinetic energy
C  Ekin=1/2*(3N-3)*kB*Temp  (yields Ekin in Ry if Temp is in Kelvin)
        tekin = eV * 0.5d0 * (3.d0 * nfree - ct) * 8.617d-5 * tt
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
      endif

C  calculate cell volume at current time
      vol = volcel( h )
C ........................

C Compute Parrinello-Rahman variables (H and scaled coordinates) ............
C Compute Inverse of H and G at current time 
      call inver(h,hi,3,3,info)
      if (info .ne. 0) stop 'anneal: INVER failed'

C Calculate scaled coordinates (referred to matrix H) at current time
      do ia = 1,natoms
        do i = 1,3
          s(i,ia) = 0.0d0
          do j = 1,3
            s(i,ia) = s(i,ia) + hi(i,j) * xa(j,ia)
          enddo
        enddo
      enddo

C Initialize variables if current time step is the first of the simulation
      if (istep .eq. 1) then
        do ia = 1,natoms
          do i = 1,3
            sold(i,ia) = 0.0d0
            do j = 1,3
              sold(i,ia) = sold(i,ia) + hi(i,j)*(xa(j,ia)-dt*va(j,ia)
     .                     + (dt2/2.0d0) * fovermp * fa(j,ia) / ma(ia))
            enddo
          enddo
        enddo
      endif
C ..................

C Compute uncorrected next positions .....................................
      do ia = 1,natoms
        do i = 1,3
          sunc(i,ia) = -sold(i,ia) + 2.0d0 * s(i,ia)
          do k = 1,3
            sunc(i,ia) = sunc(i,ia) + 
     .                   dt2 * hi(i,k) * fovermp * fa(k,ia) / ma(ia)
          enddo
        enddo
      enddo
C ...................

C Calculate uncorrected velocities at current time
      do ia = 1,natoms
        do i = 1,3
          sdot(i,ia) = (sunc(i,ia) - sold(i,ia)) / twodt
        enddo
      enddo

C Calculate pressure tensor at current time and ideal gas pressure
      do i = 1,3
        do j = 1,3
          press(i,j) = 0.0d0
        enddo
      enddo
      do ia = 1,natoms
        do i = 1,3
          hs(i) = 0.0d0
          do j = 1,3
            hs(i) = hs(i) + h(i,j) * sdot(j,ia)
          enddo
        enddo
        do j = 1,3
          do i = 1,3
            press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
          enddo
        enddo
      enddo
      pgas = 0.0d0
      do i = 1,3
        pgas = pgas + press(i,i) / vol
      enddo
      pgas = pgas / 3.0d0
      do j = 1,3
        do i = 1,3
          press(i,j) = press(i,j) / vol - stress(i,j)
        enddo
      enddo

C Compute internal pressure  (pressin = 1/3 Tr (press))   at current time
      pressin = 0.0d0
      do i = 1,3
        pressin = pressin + press(i,i)
      enddo
      pressin = pressin / 3.0d0

C Compute kinetic energy
      kin = (3.0d0 / 2.0d0) * pgas * vol

      if (ianneal .eq. 1 .or. ianneal .eq. 3) then
C Correct velocities to reach target termperature
      if (kin .eq. 0.0) then
        rfac2 = 1.0d0 + dt/taurelax
      else
        rfac2 = (1.0d0 + dt/taurelax * (tekin/kin -1.0d0))
      endif
      if (rfac2 .le. 0.0) then
        write(6,*) 'Wrong anneal parameter'
        stop
      endif
      rfac = sqrt(rfac2)

      do ia = 1,natoms
        do i = 1,3
          sdot(i,ia) = rfac * sdot(i,ia) 
        enddo
      enddo

C Compute again pressure, with corrected velocities
      do i = 1,3
        do j = 1,3
          press(i,j) = 0.0d0
        enddo
      enddo
      do ia = 1,natoms
        do i = 1,3
          hs(i) = 0.0d0
          do j = 1,3
            hs(i) = hs(i) + h(i,j) * sdot(j,ia)
          enddo
        enddo
        do j = 1,3
          do i = 1,3
            press(i,j) = press(i,j) + ma(ia) * hs(i) * hs(j) / fovermp
          enddo
        enddo
      enddo
      pgas = 0.0d0
      do i = 1,3
        pgas = pgas + press(i,i) / vol
      enddo
      pgas = pgas / 3.0d0
      do j = 1,3
        do i = 1,3
          press(i,j) = press(i,j) / vol - stress(i,j)
        enddo
      enddo

C Compute internal pressure  (pressin = 1/3 Tr (press))   at current time
      pressin = 0.0
      do i = 1,3
        pressin = pressin + press(i,i)
      enddo
      pressin = pressin / 3.0d0
      endif


C Correct new possitions according to corrected velocities
      do ia = 1,natoms
        do i = 1,3
          snew(i,ia) = sold(i,ia) + twodt * sdot(i,ia)
        enddo
      enddo

      if (ianneal .eq. 2 .or. ianneal .eq. 3) then
C Correct cell shape to reach target pressure
      do i = 1,3
        do j = 1,3
          if (i .ne. j) then
            rfac2 = 1.0 + (dt / taurelax) * press(i,j) 
     .                    / (0.5 * bulkm)
          else
            rfac2 = 1.0 + (dt / taurelax) * (press(i,i) - tp) 
     .                    / (0.5 * bulkm)
          endif
          if (rfac2 .le. 0.0) then
            write(6,*) 'Wrong anneal parameter'
            stop
          endif
          rfac = sqrt(rfac2)
          h(i,j) = rfac * h(i,j)
        enddo
      enddo
      endif

C Save current atomic positions as old ones, 
C   and next positions as current ones
      do i = 1,3
        do ia = 1,natoms
          sold(i,ia) = s(i,ia)
          s(i,ia) = snew(i,ia)
        enddo
      enddo

C Transform back to absolute coordinates 
      do ia = 1,natoms
        do i = 1,3
          xa(i,ia) = 0.0d0
          va(i,ia) = 0.0d0
          do j = 1,3
            xa(i,ia) = xa(i,ia) + h(i,j) * s(j,ia)
            va(i,ia) = va(i,ia) + h(i,j) * sdot(j,ia)
          enddo
        enddo
      enddo
C ....................

C Calculate Kinetic and potential energies ................................
C Kinetic energy of atoms 
      kin = (3.0d0 / 2.0d0) * pgas * vol

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        temp = 3.0d0*vol*pgas/(3.0d0*nfree-ct)/8.617d-5
      else
        temp = 3.0d0*vol*pgas/(3.0d0*nfree-ct)/8.617d-5/eV
      endif

C .....................

C Deallocate local memory
      call memory('D','D',size(s),'anneal')
      deallocate(s)
      call memory('D','D',size(sdot),'anneal')
      deallocate(sdot)
      call memory('D','D',size(snew),'anneal')
      deallocate(snew)
      call memory('D','D',size(sunc),'anneal')
      deallocate(sunc)

      return
      end

      subroutine verlet2(istep,iunit,iquench,natoms,fa,dt,ma,ntcon,va,
     .                   xa,kin,temp,nfree)
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
      implicit none

      integer 
     .   natoms,ntcon,istep,iquench,iunit,nfree

      double precision
     .  dt,fa(3,natoms),kin,ma(natoms),
     .  va(3,natoms),xa(3,natoms)

      external
     .  memory

C Internal variables ..........................................................
 
      integer
     .  ct,i,ia,Node

      double precision
     .  Ang,dot,dt2,dtby2,eV,fovermp,temp

      double precision, dimension(:,:), allocatable, save ::
     .  accold,vold
C ........................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
          write(6,*) 'verlet2: Wrong iunit option;  must be 1 or 2'
        stop
      endif
      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory - only done once as data must be saved. As a
C result the memory is not deallocated at the end of the routine.
      if (.not.allocated(accold)) then
        allocate(accold(3,natoms))
        call memory('A','D',3*natoms,'verlet2')
      endif
      if (.not.allocated(vold)) then
        allocate(vold(3,natoms))
        call memory('A','D',3*natoms,'verlet2')
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      dtby2 = dt/2.0d0

      Ang = 1.d0 / 0.529177d0
      eV  = 1.d0 / 13.60580d0

      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
        fovermp = 0.009579038
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
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
        temp = 2.0d0*kin/(3.0d0*nfree-ct)/8.617d-5
      else
        temp = 2.0d0*kin/(3.0d0*nfree-ct)/8.617d-5/eV
      endif

C .....................

      return
      end

      subroutine berendsen(istep,iunit,natoms,fa,dt,tau,ma,ntcon,va,xa,
     .                   kin,temp,ctemp,nfree)

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
      implicit none

      integer 
     .   natoms,ntcon,istep,iunit,nfree

      double precision
     .  dt,fa(3,natoms),kin,ma(natoms),
     .  va(3,natoms),xa(3,natoms)

      external
     .  memory

C Internal variables ..........................................................
 
      integer
     .  ct,i,ia,Node,regstep
 
      double precision
     .  Ang,dot,dt2,dtby2,eV,fovermp,temp,ctemp,tau,chi
 
      double precision, dimension(:,:), allocatable, save ::
     .  accold, vold
 
      logical bervble
      save    bervble
      data    bervble /.false./
C ........................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
          write(6,*) 'verlet2: Wrong iunit option;  must be 1 or 2'
        stop
      endif
      ct = 3 + ntcon
      if (natoms .eq. 1) ct = 0

C Allocate local memory - only done once as data must be saved. As a
C result the memory is not deallocated at the end of the routine.
      if (.not.allocated(accold)) then
        allocate(accold(3,natoms))
        call memory('A','D',3*natoms,'berendsen')
      endif
      if (.not.allocated(vold)) then
        allocate(vold(3,natoms))
        call memory('A','D',3*natoms,'berendsen')
      endif

C Define constants and conversion factors .....................................
      dt2   = dt**2
      dtby2 = dt/2.0d0

      Ang = 1.d0 / 0.529177d0
      eV  = 1.d0 / 13.60580d0

      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
        fovermp = 0.009579038
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038 * Ang**2 / eV
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
        temp = 2.0d0*kin/(3.0d0*nfree-ct)/8.617d-5
      else
        temp = 2.0d0*kin/(3.0d0*nfree-ct)/8.617d-5/eV
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
        temp = 2.0d0*kin/(3.0d0*nfree-ct)/8.617d-5
      else
        temp = 2.0d0*kin/(3.0d0*nfree-ct)/8.617d-5/eV
      endif

C .....................

      return
      end

 
