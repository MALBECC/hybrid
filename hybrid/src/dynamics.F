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
     .   natoms,ntcon,istep,iquench,iunit,nfree

      double precision
     .  dt,fa(3,natoms),kin,ma(natoms),
     .  va(3,natoms),xa(3,natoms)

C      external
C     .  memory

C Internal variables ..........................................................
 
      integer
     .  ct,i,ia,Node

      double precision
     .  dot,dt2,dtby2,fovermp,temp
c saque Ang y eV
      double precision, dimension(:,:), allocatable, save ::
     .  accold,vold
C ........................

c	 write(456,*) "istep", istep
c        write(456,*) "iunit", iunit
c        write(456,*) "iquench", iquench
c        write(456,*) "natoms", natoms
c        write(456,*) "fa", fa
c        write(456,*) "dt", dt
c        write(456,*) "ma", ma
c        write(456,*) "ntcon", ntcon
c        write(456,*) "va", va
c        write(456,*) "xa", xa
c        write(456,*) "kin", kin
c        write(456,*) "temp", temp
c        write(456,*) "nfree", nfree




      if (iunit .ne. 1 .and. iunit .ne. 2 .and. iunit .ne. 3) then
          write(6,*) 'verlet2: Wrong iunit option;  must be 1, 2 or 3'
        stop
      endif
      ct = 3 + ntcon
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

c      Ang = 1.d0 / 0.529177d0 !Dividir una distancia en Bohr por "Ang" la pasa a Angstroms
c      eV  = 1.d0 / 13.60580d0 !Dividir una energía en Ry por "eV" la pasa a eV
      
      fovermp = 0.009648534d0
      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
c        fovermp = 0.009579038d0                      !Revise la cuenta y a mi me dio 0.009648534 Jota
c         fovermp = 0.009648534d0
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
C	write(952,*) ia, ma(ia)   ! Para ver en que unidades escribe ma jota
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
C creemos que esto esta bien jota <--------------------------------------------------
C .....................



c  	write(456,*) "end"
c        write(456,*) "istep", istep
c        write(456,*) "iunit", iunit
c        write(456,*) "iquench", iquench
c        write(456,*) "natoms", natoms
c        write(456,*) "fa", fa
c        write(456,*) "dt", dt
c        write(456,*) "ma", ma
c        write(456,*) "ntcon", ntcon
c        write(456,*) "va", va
c        write(456,*) "xa", xa
c        write(456,*) "kin", kin
c        write(456,*) "temp", temp
c        write(456,*) "nfree", nfree



      return
      end

