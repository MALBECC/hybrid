c Reads the required QM information

      subroutine read_qm(na_u,nesp,isa,iza,xa,atsym, charge, spin)

      use precision
      use sys
      use fdf

      implicit          none
      character         acf*22, acf_default*22
      integer           i, j, iunit, na_u, nesp, iscale 
      integer           isa(na_u),iza(na_u),atnum(nesp)
      double precision  xa(3,na_u), spin
      character*2       atsym(nesp)
      logical           leqi
      integer           nelec,charge
      double precision  charnet,charnet_def

C Format of atomic coordinates
      acf_default = 'Ang'
      acf = fdf_string('AtomicCoordinatesFormat',acf_default)

      if (leqi(acf,'NotScaledCartesianBohr') .or.
     .    leqi(acf,'Bohr') ) then
        iscale = 0
        write(6,'(a,a)')
     .   'read: Atomic-coordinates input format  = ',
     .   '    Cartesian coordinates (in Bohr)'
      else if (leqi(acf,'NotScaledCartesianAng') .or.
     .         leqi(acf,'Ang') ) then
        iscale = 1
        write(6,'(a,a)')
     .   'read: Atomic-coordinates input format  = ',
     .   '    Cartesian coordinates (in Ang)'
      else
        write(6,"(/,'read: ',72(1h*))")
        write(6,"('read:                  INPUT ERROR')")
        write(6,'(2a)') 'read: You must use one of the following',
     .                            'coordinate options:'
        write(6,'(a)') 'read:     - NotScaledCartesianBohr (or Bohr)'
        write(6,'(a)') 'read:     - NotScaledCartesianAng (or Ang) '
        write(6,"('read: ',72(1h*))")
       call die
      endif                                                   

c Read atomic coordinates and species
      if ( fdf_block('AtomicCoordinatesAndAtomicSpecies',iunit) )
     .     then
         do i = 1, na_u 
            read(iunit,*,err=1,end=1) (xa(j,i), j=1,3), isa(i)
         enddo
      else
         call die("read: You must specify the atomic coordinates")
      endif

c change coordinates to Siesta's format
c       Scale atomic coordinates
c       Coord. option = 0 => Do nothing
c       Coord. option = 1 => Multiply by 1./0.529177 (Ang --> Bohr)

      if (iscale .eq. 1) then
        xa(1:3,1:na_u) = xa(1:3,1:na_u) / 0.529177d0
      endif

c Read the atomic labels
      if ( fdf_block('ChemicalSpeciesLabel',iunit) ) then
         do i=1,nesp
         read(iunit,*,err=2,end=2) j,atnum(i),atsym(i)
         enddo
      else 
         call die("read: You must specify the atomic labels")
      endif

c Assignates atomic number (iza)
       iza=0 
       do i=1,na_u
          do j=1,nesp
             if(isa(i).eq.j) iza(i)=atnum(j)
          enddo
             if(iza(i).eq.0) then
             call die("read: There are atomos without atomic number")
             endif
       enddo

c Calculates the total number of electrons 
      charnet_def=0.0d0
      charnet=fdf_double('NetCharge',charnet_def)
      charge=int(charnet)
      nelec=0
      do i=1,na_u
	nelec=nelec+iza(i)
      enddo
      nelec=nelec-charge
        write(6,'(a,i5)')
     .  'read: Total number of electrons        = ',nelec

C Check mumber of electrons with solute spin multiplicity
       spin = fdf_double('TotalSpin',0.0d0)
       if(mod(spin,2.0).eq.0.0.and.mod(nelec,2).ne.0) then
         call die("read:
     . Problem between the number of electrons and the total spin")
       endif
       if(mod(spin,2.0).ne.0.0.and.mod(nelec,2).eq.0) then
       call die("read:
     . Problem between the number of electrons and the total spin")
       endif

       return
 1     stop 'read: problem reading solute coordinates'
 2     stop 'read: problem reading solute chemical specie label' 
        end
c************************************************************************************
C Read the Dynamics Options

      subroutine read_md(idyn, nmove,
     .                   dt, dxmax, ftol,  
     .                   usesavecg, usesavexv , Nick_cent,
     .                   na, 
     .                   nat, nfce, wricoord, mmsteps, tempinit,
     .                   tt, tauber,mn)

C  Modules
      use precision
      use fdf
      use sys
      use scarlett, only: NEB_move_method, NEB_spring_constant,
     .   NEB_Nimages, time_steep, time_steep_max, traj_frec

      implicit none

      integer
     .  idyn,  
     .  nmove, 
     .  na, nfce, mmsteps,
     .  wricoord, nat

      double precision
     .  dt, dxmax, ftol, tempinit,
     .  tt, tauber, mn

      logical
     .    usesavecg, usesavexv, Nick_cent

C  Internal variables .................................................
      character
     .  dyntype*22

      character
     .  dyntype_default*22 

      integer 
     .  nmove_default,
     .  iunit, wrifces,
     .  NEB_Nimages_default,
     .  NEB_move_method_default

      double precision
     .  dt_default, dxmax_default,
     .  ftol_default,  
     .  dx_default,
     .  NEB_spring_constant_default,
     .  time_steep_default
      logical
     .  leqi, qnch, qnch_default

C Temperatura inicial

	tempinit = fdf_physical('MD.InitialTemperature',300.d0,'K')

C Target Temperature JOTA

        tt = fdf_physical('MD.TargetTemperature',300.d0,'K')


C Kind of dynamics
      dyntype_default='Jolie'
        dyntype = fdf_string('MD.TypeOfRun',dyntype_default)

      if (leqi(dyntype,'cg')) then
        idyn = 0
          write(6,'(a,a)') 
     .     'read: Dynamics option                  = ',
     .     '    CG coord. optimization'
          usesavecg  = fdf_boolean('MD.UseSaveCG',.false.)
          write(6,'(a,4x,l1)')
     .     'read: Use continuation files for CG    = ',
     .     usesavecg
      elseif (leqi(dyntype,'qm')) then
        idyn = 2
          write(6,'(a,a)')
     .     'read: Dynamics option                  = ',
     .     '    QM coord. optimization'
          usesavecg  = fdf_boolean('MD.UseSaveCG',.false.)
          write(6,'(a,4x,l1)')
     .     'read: Use continuation files for CG    = ',
     .     usesavecg
      elseif (leqi(dyntype,'fire')) then
        idyn = 3
          write(6,'(a,a)')
     .     'read: Dynamics option                  = ',
     .     '    QM coord. optimization'
          usesavecg  = fdf_boolean('MD.UseSaveCG',.false.)
          write(6,'(a,4x,l1)')
     .     'read: Use continuation files for CG    = ',
     .     usesavecg
      else if (leqi(dyntype,'verlet')) then
        idyn = 4
          write(6,'(a,a)')
     .     'read: Dynamics option                  = ',
     .     '    Velocity Verlet MD run'
      else if (leqi(dyntype,'berendsen')) then
        idyn = 5
          write(6,'(a,a)')
     .     'read: Dynamics option                  = ',
     .     '    Berendsen MD run'
      else if (leqi(dyntype,'nose')) then
        idyn = 6
          write(6,'(a,a)')
     .     'read: Dynamics option                  = ',
     .     '    Nose termostat MD run'
      elseif (leqi(dyntype,'neb')) then
        idyn = 1
          write(6,'(a,a)')
     .     'read: Dynamics option                  = ',
     .     '    BAND coord. optimization'
          usesavecg  = fdf_boolean('MD.UseSaveCG',.false.)
          write(6,'(a,4x,l1)')
     .     'read: Use continuation files for CG    = ',
     .     usesavecg
      else
        write(6,100) 
        write(6,101) 
        write(6,'(a)') 'read:  Wrong Dynamics Option Selected       '
        write(6,'(a)') 'read  You must choose one of the following:'
        write(6,'(a)') 'read:                                       '
        write(6,'(a)') 'read:    - CG - QM - FIRE - NEB                '
        write(6,102)
        call die
      endif 

C Berensen coupling constant
      tauber = dt
      if (idyn .eq. 5) then
      tauber = fdf_physical('MD.TauRelax',dt,'fs')
        write(6,'(a,f10.4,a)')
     .  'read: Berendsen Coupling Constant      = ',
     .   tauber,'  fs'
      endif

C Mass of Nose variable

      if (idyn.eq.6) then
        mn = fdf_physical('MD.NoseMass',1.d2,'eV*fs**2')
        if (mn .eq. 0) then
          write(6,'(/,a)') 
     .  'read: Nose mass                        =   Estimated from Ndf' 
        else
          write(6,'(a,f10.4,a)')
     .  'read: Nose mass                        = ',mn,'  eV/fs**2'
        endif 
      endif


C Read if use saved XV data
      usesavexv = fdf_boolean('MD.UseSaveXV', .false.)
          write(6,'(a,4x,l1)')
     .     'read: Use continuation files for XV    = ',
     .     usesavexv

C Maximum number of steps in CG coordinate optimization
      nmove_default = 0
        nmove = fdf_integer('MD.NumCGsteps',nmove_default)


C NEB optimization method
      NEB_move_method_default = 2
	NEB_move_method = fdf_integer('NEB.OptimizationMethod',
     .  NEB_move_method_default)
C spring constant
      NEB_spring_constant_default=0.2d0
	NEB_spring_constant = fdf_double('NEB.SpringConstant',
     .  NEB_spring_constant_default)
C NEB images
      NEB_Nimages_default = 1
        NEB_Nimages = fdf_integer('NEB.NumImages',NEB_Nimages_default)




C Maximum atomic displacement in one CG step
      dxmax_default = 0.2d0
        dxmax = fdf_physical('MD.MaxCGDispl',dxmax_default,'Bohr')

C Tolerance in the maximum atomic force (def 0.04 eV/Ang)
      ftol_default = 0.00155574d0
        ftol = fdf_physical('MD.MaxForceTol',ftol_default,'Ry/Bohr')

C re-Center system
      Nick_cent = fdf_boolean('CG.Nick_center', .false.)

      if (idyn .eq. 0) then
        write(6,'(a,i5)') 
     .  'read: Maximum number of CG moves       = ',nmove
        write(6,'(a,f10.4,a)') 
     .  'read: Max atomic displ per move        = ',dxmax,'  Bohr'
        write(6,'(a,f10.4,a)') 
     .  'read: Force tolerance                  = ',ftol,'  Ry/Bohr'
      endif
  

C Length of time step for MD
      dt_default = 0.1 
        dt = fdf_physical('MD.LengthTimeStep',dt_default,'fs')
C hay qunificar los timesteeps
      time_steep_default=1d-1
      time_steep = fdf_double('Tstep',
     .  time_steep_default)
      time_steep_max=15.d0*time_steep
C Trajectory frecuency to write coordinates and Energy 
	traj_frec = fdf_integer('MD.TrajFrec',100)

C Quench Option
      qnch_default = .false.
        qnch = fdf_boolean('MD.Quench',qnch_default)


C Sets dt for a CG run
      if(idyn.eq.0) dt = 20./1000.

C Write coordinate variable 
        wricoord = fdf_integer('WriteCoordinates',1)
        write(6,'(a,i5,a)')
     .  'read: Write coordinates each           = ',wricoord, '  steps'

C MMxQM steps
        mmsteps = fdf_integer('MD.MMxQMsteps',1)
        write(6,'(a,i5)')
     .  'read: MM x QM steps                    = ',mmsteps 

C Sets the atoms for whom the forces will be writen
      wrifces = fdf_integer('WriteForcesQMMM',1)
      if(wrifces.eq.0) then
         nfce = 0
      elseif(wrifces.eq.1) then
         nfce = na
      elseif(wrifces.eq.2) then
         nfce = nat
      else
      call die('read: WriteForcesQMMM could only be 0, 1 or 2')
      endif
        write(6,'(a,i5)')
     .  'read: Write forces QM-MM               = ',wrifces

      write(6,102)   

100   format(/,'read: ',71(1h*))
101   format('read:                  INPUT ERROR')
102   format('read: ',71(1h*))
103   format('read: ',i4,2x,3f10.5,i3) 

      return
      end

c**********************************************************************************
C Subroutine to read a crd file from an Amber run
	subroutine readcrd(na_u,nac,masst,linkatom,linkat,numlink,
     .                     rclas,vat,foundcrd,foundvat)

        use precision
        use ionew
        use sys
        use fdf
	use scarlett, only: Ang,eV

	implicit none
	integer i,j,na_u,nac,nat,numlink,linkat(numlink),
     .          natoms,iu
	double precision rclas(1:3,na_u+nac),masst(na_u+nac),
     .                   vat(1:3,na_u+nac),time
	character sname*24,fname*24,paste*24,title*30,ch*4
	logical   foundcrd,linkatom,foundvat
	external paste

        foundcrd = .false.
          sname = fdf_string( 'SystemLabel', 'siesta' )
          fname = paste( sname, '.crd' )
          inquire( file=fname, exist=foundcrd )

          if (foundcrd) then
            call io_assign( iu )
            open( iu, file=fname, status='old' )

	read(iu,err=1,end=1,fmt='(a30)') title
	ch=title(1:4)
	if(ch.eq.'File') then
	foundcrd=.false.
	else
	read(iu,*,err=1,end=1) natoms,time
	
	nat=na_u+nac-numlink
	if(nat.ne.natoms) then
	  write(*,*) 'Wrong number of atoms in ', fname
	  STOP
	end if
	
        if(.not.linkatom) then
	read(iu,err=1,end=1,fmt='(6f12.7)')
     .  ((rclas(i,j),i=1,3),j=1,natoms)
	read(iu,err=1,end=1,fmt='(6f12.7)')
     .  ((vat(i,j),i=1,3),j=1,natoms)
	else !linkatom=.true.
	if(linkat(1).eq.(na_u-numlink+1)) then
	read(iu,err=1,end=1,fmt='(6f12.7)') 
     .  ((rclas(i,j),i=1,3),j=1,na_u-numlink)
        read(iu,err=1,end=1,fmt='(6f12.7)') 
     .  ((rclas(i,j),i=1,3),j=na_u+1,na_u+nac)
        read(iu,err=1,end=1,fmt='(6f12.7)') 
     .  ((vat(i,j),i=1,3),j=1,na_u-numlink)     
        read(iu,err=1,end=1,fmt='(6f12.7)') 
     .  ((vat(i,j),i=1,3),j=na_u+1,na_u+nac)   
	else
	call die("Read CRD: link atoms must be at the input bottom")
	endif
	rclas(1:3,linkat(1:numlink))=0.0
	vat(1:3,linkat(1:numlink))=0.0    
	endif
	call io_close(iu)

        rclas(1:3,1:natoms)=rclas(1:3,1:natoms)*Ang
        vat(1:3,1:natoms)=vat(1:3,1:natoms)*Ang*20.455d0/1000d0

      write(6,'(/a)') 
     .'readcrd: Reading Coordinates and Velocities from CRD file'
	endif
	  endif !foundcrd

c change velocities if using deuterium
	if(foundcrd) then
	do i=1,na_u+nac
	if(masst(i).eq.2.0) vat(1:3,i)=vat(1:3,i)/sqrt(2.)
	enddo
	endif

c check if velocities are zero
	do j=1,3
	  do i=1,na_u+nac
	    if(vat(j,i).ne.0.0) foundvat=.true.
	  enddo
	enddo

	return
 1      stop 'read: problem reading CRD file'
	end

c**************************************************************************
