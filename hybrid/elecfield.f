c This subroutine calculates the interaction bethween the clasiscal
c atoms and a constant external electric field applied

	subroutine elecfield(nac,rclas,fdummy,pc,Eelec)

	use fdf
	implicit none

c	external
	integer nac
	double precision rclas(3,nac), fdummy(3,nac), pc(nac), Eelec 

c	internal
	integer i,iunit
	double precision efield(3),cte 
	character units
	logical ftm,found
	save ftm,found,efield
	data ftm /.true./, found /.false./
	Eelec = 0.0

c leo el campo electrico del .fdf
	if(ftm) then
	if ( fdf_block('ExternalElectricField',iunit)) then
	read(iunit,*) efield(1:3), units
	if(units.ne.'V') stop 'EField units must be V/Ang'
	write(*,*)
	write(*,'(a)') 'Running with an External Electric Field of:'
	write(*,'(3F8.4,3x,A)') efield(1:3),'V/Ang'
	write(*,'(a)') 'Warning: system must be neutral'
	cte = 6.02214D23  * 1.602177D-19 / 4184.0D0 
	efield = efield * cte
	ftm = .false.
	found = .true.
	endif
	endif

c calculo la fza del campo sobre los at clasicos
	if(found) then
	call electricfield( nac, rclas, pc, efield, fdummy, Eelec )
	endif
	end

c**************************************************************************
        subroutine electricfield( nac, r, pc, efield, force, Eelec )
        implicit none
        integer i,nac
        double precision r(3,nac), force(3,nac), efield(3),
     .                   pc(nac), Eelec 

c calculo la energia y fuerzas en el campo y se las sumo a fdummy
        do i=1,nac
          r(1:3,i)=r(1:3,i)*0.529177
          Eelec = Eelec - pc(i)*(r(1,i)*efield(1)+
     .            r(2,i)*efield(2)+r(3,i)*efield(3))
          force(1:3,i) = force(1:3,i) + pc(i) * efield(1:3)
          r(1:3,i)=r(1:3,i)/0.529177
        enddo
      end

c**************************************************************************
