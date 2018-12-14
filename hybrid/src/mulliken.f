C subroutine to read from Siesta standard output the
C mulliken charges and population analysis

      subroutine mulliken(na_u,isa,iza,nesp,
     .atsym,siestaslabel,step,nconstr)

	use ionew
	use fdf
	implicit none
	integer na_u,nspin,nesp,isa(na_u),iza(na_u),
     .          atxesp(nesp),na(100,100),nconstr,step(nconstr)
	double precision  qa(na_u,2),q(100,100),qna(na_u)
	character siestaslabel*20,slabel*30,atsym(nesp)*2,simbolo*2 
	character option*16,exp,paste*30   
	logical search
	external paste  
	integer i,j,k,l,m,n,p,iu,ico,ist,keystep(nconstr),count(nconstr)
	double precision spin
	logical sppol,spinlog,fix

	qa=0.0
	q=0.0
	qna=0.0

C Everything related to spin
        nspin = 1
        spinlog = .false.
        sppol  = fdf_boolean('SpinPolarized',.false.)
        if(sppol) nspin = 2
        if(sppol) spinlog = .true.
        fix = fdf_boolean('FixSpin',.false.)
        if (fix) spin = fdf_double('TotalSpin',0.0d0)
        if(spin.eq.0.0) spinlog=.false.

C Assigns atoms in each specie
        atxesp=0
        do i=1,nesp
          do j=1,na_u
        if(isa(j).eq.i) atxesp(i)=atxesp(i)+1
          enddo
        enddo

c calculate keystep
        keystep=0
        keystep(1)=step(1)+1
        if(nconstr.gt.1) then
        do i=2,nconstr
        keystep(i)=keystep(i-1)+step(i)
        enddo
        endif    

C Open/Close Siesta output file for ionew not to crush 
        slabel = paste( siestaslabel, '.out' )
	write(*,*) 'siestaslabel:  ',siestaslabel
        call io_assign( iu )
        open( iu, file=slabel, status='unknown' )
 1      continue
        call io_close( iu )

c re-open file 
        call io_assign( iu )
        open( iu, file=slabel, status='unknown' )

c Reads from file
	search=.true.
	do while (search)   
	read (iu,'(A16)',end=1,err=10) option
	if(option.eq.'Stopping Program') goto 2
	enddo

 2      continue
        call io_close( iu )

c loop over nconstr
	do ico=1,nconstr
	count(ico)=0

c re-open file
        call io_assign( iu )
        open( iu, file=slabel, status='unknown' )

c re-reads from file
        search=.true.
        do while (search)
 3      continue
	read (iu,'(A16)',end=5,err=10) option
	if (option.eq.'mulliken: Atomic') then
	count(ico)=count(ico)+1

	if (nspin.eq.2) then
	read(iu,*)
	read(iu,*)
	read(iu,*)
	do k=1,nspin
	do l=1,nesp
	read(iu,*,err=20,end=20) exp,simbolo
        if(simbolo.ne.'H'.and.simbolo.ne.'Na'
     .  .and.simbolo.ne.'Mg') read(iu,*)
	read(iu,*)
	read(iu,*)
	do j=1,nesp
	if(simbolo.eq.atsym(j)) then
	do m=1,atxesp(j)
	read(iu,*,err=30,end=30) na(j,m),q(j,m) 
        if(simbolo.ne.'H'.and.simbolo.ne.'Na'
     .  .and.simbolo.ne.'Mg') read(iu,*)
	enddo
	read(iu,*)
	endif
	enddo
	enddo
	read(iu,*)
	read(iu,*)
	read(iu,*)
	read(iu,*)
	do j=1,nesp
	do m=1,atxesp(j)
	qa(na(j,m),k)=q(j,m)    
	enddo
	enddo
	
	enddo

	elseif(nspin.eq.1) then
	read(iu,*)
	do l=1,nesp
	read(iu,*,err=20,end=20) exp,simbolo
        if(simbolo.ne.'H'.and.simbolo.ne.'Na'
     .  .and.simbolo.ne.'Mg') read(iu,*)
	read(iu,*)
	read(iu,*)
	do j=1,nesp
	if(simbolo.eq.atsym(j)) then
	do m=1,atxesp(j)
	read(iu,*,err=30,end=30) na(j,m),q(j,m)
        if(simbolo.ne.'H'.and.simbolo.ne.'Na'
     .  .and.simbolo.ne.'Mg') read(iu,*)
	enddo
	read(iu,*)
	endif
	enddo
	enddo
	do j=1,nesp
	do m=1,atxesp(j)
	qa(na(j,m),1)=q(j,m)
	enddo
	enddo
	endif !nspin

	if(count(ico).eq.keystep(ico)) then !write in file
C calculates Mulliken and spin populations
        do i=1,na_u
        if(iza(i).le.2) 
     .    qna(i)=iza(i)-qa(i,1)-qa(i,2)
        if(iza(i).gt.2.and.iza(i).le.10)
     .    qna(i)=iza(i)-2.0-qa(i,1)-qa(i,2)
        if(iza(i).gt.10.and.iza(i).le.18)
     .    qna(i)=iza(i)-10.0-qa(i,1)-qa(i,2)
        if(iza(i).gt.18.and.iza(i).le.36)
     .    qna(i)=iza(i)-18.0-qa(i,1)-qa(i,2)
        if(iza(i).gt.36.and.iza(i).le.54)
     .    qna(i)=iza(i)-36.0-qa(i,1)-qa(i,2)
        if(iza(i).gt.54.and.iza(i).le.71)
     .    qna(i)=iza(i)-54.0-qa(i,1)-qa(i,2)
        if(iza(i).gt.71.and.iza(i).le.86)
     .    qna(i)=iza(i)-68.0-qa(i,1)-qa(i,2)
        enddo

       if(nconstr.gt.1) then
       write(6,'(/a,i4)') 'Analysis of constrained step:  ',ico
       endif
       write(6,'(/a)') 'Mulliken Population Analisys:'
       write(6,"(i4,2x,A2,2x,F10.6)") 
     . (i,atsym(isa(i)),qna(i),i=1,na_u)   
        if(spinlog) then
       write(6,'(/a)') 'Spin Population Analisys:'
       write(6,"(i4,2x,A2,2x,F10.6)") 
     . (i,atsym(isa(i)),(qa(i,1)-qa(i,2)),i=1,na_u)
        endif

	endif !write
	else
	goto 3
	endif !option
	enddo !search
 5      continue
	call io_close (iu)
	enddo !nconstr

      return
 10   write(*,*)'mulliken: problem reading siesta output'
      stop
 20   write(*,*)'mulliken: problem reading atomic specie:  ',
     .simbolo
      stop
 30   write(*,*)'mulliken: problem reading mulliken charge:  ',
     .simbolo,m
      stop
      end
 
