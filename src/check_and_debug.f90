	subroutine print_conectivity()
	use scarlett, only: nac, atname, aaname, ng1, angexat, atange, angmxat, atangm, &
	dihexat,atdihe, dihmxat, atdihm


	implicit none
	integer :: i, j, jmax, k


	write(*,*) "------------------------------------------------------"
	do i=1,nac
	   write(*,*) "atom:", i, " name: ", atname(i), " residue: ", aaname(i)
!Bonda
	   write(*,*) "Bonds: "
	   jmax=0
	   do j=1,6
	     if (ng1(i,j) .ne. 0) jmax=j
	   end do
	   do j=1,jmax
	     write(*,*) "   atom:", ng1(i,j), " name: ", atname(ng1(i,j)), &
	     "residue: ", aaname(ng1(i,j))
	   end do
	   write(*,*) "   - - - - - - - - - - - - - - -"
!angles
	   write(*,*) "Angles: "
	   do j=1,angexat(i)
	     write(*,*) "   atom 1:",atange(i,j,1), " name: ", atname(atange(i,j,1)), &
	     "residue: ", aaname(atange(i,j,1))
	     write(*,*) "   atom 2:",atange(i,j,2), " name: ", atname(atange(i,j,2)), &
	     "residue: ", aaname(atange(i,j,2))
	   write(*,*) "   -   -   -   -   -   -   -   -"
	   end do
	   do j=1,angmxat(i)
	     write(*,*) "   atom 1:",atangm(i,j,1), " name: ", atname(atangm(i,j,1)), &
	     "residue: ", aaname(atangm(i,j,1))
	      write(*,*) "   atom 2:",atangm(i,j,2), " name: ", atname(atangm(i,j,2)), &
	     "residue: ", aaname(atangm(i,j,2))
	     write(*,*) "   -   -   -   -   -   -   -   -"
	   end do
	   write(*,*) "   - - - - - - - - - - - - - - -"
!dihedrals
	   write(*,*) "Dihedrals: "
	   do j=1, dihexat(i)
	     do k=1,3
		write(*,*) "   atom", k, ": ",atdihe(i,j,k), " name: ", atname(atdihe(i,j,k)), &
		"residue: ", aaname(atdihe(i,j,k))
	     end do
	     write(*,*) "   -   -   -   -   -   -   -   -"
	   end do
	   do j=1,dihmxat(i)
	     do k=1,3
		write(*,*) "   atom", k, ": ",atdihm(i,j,k), " name: ", atname(atdihm(i,j,k)), &
		"residue: ", aaname(atdihm(i,j,k))
	     end do
	     write(*,*) "   -   -   -   -   -   -   -   -"
	   end do
	   write(*,*) "   - - - - - - - - - - - - - - -"
	end do

	end subroutine print_conectivity


	subroutine solv_ene_fce_in(natot,na_u,nac,ng1,rclas,Em,Rm,pc, &
	Etot_amber,fcetot_amber,attype, &
	nbond,nangle,ndihe,nimp,multidihe, multiimp,kbond,bondeq, &
	kangle,angleeq,kdihe,diheeq,kimp,impeq,perdihe,perimp, &
	bondtype,angletype,dihetype,imptype, &
	bondxat,angexat,angmxat,dihexat,dihmxat,impxat,atange, &
	atangm,atdihe,atdihm,atimp, &
	ng1type,angetype,angmtype,evaldihe,evaldihm, &
	dihety,dihmty,impty,nonbonded,scale,nonbondedxat,scalexat, &
	evaldihelog,evaldihmlog,paso,nparm, &
	actualiz,listcut, &
	noat,noaa,sfc,timestep, &
	water,masst,radblommbond)     

      use scarlett, only: eV, Ang, kcal
        implicit none

       integer   i,j,k,l,na_u,natot,nac,ng1(nac,6),paso 
       double precision  pcA(nac),rclas(3,natot),ramber(3,nac), &
	EmA(nac),RmA(nac),pc(1:nac),Em(natot),Rm(natot), &
	Etot_amber,Ebond_amber,Eangle_amber,Edihe_amber,Eimp_amber, &
	Elj_amber,Eelec_amber,Elj_amber14,Eelec_amber14 
       double precision fcetot_amber(3,nac), &
	 fcebond_amber(3,nac),fceangle_amber(3,nac), &
	 fcedihe_amber(3,nac),fceimp_amber(3,nac), &
	 fcelj_amber(3,nac),fceelec_amber(3,nac), &
	 fcelj_amber14(3,nac),fceelec_amber14(3,nac), &
	 fcetotaxes_amber(3)
       character  attype(nac)*4,noat(nac)*4,noaa(nac)*4
       double precision listcut,sfc,timestep,masst(natot),ewat,fwat(3,nac)
       logical water
!c      vbles de los params de union
       integer   nbond,nangle,ndihe,nimp,nparm,multidihe(nparm), multiimp(nparm)
       double precision   kbond(nparm),bondeq(nparm),kangle(nparm), &
	 angleeq(nparm),kdihe(nparm),diheeq(nparm),kimp(nparm), &
	 impeq(nparm),perdihe(nparm),perdihe2(nparm),perimp(nparm)
	character   bondtype(nparm)*5,angletype(nparm)*8, &
	 dihetype(nparm)*11,imptype(nparm)*11
!c      vbles de bond, angle, dihe e imp
       integer   bondxat(nac),angexat(nac),angmxat(nac), &
	dihexat(nac),dihmxat(nac),impxat(nac)
       integer   atange(nac,25,2),atangm(nac,25,2), &
	 atdihe(nac,100,3),atdihm(nac,100,3),atimp(nac,25,4)
!c	vbles q faltaban
       integer  ng1type(nac,6),angetype(nac,25),angmtype(nac,25), &
	        evaldihe(nac,100,5),evaldihm(nac,100,5), &
	        dihety(nac,100),dihmty(nac,100),impty(nac,25), &
	        nonbonded(nac,100),scale(nac,100), &
	        nonbondedxat(nac),scalexat(nac) 
	logical  evaldihelog(nac,100),evaldihmlog(nac,100),actualiz
!c parche para omitir interaccion entre extremos terminales
	double precision, intent(in) :: radblommbond

!c inicializa las energias y fuerzas
      Etot_amber=0.d0
      Ebond_amber=0.d0
      Eangle_amber=0.d0
      Edihe_amber=0.d0
      Eimp_amber=0.d0
      Elj_amber=0.d0
      Eelec_amber=0.d0
      Elj_amber14=0.d0
      Eelec_amber14=0.d0
      fcetot_amber=0.d0
      fcebond_amber=0.d0
      fceangle_amber=0.d0
      fcedihe_amber=0.d0
      fceimp_amber=0.d0
      fcelj_amber=0.d0
      fceelec_amber=0.d0
      fcetotaxes_amber=0.d0   
      ewat=0.d0
      fwat=0.d0

!c cambia variables
      k=1
      do j=1,nac
      pcA(k)=pc(j)
      k=k+1
      enddo
      k=1
      do i=na_u+1,natot
      ramber(1:3,k)=rclas(1:3,i)
      k=k+1
      enddo
      k=1
      do i=na_u+1,natot
      EmA(k)=Em(i)
      RmA(k)=Rm(i)
      k=k+1
      enddo
 
!c  pasa a las unidades de Amber
      do i=1,nac
      RmA(i) = RmA(i)*(2.d0**(1.d0/6.d0))/(2.d0*Ang) ! revisar jota
      EmA(i) = EmA(i)*kcal/eV !627.5108d0 
      ramber(1:3,i)=ramber(1:3,i)/Ang   !ramber queda en angstroms
      enddo

!c  llama a subrutina q calcula la energia y fuerzas de bonds

        call amber_bonds_in(nac,ng1,ramber,Ebond_amber,attype, &
	     nbond,kbond,bondeq,bondtype,bondxat,fcebond_amber, &
	     ng1type,paso,nparm,radblommbond)

!c  llama a subrutina q calcula la energia y fuerzas de angles

        call amber_angles_in(nac,ng1,ramber,Eangle_amber,attype, &
	     nangle,kangle,angleeq,angletype,angexat,angmxat,atange, &
	     atangm,fceangle_amber,angetype,angmtype,paso,nparm)

!c  llama a subrutina q calcula la energia y fuerzas de dihedros     
        perdihe2=perdihe  
        call amber_dihes_in(nac,ng1,ramber,Edihe_amber, &
	     attype,ndihe,kdihe,diheeq,perdihe2,multidihe, &
	     dihetype,dihexat,dihmxat,atdihe,atdihm, &
	     fcedihe_amber,evaldihelog,evaldihe,dihety, &
	     evaldihmlog,evaldihm,dihmty,paso,nparm) 

!c  llama a subrutina q calcula la energia y fuerzas de impropers

!        call amber_improper_in(nac,ng1,ramber,Eimp_amber,attype, &
!	     nimp,kimp,impeq,perimp,multiimp,imptype,impxat,atimp, &
!	     fceimp_amber,impty,paso,nparm)

      return
      end
!c****************************************************************
!c subrutina q calcula la energia y fuerzas de bonds

        subroutine amber_bonds_in(nac,ng1,ramber,Ebond_amber, &
	           attype,nbond,kbond,bondeq,bondtype,bondxat, &
	           fcebond_amber,ng1type,paso,nparm,radblommbond)

       implicit none      
!c      vbles grales 
       integer   nac,ng1(nac,6),i,j,k,l,m,n,paso
       double precision   ramber(3,nac),Ebond_amber, fcebond_amber(3,nac)
       character   attype(nac)*4
!c      vbles de los params de union
       integer   nbond,nparm
       double precision   kbond(nparm),bondeq(nparm)
       character   bondtype(nparm)*5
!c      vbles de bond, angle, dihe e imp
       integer   bondxat(nac)
!c      vbles de asignacion
       character*4 ty1,ty2
       character*5 tybond
       integer ng1type(nac,6)
       double precision   rij,dist,dx1,dx2,dy1,dy2,dz1,dz2, ftotbond(3)
      logical           first                                                                  
      save              first                                                                  
      data              first /.true./    
      logical           ST
       double precision :: radblommbond

	ST=.false.

!c asignacion de tipos de union
      if(first) then
       do i=1,nac
!c barre at clasicos
        do j=1,bondxat(i)
!c barre bonds del atomo i
         do k=1,nbond
!c barre todos los bonds leidos en el amber.parm
          tybond=bondtype(k)
          ty1=tybond(1:2)
          ty2=tybond(4:5)

          if(attype(i).eq.ty1.and.attype(ng1(i,j)).eq.ty2) then
            ng1type(i,j)=k
          elseif(attype(i).eq.ty2.and.attype(ng1(i,j)).eq.ty1) then
            ng1type(i,j)=k
          endif

         enddo
        enddo
       enddo
      first=.false.
      endif !asignacion
      end
!c******************************************************
!c  subrutina q calcula la energia y fuerzas de angles
 
        subroutine  amber_angles_in(nac,ng1,ramber, &
	 Eangle_amber,attype,nangle,kangle,angleeq,angletype, &
	 angexat,angmxat,atange,atangm,fceangle_amber, &
	 angetype,angmtype,paso,nparm)

        implicit none
!c      vbles grales
       integer   nac,ng1(nac,6),i,j,k,l,m,n,paso
       double precision   ramber(3,nac),Eangle_amber, fceangle_amber(3,nac)
       character   attype(nac)*4
!c      vbles de los params de union
       integer   nangle,nparm
       double precision kangle(nparm),angleeq(nparm)
       character angletype(nparm)*8
!c      vbles de bond, angle, dihe e imp
       integer   angexat(nac),angmxat(nac)
       integer   atange(nac,25,2),atangm(nac,25,2)
!c      vbles de asignacion
       character*4 ty1,ty2,ty3
       character*8 tyangle
       integer angetype(nac,25),angmtype(nac,25)
       double precision  angulo,angle,pi,dx,dy,dz,scal,r12,r32, &
	                 scalar,ftotangle(3),dr12r32,dscalar,dist, &
	                 fesq(3,nac),fmedio(3,nac)
      logical           first                                                                  
      save              first                                                                  
      data              first /.true./    
      
       pi=DACOS(-1.d0)
       fesq=0.d0
       fmedio=0.d0
       ftotangle=0.d0

!c asignacion de los tipos de angulos 
      if(first) then
       do i=1,nac
        do j=1,angexat(i)
         do k=1,nangle 
          tyangle=angletype(k)
          ty1=tyangle(1:2)
          ty2=tyangle(4:5)
          ty3=tyangle(7:8)
          if(attype(i).eq.ty1.and.attype(atange(i,j,1)).eq.ty2.and. &
	       attype(atange(i,j,2)).eq.ty3) then
          angetype(i,j)=k
          elseif(attype(i).eq.ty3.and.attype(atange(i,j,1)).eq.ty2.and. & 
	       attype(atange(i,j,2)).eq.ty1) then
          angetype(i,j)=k
          endif
         enddo
        enddo
       enddo

       do i=1,nac
        do j=1,angmxat(i)
         do k=1,nangle
          tyangle=angletype(k)
          ty1=tyangle(1:2)
          ty2=tyangle(4:5)
          ty3=tyangle(7:8)
          if(attype(i).eq.ty2.and.attype(atangm(i,j,1)).eq.ty1.and. &
	     attype(atangm(i,j,2)).eq.ty3) then
          angmtype(i,j)=k
          elseif(attype(i).eq.ty2.and.attype(atangm(i,j,1)).eq.ty3.and. &
	     attype(atangm(i,j,2)).eq.ty1) then
          angmtype(i,j)=k
          endif
         enddo
        enddo
       enddo
      first=.false.
      endif ! asignacion

        
	end
!c****************************************************************
!c  subrutina q calcula la energia y fuerzas de dihedros
 
      subroutine  amber_dihes_in(nac,ng1,ramber,Edihe_amber, &
	          attype,ndihe,kdihe,diheeq,perdihe,multidihe, &
	          dihetype,dihexat,dihmxat,atdihe,atdihm, &
	          fcedihe_amber,evaldihelog,evaldihe,dihety, &
	          evaldihmlog,evaldihm,dihmty,paso,nparm)
	
        implicit none
 
!c      vbles grales
       integer   nac,ng1(nac,6),i,j,k,l,m,n,paso
       double precision   ramber(3,nac),Edihe_amber, fcedihe_amber(3,nac)
       character   attype(nac)*4
!c      vbles de los params de union
       integer   ndihe,nparm,multidihe(nparm)
       double precision kdihe(nparm),diheeq(nparm),perdihe(nparm)
       character dihetype(nparm)*11
!c      vbles de bond, angle, dihe e imp
       integer   dihexat(nac),dihmxat(nac)
       integer   atdihe(nac,100,3),atdihm(nac,100,3)
!c      vbles de asignacion
       character*4 ty1,ty2,ty3,ty4
       character*11 tydihe
       integer dihety(nac,100),dihmty(nac,100),evaldihe(nac,100,5), evaldihm(nac,100,5)
       double precision  pi,dihedro,dihedral,E1,dist,dx,dy,dz,ftotdihe(3), &
	                 fesq(3,nac),fmedio(3,nac), fce(12)
       logical evaldihelog(nac,100),evaldihmlog(nac,100)
       logical           first                                                                  
       save              first                                                                  
       data              first /.true./    

       pi=DACOS(-1.d0)
       ftotdihe=0.d0
       fesq=0.d0
       fmedio=0.d0

!c asignacion de los tipos de dihedros
	if(first) then
	  evaldihelog=.false.
	  evaldihmlog=.false.
	
! Para los dihedros del extremo
	  do i=1,nac            
	    do j=1,dihexat(i)    ! Barre todos los átomos, y en cada átomo, todos sus dihedros en extremo (?
	      dihety(i,j)=1        ! Dihety es una matriz de dimensión nac x 100 (100 porque es
                                   ! un número seguramente más grande que dihexat máximo). Tiene,
                                   ! para cada átomo, cuántos dihedros tiene
	      m=0
!c Comienza asignación
	      do k=1,ndihe       ! Barre TODOS los dihedros: ndihe es el número total de dihedros 
                             ! definidos en el amber.parm 
 	        tydihe=dihetype(k) ! Cada uno de los dihedros tiene asignado un tipo
	        ty1=tydihe(1:2)    !
	        ty2=tydihe(4:5)    ! Acá se fija entre qué tipo de átomos se da el dihedro
	        ty3=tydihe(7:8)    !
	        ty4=tydihe(10:11)  ! 

!c Ahora asigna tipo de átomo a cada átomo del dihedro
	        if(ty1.eq.'X ') then   ! Caso en que el dihedro empieza con X         
	          if(attype(atdihe(i,j,1)).eq.ty2.and. attype(atdihe(i,j,2)).eq.ty3)  then
	            dihety(i,j)=k
	          elseif(attype(atdihe(i,j,1)).eq.ty3.and. attype(atdihe(i,j,2)).eq.ty2)  then
	            dihety(i,j)=k
	          endif
	        elseif(ty1.ne.'X ') then ! Caso en que el dihedro no empieza con X
	          if(attype(i).eq.ty1.and.attype(atdihe(i,j,1)).eq. &
	            ty2.and.attype(atdihe(i,j,2)).eq.ty3.and. &
	            attype(atdihe(i,j,3)).eq.ty4) then 
	            dihety(i,j)=k            
	            if(perdihe(k).lt.0) then
	              evaldihelog(i,j)=.true.
	              m=m+1
 	              evaldihe(i,j,m)=k
	              evaldihe(i,j,m+1)=k+1
	            endif
	          elseif(attype(i).eq.ty4.and.attype(atdihe(i,j,1)).eq. &
	            ty3.and.attype(atdihe(i,j,2)).eq.ty2.and. &
	            attype(atdihe(i,j,3)).eq.ty1) then
	            dihety(i,j)=k
	            if(perdihe(k).lt.0) then
	              evaldihelog(i,j)=.true.
	              m=m+1
	              evaldihe(i,j,m)=k
	              evaldihe(i,j,m+1)=k+1    
	            endif
	          endif
	        endif  
	      enddo
	    enddo
	  enddo
! Para los dihedros del medio 
	  do i=1,nac
	    do j=1,dihmxat(i)
	      dihmty(i,j)=1
	      m=0
	      do k=1,ndihe
	        tydihe=dihetype(k)
	        ty1=tydihe(1:2)
	        ty2=tydihe(4:5)
	        ty3=tydihe(7:8)
	        ty4=tydihe(10:11)

	        if(ty1.eq.'X ') then
	          if(attype(i).eq.ty2.and. attype(atdihm(i,j,2)).eq.ty3)  then
	            dihmty(i,j)=k
	          elseif(attype(i).eq.ty3.and. attype(atdihm(i,j,2)).eq.ty2)  then
	            dihmty(i,j)=k
	          endif
	        elseif(ty1.ne.'X ') then

	          if(attype(atdihm(i,j,1)).eq.ty1.and.attype(i).eq. &
	            ty2.and.attype(atdihm(i,j,2)).eq.ty3.and. &
	            attype(atdihm(i,j,3)).eq.ty4) then
	            dihmty(i,j)=k
		    if(perdihe(k).lt.0.d0) then
	              evaldihmlog(i,j)=.true.
	              m=m+1
	              evaldihm(i,j,m)=k
	              evaldihm(i,j,m+1)=k+1  
	            endif
	          elseif(attype(atdihm(i,j,1)).eq.ty4.and.attype(i).eq. &
	            ty3.and.attype(atdihm(i,j,2)).eq.ty2.and. &
	            attype(atdihm(i,j,3)).eq.ty1) then
	            dihmty(i,j)=k
	            if(perdihe(k).lt.0.d0) then
	              evaldihmlog(i,j)=.true.
	              m=m+1
	              evaldihm(i,j,m)=k
	              evaldihm(i,j,m+1)=k+1
	            endif
	          endif
	        endif
	      enddo
	    enddo
	  enddo
	first=.false.
	endif !asignacion
       end
!c******************************************************************
!c  subrutina q calcula la energia y fuerzas de impropers 


