	subroutine P_fluctuation(niter,good, xnano)
	use garcha_mod, only : RMM, M, P_oscilation_analisis, P_hist
	implicit none
	integer, intent(in) :: niter
	double precision, intent(inout) :: good
	double precision, intent(in), dimension(M,M) :: xnano
	integer :: jj,kk,ll
	integer :: M2, MM
	integer :: Rposition
	double precision :: del
	M2=2*M

	if (P_oscilation_analisis) then
	
	  MM=M*(M+1)/2
	  if (niter .lt. 1) STOP "niter <0 in P_fluctuation"

	  if (niter .eq. 1) then
	    allocate (P_hist(3,MM))
	    P_hist=0.d0
	  end if
	  
	  if (niter .le. 3) then
	    do jj=1,M
	    do kk=jj,M
	      Rposition=kk+(M2-jj)*(jj-1)/2
	      P_hist(niter,Rposition)=xnano(jj,kk)
	    end do
	    end do
	  else
	    good = 0.0d0
	    do ll=1,3
	      do jj=1,M
	      do kk=jj,M
	        Rposition=kk+(M2-jj)*(jj-1)/2
	        del=xnano(jj,kk)-P_hist(ll,Rposition)
	        del=del*sqrt(2.d0)
	        good=good+del**2
	      enddo
	      enddo
	      good=sqrt(good)/float(M)
	      write(*,*) "oscilation", ll, good
	    end do
	    P_hist(1,1:MM)=P_hist(2,1:MM)
	    P_hist(2,1:MM)=P_hist(3,1:MM)
	    do jj=1,M
	    do kk=jj,M
	      Rposition=kk+(M2-jj)*(jj-1)/2
	      P_hist(3,Rposition)=xnano(jj,kk)
	    end do
	    end do
	  end if

	end if

	good = 0.0d0
	do jj=1,M
	  do kk=jj,M
	    Rposition=kk+(M2-jj)*(jj-1)/2
	    del=xnano(jj,kk)-RMM(Rposition)
	    del=del*sqrt(2.d0)
	    good=good+del**2
	    RMM(kk+(M2-jj)*(jj-1)/2)=xnano(jj,kk)
	  enddo
	enddo
	good=sqrt(good)/float(M)
	end subroutine P_fluctuation

	subroutine P_linearsearch_init(En)
	use garcha_mod, only : M, Md, OPEN, rho_lambda1, rho_lambda0, rho_lambda1_alpha, rho_lambda0_alpha, rho_lambda1_betha, rho_lambda0_betha, RMM, rhoalpha, rhobeta, Elast
	use faint_cpu77, only: int3lu
	implicit none
	integer :: MM, M1, M3, M5, M7, M9, M11, MMd
	integer :: kk
	double precision :: E1, E2, Ex
	double precision, intent(in) :: En
	write(*,*) "Doing linear search in P"
	MM=M*(M+1)/2
	MMd=Md*(Md+1)/2
	allocate(rho_lambda1(MM), rho_lambda0(MM))
	rho_lambda0(1:MM)=RMM(1:MM)
	if(OPEN) allocate(rho_lambda1_alpha(MM), rho_lambda0_alpha(MM),rho_lambda1_betha(MM), rho_lambda0_betha(MM))
	if(OPEN) rho_lambda0_alpha(1:MM)=rhoalpha(1:MM)
	if(OPEN) rho_lambda0_betha(1:MM)=rhobeta(1:MM)


      M1=1 ! first P
      M3=M1+MM ! now Pnew
      M5=M3+MM! now S, F also uses the same position after S was used
      M7=M5+MM! now G
      M9=M7+MMd ! now Gm
      M11=M9+MMd! now H
	E1=0.d0
	do kk=1,MM
            E1=E1+RMM(kk)*RMM(M11+kk-1)
        enddo

          call int3lu(E2) ! Computes Coulomb part of Fock, and energy on E2
          call g2g_solve_groups(0,Ex,0) ! XC integration / Fock elements
          Elast=E1+E2+En+Ex

	return
	end subroutine P_linearsearch_init

	subroutine P_linear_calc(niter, En, good)
!,E2, Ex,good)
	use garcha_mod, only : M, Md, OPEN, rho_lambda1, rho_lambda0, rho_lambda1_alpha, rho_lambda0_alpha, rho_lambda1_betha, rho_lambda0_betha, RMM, Pstepsize, rhoalpha, rhobeta, Elast
	use faint_cpu77, only: int3lu
	implicit none
	integer, intent(in) :: niter
	integer :: ilambda, imin
	double precision :: dlambda
	double precision :: E1, Emin
	double precision :: Blambda, del
	double precision, intent(in) :: En
	double precision, dimension(0:10) :: E_lambda
	double precision, intent(inout) :: good
	double precision :: E2, Ex
!	double precision, intent(inout) :: E2, Ex
	integer :: MM, kk, jj,  MMd
	integer :: M1, M3, M5, M7, M9, M11
	integer :: Rposition, M2
	double precision :: traza
	MM=M*(M+1)/2
	MMd=Md*(Md+1)/2
	M2=2*M

      M1=1 ! first P
      M3=M1+MM ! now Pnew
      M5=M3+MM! now S, F also uses the same position after S was used
      M7=M5+MM! now G
      M9=M7+MMd ! now Gm
      M11=M9+MMd! now H
	if (niter .eq. 1) then
	  Pstepsize=1.d0
	end if

	rho_lambda1(1:MM)=RMM(1:MM)
	
	if(OPEN) rho_lambda1_alpha(1:MM)=rhoalpha(1:MM)
	if(OPEN) rho_lambda1_betha(1:MM)=rhobeta(1:MM)
	
	E_lambda=0.d0

!skip LS
	if (.true.) then
          E1=0.D0
          do kk=1,MM
            E1=E1+RMM(kk)*RMM(M11+kk-1)
          enddo

          call int3lu(E2) ! Computes Coulomb part of Fock, and energy on E2
          call g2g_solve_groups(0,Ex,0) ! XC integration / Fock elements
  	  E_lambda(10)=E1+E2+En+Ex
	end if
!

	if (E_lambda(10) .lt. Elast) then
	  Blambda=Pstepsize
	  write(*,*) "skiping interpolation"
	else

	do ilambda=0, 9
	  dlambda=Pstepsize*dble(ilambda)/10.d0
	  RMM(1:MM)=rho_lambda0(1:MM)*(1.d0-dlambda)+rho_lambda1(1:MM)*dlambda
	  

          if(OPEN) rhoalpha(1:MM)=rho_lambda0_alpha(1:MM)*(1.d0-dlambda)+rho_lambda1_alpha(1:MM)*dlambda
          if(OPEN) rhobeta(1:MM)=rho_lambda0_betha(1:MM)*(1.d0-dlambda)+rho_lambda1_betha(1:MM)*dlambda


	  E1=0.D0
	  do kk=1,MM
	    E1=E1+RMM(kk)*RMM(M11+kk-1)
	  enddo

	  call int3lu(E2) ! Computes Coulomb part of Fock, and energy on E2
	  call g2g_solve_groups(0,Ex,0) ! XC integration / Fock elements
	  write(*,*) ilambda,E1,E2,En,Ex

	  E_lambda(ilambda)=E1+E2+En+Ex
	end do

!now find best linear combination of P
	Emin=999999999999
	imin=11
	do ilambda=0, 10
	  write(*,*) "rholinear E", ilambda, E_lambda(ilambda)
	  if (E_lambda(ilambda) .lt. Emin) then
	    Emin=E_lambda(ilambda)
	    imin=ilambda
	  end if
	end do
	if (imin .eq. 11) STOP "imin 11"
	Blambda=Pstepsize*dble(imin)/10.d0
!aca va el linear search
	write(*,*) "best lambda prediction: normal", Blambda

        call line_search(11,E_lambda, 1d0, Blambda )
	write(*,*) "best lambda prediction: 0", Blambda
	if (Blambda .ge. 1.d0) Blambda=Blambda-1.0d0
	Blambda=Blambda*Pstepsize/10.d0
	write(*,*) "best lambda prediction: 1", Blambda
	end if !Nuevo

	RMM(1:MM)=rho_lambda0(1:MM)*(1.d0-Blambda)+rho_lambda1(1:MM)*Blambda

        if(OPEN) rhoalpha(1:MM)=rho_lambda0_alpha(1:MM)*(1.d0-Blambda)+rho_lambda1_alpha(1:MM)*Blambda
        if(OPEN) rhobeta(1:MM)=rho_lambda0_betha(1:MM)*(1.d0-Blambda)+rho_lambda1_betha(1:MM)*Blambda


          E1=0.D0
          do kk=1,MM
            E1=E1+RMM(kk)*RMM(M11+kk-1)
          enddo


          call int3lu(E2) ! Computes Coulomb part of Fock, and energy on E2
          call g2g_solve_groups(0,Ex,0) ! XC integration / Fock elements

	  Elast=E1+E2+En+Ex

        good = 0.0d0
        do jj=1,M
        do kk=jj,M
	  Rposition=kk+(M2-jj)*(jj-1)/2
          del=rho_lambda0(Rposition)-(RMM(Rposition))
          del=del*sqrt(2.D0)
          good=good+del**2
        enddo
        enddo
        good=sqrt(good)/float(M)

	rho_lambda0(1:MM)=RMM(1:MM)

	if(OPEN)rho_lambda0_alpha(1:MM)=rhoalpha(1:MM)
	if(OPEN)rho_lambda0_betha(1:MM)=rhobeta(1:MM)

	write(*,*) "PLIN values 0", Blambda, Pstepsize, good

	if (Blambda .le. 1.d-1) Pstepsize=Pstepsize*0.5d0
	if (Blambda .ge. 1.0d-1) Pstepsize=Pstepsize*1.2d0
	if (Pstepsize .ge. 1.d0) Pstepsize=1.d0
	if (Blambda .le. 1.d-1 .and. Pstepsize .gt. 1d-4) good=1.d0

	write(*,*) "PLIN values 1", Blambda, Pstepsize, good
	write(*,*) "elast", niter, Elast
	return
	end subroutine P_linear_calc
	
