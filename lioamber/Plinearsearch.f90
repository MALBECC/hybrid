	subroutine P_linearsearch_init()
	use garcha_mod, only : M, OPEN, rho_lambda1, rho_lambda0, rho_lambda1_alpha, rho_lambda0_alpha, rho_lambda1_betha, rho_lambda0_betha, RMM, rhoalpha, rhobeta
	implicit none
	integer :: MM
	write(*,*) "Doing linear search in P"
	MM=M*(M+1)/2
	allocate(rho_lambda1(MM), rho_lambda0(MM))
	rho_lambda0(1:MM)=RMM(1:MM)
	if(OPEN) allocate(rho_lambda1_alpha(MM), rho_lambda0_alpha(MM),rho_lambda1_betha(MM), rho_lambda0_betha(MM))
	if(OPEN) rho_lambda0_alpha(1:MM)=rhoalpha(1:MM)
	if(OPEN) rho_lambda0_betha(1:MM)=rhobeta(1:MM)
	return
	end subroutine P_linearsearch_init

	subroutine P_linear_calc(niter, En, good)
!,E2, Ex,good)
	use garcha_mod, only : M, Md, OPEN, rho_lambda1, rho_lambda0, rho_lambda1_alpha, rho_lambda0_alpha, rho_lambda1_betha, rho_lambda0_betha, RMM, Pstepsize, rhoalpha, rhobeta
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
	if (niter .eq. 1) Pstepsize=1.d0
!        do jj=1,M
!        do kk=jj,M
!	  Rposition=kk+(M2-jj)*(jj-1)/2
!	  rho_lambda1(Rposition)=xnano(jj,kk)
!	end do
!	end do

	rho_lambda1(1:MM)=RMM(1:MM)
	
	if(OPEN) rho_lambda1_alpha(1:MM)=rhoalpha(1:MM)
	if(OPEN) rho_lambda1_betha(1:MM)=rhobeta(1:MM)
	
	E_lambda=0.d0
	do ilambda=0, 10
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

!	RMM(1:MM)=rho_lambda0(1:MM)*(1-Blambda)*Pstepsize+rho_lambda1(1:MM)*Blambda*Pstepsize
	RMM(1:MM)=rho_lambda0(1:MM)*(1.d0-Blambda)+rho_lambda1(1:MM)*Blambda
!a	traza=0.d0
!        do jj=1,M
!          kk=jj
!          Rposition=kk+(M2-jj)*(jj-1)/2
!          traza=traza+(RMM(Rposition))
!        enddo
!	write(*,*) "traza ", traza

        if(OPEN) rhoalpha(1:MM)=rho_lambda0_alpha(1:MM)*(1.d0-Blambda)+rho_lambda1_alpha(1:MM)*Blambda
        if(OPEN) rhobeta(1:MM)=rho_lambda0_betha(1:MM)*(1.d0-Blambda)+rho_lambda1_betha(1:MM)*Blambda

          call int3lu(E2) ! Computes Coulomb part of Fock, and energy on E2
          call g2g_solve_groups(0,Ex,0) ! XC integration / Fock elements


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

	if (Blambda .le. 1.d-1) Pstepsize=Pstepsize*5.d-1
	if (Blambda .ge. 1.0d-1) Pstepsize=1.d0
	if (Blambda .le. 1.d-1 .and. Pstepsize .gt. 1d-4) good=1.d0

	write(*,*) "PLIN values 1", Blambda, Pstepsize, good

	return
	end subroutine P_linear_calc
	
