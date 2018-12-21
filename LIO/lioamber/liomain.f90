!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LIOMAIN.F90  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains the liomain subroutine, which performs several common     !
! procedures before and after calling SCF either by LIO alone or when          !
! performing in tantem with AMBER/GROMACS. Currently the routine:              !
! * Allocates matrices needed for SCF calculations.                            !
! * Calls SCF or SCFOP for closed/open shell calculations.                     !
! Other interfacial routines included are:                                     !
! * do_forces        (calculates forces/gradients)                             !
! * do_dipole        (calculates dipole moment)                                !
! * do_population_analysis (performs the analysis required)                    !
! * do_fukui         (performs Fukui function calculation and printing)        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine liomain(E, dipxyz)
    use garcha_mod, only : M, Smat, RealRho, OPEN, writeforces, energy_freq,   &
                           restart_freq, npas, sqsm, mulliken, lowdin, spinpop,&
                           dipole, doing_ehrenfest, first_step,                &
                           Eorbs, fukui, print_coeffs, steep, idip, calc_propM
    use ecp_mod   , only : ecpmode, IzECP
    use ehrensubs,  only : ehrendyn_main

    implicit none
    REAL*8, intent(inout) :: dipxyz(3), E
    integer :: idip_scrach
    logical :: calc_prop

    call g2g_timer_sum_start("Total")

    if (.not.allocated(Smat))    allocate(Smat(M,M))
    if (.not.allocated(RealRho)) allocate(RealRho(M,M))
    if (.not.allocated(sqsm))    allocate(sqsm(M,M))
    if (.not.allocated(Eorbs))   allocate(Eorbs(M))


    if (steep) then
      idip_scrach=idip
      idip=0 !skip dipole calculation in geometry optimization
      call do_steep(E)
      idip=idip_scrach
    end if


!------------------------------------------------------------------------------!
! FFR - Option to do ehrenfest
    if ( doing_ehrenfest ) then
       if ( first_step ) call SCF( E, dipxyz )
       call ehrendyn_main( E, dipxyz )
    else
       call SCF(E)
    endif


    calc_prop=.false.
    if (MOD(npas, energy_freq).eq.0) calc_prop=.true.
    if (calc_propM) calc_prop=.true.
    if ((restart_freq.gt.0).and.(MOD(npas, restart_freq).eq.0)) call do_restart(88)


    ! Perform Mulliken and Lowdin analysis, get fukui functions and dipole.
    if (calc_prop) then

        if (mulliken.or.lowdin.or.spinpop) call do_population_analysis()

        if (dipole) call do_dipole(dipxyz, 69)

        if (fukui) call do_fukui()

        if (writeforces) then
            if (ecpmode) stop "ECP does not feature forces calculation."
            call do_forces(123)
        endif

        if (print_coeffs) call write_orbitals(29)
    endif

    call g2g_timer_sum_pause("Total")

    return
end subroutine liomain


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates forces for QM and MM regions and writes them to output.           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_forces(uid)

    use garcha_mod, only : natom, nsol

    implicit none
    integer, intent(in) :: uid
    integer             :: k
    real*8, allocatable :: dxyzqm(:,:), dxyzcl(:,:)

    call g2g_timer_start('Forces')
    open(unit=uid, file='forces')

    allocate ( dxyzqm(3, natom) )
    dxyzqm = 0.0

    call dft_get_qm_forces(dxyzqm)
    if (nsol.gt.0) then
        allocate ( dxyzcl(3, natom+nsol) )
        dxyzcl = 0.0
        call dft_get_mm_forces(dxyzcl, dxyzqm)
    endif

    call write_forces(dxyzqm, natom, 0, uid)

    deallocate (dxyzqm)



    if(nsol.gt.0) then
        call write_forces(dxyzcl, nsol, 0, uid)
        deallocate (dxyzcl)
    endif
    call g2g_timer_stop('Forces')

    return
end subroutine do_forces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_DIPOLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Sets variables up and calls dipole calculation.                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_dipole(dipxyz, uid)
    implicit none
    integer, intent(in)    :: uid
    real*8 , intent(inout) :: dipxyz(3)
    real*8                 :: u

    call g2g_timer_start('Dipole')
    call dip(dipxyz)
    u = sqrt(dipxyz(1)**2 + dipxyz(2)**2 + dipxyz(3)**2)

    call write_dipole(dipxyz, u, uid, .true.)
    call write_dipole(dipxyz, u, uid, .false.)
    call g2g_timer_stop('Dipole')

    return
end subroutine do_dipole
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_POPULATION_ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs the different population analyisis available.                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_population_analysis()
   use garcha_mod, only : RMM, Smat, RealRho, M, Enucl, Nuc, Iz, natom, &
                          mulliken, lowdin, sqsm, a, c, d, r, Iz, ncont, NORM,&
                          M, Md, spinpop,OPEN, rhoalpha,rhobeta

   use ECP_mod   , only : ecpmode, IzECP
   use faint_cpu, only: int1

   implicit none
   integer :: M1, M5, IzUsed(natom), kk
   real*8, allocatable, dimension(:,:)  :: RealRho_alpha, RealRho_betha
   real*8  :: q(natom)

   ! Needed until we dispose of RMM.
   M1=1 ; M5=1+M*(M+1)

   ! Iz used to write the population file.
   IzUsed = Iz
   if (ecpmode) IzUsed = IzECP

   ! Decompresses and fixes S and RealRho matrixes, which are needed for
   ! population analysis.
   call int1(Enucl,RMM,Smat,Nuc,a,c,d,r,Iz,ncont,NORM,natom,M,Md)
   call spunpack('L',M,RMM(M5),Smat)
   call spunpack('L',M,RMM(M1),RealRho)
   call fixrho(M,RealRho)

   do kk=1,natom
       q(kk) = real(IzUsed(kk))
   enddo

   ! Performs Mulliken Population Analysis if required.
   if (mulliken) then
       call g2g_timer_start('Mulliken')
       call mulliken_calc(natom, M, RealRho, Smat, Nuc, q)
       call write_population(85, natom, Iz, q, 0)
       call g2g_timer_stop('Mulliken')
   endif

   ! Performs Löwdin Population Analysis if required.
   if (lowdin) then
       call g2g_timer_start('Lowdin')
       call lowdin_calc(M, natom, RealRho, sqsm, Nuc, q)
       call write_population(85, natom, Iz, q, 1)
       call g2g_timer_stop('Lowdin')
   endif

   if (spinpop) then
       if (.not. OPEN) then
         write(*,*) "cant perform a spin population analysis in a close shell calculation"
         return
       end if
       allocate (RealRho_alpha(M,M), RealRho_betha(M,M))
       call spunpack('L',M,rhoalpha(1),RealRho_alpha) !pasa vector a matriz
       call spunpack('L',M,rhobeta(1),RealRho_betha) !pasa vector a matriz
       q=0
       call spin_pop_calc(natom, M, RealRho_alpha, RealRho_betha, Smat, Nuc, q)
       call write_population(natom, IzUsed, q, 2, 85)
   end if

   return
endsubroutine do_population_analysis
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_FUKUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs Fukui function calls and printing.                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_fukui()
    use garcha_mod, only : X, NCO, M, natom, Nuc, Smat, Eorbs, Iz, OPEN

    implicit none
    real*8  :: fukuim(natom), fukuin(natom), fukuip(natom), softness


    call g2g_timer_start("Fukui")
    if (OPEN) then
        ! TO-DO
        ! Add call fukui_calc_OS when openshell is available.
    else
        call fukui_calc(X(1,M*2+1), NCO, M, natom, Nuc, Smat, fukuim, fukuip, &
                        fukuin, Eorbs)
        call get_softness(Eorbs(NCO-1), Eorbs(NCO), Eorbs(NCO-1), Eorbs(NCO), &
                          softness)
        call write_fukui(fukuim, fukuip, fukuin, natom, Iz, softness)
    endif
    call g2g_timer_stop("Fukui")

    return
end subroutine do_fukui
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_FUKUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs Fukui function calls and printing.                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_restart(UID)
   use garcha_mod, only : OPEN, NCO, NUNP, M, MO_coef_at, MO_coef_at_b, indexii
   use fileio    , only : write_coef_restart
   implicit none
   integer, intent(in) :: UID
   integer             :: NCOb, icount, jcount, coef_ind
   real*8, allocatable :: coef(:,:), coef_b(:,:)

   allocate(coef(M, NCO))
   do icount=1, M
   do jcount=1, NCO
      coef_ind = icount + M*(jcount-1)
      coef(indexii(icount), jcount) = MO_coef_at(coef_ind)
   enddo
   enddo


   if (OPEN) then
      NCOb = NCO + NUNP
      allocate(coef_b(M, NCOb))

      do icount=1, M
      do jcount=1, NCOb
         coef_ind = icount + M*(jcount-1)
         coef_b(indexii(icount), jcount) = MO_coef_at_b(coef_ind)
      enddo
      enddo

      call write_coef_restart(coef, coef_b, M, NCO, NCOb, UID)
      deallocate(coef_b)
   else
      call write_coef_restart(coef, M, NCO, UID)
   endif

   deallocate(coef)
   return
end subroutine do_restart
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!