      subroutine fe_opt(rcorteqmmm, radbloqmmm, Etot,&
       do_SCF, do_QM_forces, do_properties, istp, step,&
       nbond, nangle, ndihe, nimp, Etot_amber, Elj,&
       Etots, constropt,nconstr, nstepconstr, typeconstr, kforce, ro,&
       rt, coef, atmsconstr, ndists, istepconstr, rcortemm,&
       radblommbond, optimization_lvl, dt, sfc, water,&
       imm,rini,rfin,maxforce,maxforceatom,rconverged,ntcon,&
       nfree,cmcf)

      use ionew
      use scarlett, only: natmsconstr, natot, Ang, eV, kcal, rshiftm,&
        rshiftm2, rref, fef, fdummy, cfdummy, tt, masst, kn,vn, mn, &
        rclas, vat, Ekinion, tempion, tt, tempqm, tempinit, blockall, &
        cov_matrix, rshxrshm, rshiftsd, fedynamic, tauber, innermax
 
      implicit none
      integer i, j, k, inneri, at1, at2, k1, k2
      

!------------------------------------------------ Variables required for free energy gradient calculations

      double precision, intent(inout) :: maxforce  
      integer, intent (out) :: maxforceatom
      double precision :: tempforce, kf
      logical :: rconverged
      double precision :: Fmax,F_i

!------------------------------------------------ Variables required for inner MD

      integer :: ntcon !number of restrained degrees of freedom JOTA 
      integer :: nfree !number of atoms without a contrain of movement
      integer :: cmcf !Center of Mass Coordinates Fixed JOTA

!------------------------------------------------ Variables required for energy and forces calculations

      double precision, intent(in) :: rcortemm ! distance for LJ & Coulomb MM interaction
      double precision, intent(in) :: radblommbond !parche para omitir bonds en extremos terminales, no se computan bonds con distancias mayores a radblommbond
      double precision :: rcorteqmmm !Distance for QM-MM interaction
      double precision, intent(in) :: radbloqmmm !Distance that allow to move MM atoms from QM sub-system
      double precision, intent(inout) :: Etot ! QM + QM-MM interaction energy
      double precision, intent(in) :: dt !time step
      logical, intent(in) :: do_SCF, do_QM_forces !control for make new calculation of rho, forces in actual step
      logical, intent(in) :: do_properties !control for lio properties calculation
      integer, intent(inout) :: step !total number of steps in a MMxQM movement (not tested in this version)
      integer, intent(in) :: istp !number of move step for each restrain starting on 1
      integer, intent(in) :: imm !MM step by QM each step
      integer, intent(in) :: nbond, nangle, ndihe, nimp !number of bonds, angles, dihedrals and impropers defined in amber.parm
      double precision, intent(out) :: Etot_amber !total MM energy
      double precision, intent(out) :: Elj !LJ interaction (only QMMM)
      double precision, intent(out) :: Etots !QM+QMMM+MM energy
      logical, intent(in) :: constropt !activate restrain optimizaion
      integer, intent(in) :: nconstr !number of constrains
      integer, intent(in) :: nstepconstr !numero de pasos en los que barrera la coordenada de reaccion (limitado a 1-100 hay q cambiar esto, Nick)
      integer, dimension(20), intent(in) :: typeconstr !type of cosntrain (1 to 8)
      double precision, dimension(20), intent(in) :: kforce !force constant of constrain i
      double precision :: rini,rfin  !initial and end value of reaction coordinate
      double precision, dimension(20), intent(in) :: ro ! fixed value of constrain in case nconstr > 1 for contrains 2+
      double precision, dimension(20), intent(inout) :: rt ! value of reaction coordinate in constrain i
      double precision, dimension(20,10), intent(in) :: coef ! coeficients for typeconstr=8
      integer, dimension(20,20), intent(in) :: atmsconstr
      integer, dimension(20), intent(in) :: ndists !atomos incluidos en la coordenada de reaccion
      integer, intent(in) :: istepconstr !step of restraint 
      integer, intent(in) :: optimization_lvl ! level of movement in optimization scheme (1 only QM atoms with restrain,2 only MM atoms, 3 all)
      double precision, intent(in) :: sfc
      logical, intent(in) :: water

!------------------------------------------------- 

      rref=rclas
      rshiftm=0.d0
      rshiftm2=0.d0
      inneri=1
      rconverged=.false.       

!         rshiftsd=0.d0
!         if (.not. relaxd) then
!           inneri=1
!           rconverged=.false.
       
      call vmb(natot,tempinit,masst,vat,cmcf,blockall,ntcon)

      if (fedynamic .eq. 1) then
        mn=dble(3*natot-ntcon-cmcf)*tt*8.617d-5*(50.d0*dt)**2
        write(6,'(/,a)') 'Calculating Nose mass as Ndf*Tt*KB*(50dt)**2'
        write(6,'(a,2x,F30.18)') "mn =", mn
      endif

      do while ((.not. rconverged) .and. (inneri .le. innermax))  ! <<<<<<<<<<<<<<< DM in FE Calculations para MB con ts de 0.1 fs

        call do_energy_forces(rcorteqmmm, radbloqmmm, Etot,&
        do_SCF, do_QM_forces, do_properties, istp, step,&
        nbond, nangle, ndihe, nimp, Etot_amber, Elj,&
        Etots, constropt,nconstr, nstepconstr, typeconstr, kforce, ro,&
        rt, coef, atmsconstr, ndists, istepconstr, rcortemm,& 
        radblommbond, optimization_lvl, dt, sfc, water,&
        imm,rini,rfin)

       if (fedynamic .eq. 0) then !berendsen
         call berendsen(inneri,3,natot,cfdummy,dt,tauber,masst, &
         ntcon,vat,rclas,Ekinion,tempion,tt,nfree,cmcf)
       elseif (fedynamic .eq. 1) then !nose
         call nose(inneri,natot,cfdummy,tt,dt,masst,mn,ntcon,vat,&
         rclas,Ekinion,kn,vn,tempion,nfree,cmcf)
       endif


! Save summ of rshifm and rshiftm2
        do i=1,natmsconstr
          do j=1,3
            at1=atmsconstr(1,i)
            rshiftm(j,at1)=rshiftm(j,at1)+(rclas(j,at1)-rref(j,at1))
            rshiftm2(j,at1)=rshiftm2(j,at1)+((rclas(j,at1)-rref(j,at1))**2)
          enddo
        enddo

! Save rshxrshm elements
        do i=1,natmsconstr
          at1=atmsconstr(1,i)
          do j=1,natmsconstr
            at2=atmsconstr(1,j)
              do k1=1,3
                do k2=1,3
                  rshxrshm((3*j-2)+k1-1,(3*i-2)+k2-1)=&
                  rshxrshm((3*j-2)+k1-1,(3*i-2)+k2-1)+&
                  (rclas(k2,at1)-rref(k2,at1))*&
                  (rclas(k1,at2)-rref(k1,at2))
                enddo
              enddo
          enddo
        enddo

! Check rshiftm convergency every 1000 steps, after the first 5000 steps
      if (MOD(inneri,1000) .eq. 0 .and. inneri .gt. 5000) then
        do i=1,natmsconstr
          at1=atmsconstr(1,i)
          do j=1,3
            rshiftsd(j,at1)=&
            sqrt(rshiftm2(j,at1)-(rshiftm(j,at1)**2))/inneri
          enddo
        enddo
        rconverged=.false.
        k=0
        do i=1,natmsconstr
          at1=atmsconstr(1,i)
          do j=1,3
            if (abs((rshiftsd(j,at1)*inneri/rshiftm(j,at1))) .le. 0.1d0) k=k+1
          enddo
        enddo
        rconverged=(k .eq. 3*natmsconstr)
        if(rconverged) write(*,*) "FEG converged in ",inneri," steps"
       endif

!Escribe cosas 

      call calculateTemp(Ekinion,tempion,tempqm,vat,ntcon,&
      nfree,cmcf)


      if (MOD(inneri,1000) .eq. 0 .and. inneri .gt. 5000) then
      do i=1,natmsconstr
        at1=atmsconstr(1,i)
        do j=1,natmsconstr
          at2=atmsconstr(1,j)
          do k1=1,3
            do k2=1,3
              cov_matrix((3*j-2)+k1-1,(3*i-2)+k2-1)=&
              (rshxrshm((3*j-2)+k1-1,(3*i-2)+k2-1)-&
              rshiftm(k2,at1)*rshiftm(k1,at2)/dble(inneri))/dble(inneri)
            enddo
          enddo
        enddo
      enddo
      endif

! sets variables for next cycle
          fdummy = 0.d0
          inneri = inneri + 1
      enddo! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< DM in FE Calculations

      rshiftm=rshiftm/dble(inneri-1)
      rshxrshm=rshxrshm/dble(inneri-1)

      do i=1,natmsconstr
        at1=atmsconstr(1,i)
        do j=1,natmsconstr
          at2=atmsconstr(1,j)
          do k1=1,3
            do k2=1,3
              cov_matrix((3*j-2)+k1-1,(3*i-2)+k2-1)=&
              rshxrshm((3*j-2)+k1-1,(3*i-2)+k2-1)-&
              rshiftm(k2,at1)*rshiftm(k1,at2)
            enddo
          enddo
        enddo
      enddo

! Once rshift convergence is achieved (or inneri .gt. 25000), calculate fef

      kf=kforce(1)
      do i=1,natmsconstr
        do j=1,3
          at1=atmsconstr(1,i)
          fef(j,at1)=2.d0*kf*rshiftm(j,at1)/Ang !Son fuerzas, no gradientes
        enddo
      enddo


! fef queda en Hartree/Bohr

      fef=fef*eV/(Ang*kcal)

      maxforce=0.d0
      maxforceatom=0

      do i=1,natmsconstr
        at1=atmsconstr(1,i)
        tempforce=0.d0
        do j=1,3
          tempforce=tempforce+fef(j,at1)**2
        enddo
      if(tempforce .gt. maxforce) then
        maxforce = tempforce
        maxforceatom = at1
      endif
      enddo

!Pone en rclas las posiciones originales de todos los atomos (rever si
!esto es necesario
!     rclas=rref

!Pisa en rclas las posiciones los atomos con restraint con sus posiciones medias

      do i=1,natmsconstr
        at1=atmsconstr(1,i)
          do j=1,3
            rclas(j,at1)=rshiftm(j,at1)+rref(j,at1)
          enddo
      enddo

!Pone fuerzas medias en cfdummy

      cfdummy=fef

!Mata velocidades

!      vat=0

!Check max force
        Fmax=0.d0
        do i=1, natot
          F_i=cfdummy(1,i)**2 + cfdummy(2,i)**2 + cfdummy(3,i)**2
          F_i=sqrt(F_i)
          if (F_i .gt. Fmax) Fmax=F_i
        end do
!       write (456456,*) istp, Fmax
      return  
 
      end subroutine fe_opt
