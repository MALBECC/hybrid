      subroutine calculate_fef(atmsconstr,kforce,maxforce,maxforceatom)
      use ionew
      use scarlett, only: natmsconstr, natot, Ang, eV, kcal, rshiftm,& 
        rref,fef

      implicit none
      integer i,j, at1
      double precision, dimension (25), intent (in) :: kforce
      double precision :: maxforce, tempforce
      integer, intent (out) :: maxforceatom
      integer, dimension(25,25), intent(in) :: atmsconstr
      double precision kf

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
      return
      end subroutine calculate_fef
! -------------------------------------------------------------------
