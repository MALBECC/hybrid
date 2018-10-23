      subroutine compute_cutsqmmm(at_MM_cut_QMMM, &
      istepconstr,radbloqmmm,rcorteqmmm,nroaa,atxres)
      
      use scarlett, only: r_cut_list_QMMM,MM_freeze_list,na_u, &
      MM_freeze_list,rclas,nac,listqmmm

      implicit none
      double precision, intent(in) :: rcorteqmmm, radbloqmmm
!      double precision, allocatable, dimension(:,:) , intent(out):: &
!      r_cut_QMMM, F_cut_QMMM
!      double precision, allocatable, dimension(:), intent(out) :: &
!      Iz_cut_QMMM
      integer, intent(in) :: nroaa, atxres(20000)
      integer, intent(inout) :: at_MM_cut_QMMM

! internal variables
      integer :: r_cut_pos
      double precision :: r12, cm(3,20000), dist, dist2
      integer :: i, j, k, l, i_qm, i_mm, istepconstr 
      logical :: done, done_freeze, done_QMMM !control variables


!	    if (allocated(r_cut_QMMM)) deallocate(r_cut_QMMM)
!	    if (allocated(F_cut_QMMM)) deallocate(F_cut_QMMM)
!	    if (allocated(Iz_cut_QMMM)) deallocate(Iz_cut_QMMM)
	    r_cut_list_QMMM=0
	    r_cut_pos=0
	    at_MM_cut_QMMM=0

	    if (istepconstr.eq.1) then
		MM_freeze_list=.true.
		do i_qm=1,na_u
		  MM_freeze_list(i_qm)=.false.
		end do
	    end if

	    do i_mm=1, nac !MM atoms
	      i_qm=0
	      done=.false.
              done_freeze=.false.
              done_QMMM=.false.
	      do while (i_qm .lt. na_u .and. .not. done) !QM atoms
	        i_qm=i_qm+1
                r12=(rclas(1,i_qm)-rclas(1,i_mm+na_u))**2.d0 + &
                    (rclas(2,i_qm)-rclas(2,i_mm+na_u))**2.d0 + &
                    (rclas(3,i_qm)-rclas(3,i_mm+na_u))**2.d0

	        if(r12 .lt. rcorteqmmm .and. .not. done_QMMM) then
	          done_QMMM=.true.
	          at_MM_cut_QMMM=at_MM_cut_QMMM+1
	          r_cut_pos=r_cut_pos+1
	          r_cut_list_QMMM(i_mm)=r_cut_pos
	        end if



		if (istepconstr.eq.1) then !define lista de movimiento para la 1er foto
	          if(r12 .lt. radbloqmmm .and. .not. done_freeze) then
	             MM_freeze_list(i_mm+na_u)=.false.
	             done_freeze=.true.
	          end if
		end if



		done=done_QMMM .and. done_freeze
	      end do
	    end do


!	  allocate (r_cut_QMMM(3,at_MM_cut_QMMM+na_u), &
!          F_cut_QMMM(3,at_MM_cut_QMMM+na_u), &
!          Iz_cut_QMMM(at_MM_cut_QMMM+na_u))

!	  r_cut_QMMM=0.d0
!	  F_cut_QMMM=0.d0
!	  Iz_cut_QMMM=0

!Recalcula list_qmmm

!Calcula centro de masa de todos los residuos        
        cm=0.d0
        k=na_u+1
        do i=1,nroaa
         do j=1,atxres(i)
        cm(1:3,i)=cm(1:3,i)+rclas(1:3,k)
        k=k+1
         enddo
        cm(1:3,i)=cm(1:3,i)/atxres(i)
        enddo

!Se fija si tiene que hacer lj o no

        dist=0.0
        dist2=rcorteqmmm**2
        listqmmm=1
        do l=1,na_u
        k=1
         do i=1,nroaa
          dist=(rclas(1,l)-cm(1,i))**2+ &
               (rclas(2,l)-cm(2,i))**2+ &
               (rclas(3,l)-cm(3,i))**2
          do j=1,atxres(i)
           if(dist.le.dist2) listqmmm(k)=0 !0 interacción ON
           k=k+1
          enddo
         enddo
        enddo

      end subroutine compute_cutsqmmm
