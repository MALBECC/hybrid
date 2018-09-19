!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Gets_orthog_4m( this, method_id, maxval_ld, Xmat, Ymat, Xtrp, Ytrp )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   implicit none
   class(sop), intent(in)  :: this
   integer   , intent(in)  :: method_id
   real*8    , intent(in)  :: maxval_ld
   real*8    , intent(out) :: Xmat(:,:)
   real*8    , intent(out) :: Ymat(:,:)
   real*8    , intent(out) :: Xtrp(:,:)
   real*8    , intent(out) :: Ytrp(:,:)

   logical                 :: error_found


!  Checks and preps
!------------------------------------------------------------------------------!
   if ( this%Nbasis <= 0 ) then
      print*, "ERROR INSIDE Gets_orthog_2m: overlap matrix was never set"
      print*, "ABORTING RUN"; stop
   endif

   error_found = .false.
   error_found = (error_found).or.( this%Nbasis /= size(Xmat,1) )
   error_found = (error_found).or.( this%Nbasis /= size(Xmat,2) )
   error_found = (error_found).or.( this%Nbasis /= size(Ymat,1) )
   error_found = (error_found).or.( this%Nbasis /= size(Ymat,2) )
   error_found = (error_found).or.( this%Nbasis /= size(Xtrp,1) )
   error_found = (error_found).or.( this%Nbasis /= size(Xtrp,2) )
   error_found = (error_found).or.( this%Nbasis /= size(Ytrp,1) )
   error_found = (error_found).or.( this%Nbasis /= size(Ytrp,2) )


!  Returns the appropriate basis matrix (note: Ytrp = Xinv)
!------------------------------------------------------------------------------!

	write(*,*) "flag 1 Gets_orthog_4m"
   call this%Gets_orthog_2m( method_id, maxval_ld, Xmat, Ymat )
	write(*,*) "flag 2 Gets_orthog_4m"
   select case (method_id)
      case (0)
!        Use last method saved inside the object
	write(*,*) "flag 2.1 Gets_orthog_4m"

         Xtrp = this%Xtrp
         Ytrp = this%Ytrp

      case (1,3)
	write(*,*) "flag 2.2 Gets_orthog_4m"

!        Cholesky Decomposition
         Xtrp = transpose( Xmat )
         Ytrp = transpose( Ymat )

      case (2)
	write(*,*) "flag 2.3 Gets_orthog_4m"

!        Symetric/Lowdin Orthogonalization
         Xtrp = Xmat
         Ytrp = Ymat

      case default
         print*,"ERROR INSIDE Gets_orthog_4m: wrong method_id=", method_id
         print*,"ABORTING RUN"; stop

   end select

end subroutine Gets_orthog_4m
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
