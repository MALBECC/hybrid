!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine Gets_orthog_2m( this, method_id, maxval_ld, Xmat, Ymat )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   use liosubs_math, only: purge_zeros_m, matmul3
   implicit none
   class(sop), intent(in)  :: this
   integer   , intent(in)  :: method_id
   real*8    , intent(in)  :: maxval_ld
   real*8    , intent(out) :: Xmat(:,:)
   real*8    , intent(out) :: Ymat(:,:)

   real*8    , allocatable :: Gmat_li(:,:), Ginv_li(:,:)
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

   allocate( Gmat_li( this%Nbasis, this%Nbasis ) )
   allocate( Ginv_li( this%Nbasis, this%Nbasis ) )

   Gmat_li = this%Gmat
   Ginv_li = this%Ginv
   call this%Drop_ldvals( maxval_ld, Gmat_li, Ginv_li )


!  Returns the appropriate basis matrix (note: Ytrp = Xinv)
!------------------------------------------------------------------------------!
   select case (method_id)
      case (0)
!        Use last method saved inside the object
         Xmat = this%Xmat
         Ymat = this%Ymat

      case (1)
!        Cholesky Decomposition
!        S = Y * Yt = L * Lt
         call matmul3( this%Umat, Ginv_li, this%Vtrp, Xmat )
         call matmul3( this%Umat, Gmat_li, this%Vtrp, Ymat )

      case (2)
!        Symetric/Lowdin Orthogonalization
!        S = Y * Yt = (U*sq*Ut)*(U*sq*Ut)
         call matmul3( this%Umat, Ginv_li, this%Utrp, Xmat )
         call matmul3( this%Umat, Gmat_li, this%Utrp, Ymat )

      case (3)
!        Canonical Orthogonalization
!        S = Y * Yt = (U*sq)*(sq*Ut)
         Xmat = matmul( this%Umat, Ginv_li )
         Ymat = matmul( this%Umat, Gmat_li )

      case default
         print*,"ERROR INSIDE Gets_orthog_2m: wrong method_id=", method_id
         print*,"ABORTING RUN"; stop

   end select

   call purge_zeros_m( 1.0d-10, Xmat )
   call purge_zeros_m( 1.0d-10, Ymat )

end subroutine Gets_orthog_2m
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
