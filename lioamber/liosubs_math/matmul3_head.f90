!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine matmul3_ddd( Amat, Bmat, Cmat, Dmat )
   implicit none
   real*8    , intent(in) :: Amat(:,:)
   real*8    , intent(in) :: Bmat(:,:)
   real*8    , intent(in) :: Cmat(:,:)
   real*8                 :: Dmat(size(Bmat,1),size(Bmat,2))
#  include "matmul3_body.f90"
end subroutine matmul3_ddd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine  matmul3_dcd( Amat, Bmat, Cmat, Dmat )
   implicit none
   real*8    , intent(in) :: Amat(:,:)
   complex*16, intent(in) :: Bmat(:,:)
   real*8    , intent(in) :: Cmat(:,:)
   complex*16             :: Dmat(size(Bmat,1),size(Bmat,2))
#  include "matmul3_body.f90"
end subroutine matmul3_dcd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
