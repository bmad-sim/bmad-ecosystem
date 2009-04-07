!+
! Subroutine mat_symp_conj (mat1, mat2)
!
! Subroutine to take the symplectic conjugate of a square matrix.
!
! Modules needed:
!   use dcslib
!
! Input:
!   mat1(:, :) -- Real(rp): Input matrix.
!
! Output:
!   mat2(:, :) -- Real(rp): Symplectic conjugate of mat1.
!-

#include "CESR_platform.inc"

subroutine mat_symp_conj(mat1, mat2)

  use precision_def

  implicit none

  integer i, j, nn

  real(rp) mat1(:,:), mat2(:,:)
  real(rp), allocatable, save :: mat22(:,:)

! check bounds

  nn = size(mat1, 1)

  if (mod(nn, 2) /= 0) then
    print *, 'ERROR IN MAT_SYMP_CONJ SUBROUTINE: ARRAY SIZE IS NOT EVEN'
    call err_exit
  endif

  if (allocated(mat22)) then
    if (size(mat22, 1) < nn) then
      deallocate(mat22)
      allocate(mat22(nn,nn))
    endif
  else
    allocate (mat22(nn,nn))
  endif

! compute conjugate
! mat22 is used in case mat1 and mat2 are the same.

  do i = 1, nn, 2
    do j = 1, nn, 2
      mat22(i,   j)   =  mat1(j+1, i+1)
      mat22(i,   j+1) = -mat1(j,   i+1)
      mat22(i+1, j)   = -mat1(j+1, i)
      mat22(i+1, j+1) =  mat1(j,   i)
    enddo
  enddo

  mat2 = mat22(1:nn,1:nn)

end subroutine
