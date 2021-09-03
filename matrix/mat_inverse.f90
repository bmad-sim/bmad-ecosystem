!+
! Subroutine mat_inverse (mat, mat_inv, ok, print_error)
!
! Takes the inverse of a square matrix using LU Decomposition from
! Numerical Recipes.
!
! Input:
!   mat(:,:)     -- Real(rp): Input matrix array
!   print_err    -- Logical, optional: If True then the subroutine will type out
!                         a warning message. Default is False.
! Output:
!   mat_inv(:,:) -- Real(rp): inverse of mat1
!   ok           -- Logical, optional: Set False for a singular input matrix.
!-

subroutine mat_inverse (mat, mat_inv, ok, print_err)

use output_mod, except => mat_inverse
use nr, only: ludcmp, lubksb

implicit none

real(rp) :: mat(:,:)
real(rp) :: mat_inv(:,:)

real(rp) :: mat2(size(mat, 1), size(mat, 1)), vec(size(mat, 1))
integer :: indx(size(mat, 1))
real(rp) d

integer n, i
logical, optional :: ok, print_err
character(16) :: r_name = 'mat_inverse'

!

n = size (mat, 1)

mat2 = mat  ! use temp mat so as to not change mat

vec(1:n) = maxval(abs(mat), dim = 2)
if (any(vec(1:n) == 0)) then
  if (logic_option(.false., print_err)) call out_io (s_error$, r_name, 'SINGULAR MATRIX.')
  if (present(ok)) ok = .false.
  return
endif
if (present(ok)) ok = .true.

call ludcmp (mat2, indx, d)

mat_inv(1:n,1:n) = 0
forall (i = 1:n) mat_inv(i,i) = 1

do i = 1, n
  call lubksb (mat2, indx, mat_inv(1:n,i))
enddo

end subroutine

