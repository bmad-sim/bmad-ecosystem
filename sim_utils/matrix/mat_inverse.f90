!+
! Subroutine mat_inverse (mat, mat_inv, vec0, vec0_inv, ok, print_error)
!
! Takes the inverse, using LU Decomposition, of a map
!   r_out = mat * r_in + vec0
! whee mat is a square matrix.
!
! Input:
!   mat(:,:)      -- real(rp): Input matrix 
!   vec0(:)       -- real(rp), optional: O^th order part of the map.
!   print_err     -- logical, optional: If True then the subroutine will type out
!                         a warning message. Default is False.
! Output:
!   mat_inv(:,:)  -- real(rp): inverse of mat1
!   vec0_inv(:)   -- real(rp): 0th order part of the inverse map.
!   ok            -- logical, optional: Set False for a singular input matrix.
!-

subroutine mat_inverse (mat, mat_inv, vec0, vec0_inv, ok, print_err)

use output_mod, except => mat_inverse
use f95_lapack, only: la_getrf, la_getri

implicit none

real(rp) :: mat(:,:), mat_inv(:,:)
real(rp), optional :: vec0(:), vec0_inv(:)
real(rp) rcond

integer :: indx(size(mat, 1))

integer info1, info2
logical, optional :: ok, print_err
character(*), parameter :: r_name = 'mat_inverse'

!

mat_inv = mat

call la_getrf (mat_inv, indx, rcond, info = info1)
call la_getri (mat_inv, indx, info2)

if (info1 /= 0 .or. info2 /= 0) then
  if (logic_option(.false., print_err)) call out_io (s_error$, r_name, 'SINGULAR MATRIX.')
  if (present(ok)) ok = .false.
  return
else
  if (present(ok)) ok = .true.
  if (rcond < 1e-13 .and. logic_option(.false., print_err)) call out_io (s_warn$, r_name, 'NEARLY SINGULAR MATRIX.')
endif

if (present(vec0)) then
  vec0_inv = -matmul(mat_inv, vec0)
endif

end subroutine
