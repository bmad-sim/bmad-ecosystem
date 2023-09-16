!+
! Subroutine cplx_mat_inverse (mat, mat_inv, ok, print_error)
!
! Takes the inverse of a square matrix using LU Decomposition.
!
! Input:
!   mat(:,:)     -- complex(rp): Input matrix array
!   print_err    -- Logical, optional: If True then the subroutine will type out
!                         a warning message. Default is False.
! Output:
!   mat_inv(:,:) -- complex(rp): inverse of mat1
!   ok           -- Logical, optional: Set False for a singular input matrix.
!-

subroutine cplx_mat_inverse (mat, mat_inv, ok, print_err)

use output_mod, except => cplx_mat_inverse
use f95_lapack, only: la_getrf, la_getri

implicit none

complex(rp) :: mat(:,:), mat_inv(:,:)

real(rp) rcond

integer :: indx(size(mat, 1))

integer info1, info2
logical, optional :: ok, print_err
character(*), parameter :: r_name = 'cplx_mat_inverse'

!

mat_inv = mat

call la_getrf (mat_inv, indx, rcond, info = info1)
call la_getri (mat_inv, indx, info2)

if (info1 /= 0 .or. info2 /= 0) then
  if (logic_option(.false., print_err)) call out_io (s_error$, r_name, 'SINGULAR MATRIX.')
  if (present(ok)) ok = .false.
else
  if (present(ok)) ok = .true.
  if (rcond < 1e-13 .and. logic_option(.false., print_err)) call out_io (s_warn$, r_name, 'NEARLY SINGULAR MATRIX.')
endif

end subroutine

