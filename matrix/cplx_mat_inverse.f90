!+
! Subroutine cplx_mat_inverse (mat, mat_inv, ok, print_error)
!
! Takes the inverse of a square matrix using LU Decomposition from
! Numerical Recipes.
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

implicit none

complex(rp) :: mat(:,:)
complex(rp) :: mat_inv(:,:)

complex(rp) :: mat2(size(mat, 1), size(mat, 2))
integer :: indx(size(mat, 1))

real(rp) d

integer n, i
logical, optional :: ok, print_err
character(16) :: r_name = 'mat_inverse'

!

n = size (mat, 1)
mat2 = mat  ! use temp mat so as to not change mat

if (any(maxval(abs(mat), dim = 2) == 0)) then
  if (logic_option(.false., print_err)) call out_io (s_error$, r_name, 'SINGULAR MATRIX.')
  if (present(ok)) ok = .false.
  return
endif
if (present(ok)) ok = .true.

call cplx_ludcmp (mat2, indx, d)

mat_inv(1:n,1:n) = (0,0) 
forall (i = 1:n) mat_inv(i,i) = (1,0)

do i = 1, n
  call cplx_lubksb (mat2, indx, mat_inv(1:n,i))
enddo

end subroutine

