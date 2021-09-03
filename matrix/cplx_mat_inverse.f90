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

subroutine cplx_mat_inverse (mat_r, mat_i, inv_r, inv_i, ok, print_err)

use output_mod, except => cplx_mat_inverse

implicit none

real(rp) :: mat_r(:,:)
real(rp) :: mat_i(:,:)
real(rp) :: inv_r(:,:)
real(rp) :: inv_i(:,:)

complex(rp) :: mat(size(mat_r,1),size(mat_r,1))
complex(rp) :: mat_inv(size(mat_r,1),size(mat_r,1))

complex(rp) :: mat2(size(mat_r, 1), size(mat_r, 1))
real(rp) :: vec(size(mat_r, 1))
integer :: indx(size(mat_r, 1))
real(rp) d

integer n, i
logical, optional :: ok, print_err
character(16) :: r_name = 'mat_inverse'

!

mat = CMPLX(mat_r,mat_i)

n = size (mat, 1)

mat2 = mat  ! use temp mat so as to not change mat

vec(1:n) = maxval(abs(mat), dim = 2)
if (any(vec(1:n) == 0)) then
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

inv_r = REAL(mat_inv)
inv_i = AIMAG(mat_inv)

end subroutine

