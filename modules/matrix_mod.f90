module matrix_mod

use precision_def
use sim_utils
use basic_bmad_mod

contains
    
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine mat_inverse (mat, mat_inv, ok, print_error)
!
! Takes the inverse of a square matrix using LU Decomposition from
! Numerical Recipes.
!
! Modules needed:
!   bmad_interface
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

use nr

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

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Function mat_symp_error (mat) result (error)
!
! Routine to check the symplecticity of a square matrix. The error is
! defined to via:
!     error = maxval (abs (Mat_transpose * S * Mat - S))
!
! Modules Needed:
!
! Input:
!   mat(:,:) -- Real(rp): Matrix to check
!
! Output:
!   error -- Real(rp): difference from symplecticity:
!             = 0    --> perfect.
!             = 1e-4 --> Fair, but I wouldn't use for long term tracking.
!             = 1    --> Terrible.
!-

function mat_symp_error (mat) result (error)

implicit none

integer i, j, n

real(rp), intent(in) :: mat(:,:)
real(rp) :: m2(size(mat, 1), size(mat, 1))
real(rp) error

logical :: debug = .false.

!

n = ubound(mat, 1)

if (mod(n, 2) /= 0) then
  print *, 'ERROR IN MAT_SYMP_ERROR: MATRIX DOES NOT HAVE EVEN SIZE'
  if (bmad_status%exit_on_error) call err_exit
endif

!

do j = 2, n, 2                  ! m2 = mat_transpose * S
  m2(1:n,j)   = -mat(j-1,1:n)
  m2(1:n,j-1) =  mat(j,  1:n)
enddo

m2 = matmul (m2, mat)

do j = 2, n, 2
  m2(j, j-1) = m2(j, j-1) - 1
  m2(j-1, j) = m2(j-1, j) + 1
enddo

error = maxval(abs(m2))

if (debug) then
  print *, 'MAT_SYMP_ERROR:', error
  do i = 1, n
    print '(5x, (6f11.5))', &
                      (max(min(m2(i, j), 999.0_rp), -999.0_rp), j = 1, n)
  enddo
endif

end function mat_symp_error

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine mat_symplectify (mat_in, mat_symp, f_scale)
!
! Subroutine to form the symplectic matrix that is as "close" as possible to the input matrix.
! If mat_in is symplectic, and f_length is 1, then mat_sym will be equal to mat_in.
! 
! The f_scale argument is used to "scale" mat_symp such that if mat_in is symplectic the
! following is true:
!   f_scale = 0   --> mat_symp = I (unit matrix, always true independent of mat_in)
!   f_scale = 1/N --> mat_symp^N = mat_in
!   f_scale = 1   --> mat_symp = mat_in
!   
! Modules needed:
!   matrix_mod
!
! Input:
!   mat_in(:,:) -- Real(rp): Input matrix to symplectify.
!   f_scale     -- Real(rp), optional: Scaling factor. Default is 1.
!
! Output:
!   mat_symp(:,:) -- Real(rp): Symplectic output matrix
!-

subroutine mat_symplectify (mat_in, mat_symp, f_scale)

implicit none

real(rp), intent(in)  :: mat_in(:,:)
real(rp), intent(out) :: mat_symp(:,:)
real(rp), intent(in), optional :: f_scale

real(rp), dimension(size(mat_in, 1), size(mat_in, 2)) :: m1, m2, m3, m_symmetric

integer n,i,j

!

n = ubound(mat_in, 1)

if (ubound(mat_symp, 1) /= n) then
  print *, 'ERROR IN MAT_SYMPLECTIFY: UNEQUAL MATRIX SIZES.'
  if (bmad_status%exit_on_error) call err_exit
endif

! Form the symmetrix matrix:
!         m1 = S * (1 - mat_in) * (1 + mat_in)^-1
!         m_symmetric = (m1 + m1_transpose) / 2

m1 =-mat_in
forall (i = 1:n) m1(i,i) = m1(i,i) + 1

m2 = mat_in
forall (i = 1:n) m2(i,i) = m2(i,i) + 1
call mat_inverse (m2, m2)

m3 = matmul(m1, m2)

do i = 2, n, 2
  m1(i-1,1:n) =  m3(i,  1:n)
  m1(i,  1:n) = -m3(i-1,1:n)
enddo

forall (i = 1:n, j = 1:n) m_symmetric(i,j) = (m1(i,j) + m1(j,i)) / 2

! Scale by f_scale

if (present(f_scale)) m_symmetric = f_scale * m_symmetric

! Form the symplectic matrix
!         m_symp = (1 + S * m_symmetric) * (1 - S * m_symmetric)^-1

do i = 2, n, 2
  m1(i-1,1:n) =  m_symmetric(i,  1:n)
  m1(i,  1:n) = -m_symmetric(i-1,1:n)
enddo

m2 = m1
forall (i = 1:n) m2(i,i) = m2(i,i) + 1

m3 = -m1
forall (i = 1:n) m3(i,i) = m3(i,i) + 1
call mat_inverse (m3, m3)

mat_symp = matmul (m2, m3)

end subroutine

end module
