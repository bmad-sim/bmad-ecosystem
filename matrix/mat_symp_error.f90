!+
! Function mat_symp_error (mat) result (error)
!
! Routine to check the symplecticity of a square matrix. The error is
! defined to via:
!     error = maxval (abs (Mat_transpose * S * Mat - S))
!
! Module Needed:
!   use sim_utils_interface
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

use precision_def
use output_mod

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
  call err_exit
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
    print '(5x, (6f11.5))', (max(min(m2(i, j), 999.0_rp), -999.0_rp), j = 1, n)
  enddo
endif

end function mat_symp_error
