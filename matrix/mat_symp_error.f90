!+
! Function mat_symp_error (mat, p0_ratio, err_mat) result (error)
!
! Routine to check the symplecticity of a square matrix. The error is
! defined to via:
!     error = maxval (abs (err_mat))
! where
!     err_mat = Mat_transpose * S * Mat - S
!
! When using coordinates like (x, Px/P0, y, Py/P0, z, dP/P0), where P0 is
! the reference momentum, then:
!   p0_ratio = P0(exit coords) / P0(entrance coords). 
! Normally this is 1 except when dealing with something like an accelerating RF cavity.
! When the p0_ratio is not 1, mat is scalled to coordinates p0_ratio is 1 before the
! symplectic error is computed.
!
! Input:
!   mat(:,:)      -- Real(rp): Matrix to check
!   p0_ratio      -- Real(rp), optional: Ratio of p0_exit / p0_entrance. Default is 1.
!                              If present, the size of mat_in must be 6.
!
! Output:
!   error         -- Real(rp): difference from symplecticity:
!                     = 0    --> perfect.
!                     = 1d-4 --> Fair, but I wouldn't use mat for long term tracking.
!                     = 1    --> Terrible.
!   err_mat(:,:)  -- Real(rp), optional: Error matrix. Examination of this matrix can give clues 
!                     as to what part of phase space is contributing most to non-symplectivity.
!-

function mat_symp_error (mat, p0_ratio, err_mat) result (error)

use output_mod, except => mat_symp_error

implicit none

integer i, j, n

real(rp), intent(in) :: mat(:,:)
real(rp), optional :: p0_ratio, err_mat(:,:)
real(rp) :: m1(size(mat, 1), size(mat, 1)), m2(size(mat, 1), size(mat, 1))
real(rp) error

logical :: debug = .false.

!

n = ubound(mat, 1)

if (mod(n, 2) /= 0) then
  print *, 'ERROR IN MAT_SYMP_ERROR: MATRIX DOES NOT HAVE EVEN SIZE'
  if (global_com%exit_on_error) call err_exit
endif

!

m1 = mat_scale_p0 (mat, p0_ratio, .true.)

do j = 2, n, 2                  ! m2 = mat_transpose * S
  m2(1:n,j)   = -m1(j-1,1:n)
  m2(1:n,j-1) =  m1(j,  1:n)
enddo

m2 = matmul (m2, m1)

do j = 2, n, 2
  m2(j, j-1) = m2(j, j-1) - 1
  m2(j-1, j) = m2(j-1, j) + 1
enddo

error = maxval(abs(m2))

if (present(err_mat)) err_mat = m2

if (debug) then
  print *, 'MAT_SYMP_ERROR:', error
  do i = 1, n
    print '(5x, (6f11.5))', (max(min(m2(i, j), 999.0_rp), -999.0_rp), j = 1, n)
  enddo
endif

end function mat_symp_error
