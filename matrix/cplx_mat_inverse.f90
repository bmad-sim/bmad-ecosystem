#include "CESR_platform.inc"

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine cplx_mat_inverse (mat_r, mat_i, inv_r, inv_i, ok, print_error)
!
! Computes the inverse of a complex matrix by computing inverse using mat_inverse, which
! computes the inverse of real matrixes.
! Based on : "Method to Calculate the Inverse of a Complex Matrix using Real Matrix Inversion"
!            by Andreas Falkenberg
!
! Input:
!   mat_r(:,:)     -- Real(rp): Real part of input matrix
!   mat_i(:,:)     -- Real(rp): Imaginary part of input matrix
!   print_err    -- Logical, optional: If True then the subroutine will type out
!                         a warning message. Default is False.
! Output:
!   inv_r(:,:)     -- Real(rp): Real part of inverted matrix
!   inv_i(:,:)     -- Real(rp): Imaginary part of inverted matrix
!   ok           -- Logical, optional: Set False for a singular input matrix.
!-
subroutine cplx_mat_inverse(mat_r, mat_i, inv_r, inv_i, ok, print_err)

use precision_def
use sim_utils_interface

implicit none

real(rp) :: mat_r(:,:)
real(rp) :: mat_i(:,:)
real(rp) :: inv_r(:,:)
real(rp) :: inv_i(:,:)

real(rp) :: r0(1:size(mat_r,1),1:size(mat_r,1))
real(rp) :: a_inv(1:size(mat_r,1),1:size(mat_r,1))
real(rp) :: y00(1:size(mat_r,1),1:size(mat_r,1))

logical, optional :: ok, print_err
character(16) :: r_name = 'cplx_mat_inverse'

call mat_inverse(mat_r, a_inv, ok, print_err)
r0 = matmul(a_inv,mat_i)
call mat_inverse (matmul(mat_i,r0)+mat_r, inv_r, ok, print_err)
inv_i = -1.0 * matmul(r0,inv_r)

end subroutine cplx_mat_inverse


