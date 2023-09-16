!+
! Subroutine mat_scale_p0 (mat_in, p0_ratio, invert) result (mat_out)
!
! Routine to scale a matrix so that the ratio between the input coordinates
! and the output coordinates is changed by p0_ratio
!
! When using coordinates like (x, Px/P0, y, Py/P0, z, dP/P0), where P0 is 
! the reference momentum, then the ratio
!   P0(exit coords) / P0(entrance coords). 
! is Normally 1 except when dealing with something like an accelerating RF cavity.
! When this ratio is not 1, a matrix must be scalled using this routine before 
! symplectic manipulations are done. 
!
! Note: This routine scales the exit P0 keeping the entrance P0 fixed.
!
! Input:
!   mat_in(:,:) -- Real(rp): Matrix
!   p0_ratio    -- Real(rp): optional: p0_exit / p0_entrance. is scaled by p0_ratio.
!                            If not present or if p0_ratio = 1 then mat_out = mat_in.
!                            If present, the size of mat_in must be 6.
!   invert      -- Logical, optional: If present and True then scale by 1/p0_ratio.
!
! Output:
!   mat_out(:,:) -- Real(rp): Scalled matrix.
!-

function mat_scale_p0 (mat_in, p0_ratio, invert) result (mat_out)

use output_mod, except => mat_scale_p0

implicit none

real(rp), intent(in) :: mat_in(:,:)
real(rp), optional :: p0_ratio
real(rp) :: mat_out(size(mat_in, 1), size(mat_in, 2)), r_p0
logical, optional :: invert

!

mat_out = mat_in

if (.not. present(p0_ratio)) return
if (p0_ratio == 1) return

if (size(mat_in) /= 36) then
  print *, 'ERROR IN MAT_SCALE_P0: MATRIX SIZE NOT 6X6'
  if (global_com%exit_on_error) call err_exit
  mat_out = 0
endif

!

if (logic_option(.false., invert)) then
  r_p0 = 1 / p0_ratio
else
  r_p0 = p0_ratio
endif

mat_out(2,:) = mat_out(2,:) / r_p0
mat_out(4,:) = mat_out(4,:) / r_p0
mat_out(6,:) = mat_out(6,:) / r_p0

! mat_out(:,2) = mat_out(:,2) * r_p0
! mat_out(:,4) = mat_out(:,4) * r_p0
! mat_out(:,6) = mat_out(:,6) * r_p0

end function mat_scale_p0
