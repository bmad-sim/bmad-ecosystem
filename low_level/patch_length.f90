!+
! Function patch_length (patch, ref_coords) result (length)
!
! Routine to calculate the length of a patch element.
! See the Bmad manual for further information.
!
! Input:
!   patch       -- ele_struct: Patch element.
!   ref_coords  -- integer, optional: Reference coords to use. entrance_end$, exit_end$
!                    Default is nint(patch%value(ref_coords$)).
!
! Output:
!   length      -- real(rp): Length of patch.
!-

function patch_length (patch, ref_coords) result (length)

use bmad_routine_interface, dummy => patch_length

implicit none

type (ele_struct) patch
real(rp) length, ww(3,3)
integer, optional :: ref_coords

!

if (is_true(patch%value(user_sets_length$))) then
  length = patch%value(l$)
  return
endif

select case (integer_option(nint(patch%value(ref_coords$)), ref_coords))
case (entrance_end$)
  length = patch%value(z_offset$)

case (exit_end$)
  call floor_angles_to_w_mat (patch%value(x_pitch$), patch%value(y_pitch$), patch%value(tilt$), w_mat_inv = ww)
  length = (ww(3,1) * patch%value(x_offset$) + ww(3,2) * patch%value(y_offset$) + ww(3,3) * patch%value(z_offset$))

case default
  call err_exit  ! Should not be here
end select

end function
