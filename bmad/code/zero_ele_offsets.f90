!+
! Subroutine zero_ele_offsets (ele)
!
! Subroutine to zero the offsets, pitches, tilt and ref_tilt of an element.
! Also see: ele_has_nonzero_offset, zero_ele_kicks, ele_has_nonzero_kick
!
! Input
!   ele -- Ele_struct: Element with possible nonzero offsets, etc.
!
! Output:
!   ele -- Ele_struct: Element with no (mis)orientation.
!-

subroutine zero_ele_offsets (ele)

use attribute_mod, dummy => zero_ele_offsets

implicit none

type (ele_struct) ele

!

if (.not. has_orientation_attributes(ele)) return

select case (ele%key)
case (sbend$, rf_bend$)
  ele%value(roll$) = 0
  ele%value(roll_tot$) = 0
  ele%value(ref_tilt$) = 0
  ele%value(ref_tilt_tot$) = 0
case (mirror$, multilayer_mirror$, crystal$)
  ele%value(tilt$) = 0
  ele%value(tilt_tot$) = 0
  ele%value(ref_tilt$) = 0
  ele%value(ref_tilt_tot$) = 0
case default
  ele%value(tilt$) = 0
  ele%value(tilt_tot$) = 0
end select

ele%value(x_pitch$) = 0
ele%value(y_pitch$) = 0
ele%value(x_offset$) = 0
ele%value(y_offset$) = 0
ele%value(z_offset$) = 0

ele%value(x_pitch_tot$) = 0
ele%value(y_pitch_tot$) = 0
ele%value(x_offset_tot$) = 0
ele%value(y_offset_tot$) = 0
ele%value(z_offset_tot$) = 0

end subroutine zero_ele_offsets
