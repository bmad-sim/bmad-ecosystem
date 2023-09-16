!+
! Function ele_has_nonzero_offset (ele) result (has_offset)
!
! Function to tell if an element has a non-zero offset, pitch or tilt.
! Also see: zero_ele_offsets, zero_ele_kicks, ele_has_nonzero_kick
!
! Input
!   ele -- Ele_struct: Element with possible nonzero offsets.
!
! Output:
!   has_offset -- Logical: Set true is element has a non-zero offset.
!-

function ele_has_nonzero_offset (ele) result (has_offset)

use bmad_interface, dummy => ele_has_nonzero_offset
implicit none

type (ele_struct) ele
logical has_offset

!

has_offset = .false.
if (.not. has_orientation_attributes(ele)) return

select case (ele%key)
case (sbend$, rf_bend$)
  if (ele%value(roll_tot$) /= 0) has_offset = .true.
  if (ele%value(ref_tilt_tot$) /= 0) has_offset = .true.
case (mirror$, multilayer_mirror$, crystal$)
  if (ele%value(tilt_tot$) /= 0) has_offset = .true.
  if (ele%value(ref_tilt_tot$) /= 0) has_offset = .true.
case default
  if (ele%value(tilt_tot$) /= 0) has_offset = .true.
end select

if (ele%value(x_pitch_tot$) /= 0) has_offset = .true.
if (ele%value(y_pitch_tot$) /= 0) has_offset = .true.
if (ele%value(x_offset_tot$) /= 0) has_offset = .true.
if (ele%value(y_offset_tot$) /= 0) has_offset = .true.
if (ele%value(z_offset_tot$) /= 0) has_offset = .true.

end function ele_has_nonzero_offset

