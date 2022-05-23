!+
! Function ele_has_nonzero_kick (ele) result (has_kick)
!
! Function to determine if an element has nonzero kick values.
! Kicks are something hkick$, bl_vkick$, etc.
! See also: zero_ele_kicks, ele_has_nonzero_offset, zero_ele_offsets.
!
! Input
!   ele -- Ele_struct: Element with possible nonzero kicks.
!
! Output:
!   ele -- Ele_struct: Element with no kicks.
!-

function ele_has_nonzero_kick (ele) result (has_kick)

use bmad_struct

implicit none

type (ele_struct) ele
logical has_kick

!

has_kick = .false.

if (has_hkick_attributes(ele%key)) then
  if (ele%value(bl_hkick$) /= 0) has_kick = .true.
  if (ele%value(bl_vkick$) /= 0) has_kick = .true.

elseif (has_kick_attributes(ele%key)) then
  if (ele%value(bl_kick$) /= 0) has_kick = .true.
endif

if (ele%key == lcavity$ .or. ele%key == rfcavity$) then
  if (ele%value(coupler_strength$) /= 0) has_kick = .true.
endif

end function ele_has_nonzero_kick
