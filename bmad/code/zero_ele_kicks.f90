!+
! Subroutine zero_ele_kicks (ele)
!
! Subroutine to zero any kick attributes like hkick$, bl_vkick$, etc.
! See also: ele_has_nonzero_kick, ele_has_nonzero_offset, zero_ele_offsets.
!
! Input
!   ele -- Ele_struct: Element with possible nonzero kicks.
!
! Output:
!   ele -- Ele_struct: Element with no kicks.
!-

subroutine zero_ele_kicks (ele)

use bmad_struct
implicit none

type (ele_struct) ele

!

if (has_hkick_attributes(ele%key)) then
  ele%value(hkick$) = 0
  ele%value(vkick$) = 0
  ele%value(bl_hkick$) = 0
  ele%value(bl_vkick$) = 0

elseif (has_kick_attributes(ele%key)) then
  ele%value(kick$) = 0
  ele%value(bl_kick$) = 0
endif

if (ele%key == lcavity$ .or. ele%key == rfcavity$) then
  ele%value(coupler_strength$) = 0
endif

end subroutine zero_ele_kicks
