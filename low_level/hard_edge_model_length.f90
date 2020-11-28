!+
! Function hard_edge_model_length (ele) result (l_hard)
!
! Input:
!   ele -- ele_struct: Element
!
! Output:
!   l_hard -- real(rp): Length of the hard edge model.
!-

function hard_edge_model_length (ele) result (l_hard)

use attribute_mod, except_dummy => hard_edge_model_length

implicit none

type (ele_struct) ele
real(rp) l_hard

!

if (has_attribute(ele, 'L_HARD_EDGE')) then
  l_hard = ele%value(l_hard_edge$)
else
  l_hard = ele%value(l$)
endif

end function hard_edge_model_length

