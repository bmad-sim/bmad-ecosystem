!+
! Function num_field_eles (ele) result (n_field_ele)
!
! Routine to return the number of field elements associated with a given element.
! n_field_ele will be zero for groups, overlays, and girders.
!
! See the routine pointer_to_field_ele for more info.
!
! Input:
!   ele       -- ele_struct: Element with sum number of associated field elements.
!
! Output:
!   n_field_ele -- integer: Number of associated field elements.
!-

function num_field_eles (ele) result (n_field_ele)

use bmad_routine_interface, dummy => num_field_eles

implicit none

type (ele_struct) ele
type (ele_struct), pointer :: f_ele
integer n_field_ele

!

n_field_ele = 0

do
  f_ele => pointer_to_field_ele(ele, n_field_ele+1)
  if (.not. associated(f_ele)) exit
  n_field_ele = n_field_ele + 1
enddo

end function num_field_eles
