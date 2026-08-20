!+
! Function aperture_at_is_wall_transition (ele) result (is_wall_transition)
!
! Routine to determine if a particle crossing the aperture of an element is crossing a wall
! transition. That is, is ele%aperture_at = wall_transition$.
!
! Super_slave elements with multiple lords have ele%aperture_at = lord_defined$ and in this case
! the lord(s) must be queried. Note: If the slave has multiple lords with apertures, True is
! returned if any one of the lords has a wall transition aperture.
!
! Input:
!   ele   -- ele_struct: Element under consideration.
!
! Output:
!   is_wall_transition  -- logical: True if the element, or a super_lord of the element, has
!                            aperture_at = wall_transition$.
!-

function aperture_at_is_wall_transition (ele) result (is_wall_transition)

use bmad_routine_interface, dummy => aperture_at_is_wall_transition

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord

integer i
logical is_wall_transition

!

if (ele%aperture_at /= lord_defined$) then
  is_wall_transition = (ele%aperture_at == wall_transition$)
  return
endif

is_wall_transition = .true.
do i = 1, ele%n_lord
  lord => pointer_to_lord(ele, i)
  if (lord%lord_status /= super_lord$) cycle
  if (lord%aperture_at == wall_transition$) return
enddo

is_wall_transition = .false.

end function aperture_at_is_wall_transition
