!+
! Function on_a_girder(ele) result (is_on_girder)
!
! Routine to determine if an element is being supported by a girder element.
!
! Input:
!   ele     -- ele_struct: Element to check.
!
! Output:
!   is_on_girder -- Logical: True if supported. False otherwise.
!- 

function on_a_girder (ele) result (is_on_girder)

use bmad_routine_interface, dummy => on_a_girder

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
logical is_on_girder
integer i

!

is_on_girder = .false.

do i = 1, ele%n_lord
  lord => pointer_to_lord(ele, i)
  if (lord%key /= girder$) cycle
  is_on_girder = .true.
  return
enddo

end function on_a_girder
