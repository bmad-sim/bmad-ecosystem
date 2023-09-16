!+
! Function tao_pointer_to_building_wall_shape (wall_name) result (e_shape)
!
! Routine to return the shape associated with a building_wall section.
!
! Input:
!   wall_name   -- character(*): Name of the wall.
!
! Output:
!   e_shape     -- tao_ele_shape_struct, pointer: Associated shape. 
!                       Nullified if there is no associated shape.
!-

function tao_pointer_to_building_wall_shape (wall_name) result (e_shape)

use tao_interface, dummy => tao_pointer_to_building_wall_shape

implicit none

type (tao_ele_shape_struct), pointer :: e_shape
character(*) wall_name
integer is

!

do is = 1, size(s%plot_page%floor_plan%ele_shape)
  e_shape => s%plot_page%floor_plan%ele_shape(is)
  if (e_shape%ele_id(1:15) /= 'building_wall::') cycle
  if (.not. match_wild(wall_name, e_shape%ele_id(16:))) cycle
  if (.not. e_shape%draw) cycle
  return
enddo

nullify(e_shape)

end function tao_pointer_to_building_wall_shape 
