!+
! Subroutine track_a_pickup (orbit, ele, param, err_flag, mat6, make_matrix)
!
! Bmad_standard tracking through a pickup element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Pickup element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_pickup (orbit, ele, param, err_flag, mat6, make_matrix)

use bmad_routine_interface, except_dummy => track_a_pickup

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)

logical, optional :: make_matrix
logical, optional :: err_flag

!

end subroutine
