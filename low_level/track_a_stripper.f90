!+
! Subroutine track_a_stripper (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through an stripper element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: stripper element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is False.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_stripper (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_stripper
use random_mod
use super_recipes_mod
use spline_mod

implicit none

type (coord_struct) :: orbit, orb0
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)


logical, optional :: make_matrix

character(*), parameter :: r_name = 'track_a_stripper'

!

end subroutine track_a_stripper
