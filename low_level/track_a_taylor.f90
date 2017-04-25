!+
! Subroutine track_a_taylor (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through a taylor element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Taylor element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_taylor (orbit, ele, param, mat6, make_matrix)

use ptc_interface_mod, except_dummy => track_a_taylor

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (taylor_struct) taylor1(6), taylor2(6)

real(rp), optional :: mat6(6,6)

logical, optional :: make_matrix

!

if (ele%orientation * orbit%direction == 1) then
  call track1_taylor (orbit, ele, param, orbit)

else
  call taylor_inverse (ele%taylor, taylor1)
  taylor2 = ele%taylor
  ele%taylor = taylor1
  call track1_taylor (orbit, ele, param, orbit)
  ele%taylor = taylor2
  call kill_taylor(taylor1)
endif

end subroutine
