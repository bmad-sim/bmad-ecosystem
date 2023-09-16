!+
! Subroutine track_a_gkicker (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through a gkicker element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Gkicker
!   param       -- lat_param_struct: Lattice parameters.
!   mat6(6,6)   -- Real(rp), optional: Transfer matrix before the element.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_gkicker (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_gkicker

implicit none

type (coord_struct) :: orbit, start_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) zb_old
logical, optional :: make_matrix

!

zb_old = orbit%vec(5) / (c_light * orbit%beta)
orbit%vec = orbit%vec + ele%value(x_kick$:pz_kick$)

call convert_pc_to(orbit%p0c*(1.0_rp+orbit%vec(6)), orbit%species, beta = orbit%beta)
orbit%t = orbit%t + zb_old - orbit%vec(5) / (c_light * orbit%beta)

if (logic_option(.false., make_matrix)) call mat_make_unit(mat6)

end subroutine
