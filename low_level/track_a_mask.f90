!+
! Subroutine track_a_mask (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through a mask element. 
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: Mask element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_mask (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_mask

implicit none

type (coord_struct) :: orbit, temp_orb
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)

integer ix_sec

logical, optional :: make_matrix

! If the plate/mask is turned off then all particles are simply transmitted through.

if (.not. ele%is_on) return

! Particle is lost if in an opaque section

temp_orb = orbit
call offset_particle (ele, set$, temp_orb, set_hvkicks = .false.)

ix_sec = diffraction_plate_or_mask_hit_spot (ele, temp_orb)
if (ix_sec == 0) orbit%state = lost$

end subroutine
