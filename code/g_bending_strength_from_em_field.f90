!+
! Subroutine g_bending_strength_from_em_field (ele, param, s_rel, orbit, local_ref_frame, g, dg)
!
! Subroutine to calculate the g bending strength (1/bend_radius) felt by a particle in a element.
! 
! Example: A simple bend has g = (1/rho, 0, 0) where rho is the design bending radius.
!
! Input:
!   ele              -- Ele_struct: Element being tracked thorugh.
!   param            -- lat_param_struct: Lattice parameters.
!   s_rel            -- Real(rp): Distance from the start of the element to the particle.
!   orbit            -- Coord_struct: Particle position in lab (not element) frame.
!   local_ref_frame  -- Logical, If True then take the input coordinates and output g
!                         as being with respect to the frame of referene of the element (ignore misalignments).
!
! Output:
!   g(3)              -- Real(rp): g = (g_x, g_y, g_s) bending strength vector (|g| = 1/bend_radius).
!   dg(3,3)           -- Real(rp), optional: dg(:)/dr gradient. 
!                         Takes into account dg_x/dx in a bend due to curvilinear coords.
!-

subroutine g_bending_strength_from_em_field (ele, param, s_rel, orbit, local_ref_frame, g, dg)

use bmad_routine_interface, dummy => g_bending_strength_from_em_field

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (em_field_struct) field
type (coord_struct) orbit

real(rp), intent(in) :: s_rel
real(rp), intent(out) :: g(3)
real(rp), optional :: dg(3,3)
real(rp) vel_unit(3), fact
real(rp) f

logical local_ref_frame

! calculate the field

call em_field_calc (ele, param, s_rel, orbit, local_ref_frame, field, present(dg))

! vel_unit is the velocity normalized to unit length

vel_unit(1:2) = [orbit%vec(2), orbit%vec(4)] / (1 + orbit%vec(6))
vel_unit(3) = sqrt(1 - vel_unit(1)**2 - vel_unit(2)**2)
fact = 1 / (orbit%beta * ele%value(p0c$) * (1 + orbit%vec(6)))
g = g_from_field (field%B, field%E, orbit, vel_unit, param, fact)

! Derivative

if (present(dg)) then
  dg(:,1) = g_from_field (field%dB(:,1), field%dE(:,1), orbit, vel_unit, param, fact)
  dg(:,2) = g_from_field (field%dB(:,2), field%dE(:,2), orbit, vel_unit, param, fact)
  dg(:,3) = g_from_field (field%dB(:,3), field%dE(:,3), orbit, vel_unit, param, fact)
endif

!---------------------------------------------------------------
contains

function g_from_field (B, E, orbit, vel_unit, param, fact) result (g_bend)

type (coord_struct) orbit
type (lat_param_struct) param

real(rp) B(3), E(3), g_bend(3)
real(rp) force(3), force_perp(3)
real(rp) vel_unit(3), fact

! force_perp is the perpendicular component of the force.

force = (E + cross_product(vel_unit, B) * orbit%beta * c_light) * charge_of(param%particle)
force_perp = force - vel_unit * (dot_product(force, vel_unit))
g_bend = -force_perp * fact

end function g_from_field

end subroutine g_bending_strength_from_em_field 

