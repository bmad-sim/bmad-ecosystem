!+
! Subroutine rotate_spin_a_step (orbit, field, ele, ds)
!
! Routine to rotate the spin through an integration step.
! Note: It is assumed that the orbit coords are in the element body frame and not the lab frame.
!
! Input:
!   orbit   -- coord_struct: Initial orbit.
!   field   -- em_field_struct: EM Field 
!   ele     -- ele_struct, Element being tracked through. 
!   ds      -- real(rp): Longitudinal step in element body frame.
!
! Output:
!   orbit   -- coord_struct: Orbit with rotated spin
!-

subroutine rotate_spin_a_step (orbit, field, ele, ds)

use equal_mod, dummy_except => rotate_spin_a_step

implicit none

type (ele_struct) ele
type (coord_struct) orbit
type (em_field_struct) field

real(rp) ds, omega(3)
integer sign_z_vel

!

sign_z_vel = ele%orientation * orbit%direction
omega = spin_omega (field, orbit, sign_z_vel)

if (ele%key == sbend$ .or. ele%key == rf_bend$) then
  omega = (1 + ele%value(g$) * orbit%vec(1)) * omega + [0.0_rp, ele%value(g$)*sign_z_vel, 0.0_rp]
endif

call rotate_spin (ds * omega, orbit%spin)

end subroutine rotate_spin_a_step


