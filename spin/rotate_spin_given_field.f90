!+
! Subroutine rotate_spin_given_field (orbit, sign_z_vel, BL, EL, qrot)
!
! Routine to rotate a spin given the integrated magnetic and/or electric field strengths.
!
! Integrated field is the field * length which is independent of particle direction of travel.
!
! Input:
!   orbit       -- coord_struct: Initial orbit.
!   sign_z_vel  -- integer: +/- 1. Sign of direction of travel relative to the element.
!   BL(3)       -- real(rp), optional: Integrated field strength. Assumed zero if not present.
!   EL(3)       -- real(rp), optional: Integrated field strength. Assumed zero if not present.
!   qrot(0:3)   -- real(rp), optional: Initial rotation quaternion.
!
! Output:
!   orbit       -- coord_struct: Orbit with rotated spin
!   qrot(0:3)   -- real(rp), optional: Rotation quaternion with rotation due to the field added in.
!-

subroutine rotate_spin_given_field (orbit, sign_z_vel, BL, EL, qrot)

use equal_mod, dummy_except => rotate_spin_given_field

implicit none

type (coord_struct) orbit
type (em_field_struct) field

real(rp), optional :: BL(3), EL(3), qrot(0:3)
real(rp)  omega(3)

integer sign_z_vel

!

if (present(BL)) then
  field%B = orbit%time_dir * BL
else
  field%B = 0
endif

if (present(EL)) then
  field%E = orbit%time_dir * EL
else
  field%E = 0
endif

omega = spin_omega (field, orbit, sign_z_vel)
call rotate_spin (omega, orbit%spin)

if (present(qrot)) qrot = quat_mul(omega_to_quat(omega), qrot)

end subroutine rotate_spin_given_field

