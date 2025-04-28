!+
! Subroutine rotate_spin (rot_vec, spin, qrot)
!
! Routine to rotate a spin.
!
! Input:
!   rot_vec(3)  -- real(rp): Rotation axis. Magnitude of rot_vec is the rotation angle.
!   spin(3)     -- real(rp): Initial coords.
!
! Output:
!   spin(3)     -- real(rp): Final coords.
!   qrot(0:3)   -- real(rp), optional :: rotation quaternion.
!-

subroutine rotate_spin (rot_vec, spin, qrot)

use bmad_routine_interface, dummy_except => rotate_spin

implicit none

real(rp), optional :: qrot(0:3)
real(rp) :: spin(3), rot_vec(3), axis(3), angle

!

angle = norm2(rot_vec)
if (angle == 0) return

axis = rot_vec / angle

spin = rotate_vec_given_axis_angle (spin, axis, angle)

if (present(qrot)) qrot = axis_angle_to_quat(axis, angle)

end subroutine rotate_spin


