module rotation_3d_mod

use sim_utils

implicit none

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine w_mat_to_axis_angle (w_mat, axis, angle)
!
! Routine to find the rotation axis and rotation angle corresponding  to a given
! 3D rotation matrix.
!
! The rotation axis is choisen to have a non-negative dot production with the
! vector (1, 1, 1).
!
! The rotation angle is chosen in the range [-pi, pi].
!
! Module needed:
!   use rotation_3d_mod
!
! Input:
!   w_mat(3,3) -- real(rp): Rotation matrix
!
! Output:
!   axis(3)    -- real(rp): Rotation axis. Normalized to 1.
!   angle      -- real(rp): Rotation angle in the range [-pi, pi].
!-

subroutine w_mat_to_axis_angle (w_mat, axis, angle)

real(rp) w_mat(3,3), axis(3), angle
real(rp) sin_ang, cos_ang

!

axis(1) = w_mat(3,2) - w_mat(2,3)
axis(2) = w_mat(1,3) - w_mat(3,1)
axis(3) = w_mat(2,1) - w_mat(1,2)

sin_ang = norm2(axis) / 2
if (sin_ang == 0) then
  axis = [1, 0, 0]
  angle = 0
  return
endif

axis = axis / (2 * sin_ang)

!

cos_ang = (w_mat(1,1) + w_mat(2,2) + w_mat(3,3) - 1) / 2
angle = atan2(sin_ang, cos_ang)

! Align to axis to point in the general direction of (1,1,1)

if (sum(axis) < 0) then
  axis = -axis 
  angle = -angle
endif

end subroutine w_mat_to_axis_angle

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine axis_angle_to_w_mat (axis, angle, w_mat)
!
! Routine to construct the 3D rotation matrix w_mat given an axis of rotation
! and a rotation angle.
!
! Module needed:
!   use rotation_3d_mod
!
! Input:
!   axis(3)    -- real(rp): Rotation axis. Does not have to be normalized.
!   angle      -- real(rp): Rotation angle in the range [-pi, pi].
!
! Output:
!   w_mat(3,3) -- real(rp): Rotation matrix
!-

subroutine axis_angle_to_w_mat (axis, angle, w_mat)

real(rp) w_mat(3,3), axis(3), angle
real(rp) sin_a, cos_a, norm, x, y, z
character(*), parameter :: r_name = 'axis_angle_to_w_mat'

!

if (angle == 0) then
  call mat_make_unit (w_mat)
  return
endif

!

norm = norm2(axis)
if (norm == 0) then
  w_mat = 0
  call out_io (s_fatal$, r_name, 'ZERO AXIS LENGTH WITH NON-ZERO ROTATION!')
  if (global_com%exit_on_error) call err_exit
  return
endif

x = axis(1) / norm; y = axis(2) / norm; z = axis(3) / norm
cos_a = cos(angle); sin_a = sin(angle)

w_mat(1,1:3) = [x*x + (1 - x*x) * cos_a,       x*y * (1 - cos_a) - z * sin_a, x*z * (1 - cos_a) + y * sin_a]
w_mat(2,1:3) = [x*y * (1 - cos_a) + z * sin_a, y*y + (1 - y*y) * cos_a,       y*z * (1 - cos_a) - x * sin_a]
w_mat(3,1:3) = [x*z * (1 - cos_a) - y * sin_a, y*z * (1 - cos_a) + x * sin_a, z*z + (1 - z*z) * cos_a]

end subroutine axis_angle_to_w_mat

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function quat_to_omega (quat) result (omega)
! 
! Routine to convert rotation from quaternion representation to omega (axis + angle).
!
! Input:
!   quat(0:3)   -- real(rp): Rotation quaternion. Assumed normalized.
!
! Output:
!   omega(3)    -- real(rp): Axis of rotation + magnitude = rotation angle.
!-

function quat_to_omega (quat) result (omega)

real(rp) omega(3), quat(0:3)
real(rp) anorm, angle

!

anorm = norm2(quat(1:3))

if (anorm == 0) then
  omega = 0
else
  angle = 2 * atan2(anorm, quat(0))
  omega = angle * quat(1:3) / anorm
endif

end function quat_to_omega 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function omega_to_quat (omega) result (quat)
! 
! Routine to convert from omega + angle representation to a quaternion.
!
! Input:
!   omega(3)    -- real(rp): Axis of rotation + magnitude = rotation angle.
!
! Output:
!   quat(0:3)   -- real(rp): Rotation quaternion.
!-

function omega_to_quat (omega) result (quat)

real(rp) omega(3), angle, quat(0:3)

!

angle = norm2(omega)

if (angle == 0) then
  quat = [1, 0, 0, 0]
else
  quat = [cos(angle/2), omega * (sin(angle/2) / angle)]
endif

end function omega_to_quat 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine quat_to_axis_angle (quat, axis, angle)
! 
! Routine to convert from axis + angle representation to a quaternion.
!
! Input:
!   quat(0:3)   -- real(rp): Rotation quaternion. Assumed normalized.
!
! Output:
!   axis(3)     -- real(rp): Axis of rotation.
!   angle       -- real(rp): angle of rotation.
!-

subroutine quat_to_axis_angle (quat, axis, angle)

real(rp) axis(3), angle, quat(0:3)
real(rp) anorm

!

anorm = norm2(quat(1:3))

if (anorm == 0) then
  angle = 0
  axis = [0, 0, 1]  ! Arbitrary.
else
  angle = 2 * atan2(anorm, quat(0))
  axis = quat(1:3) / anorm
endif

end subroutine quat_to_axis_angle 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function axis_angle_to_quat (axis, angle) result (quat)
! 
! Routine to convert from axis + angle representation to a quaternion.
!
! Input:
!   axis(3)     -- real(rp): Axis of rotation.
!   angle       -- real(rp): angle of rotation.
!
! Output:
!   quat(0:3)   -- real(rp): Rotation quaternion.
!-

function axis_angle_to_quat (axis, angle) result (quat)

real(rp) axis(3), angle, quat(0:3)

!

quat = [cos(angle/2), axis * sin(angle/2)]

end function axis_angle_to_quat 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function quat_inverse (q_in) result (q_out)
!
! Routine to create the inverse of a quaternian. 
!
! Input:
!   q_in(0:3)     -- real(rp): Quaternion input.
!
! Output:
!   q_out(0:3)    -- real(rp): Inverse quaternion.
!-

function quat_inverse (q_in) result (q_out)

real(rp) q_in(0:3), q_out(0:3)

q_out = [q_in(0), -q_in(1:3)] / dot_product(q_in, q_in)

end function quat_inverse

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function quat_mul (q1, q2) result (q_out)
!
! Routine to multiply two quaternions q_out = q1 * q2.
! Note: q_out = q1 * q2 represents a rotation of q2 first followed by q1.
!
! Input:
!   q1(0:3), q2(0:3)  -- real(rp): Quaternions.
!
! Output:
!   q_out(0:3)        -- real(rp): Resultant q1 * q2
!-

function quat_mul (q1, q2) result (q_out)

real(rp) q1(0:3), q2(0:3), q_out(0:3)
real(rp) a, b, c, d, e, f, g, h

!

A = (q1(0) + q1(1)) * (q2(0) + q2(1))
B = (q1(3) - q1(2)) * (q2(2) - q2(3))
C = (q1(0) - q1(1)) * (q2(2) + q2(3)) 
D = (q1(2) + q1(3)) * (q2(0) - q2(1))
E = (q1(1) + q1(3)) * (q2(1) + q2(2))
F = (q1(1) - q1(3)) * (q2(1) - q2(2))
G = (q1(0) + q1(2)) * (q2(0) - q2(3))
H = (q1(0) - q1(2)) * (q2(0) + q2(3))

q_out(0) = B + (-E - F + G + H)/2
q_out(1) = A -  (E + F + G + H)/2 
q_out(2) = C +  (E - F + G - H)/2 
q_out(3) = D +  (E - F - G + H)/2

end function quat_mul

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine rotate_vec_given_quat (vec, quat)
!
! Routine to rotate a vector using a quaternion..
!
! Input:
!   vec(3)      -- real(rp): Initial vector.
!   q_out(0:3)  -- real(rp): Quaternion to rotate with.
!
! Output:
!   vec(3)    -- real(rp): Final vector.
!-

subroutine rotate_vec_given_quat (vec, quat)

real(rp) :: vec(3), quat(0:3)
real(rp) :: q0_inv

! Use quaternion rotation formula: vec -> q * vec * q^-1

q0_inv = -(quat(1)*vec(1) + quat(2)*vec(2) + quat(3)*vec(3))

vec = [quat(0)*vec(1) + quat(2)*vec(3) - quat(3)*vec(2), &
       quat(0)*vec(2) + quat(3)*vec(1) - quat(1)*vec(3), &
       quat(0)*vec(3) + quat(1)*vec(2) - quat(2)*vec(1)]


vec = [quat(0)*vec(1) + quat(2)*vec(3) - quat(3)*vec(2) - q0_inv*quat(1), &
       quat(0)*vec(2) + quat(3)*vec(1) - quat(1)*vec(3) - q0_inv*quat(2), &
       quat(0)*vec(3) + quat(1)*vec(2) - quat(2)*vec(1) - q0_inv*quat(3)]

end subroutine rotate_vec_given_quat

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine rotate_vec_given_axis_angle (vec, axis, angle)
!
! Routine to rotate a vector.
!
! Input:
!   vec(3)    -- real(rp): Initial vector.
!   axis(3)   -- real(rp): Axis of rotation. Must be normalized to 1.
!   angle     -- real(rp): Angle to rotate by
!
! Output:
!   vec(3)    -- real(rp): Final vector.
!-

subroutine rotate_vec_given_axis_angle (vec, axis, angle)

real(rp) :: vec(:), angle, ca, sa, q0_inv
real(rp) :: axis(:)

! Use quaternion rotation formula: vec -> q * vec * q^-1
! q = cos(angle/2) + sin(angle/2) * (axis(1), axis(2), axis(3))

if (angle == 0) return ! Simple case 

ca = cos(angle/2)
sa = sin(angle/2)

q0_inv = -sa * (axis(1)*vec(1) + axis(2)*vec(2) + axis(3)*vec(3))

vec = [ca*vec(1) + sa * (axis(2)*vec(3) - axis(3)*vec(2)), &
       ca*vec(2) + sa * (axis(3)*vec(1) - axis(1)*vec(3)), &
       ca*vec(3) + sa * (axis(1)*vec(2) - axis(2)*vec(1))]


vec = [ca*vec(1) + sa * (axis(2)*vec(3) - axis(3)*vec(2) - q0_inv*axis(1)), &
       ca*vec(2) + sa * (axis(3)*vec(1) - axis(1)*vec(3) - q0_inv*axis(2)), &
       ca*vec(3) + sa * (axis(1)*vec(2) - axis(2)*vec(1) - q0_inv*axis(3))]

end subroutine rotate_vec_given_axis_angle

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine rotate_vec (vec, axis, angle)
!            
! Basic routine to rotate vector components around the x, y, or z axis.
! 
! Module needed:
!   use rotation_3d_mod
!
! Input:
!   vec(3)       -- real(rp): vector
!   axis         -- integer: x_axis$, y_axis$, or z_axis$
!   angle        -- real(rp): cosine of the angle to rotate
!
! Output:
!   vec(3)       -- real(rp): Rotated vector.
!-

subroutine rotate_vec (vec, axis, angle)

real(rp) :: vec(:), angle, ca, sa
integer :: axis

!

if (angle == 0) return ! Simple case 

ca = cos(angle)
sa = sin(angle) 

select case (axis)
case (x_axis$)
  vec(2:3) =   [ca*vec(2) - sa*vec(3), sa*vec(2) + ca*vec(3)]

case (y_axis$)
  vec(1:3:2) = [sa*vec(3) + ca*vec(1), ca*vec(3) - sa*vec(1)]

case (z_axis$)
  vec(1:2) =   [ca*vec(1) - sa*vec(2), sa*vec(1) + ca*vec(2)]
end select
  
end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine rotate_mat (mat, axis, angle)
!            
! Basic routine to apply a rotation matrix to mat that rotates around the x, y, or z axis.
! 
! Module needed:
!   use rotation_3d_mod
!
! Input:
!   mat(3,3)        -- real(rp): matrix
!   axis            -- integer: x_axis$, y_axis$, or z_axis$
!   angle           -- real(rp): cosine of the angle to rotate
!
! Output:
!   mat(3,3)        -- real(rp): Rotated matrix
!-

subroutine rotate_mat(mat, axis, angle)

real(rp) :: mat(:,:),  angle, temp, ca, sa
integer :: axis

!

if (angle == 0) return ! Simple case

ca = cos(angle)
sa = sin(angle) 

select case (axis)
case (x_axis$)
  mat(2:3, 1) =   [ca*mat(2,1) - sa*mat(3,1), sa*mat(2,1) + ca*mat(3,1)]
  mat(2:3, 2) =   [ca*mat(2,2) - sa*mat(3,2), sa*mat(2,2) + ca*mat(3,2)]
  mat(2:3, 3) =   [ca*mat(2,3) - sa*mat(3,3), sa*mat(2,3) + ca*mat(3,3)]

case (y_axis$)
  mat(1:3:2, 1) = [sa*mat(3,1) + ca*mat(1,1), ca*mat(3,1) - sa*mat(1,1)]
  mat(1:3:2, 2) = [sa*mat(3,2) + ca*mat(1,2), ca*mat(3,2) - sa*mat(1,2)]
  mat(1:3:2, 3) = [sa*mat(3,3) + ca*mat(1,3), ca*mat(3,3) - sa*mat(1,3)]

case (z_axis$)
  mat(1:2, 1) =   [ca*mat(1,1) - sa*mat(2,1), sa*mat(1,1) + ca*mat(2,1)]
  mat(1:2, 2) =   [ca*mat(1,2) - sa*mat(2,2), sa*mat(1,2) + ca*mat(2,2)]
  mat(1:2, 3) =   [ca*mat(1,3) - sa*mat(2,3), sa*mat(1,3) + ca*mat(2,3)]
end select
  
end subroutine

end module
