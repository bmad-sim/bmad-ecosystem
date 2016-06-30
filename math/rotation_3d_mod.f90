module rotation_3d_mod

use sim_utils

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

implicit none

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

implicit none

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

implicit none

real(rp) :: vec(:), angle, ca, sa, q0
real(rp) :: axis(:)

! Use quaternion rotation formula: vec -> q * vec * q^-1
! q = cos(angle/2) + sin(angle/2) * (axis(1), axis(2), axis(3))

if (angle == 0) return ! Simple case 

ca = cos(angle/2)
sa = sin(angle/2)

q0 = -sa * (axis(1)*vec(1) + axis(2)*vec(2) + axis(3)*vec(3))
vec = [ca*vec(1) + sa * (axis(2)*vec(3) - axis(3)*vec(2)), &
       ca*vec(2) + sa * (axis(3)*vec(1) - axis(1)*vec(3)), &
       ca*vec(3) + sa * (axis(1)*vec(2) - axis(2)*vec(1))]


vec = [ca*vec(1) + sa * (axis(2)*vec(3) - axis(3)*vec(2) - q0*axis(1)), &
       ca*vec(2) + sa * (axis(3)*vec(1) - axis(1)*vec(3) - q0*axis(2)), &
       ca*vec(3) + sa * (axis(1)*vec(2) - axis(2)*vec(1) - q0*axis(3))]

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

implicit none

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

implicit none

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
