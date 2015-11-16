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
! Subroutine rotate_vec(vec, i, j, angle)
!            rotate_vec(vec, i, j, cos_angle, sin_angle)
!            
! Basic routine to rotate vector components i and j in place by angle in the ij plane.
! If sin_angle is not present, the cosines and sines will be computed. 
! 
! In 3D, i and j are: 2, 3 about x-axis
!                     3, 1 about y-axis
!                     1, 2 about z-axis
!
! Module needed:
!   use rotation_3d_mod
!
! Input:
!   vec(:)             -- real(rp): vector
!   i, j               -- integer: indices to rotate
!   cos_angle or angle -- real(rp): cosine of the angle to rotate
!   sin_angle          -- real(rp): sine of the angle to rotate
!
! Output:
!   vec(:)   -- real(rp): vector 
!-

subroutine rotate_vec(vec, i, j, cos_angle, sin_angle)
implicit none
real(rp) :: vec(:),  cos_angle, temp, ca, sa
real(rp), optional :: sin_angle
integer :: i, j
if (present(sin_angle)) then
  ! cos, sin given
  temp   =  cos_angle*vec(i) - sin_angle*vec(j)
  vec(j) =  sin_angle*vec(i) + cos_angle*vec(j)
  vec(i) = temp
else
  ! angle given
  if (cos_angle == 0) return ! Simple case 
  ca = cos(cos_angle)
  sa = sin(cos_angle) 
  temp   =  ca*vec(i) - sa*vec(j)
  vec(j) =  sa*vec(i) + ca*vec(j)
  vec(i) = temp  
endif
  
end subroutine

!------------------------------------------------------------------------------
!+
! Subroutine rotate_vec_x(vec, angle)
!   wrapper for rotate_vec(vec, 2, 3, angle)
! -
subroutine rotate_vec_x(vec, angle)
real(rp) :: vec(3), angle
call rotate_vec(vec, 2, 3, angle)
end subroutine 

!------------------------------------------------------------------------------
!+
! Subroutine rotate_vec_y(vec, angle)
!   wrapper for rotate_vec(vec, 3, 1, angle)
!-
subroutine rotate_vec_y(vec, angle)
real(rp) :: vec(3), angle
call rotate_vec(vec, 3, 1, angle)
end subroutine 

!------------------------------------------------------------------------------
!+
! Subroutine rotate_vec_y(vec, angle)
!   wrapper for rotate_vec(vec, 1, 2, angle)
!-
subroutine rotate_vec_z(vec, angle)
real(rp) :: vec(3), angle
call rotate_vec(vec, 1, 2, angle)
end subroutine


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine rotate_mat(mat, i, j, angle)
!            rotate_mat(mat, i, j, cos_angle, sin_angle)
!            
! Basic routine to apply a rotation matrix to mat that rotates components i and j 
! in place by angle in the ij plane.
!
! mat = rotation(angle) . mat
!
! If sin_angle is not present, the cosines and sines will be computed. 
! 
! Module needed:
!   use rotation_3d_mod
!
! Input:
!   mat(3,3)           -- real(rp): matrix
!   i, j               -- integer: indices to rotate
!   cos_angle or angle -- real(rp): cosine of the angle to rotate
!   sin_angle          -- real(rp): sine of the angle to rotate
!
! Output:
!   vec(3)            -- real(rp): vector 
!-

subroutine rotate_mat(mat, i, j, cos_angle, sin_angle)
implicit none
real(rp) :: mat(:,:),  cos_angle, temp, ca, sa
real(rp), optional :: sin_angle
integer :: i, j, col

!

if (present(sin_angle)) then
  do col = lbound(mat,2), ubound(mat,2)
    call rotate_vec(mat(:, col), i, j, cos_angle, sin_angle)
  end do 
else
  ! angle only given
  if (cos_angle == 0) return ! Simple case
  ca = cos(cos_angle)
  sa = sin(cos_angle) 
  do col = lbound(mat,2), ubound(mat,2)
    call rotate_vec(mat(:, col), i, j, ca, sa)
  end do 
endif
  
end subroutine


!------------------------------------------------------------------------------
!+
! Subroutine rotate_mat_x(mat, angle)
!   wrapper for rotate_mat(vec, 2, 3, angle)
!-

subroutine rotate_mat_x(mat, angle)
real(rp) :: mat(:,:), angle
call rotate_mat(mat, 2, 3, angle)
end subroutine 

!------------------------------------------------------------------------------
!+
! Subroutine rotate_mat_y(mat, angle)
!   wrapper for rotate_mat(vec, 3, 1, angle)
!-

subroutine rotate_mat_y(mat, angle)
real(rp) :: mat(:,:), angle
call rotate_mat(mat, 3, 1, angle)
end subroutine 

!------------------------------------------------------------------------------
!+
! Subroutine rotate_mat_z(mat, angle)
!   wrapper for rotate_mat(vec, 1, 2, angle)
!-

subroutine rotate_mat_z(mat, angle)
real(rp) :: mat(:,:), angle
call rotate_mat(mat, 1, 2, angle)
end subroutine 

end module
