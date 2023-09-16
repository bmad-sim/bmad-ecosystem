module rotation_3d_mod

use sim_utils_interface
use output_mod
use modulo2_mod

implicit none

interface quat_conj
  module procedure quat_conj_real
  module procedure quat_conj_complex
end interface

interface quat_mul
  module procedure quat_mul_real
  module procedure quat_mul_complex
end interface

interface quat_rotate
  module procedure quat_rotate_real
  module procedure quat_rotate_complex
end interface

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine w_mat_to_axis_angle (w_mat, axis, angle, ref_axis)
!
! Routine to find the rotation axis and rotation angle corresponding to a given
! 3D rotation matrix.
!
! The rotation axis is chosen to have a non-negative dot product with the
! reference axis ref_axis.
!
! The rotation angle is chosen in the range [-pi, pi].
!
! Input:
!   w_mat(3,3)    -- real(rp): Rotation matrix.
!   ref_axis(3)   -- real(rp), optional: Reference axis. Default is [1, 1, 1].
!
! Output:
!   axis(3)       -- real(rp): Rotation axis. Normalized to 1.
!   angle         -- real(rp): Rotation angle in the range [-pi, pi].
!-

subroutine w_mat_to_axis_angle (w_mat, axis, angle, ref_axis)

real(rp) w_mat(3,3), axis(3), angle
real(rp), optional :: ref_axis(3)
real(rp) sin2_ang, quat(0:3)

!

quat = w_mat_to_quat(w_mat)
if (all(quat(1:3) == 0)) then
  axis = [1, 0, 0]
  angle = 0
  return
endif

!

sin2_ang = norm2(quat(1:3))
axis = quat(1:3) / sin2_ang
angle = modulo2(2 * atan2(sin2_ang, quat(0)), pi)

! Align to axis to point in the general direction of (1,1,1)

if (present(ref_axis)) then
  if (dot_product(ref_axis, axis) < 0) then
    axis = -axis 
    angle = -angle
  endif
elseif (sum(axis) < 0) then
  axis = -axis 
  angle = -angle
endif

end subroutine w_mat_to_axis_angle

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function w_mat_to_quat (w_mat) result (quat)
!
! Routine to find the quaternion corresponding to a given 3D rotation matrix.
!
! Input:
!   w_mat(3,3) -- real(rp): Rotation matrix
!
! Output:
!   quat(0:3)  -- real(rp): Quaternion.
!-

function w_mat_to_quat (w_mat) result (quat)

real(rp) w_mat(3,3), quat(0:3)
real(rp) trace, s

!

trace = w_mat(1,1) + w_mat(2,2) + w_mat(3,3)

if (trace > 0) then
  s = 0.5_rp / sqrt(trace + 1)
  quat(0) = 1 / (4 * s)
  quat(1) = (w_mat(3,2) - w_mat(2,3)) * s
  quat(2) = (w_mat(1,3) - w_mat(3,1)) * s
  quat(3) = (w_mat(2,1) - w_mat(1,2)) * s

elseif (w_mat(1,1) > w_mat(2,2) .and. w_mat(1,1) > w_mat(3,3)) then
  s = 2 * sqrt(1 + w_mat(1,1) - w_mat(2,2) - w_mat(3,3))
  quat(0) = (w_mat(3,2) - w_mat(2,3)) / s
  quat(1) = 0.25_rp * s
  quat(2) = (w_mat(1,2) + w_mat(2,1)) / s
  quat(3) = (w_mat(1,3) + w_mat(3,1)) / s

elseif (w_mat(2,2) > w_mat(3,3)) then
  s = 2 * sqrt(1 + w_mat(2,2) - w_mat(1,1) - w_mat(3,3))
  quat(0) = (w_mat(1,3) - w_mat(3,1)) / s
  quat(1) = (w_mat(1,2) + w_mat(2,1)) / s
  quat(2) = 0.25_rp * s
  quat(3) = (w_mat(2,3) + w_mat(3,2)) / s

else
  s = 2 * sqrt(1 + w_mat(3,3) - w_mat(1,1) - w_mat(2,2))
  quat(0) = (w_mat(2,1) - w_mat(1,2)) / s
  quat(1) = (w_mat(1,3) + w_mat(3,1)) / s
  quat(2) = (w_mat(2,3) + w_mat(3,2)) / s
  quat(3) = 0.25_rp * s
endif

end function w_mat_to_quat

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function quat_to_w_mat (quat) result (w_mat)
!
! Routine to construct the 3D rotation matrix w_mat given a rotation quaternion
!
! Input:
!   quat(0:3)  -- real(rp): Quaternion.
!
! Output:
!   w_mat(3,3) -- real(rp): Rotation matrix
!-

function quat_to_w_mat (quat) result (w_mat)

real(rp) quat(0:3), w_mat(3,3)
real(rp) sq1, sqx, sqy, sqz, invs, tmp1, tmp2

sq1 = quat(0) * quat(0)
sqx = quat(1) * quat(1)
sqy = quat(2) * quat(2)
sqz = quat(3) * quat(3)

! invs (inverse square length) is only required if quaternion is not already normalised

invs = 1 / (sqx + sqy + sqz + sq1)
w_mat(1,1) = ( sqx - sqy - sqz + sq1) * invs ! since sq1 + sqx + sqy + sqz =1/invs * invs
w_mat(2,2) = (-sqx + sqy - sqz + sq1) * invs
w_mat(3,3) = (-sqx - sqy + sqz + sq1) * invs

tmp1 = quat(1) * quat(2)
tmp2 = quat(3) * quat(0)
w_mat(2,1) = 2.0_rp * (tmp1 + tmp2) * invs
w_mat(1,2) = 2.0_rp * (tmp1 - tmp2) * invs

tmp1 = quat(1) * quat(3)
tmp2 = quat(2) * quat(0)
w_mat(3,1) = 2.0_rp * (tmp1 - tmp2) * invs
w_mat(1,3) = 2.0_rp * (tmp1 + tmp2) * invs
tmp1 = quat(2) * quat(3)
tmp2 = quat(1) * quat(0)
w_mat(3,2) = 2.0_rp * (tmp1 + tmp2) * invs
w_mat(2,3) = 2.0_rp * (tmp1 - tmp2) * invs    

end function quat_to_w_mat

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine axis_angle_to_w_mat (axis, angle, w_mat)
!
! Routine to construct the 3D rotation matrix w_mat given an axis of rotation
! and a rotation angle.
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
  quat = [cos(angle*0.5_rp), omega * (sin(angle*0.5_rp) / angle)]
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

quat = [cos(angle*0.5_rp), axis * sin(angle*0.5_rp)]

end function axis_angle_to_quat 

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function quat_conj_real (q_in) result (q_out)
!
! Routine to create the conjugate of a quaternian. 
! Overloaded by quat_conj.
!
! Input:
!   q_in(0:3)     -- real(rp): Quaternion input.
!
! Output:
!   q_out(0:3)    -- real(rp): Conjugate quaternion.
!-

function quat_conj_real (q_in) result (q_out)

real(rp) q_in(0:3), q_out(0:3)

q_out = [q_in(0), -q_in(1:3)]

end function quat_conj_real

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function quat_conj_complex (q_in) result (q_out)
!
! Routine to create the conjugate of a quaternian. 
! Overloaded by quat_conj.
!
! Input:
!   q_in(0:3)     -- real(rp): Quaternion input.
!
! Output:
!   q_out(0:3)    -- real(rp): Conjugate quaternion.
!-

function quat_conj_complex (q_in) result (q_out)

complex(rp) q_in(0:3), q_out(0:3)

q_out = [q_in(0), -q_in(1:3)]

end function quat_conj_complex

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
! Function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)
!
! Routine to multiply quaternions q_out = q1 * q2 * q3 * q4 * ...
! Overloaded by quat_mul.
! Note: q_out = q1 * q2 represents a rotation of q2 first followed by q1.
!
! Input:
!   q1(0:3), q2(0:3)        -- real(rp): Quaternions.
!   q3(0:3), ..., q9(0:3)   -- real(rp), optional: More quaternions.
!
! Output:
!   q_out(0:3)              -- real(rp): Resultant q1 * q2
!-

function quat_mul_real (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)

real(rp) q1(0:3), q2(0:3), q_out(0:3)
real(rp), optional :: q3(0:3), q4(0:3), q5(0:3), q6(0:3), q7(0:3), q8(0:3), q9(0:3)

!

q_out = q_mul(q1, q2)
if (.not. present(q3)) return
q_out = q_mul(q_out, q3)
if (.not. present(q4)) return
q_out = q_mul(q_out, q4)
if (.not. present(q5)) return
q_out = q_mul(q_out, q5)
if (.not. present(q6)) return
q_out = q_mul(q_out, q6)
if (.not. present(q7)) return
q_out = q_mul(q_out, q7)
if (.not. present(q8)) return
q_out = q_mul(q_out, q8)
if (.not. present(q9)) return
q_out = q_mul(q_out, q9)

!----------------------------------------------------
contains

function q_mul(q1, q2) result (q_out)

real(rp) q1(0:3), q2(0:3), q_out(0:3)
real(rp) a, b, c, d, e, f, g, h

A = (q1(0) + q1(1)) * (q2(0) + q2(1))
B = (q1(3) - q1(2)) * (q2(2) - q2(3))
C = (q1(0) - q1(1)) * (q2(2) + q2(3)) 
D = (q1(2) + q1(3)) * (q2(0) - q2(1))
E = (q1(1) + q1(3)) * (q2(1) + q2(2))
F = (q1(1) - q1(3)) * (q2(1) - q2(2))
G = (q1(0) + q1(2)) * (q2(0) - q2(3))
H = (q1(0) - q1(2)) * (q2(0) + q2(3))

q_out(0) = B + (-E - F + G + H)*0.5_rp
q_out(1) = A -  (E + F + G + H)*0.5_rp 
q_out(2) = C +  (E - F + G - H)*0.5_rp 
q_out(3) = D +  (E - F - G + H)*0.5_rp

end function q_mul

end function quat_mul_real

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function quat_mul_complex (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)
!
! Routine to multiply quaternions q_out = q1 * q2 * q3 * q4 * ...
! Overloaded by quat_mul.
! Note: q_out = q1 * q2 represents a rotation of q2 first followed by q1.
!
! Input:
!   q1(0:3), q2(0:3)        -- complex(rp): Quaternions.
!   q3(0:3), ..., q9(0:3)   -- complex(rp), optional: More quaternions.
!
! Output:
!   q_out(0:3)              -- complex(rp): Resultant q1 * q2
!-

function quat_mul_complex (q1, q2, q3, q4, q5, q6, q7, q8, q9) result (q_out)

complex(rp) q1(0:3), q2(0:3), q_out(0:3)
complex(rp), optional :: q3(0:3), q4(0:3), q5(0:3), q6(0:3), q7(0:3), q8(0:3), q9(0:3)

!

q_out = q_mul(q1, q2)
if (.not. present(q3)) return
q_out = q_mul(q_out, q3)
if (.not. present(q4)) return
q_out = q_mul(q_out, q4)
if (.not. present(q5)) return
q_out = q_mul(q_out, q5)
if (.not. present(q6)) return
q_out = q_mul(q_out, q6)
if (.not. present(q7)) return
q_out = q_mul(q_out, q7)
if (.not. present(q8)) return
q_out = q_mul(q_out, q8)
if (.not. present(q9)) return
q_out = q_mul(q_out, q9)

!----------------------------------------------------
contains

function q_mul(q1, q2) result (q_out)

complex(rp) q1(0:3), q2(0:3), q_out(0:3)
complex(rp) a, b, c, d, e, f, g, h

A = (q1(0) + q1(1)) * (q2(0) + q2(1))
B = (q1(3) - q1(2)) * (q2(2) - q2(3))
C = (q1(0) - q1(1)) * (q2(2) + q2(3)) 
D = (q1(2) + q1(3)) * (q2(0) - q2(1))
E = (q1(1) + q1(3)) * (q2(1) + q2(2))
F = (q1(1) - q1(3)) * (q2(1) - q2(2))
G = (q1(0) + q1(2)) * (q2(0) - q2(3))
H = (q1(0) - q1(2)) * (q2(0) + q2(3))

q_out(0) = B + (-E - F + G + H)*0.5_rp
q_out(1) = A -  (E + F + G + H)*0.5_rp 
q_out(2) = C +  (E - F + G - H)*0.5_rp 
q_out(3) = D +  (E - F - G + H)*0.5_rp

end function q_mul

end function quat_mul_complex

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function quat_rotate_real (quat, vec_in) result (vec_out)
!
! Routine to rotate a vector using a quaternion..
!
! Input:
!   quat(0:3 )  -- real(rp): Quaternion to rotate with. Does not have to be normalized.
!   vec_in(3)   -- real(rp): Initial vector.
!
! Output:
!   vec_out(3)  -- real(rp): Final vector.
!-

function quat_rotate_real (quat, vec_in) result (vec_out)

real(rp) :: vec_in(3), vec_out(3), quat(0:3)
real(rp) :: q0_inv

! Use quaternion rotation formula: vec -> q * vec * q^t

q0_inv = -(quat(1)*vec_in(1) + quat(2)*vec_in(2) + quat(3)*vec_in(3))

vec_out = [quat(0)*vec_in(1) + quat(2)*vec_in(3) - quat(3)*vec_in(2), &
           quat(0)*vec_in(2) + quat(3)*vec_in(1) - quat(1)*vec_in(3), &
           quat(0)*vec_in(3) + quat(1)*vec_in(2) - quat(2)*vec_in(1)]

vec_out = [quat(0)*vec_out(1) + quat(2)*vec_out(3) - quat(3)*vec_out(2) - q0_inv*quat(1), &
           quat(0)*vec_out(2) + quat(3)*vec_out(1) - quat(1)*vec_out(3) - q0_inv*quat(2), &
           quat(0)*vec_out(3) + quat(1)*vec_out(2) - quat(2)*vec_out(1) - q0_inv*quat(3)] * &
                                                                     (1.0_rp / dot_product(quat, quat))

end function quat_rotate_real

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function quat_rotate_complex (quat, vec_in) result (vec_out)
!
! Routine to rotate a vector using a quaternion..
!
! Input:
!   quat(0:3)   -- complex(rp): Quaternion to rotate with. Does not have to be normalized.
!   vec_in(3)   -- complex(rp): Initial vector.
!
! Output:
!   vec_out(3)  -- complex(rp): Final vector.
!-

function quat_rotate_complex (quat, vec_in) result (vec_out)

complex(rp) :: vec_in(3), vec_out(3), quat(0:3)
complex(rp) :: q0_inv

! Use quaternion rotation formula: vec -> q * vec * q^t

q0_inv = -(quat(1)*vec_in(1) + quat(2)*vec_in(2) + quat(3)*vec_in(3))

vec_out = [quat(0)*vec_in(1) + quat(2)*vec_in(3) - quat(3)*vec_in(2), &
           quat(0)*vec_in(2) + quat(3)*vec_in(1) - quat(1)*vec_in(3), &
           quat(0)*vec_in(3) + quat(1)*vec_in(2) - quat(2)*vec_in(1)]

vec_out = [quat(0)*vec_out(1) + quat(2)*vec_out(3) - quat(3)*vec_out(2) - q0_inv*quat(1), &
           quat(0)*vec_out(2) + quat(3)*vec_out(1) - quat(1)*vec_out(3) - q0_inv*quat(2), &
           quat(0)*vec_out(3) + quat(1)*vec_out(2) - quat(2)*vec_out(1) - q0_inv*quat(3)] * &
                                                                     (1.0_rp / dot_product(quat, quat))

end function quat_rotate_complex

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function rotate_vec_given_axis_angle (vec_in, axis, angle) result (vec_out)
!
! Routine to rotate a vector.
!
! Input:
!   vec_in(3)   -- real(rp): Initial vector.
!   axis(3)     -- real(rp): Axis of rotation. Must be normalized to 1.
!   angle       -- real(rp): Angle to rotate by
!
! Output:
!   vec_out(3)  -- real(rp): Final vector.
!-

function rotate_vec_given_axis_angle (vec_in, axis, angle) result (vec_out)

real(rp) :: vec_in(3), vec_out(3), angle, ca, sa, q0_inv
real(rp) :: axis(:)

! Use quaternion rotation formula: vec_out = q * vec_in * q^-1
! q = cos(angle/2) + sin(angle/2) * (axis(1), axis(2), axis(3))

if (angle == 0) then  ! Simple case 
  vec_out = vec_in
  return 
endif

!

ca = cos(angle*0.5_rp)
sa = sin(angle*0.5_rp)

q0_inv = -sa * (axis(1)*vec_in(1) + axis(2)*vec_in(2) + axis(3)*vec_in(3))

vec_out = [ca*vec_in(1) + sa * (axis(2)*vec_in(3) - axis(3)*vec_in(2)), &
           ca*vec_in(2) + sa * (axis(3)*vec_in(1) - axis(1)*vec_in(3)), &
           ca*vec_in(3) + sa * (axis(1)*vec_in(2) - axis(2)*vec_in(1))]


vec_out = [ca*vec_out(1) + sa * (axis(2)*vec_out(3) - axis(3)*vec_out(2) - q0_inv*axis(1)), &
           ca*vec_out(2) + sa * (axis(3)*vec_out(1) - axis(1)*vec_out(3) - q0_inv*axis(2)), &
           ca*vec_out(3) + sa * (axis(1)*vec_out(2) - axis(2)*vec_out(1) - q0_inv*axis(3))]

end function rotate_vec_given_axis_angle

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine rotate_vec (vec, axis, angle)
!            
! Basic routine to rotate vector components around the x, y, or z axis.
! 
! Input:
!   vec(3)          -- real(rp): vector
!   axis            -- integer: x_axis$, y_axis$, or z_axis$
!   angle           -- real(rp): angle to rotate.
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
! Subroutine rotate_mat (mat, axis, angle, right_multiply)
!            
! Basic routine to apply a rotation matrix to mat that rotates around the x, y, or z axis.
! 
! Input:
!   mat(3,3)        -- real(rp): matrix
!   axis            -- integer: x_axis$, y_axis$, or z_axis$
!   angle           -- real(rp): cosine of the angle to rotate
!   right_multiply  -- logical, optional: If present and True then muliply rotation
!                       matrix on the right: mat = mat . Rot. Default is False.
!
! Output:
!   mat(3,3)        -- real(rp): Rotated matrix
!-

subroutine rotate_mat(mat, axis, angle, right_multiply)

real(rp) :: mat(:,:),  angle, temp, ca, sa
integer :: axis
logical, optional :: right_multiply

!

if (angle == 0) return ! Simple case

ca = cos(angle)
sa = sin(angle) 

if (logic_option(.false., right_multiply)) then
  select case (axis)
  case (x_axis$)
    mat(1, 2:3) =   [mat(1,2)*ca + mat(1,3)*sa, -mat(1,2)*sa + mat(1,3)*ca]
    mat(2, 2:3) =   [mat(2,2)*ca + mat(2,3)*sa, -mat(2,2)*sa + mat(2,3)*ca]
    mat(3, 2:3) =   [mat(3,2)*ca + mat(3,3)*sa, -mat(3,2)*sa + mat(3,3)*ca]

  case (y_axis$)
    mat(1, 1:3:2) = [mat(1,1)*ca - mat(1,3)*sa, mat(1,1)*sa + mat(1,3)*ca]
    mat(2, 1:3:2) = [mat(2,1)*ca - mat(2,3)*sa, mat(2,1)*sa + mat(2,3)*ca]
    mat(3, 1:3:2) = [mat(3,1)*ca - mat(3,3)*sa, mat(3,1)*sa + mat(3,3)*ca]

  case (z_axis$)
    mat(1, 1:2) =   [mat(1,1)*ca + mat(1,2)*sa, -mat(1,1)*sa + mat(1,2)*ca]
    mat(2, 1:2) =   [mat(2,1)*ca + mat(2,2)*sa, -mat(2,1)*sa + mat(2,2)*ca]
    mat(3, 1:2) =   [mat(3,1)*ca + mat(3,2)*sa, -mat(3,1)*sa + mat(3,2)*ca]
  end select

else
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
endif
  
end subroutine

end module
