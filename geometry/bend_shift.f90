!+
! Function bend_shift(position1, g, delta_s, w_mat, ref_tilt) result(position2)
!
! Function to shift a particle's coordinates from one coordinate frame of reference 
! to another within a bend with curvature g and reference tilt ref_tilt.
! Note: position2%theta, %phi, and %psi are not calculated. 
! 
! Input:
!   position1    -- floor_position_struct: Position of particle in inital coordinate frame.
!   g            -- real(rp): Curvature (1/rho)
!   delta_s      -- real(rp): S-position of final frame relative to the initial frame.
!   ref_tilt     -- real(rp), Optional: ref_tilt. Default: 0
!
! Output:
!   position2    -- floor_position_struct: particle coordinates relative to the final frame.
!   w_mat(3,3)   -- real(rp), Optional: W matrix used in the transformation   
!-

function bend_shift (position1, g, delta_s, w_mat, ref_tilt) result(position2)

use bmad_struct

implicit none

type (floor_position_struct) :: position1, position2
real(rp) :: g, delta_s, S_mat(3,3), L_vec(3), tlt, angle
real(rp), optional :: w_mat(3,3), ref_tilt

!

angle = delta_s * g

if (angle == 0) then
  position2 = position1
  position2%r(3) = position2%r(3) - delta_s
  if (present(w_mat)) call mat_make_unit(w_mat)
  return
endif

!

tlt = real_option(0.0_rp, ref_tilt)
call mat_make_unit(S_mat)

if (tlt /= 0) then
  call rotate_mat(S_mat, z_axis$, -ref_tilt)
  call rotate_mat(S_mat, y_axis$,  angle)
  call rotate_mat(S_mat, z_axis$,  ref_tilt)
else
  call rotate_mat(S_mat, y_axis$, angle)
endif

L_vec = [cos_one(angle), 0.0_rp, -sin(angle)]/g
if (present(ref_tilt)) call rotate_vec(L_vec, z_axis$, ref_tilt)

position2%r = matmul(S_mat, position1%r) + L_vec
position2%w = matmul(S_mat, position1%w)

if (present(w_mat)) w_mat = s_mat

end function bend_shift
