!+
! Function bend_shift(position1, g, delta_s, w_mat, tilt) result(position2)
!
! Function to shift frame of reference within a bend with curvature g and tilt.
! Note: position2%theta, %phi, and %psi are not calculated. 
! 
! Input:
!   position1    -- floor_position_struct: Position of particle in frame 1 coordinates (Caretesian).
!   g            -- real(rp): curvature (1/rho)
!   delta_s      -- real(rp): relative s-position of frame 2 relative to frame 1
!   tilt         -- real(rp), optional: tilt. Default: 0
!
! Output:
!   position2    -- floor_position_struct: Coordinates relative to frame 2
!   w_mat(3,3)   -- real(rp), optional: W matrix used in the transformation   
!-

function bend_shift (position1, g, delta_s, w_mat, tilt) result(position2)

use bmad_struct

implicit none

type (floor_position_struct) :: position1, position2
real(rp) :: g, delta_s, S_mat(3,3), L_vec(3), tlt, angle
real(rp), optional :: w_mat(3,3), tilt

!

angle = delta_s * g

if (angle == 0) then
  position2 = position1
  position2%r(3) = position2%r(3) - delta_s
  if (present(w_mat)) call mat_make_unit(w_mat)
  return
endif

!

tlt = real_option(0.0_rp, tilt)
call mat_make_unit(S_mat)

if (tlt /= 0) then
  call rotate_mat(S_mat, z_axis$, -tilt)
  call rotate_mat(S_mat, y_axis$,  angle)
  call rotate_mat(S_mat, z_axis$,  tilt)
else
  call rotate_mat(S_mat, y_axis$, angle)
endif

L_vec = [cos_one(angle), 0.0_rp, -sin(angle)]/g
if (present(tilt)) call rotate_vec(L_vec, z_axis$, tilt)

position2%r = matmul(S_mat, position1%r) + L_vec
position2%w = matmul(S_mat, position1%w)

if (present(w_mat)) w_mat = s_mat

end function bend_shift
