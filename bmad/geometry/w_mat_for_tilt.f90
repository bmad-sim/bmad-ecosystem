!+
! Function w_mat_for_tilt (tilt, return_inverse) result (w_mat)
! 
! Routine to return the transformation matrix for an tilt.
!
! Input:
!   tilt           -- real(rp): pitch angle
!   return_inverse -- logical, optional: If True, return the inverse matrix. Default is False.
!
! Output:
!   w_mat(3,3)     -- real(rp): Transformation matrix.
!-   

function w_mat_for_tilt (tilt, return_inverse) result (w_mat)

use bmad_struct

implicit none

real(rp) tilt, c_ang, s_ang
real(rp) :: w_mat(3,3)
logical, optional :: return_inverse

! An tilt corresponds to a rotation around the y axis.

c_ang = cos(tilt); s_ang = sin(tilt)

if (logic_option(.false., return_inverse)) then
  w_mat(1,:) = [ c_ang,  s_ang,  0.0_dp ]
  w_mat(2,:) = [-s_ang,  c_ang,  0.0_dp ]
  w_mat(3,:) = [0.0_dp,  0.0_dp, 1.0_dp ]
else
  w_mat(1,:) = [c_ang,  -s_ang,  0.0_dp ]
  w_mat(2,:) = [s_ang,   c_ang,  0.0_dp ]
  w_mat(3,:) = [0.0_dp,  0.0_dp, 1.0_dp ]
endif

end function w_mat_for_tilt
