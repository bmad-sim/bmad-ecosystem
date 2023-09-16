!+
! Function w_mat_for_x_pitch (x_pitch, return_inverse) result (w_mat)
! 
! Routine to return the transformation matrix for an x_pitch.
!
! Input:
!   x_pitch        -- real(rp): pitch angle
!   return_inverse -- logical, optional: If True, return the inverse matrix. Default is False.
!
! Output:
!   w_mat(3,3)     -- real(rp): Transformation matrix.
!-   

function w_mat_for_x_pitch (x_pitch, return_inverse) result (w_mat)

use bmad_struct

implicit none

real(rp) x_pitch, c_ang, s_ang
real(rp) :: w_mat(3,3)
logical, optional :: return_inverse

! An x_pitch corresponds to a rotation around the y axis.

c_ang = cos(x_pitch); s_ang = sin(x_pitch)

if (logic_option(.false., return_inverse)) then
  w_mat(1,:) = [ c_ang, 0.0_rp,  -s_ang]
  w_mat(2,:) = [0.0_rp, 1.0_rp,  0.0_rp]
  w_mat(3,:) = [ s_ang, 0.0_rp,   c_ang]
else
  w_mat(1,:) = [ c_ang, 0.0_rp,   s_ang]
  w_mat(2,:) = [0.0_rp, 1.0_rp,  0.0_rp]
  w_mat(3,:) = [-s_ang, 0.0_rp,   c_ang]
endif

end function w_mat_for_x_pitch
