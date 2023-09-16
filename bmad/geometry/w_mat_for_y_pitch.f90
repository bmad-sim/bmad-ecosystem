
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function w_mat_for_y_pitch (y_pitch, return_inverse) result (w_mat)
! 
! Routine to return the transformation matrix for an y_pitch.
!
! Input:
!   y_pitch        -- real(rp): pitch angle
!   return_inverse -- logical, optional: If True, return the inverse matrix. Default is False.
!
! Output:
!   w_mat(3,3)     -- real(rp): Transformation matrix.
!-   

function w_mat_for_y_pitch (y_pitch, return_inverse) result (w_mat)

use bmad_struct

implicit none

real(rp) y_pitch, c_ang, s_ang
real(rp) :: w_mat(3,3)
logical, optional :: return_inverse

! An y_pitch corresponds to a rotation around the y axis.

c_ang = cos(y_pitch); s_ang = sin(y_pitch)

if (logic_option(.false., return_inverse)) then
  w_mat(1,:) = [1.0_rp,  0.0_rp, 0.0_rp]
  w_mat(2,:) = [0.0_rp,  c_ang,  -s_ang]
  w_mat(3,:) = [0.0_rp,  s_ang,   c_ang]
else
  w_mat(1,:) = [1.0_rp,  0.0_rp, 0.0_rp]
  w_mat(2,:) = [0.0_rp,  c_ang,   s_ang]
  w_mat(3,:) = [0.0_rp, -s_ang,   c_ang]
endif

end function w_mat_for_y_pitch
