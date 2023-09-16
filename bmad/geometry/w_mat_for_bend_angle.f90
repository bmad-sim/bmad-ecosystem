!+
! Function w_mat_for_bend_angle (angle, ref_tilt, r_vec) result (w_mat)
!
! Routine to compute the W matrix for the angle transformation in a bend.
! Using the notation in the Bmad manual:
!   w_mat = R_z(ref_tilt) . R_y(-angle) . R_z(-ref_tilt)
!
! Input:
!   angle       -- real(rp): Bending angle.
!   ref_tilt    -- real(rp): Reference tilt.
!   r_vec(3)    -- real(rp), optional: Starting position.
!
! Output:
!   w_mat(3,3)  -- real(rp): W matrix
!   r_vec(3)    -- real(rp), optional: position with ref_tilt transformation
!-

function w_mat_for_bend_angle (angle, ref_tilt, r_vec) result (w_mat)

use bmad_interface, dummy => w_mat_for_bend_angle

implicit none

real(rp) angle, ref_tilt, w_mat(3,3), t_mat(3,3)
real(rp), optional :: r_vec(3)

! By definition, positive angle is equivalent to negative x_pitch

w_mat = w_mat_for_x_pitch(-angle)

if (ref_tilt == 0) return

t_mat = w_mat_for_tilt (ref_tilt)

if (present(r_vec)) r_vec = matmul (t_mat, r_vec)

w_mat = matmul (t_mat, w_mat)
t_mat(1,2) = -t_mat(1,2); t_mat(2,1) = -t_mat(2,1) ! form inverse
w_mat = matmul (w_mat, t_mat)

end function w_mat_for_bend_angle
