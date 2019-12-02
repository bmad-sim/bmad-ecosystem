!+
! Subroutine tilt_coords (tilt_val, coord, mat6, make_matrix)
!
! Subroutine to effectively tilt (rotate in the x-y plane) an element
! by instead rotating the particle position with negative the angle.
!
! Input:
!   tilt_val       -- Real(rp): Tilt value (could be the roll value for a bend)
!   coord(:)       -- Real(rp): Coordinates of particle before rotation.
!   mat6(6,6)      -- Real(rp), optional: Transfer matrix before tilt.
!   make_matrix    -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   coord(:)   -- Real(rp): Coordinates of particle after rotation.
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix transfer matrix after tilt applied.
!-

subroutine tilt_coords (tilt_val, coord, mat6, make_matrix)

use sim_utils

implicit none

real(rp) :: cos_ang, sin_ang
real(rp) tilt_val, coord(:), rot_mat(2,2), mm(6,6)
real(rp), optional :: mat6(6,6)

logical, optional :: make_matrix

!

if (tilt_val == 0) return

sin_ang = sin(tilt_val)
cos_ang = cos(tilt_val)

rot_mat(1,1) =  cos_ang
rot_mat(1,2) =  sin_ang
rot_mat(2,1) = -sin_ang
rot_mat(2,2) =  cos_ang

coord(1:3:2) = matmul(rot_mat, coord(1:3:2))
coord(2:4:2) = matmul(rot_mat, coord(2:4:2))

if (logic_option(.false., make_matrix)) then
  mm = mat6
  mat6(1,:) = cos_ang * mm(1,:) + sin_ang * mm(3,:)
  mat6(2,:) = cos_ang * mm(2,:) + sin_ang * mm(4,:)
  mat6(3,:) = cos_ang * mm(3,:) - sin_ang * mm(1,:)
  mat6(4,:) = cos_ang * mm(4,:) - sin_ang * mm(2,:)
endif

end subroutine
