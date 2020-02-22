!+
! Subroutine tilt_coords_photon (tilt_val, coord, w_mat)
!
! Subroutine to effectively tilt (rotate in the x-y plane) an element
! by instead rotating the particle position with negative the angle.
!
! Input:
!   tilt_val       -- Real(rp): Tilt value (could be the roll value for a bend)
!   coord(:)       -- Real(rp): Coordinates of particle before rotation.
!   w_mat(3,3)     -- Real(rp), optional: Rotation matrix before tilt.
!
! Output:
!   coord(:)       -- Real(rp): Coordinates of particle after rotation.
!   w_mat(3,3)     -- Real(rp), optional: Rotation matrix after tilt.
!-

subroutine tilt_coords_photon (tilt_val, coord, w_mat)

use bmad_interface, dummy => tilt_coords_photon

implicit none

real(rp) tilt_val, coord(:), wm(3,3)
real(rp), optional :: w_mat(3,3)

!

if (tilt_val == 0) return
wm =  w_mat_for_tilt (-tilt_val)

coord(1:3:2) = matmul(wm(1:2,1:2), coord(1:3:2))
coord(2:4:2) = matmul(wm(1:2,1:2), coord(2:4:2))

if (present(w_mat)) w_mat = matmul(w_mat, wm)

end subroutine
