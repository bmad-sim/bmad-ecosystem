!+
! Subroutine floor_angles_to_w_mat (theta, phi, psi, w_mat, w_mat_inv)
!
! Routine to construct the W matrix that specifies the orientation of an element
! in the global "floor" coordinates. See the Bmad manual for more details.
!
! Input:
!   theta -- Real(rp): Azimuth angle.
!   phi   -- Real(rp): Pitch angle.
!   psi   -- Real(rp): Roll angle.
!
! Output:
!   w_mat(3,3)     -- Real(rp), optional: Orientation matrix.
!   w_mat_inv(3,3) -- Real(rp), optional: Inverse Orientation matrix.
!-

subroutine floor_angles_to_w_mat (theta, phi, psi, w_mat, w_mat_inv)

use bmad_struct

implicit none

real(rp), optional :: w_mat(3,3), w_mat_inv(3,3)
real(rp) theta, phi, psi
real(rp) s_the, c_the, s_phi, c_phi, s_psi, c_psi

!

s_the = sin(theta); c_the = cos(theta)
s_phi = sin(phi);   c_phi = cos(phi)
s_psi = sin(psi);   c_psi = cos(psi)

if (present(w_mat)) then
  w_mat(1,1) =  c_the * c_psi - s_the * s_phi * s_psi
  w_mat(1,2) = -c_the * s_psi - s_the * s_phi * c_psi
  w_mat(1,3) =  s_the * c_phi
  w_mat(2,1) =  c_phi * s_psi
  w_mat(2,2) =  c_phi * c_psi
  w_mat(2,3) =  s_phi 
  w_mat(3,1) = -s_the * c_psi - c_the * s_phi * s_psi
  w_mat(3,2) =  s_the * s_psi - c_the * s_phi * c_psi 
  w_mat(3,3) =  c_the * c_phi
endif

if (present(w_mat_inv)) then
  w_mat_inv(1,1) =  c_the * c_psi - s_the * s_phi * s_psi
  w_mat_inv(1,2) =  c_phi * s_psi 
  w_mat_inv(1,3) = -s_the * c_psi - c_the * s_phi * s_psi 
  w_mat_inv(2,1) = -c_the * s_psi - s_the * s_phi * c_psi 
  w_mat_inv(2,2) =  c_phi * c_psi
  w_mat_inv(2,3) =  s_the * s_psi - c_the * s_phi * c_psi
  w_mat_inv(3,1) =  s_the * c_phi
  w_mat_inv(3,2) =  s_phi
  w_mat_inv(3,3) =  c_the * c_phi
endif

end subroutine floor_angles_to_w_mat 
