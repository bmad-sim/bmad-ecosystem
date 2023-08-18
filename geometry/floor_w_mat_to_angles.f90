!+
! Subroutine floor_w_mat_to_angles (w_mat, theta, phi, psi, floor0)
!
! Routine to construct the angles that define the orientation of an element
! in the global "floor" coordinates from the W matrix. See the Bmad manual for more details.
!
! Input:
!   w_mat(3,3) -- Real(rp): Orientation matrix.
!   floor0     -- floor_position_struct, optional: There are two solutions related by:
!                   [theta, phi, psi] & [pi+theta, pi-phi, pi+psi]
!                 If floor0 is present, choose the solution "nearest" the angles in floor0.
!
! Output:
!   theta -- Real(rp): Azimuth angle.
!   phi   -- Real(rp): Pitch angle.
!   psi   -- Real(rp): Roll angle.
!-

subroutine floor_w_mat_to_angles (w_mat, theta, phi, psi, floor0)

use bmad_struct

implicit none

type (floor_position_struct), optional :: floor0
type (floor_position_struct) f0
real(rp) theta, phi, psi, w_mat(3,3)
real(rp) diff1(3), diff2(3)

! special degenerate case

if (abs(w_mat(1,3)) + abs(w_mat(3,3)) < 1d-12) then 
  ! Note: Only theta +/- psi is well defined here so this is rather arbitrary.
  if (present(floor0)) then
    theta = floor0%theta
  else
    theta = 0
  endif

  if (w_mat(2,3) > 0) then
    phi = pi/2
    psi = atan2(-w_mat(3,1), w_mat(1,1)) - theta
  else
    phi = -pi/2
    psi = atan2(w_mat(3,1), w_mat(1,1)) + theta
  endif

! normal case

else 
  theta = atan2 (w_mat(1,3), w_mat(3,3))
  phi = atan2 (w_mat(2,3), sqrt(w_mat(1,3)**2 + w_mat(3,3)**2))
  psi = atan2 (w_mat(2,1), w_mat(2,2))

  if (present(floor0)) then
    f0 = floor0
    diff1 = [modulo2(theta-f0%theta, pi), modulo2(phi-f0%phi, pi), modulo2(psi-f0%psi, pi)]
    diff2 = [modulo2(pi+theta-f0%theta, pi), modulo2(pi-phi-f0%phi, pi), modulo2(pi+psi-f0%psi, pi)]
    if (sum(abs(diff2)) < sum(abs(diff1))) diff1 = diff2
    theta = diff1(1) + f0%theta
    phi   = diff1(2) + f0%phi
    psi   = diff1(3) + f0%psi
  else
    theta = theta - twopi * nint((theta ) / twopi)
  endif

endif

end subroutine floor_w_mat_to_angles 
