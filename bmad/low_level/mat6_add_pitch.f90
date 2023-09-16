!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat6_add_pitch (x_pitch_tot, y_pitch_tot, orientation, mat6)
!
! Subroutine to modify a first order transfer matrix to include the affect
! of an element pitch. Note that this routine does not correct the 0th order
! part of the map. It is assumed that on input the transfer map
! does not include the affect of any pitches.
!
! Input:
!   x_pitch_tot -- Real(rp): Horizontal pitch
!   y_pitch_tot -- Real(rp): Vertical pitch
!   orientation -- integer: Element longitudinal orientation. +1 or -1.
!   mat6(6,6)   -- Real(rp): 1st order part of the transfer map (Jacobian).
!
! Output:
!   mat6(6,6) -- Real(rp): 1st order xfer map with pitches.
!-

subroutine mat6_add_pitch (x_pitch_tot, y_pitch_tot, orientation, mat6)

use bmad_struct

implicit none

real(rp) mat6(6,6), x_pitch_tot, y_pitch_tot
integer orientation

!

if (x_pitch_tot == 0 .and. y_pitch_tot == 0) return

! The equations below are performing matrix multiplication. The original matrix
! is being multiplied from left and right by matrices that correspond to the pitches. 
! The pitch matrices are obtained by differentiating the corresponding equations in   
! the offset_particle subroutine. The (i,j) numbers mentioned as comments refer to  
! the non-zero elements present in the pitch matrices. 

mat6(:,6) = mat6(:,6) - mat6(:,2) * orientation * x_pitch_tot ! (2,6)
mat6(:,1) = mat6(:,1) + mat6(:,5) * orientation * x_pitch_tot ! (5,1)

mat6(:,6) = mat6(:,6) - mat6(:,4) * orientation * y_pitch_tot ! (4,6)
mat6(:,3) = mat6(:,3) + mat6(:,5) * orientation * y_pitch_tot ! (5,3)

mat6(2,:) = mat6(2,:) + orientation * x_pitch_tot * mat6(6,:) ! (2,6)
mat6(5,:) = mat6(5,:) - orientation * x_pitch_tot * mat6(1,:) ! (5,1)

mat6(4,:) = mat6(4,:) + orientation * y_pitch_tot * mat6(6,:) ! (4,6)
mat6(5,:) = mat6(5,:) - orientation * y_pitch_tot * mat6(3,:) ! (5,3)

end subroutine mat6_add_pitch

