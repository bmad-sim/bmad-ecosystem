!+
! Subroutine twiss_to_1_turn_mat (twiss, phi, mat2)
!
! Subroutine to form the 2x2 1-turn transfer matrix from the twiss parameters.
!
! Input:
!   twiss -- Twiss_struct: Structure holding the Twiss parameters.
!     %beta
!     %alpha
!   phi   -- Real(rp): Tune in radians.
!
! Output:
!   mat2(2,2) -- Real(rp): 1-turn matrix.
!-

subroutine twiss_to_1_turn_mat (twiss, phi, mat2)

use bmad_struct

implicit none

type (twiss_struct) twiss

real(rp) phi, mat2(2,2), c, s

!

c = cos(phi)
s = sin(phi)

mat2(1,1) =  c + s * twiss%alpha
mat2(1,2) =  s * twiss%beta
mat2(2,1) = -s * (1 + twiss%alpha**2) / twiss%beta
mat2(2,2) =  c - s * twiss%alpha

end subroutine twiss_to_1_turn_mat

