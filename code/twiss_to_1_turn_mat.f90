!+
! Subroutine twiss_to_1_turn_mat (twiss, phi, mat2)
!
! Subroutine to form the 2x2 1-turn transfer matrix from the twiss parameters.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   twiss -- Twiss_struct: Structure holding the Twiss parameters.
!     %beta
!     %alpha
!     %gamma
!   phi   -- Real: Tune in radians.
!
! Output:
!   mat2(2,2) -- Real: 1-turn matrix.
!-

subroutine twiss_to_1_turn_mat (twiss, phi, mat2)

  use bmad_struct

  type (twiss_struct) twiss

  real phi, mat2(2,2), c, s

!

  c = cos(phi)
  s = sin(phi)

  mat2(1,1) =  c + s * twiss%alpha
  mat2(1,2) =  s * twiss%beta
  mat2(2,1) = -s * twiss%gamma
  mat2(2,2) =  c - s * twiss%alpha

end subroutine
