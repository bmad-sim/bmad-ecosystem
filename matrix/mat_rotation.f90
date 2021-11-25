!+
! Subroutine mat_rotation (mat, angle, bet_1, bet_2, alph_1, alph_2)
!
! Subroutine to construct a 2x2 rotation matrix for translation from
! point 1 to point 2
!
! Input:
!   angle          -- Real(rp):  Rotation angle in radians.
!   bet_1, bet_2   -- Real(rp):  Beta values at points 1 and 2.
!   alph_1, alph_2 -- Real(rp):  Alpha values at points 1 and 2.
!
! Output: 
!   mat(:,:)       -- Real(rp): 2x2 rotation matrix.
!-

subroutine mat_rotation (mat, angle, bet_1, bet_2, alph_1, alph_2)

  use precision_def

  implicit none

  real(rp) angle, bet_1, bet_2, alph_1, alph_2, mat(:,:)
  real(rp) c, s

!

  c = cos(angle)
  s = sin(angle)

  mat(1, 1) = sqrt(bet_2/bet_1) * (c + alph_1 * s)
  mat(1, 2) = sqrt(bet_1*bet_2) * s
  mat(2, 1) = -((1 + alph_1*alph_2)*s + (alph_2-alph_1)*c) / sqrt(bet_1*bet_2)
  mat(2, 2) = sqrt(bet_1/bet_2) * (c - alph_2*s)

end subroutine
