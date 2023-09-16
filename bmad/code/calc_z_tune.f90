!+
! Subroutine calc_z_tune (branch)
!
! Routine to calculate the synchrotron tune from the full 6X6 1-turn matrix.
!
! Note: The tune will be negative above transition which corresponds to
! counter-clockwise rotation in z-pz space.
!
! Input:
!   branch   -- branch_struct: Lattice branch
!
! Output:
!   branch -- branch_struct
!     branch(ix_branch)%z%tune            -- Synchrotron tune (radians). If unstable tune = 0.
!     branch(ix_branch)%z%stable          -- Is the mode stable? If no rf then tune is zero but is stable.
!     branch(ix_branch)%param%t1_with_RF  -- 6x6 1-turn matrix.
!-

subroutine calc_z_tune (branch)

use bmad_interface, except_dummy => calc_z_tune

implicit none

type (branch_struct), target :: branch

real(rp) a(6,6), cos_z, denom
complex(rp) eval(6), evec(6,6)

integer i, sgn

logical err

!

call transfer_matrix_calc (branch%lat, a, ix_branch = branch%ix_branch)
branch%param%t1_with_RF = a

denom = 2 * (a(5,5)*a(6,6) - a(5,6)*a(6,5))
if (denom == 0) then
  branch%z%tune = 0
  branch%z%stable = .false.
  return
endif

cos_z = (a(5,5) + a(6,6)) / denom
sgn = sign_of(a(5,6))

! Previous cutoff of 1d-7 was too restrictive for IOTA lattices.

if (cos_z - 1 > -1d-9) then
  branch%z%tune = 0
  branch%z%stable = (abs(cos_z-1) < 1d-9)
  return
endif

call mat_eigen(a, eval, evec, err)

branch%z%tune = sgn * abs(atan2(aimag(eval(5)), real(eval(5), rp)))
branch%z%stable = .true.

end subroutine
