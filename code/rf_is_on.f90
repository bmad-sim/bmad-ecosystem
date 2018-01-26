!+
! Function rf_is_on (branch) result (is_on)
!
! Routine to check if any rfcavity is powered in a branch.
!
! Input:
!   branch -- branch_struct: Lattice branch to check.
!
! Output:
!   is_on  -- Logical: True if any rfcavity is powered. False otherwise.
!-

function rf_is_on (branch) result (is_on)

use bmad_struct
implicit none

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele

integer i
logical is_on

!

is_on = .false.

do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  if (ele%key == rfcavity$ .and. ele%is_on .and. ele%value(voltage$) /= 0) then
    is_on = .true.
    return
  endif
enddo

end function rf_is_on
