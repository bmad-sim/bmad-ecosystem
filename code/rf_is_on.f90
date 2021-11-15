!+
! Function rf_is_on (branch, ix_ele1, ix_ele2) result (is_on)
!
! Routine to check if any rfcavity element is powered in a branch in the range (ix_ele1, ix_ele2].
! Note: the range is not inclusive of ix_ele1.
!
! If ix_ele2 <= ix_ele1 then the range is wrapped around the lattice ends.
!
! Input:
!   branch  -- branch_struct: Lattice branch to check.
!   ix_ele1 -- integer, optional: Start of range of elements to check. Default is 0.
!   ix_ele2 -- integer, optional: End of range of elements to check. Default is branch%n_ele_track.
!
! Output:
!   is_on   -- Logical: True if any rfcavity is powered. False otherwise.
!-

function rf_is_on (branch, ix_ele1, ix_ele2) result (is_on)

use bmad_struct
implicit none

type (branch_struct), target :: branch
type (ele_struct), pointer :: ele

integer, optional :: ix_ele1, ix_ele2
integer i, ie1, ie2
logical is_on

!

is_on = .false.

ie1 = integer_option(0, ix_ele1)
ie2 = integer_option(branch%n_ele_track, ix_ele2)

if (ie2 > ie1) then
  do i = ie1+1, ie2
    ele => branch%ele(i)
    if (ele%key /= rfcavity$ .or. .not. ele%is_on .or. ele%value(voltage$) == 0) cycle
    is_on = .true.
    return
  enddo

else
  do i = ie1+1, branch%n_ele_track
    ele => branch%ele(i)
    if (ele%key /= rfcavity$ .or. .not. ele%is_on .or. ele%value(voltage$) == 0) cycle
    is_on = .true.
    return
  enddo

  do i = 1, ie2
    ele => branch%ele(i)
    if (ele%key /= rfcavity$ .or. .not. ele%is_on .or. ele%value(voltage$) == 0) cycle
    is_on = .true.
    return
  enddo
endif

end function rf_is_on
