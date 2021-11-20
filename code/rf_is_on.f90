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

use bmad_routine_interface, dummy => rf_is_on
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
    if (ele%key /= rfcavity$) cycle
    call is_this_cavity_on(ele, is_on)
    if (is_on) return
  enddo

else
  do i = ie1+1, branch%n_ele_track
    ele => branch%ele(i)
    if (ele%key /= rfcavity$) cycle
    call is_this_cavity_on(ele, is_on)
    if (is_on) return
  enddo

  do i = 1, ie2
    ele => branch%ele(i)
    if (ele%key /= rfcavity$) cycle 
    call is_this_cavity_on(ele, is_on)
    if (is_on) return
  enddo
endif

!-----------------------------------------------
contains

subroutine is_this_cavity_on(ele, is_on)

type (ele_struct) ele
logical is_on

! It can happen that the User toggles the %is_on attribute but the program has not yet 
! called make_mat6. So call make_mat6 if needed.

if ((ele%is_on .and. ele%value(voltage$) /= 0 .and. ele%mat6(6,5) == 0) .or. &
        (.not. ele%is_on .and. ele%mat6(6,5) /= 0)) call make_mat6 (ele, branch%param, ele%map_ref_orb_in)

if (ele%mat6(6,5) /= 0) is_on = .true.

end subroutine is_this_cavity_on

end function rf_is_on
