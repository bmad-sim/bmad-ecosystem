!+
! Subroutine set_lords_status_stale (ele, stat_group, control_bookkeeping, flag)
!
! Routine to recursively set the status flag of all lords of an element.
!
! Input:
!   ele        -- ele_struct: Element
!   stat_group -- Integer: which status group to set. floor_position_group$, etc.
!                   See set_ele_status_stale for more details.
!   control_bookkeeping
!              -- logical, optional: Call control_bookkeeper for each lord if needed? 
!                   Default if False.
!   flag       -- integer, optional: Do not use. For coordinating recursion.
!
! Output:
!   ele%lat    -- Lat_struct: Lattice with status flags of lords of ele set.
!-

recursive subroutine set_lords_status_stale (ele, stat_group, control_bookkeeping, flag)

use bookkeeper_mod, dummy => set_lords_status_stale

implicit none

type (ele_struct) ele
type (ele_struct), pointer :: lord
integer stat_group, i
logical, optional :: control_bookkeeping
integer, optional :: flag

! First time through the flag argument will not be present.
! Do not set status first time through since this is the original element.
! That is, only want to set the flags of the lords.

if (present(flag)) call set_ele_status_stale (ele, stat_group, .false.)
if (logic_option(.false., control_bookkeeping) .and. &
      (ele%bookkeeping_state%control /= ok$ .or. ele%bookkeeping_state%attributes /= ok$)) then
  call control_bookkeeper (ele%branch%lat, ele)
endif

do i = 1, ele%n_lord
  lord => pointer_to_lord (ele, i)
  call set_lords_status_stale (lord, stat_group, control_bookkeeping, 1)
enddo

end subroutine set_lords_status_stale

