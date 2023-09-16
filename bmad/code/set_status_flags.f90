!+
! Subroutine set_status_flags (bookkeeping_state, stat)
!
! Routine to set the bookkeeping status block.
!
! Input:
!   stat          -- Integer: bookkeeping status. ok$, stale$, etc.
!
! Output:
!   bookkeeping_state -- bookkeeping_state_struct: 
!-

subroutine set_status_flags (bookkeeping_state, stat)

use bmad_struct

implicit none

type (bookkeeping_state_struct) bookkeeping_state
integer stat

!

bookkeeping_state%control        = stat
bookkeeping_state%s_position     = stat
bookkeeping_state%floor_position = stat
bookkeeping_state%ref_energy     = stat
bookkeeping_state%attributes     = stat
bookkeeping_state%mat6           = stat
bookkeeping_state%rad_int        = stat
bookkeeping_state%ptc            = stat

end subroutine set_status_flags

