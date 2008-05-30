!+
! Subroutine tao_view_cmd (i_universe)
!
! Routine to set the default universe.
! 
! Input:
!   i_universe -- Integer: Universe to view.
!
!  Output:
!-

subroutine tao_view_cmd (i_universe)

use tao_mod

implicit none


integer i_universe

character(20) :: r_name = 'tao_view_cmd'

! Check range

if (i_universe < lbound(s%u, 1) .or. ubound(s%u, 1) < i_universe) then
  call out_io (s_error$, r_name, 'UNIVERSE NUMBER OUT OF RANGE.')
  return
endif

s%global%u_view = i_universe

end subroutine
