!+
! Subroutine tao_view_cmd (s, i_universe)
!
! Routine to set the default universe.
! 
! Input:
!   i_universe -- Integer: Universe to view.
!
!  Output:
!   s     -- Tao_super_universe_struct
!-

subroutine tao_view_cmd (s, i_universe)

use tao_mod

implicit none

type (tao_super_universe_struct) s

integer i_universe

character(20) :: r_name = 'tao_view_cmd'

! Check range

if (i_universe < 1 .or. size(s%u) < i_universe) then
  call out_io (s_error$, r_name, 'UNIVERSE NUMBER OUT OF RANGE.')
  return
endif

s%global%u_view = i_universe

end subroutine
