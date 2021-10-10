!+
! Subroutine tao_clear_cmd (cmd_line)
!
! Routine to clear Taylor maps from lattice elements.
!
! Input:
!   cmd_line    -- Character(*): Should be set to 'maps'.
!-

subroutine tao_clear_cmd (cmd_line)

use tao_struct

implicit none

integer iu
character(*) cmd_line
character(*), parameter :: r_name = 'cmd_line'

!

do iu = lbound(s%u,1), ubound(s%u,1)
  call clear_taylor_maps_from_elements(s%u(iu)%model%lat)
enddo

end subroutine
