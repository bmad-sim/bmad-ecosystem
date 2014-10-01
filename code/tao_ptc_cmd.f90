!+
! Subroutine tao_ptc_cmd (what, input_str)
!
! Interface to ptc layout.
!
! Input:
!   input_str  -- Character(*): What to show.
!-

subroutine tao_ptc_cmd (what, input_str)

use tao_mod
use ptc_layout_mod

implicit none

integer ix_cmd

character(*) what, input_str
character(16) :: command(2) = ['init  ', 'script']
character(16) cmd_name
character(*), parameter :: r_name = 'tao_ptc_cmd'

!

call match_word (what, command, ix_cmd, .true., matched_name = cmd_name)
if (ix_cmd == 0) then
  call out_io (s_error$, r_name, 'BAD COMMAND (SHOULD BE "init", or "script")')
  return
endif

select case (cmd_name)

! init

case ('init')
  call lat_to_ptc_layout (s%u(1)%model%lat)

! script

case ('call')
  call read_ptc_command77 (input_str)

end select

end subroutine
