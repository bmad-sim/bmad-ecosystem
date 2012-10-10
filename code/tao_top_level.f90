!+
! Subroutine tao_top_level ()
!
! Top level tao routine.
!
! Modules needed:
!   use tao_mod
!-

subroutine tao_top_level ()

use tao_command_mod, dummy => tao_top_level

implicit none

type (tao_super_universe_struct), pointer :: s_ptr  ! For debug purposes
type (tao_common_struct), pointer :: t_ptr          ! For debug purposes
type (tao_universe_struct), pointer :: u

character(200) cmd_line
character(16) :: r_name = 'tao_top_level'

logical found, err

! init

s_ptr => s
t_ptr => tao_com

! Read command line arguments.

call tao_parse_command_args (err)
if (err) stop

! And init everything.

call tao_init (err)
if (err) then
  call out_io (s_fatal$, r_name, 'TAO INIT FILE NOT FOUND. STOPPING.')
  stop
endif

u => s%u(1)

! Command loop

do
  err = .false.
  call tao_get_user_input (cmd_line)
  if (tao_com%single_mode) then
    ! single mode
    call tao_single_mode (cmd_line(1:1))
    ! Do the standard calculations and plotting after command execution.
    call tao_cmd_end_calc
  else
    ! command line mode
    call tao_hook_command (cmd_line, found)
    if (.not. found) call tao_command (cmd_line, err)
  endif
  if (.not. err) call tao_cmd_history_record (cmd_line)
enddo

end subroutine
