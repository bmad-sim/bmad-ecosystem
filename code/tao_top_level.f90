!+
! Subroutine tao_top_level (command, errcode)
!
! Top level tao routine.
!
! Modules needed:
!   use tao_mod
!
! Input:
!   command    -- character(*), optional: Tao command string. 
!                                         If present, getting user input from the terminal is bypassed. 
!                                          
! Output:
!   errcode    -- integer, optional: Return error code
!-

subroutine tao_top_level (command, errcode)

use tao_command_mod, dummy => tao_top_level
!use tao_mpi_mod

implicit none

type (tao_super_universe_struct), pointer :: s_ptr  ! For debug purposes
type (tao_universe_struct), pointer :: u

integer, optional :: errcode

character(*), optional :: command
character(1000) :: cmd_line
character(16) :: r_name = 'tao_top_level'

logical found, err, will_need_cmd_line_input, interactive

! init

s_ptr => s       ! Used for debugging

! Set interactive flags
if (present(command)) then
  cmd_line = command
  interactive = .false. 
  s%com%shell_interactive = .false.
else
  interactive = .true. 
  s%com%shell_interactive = .true.
endif
will_need_cmd_line_input = .false.


! Read command line arguments.

if (interactive) then
  call tao_parse_command_args (err)
  if (err) stop
endif

! Turn off plotting for slaves
!if (.not. s%mpi%master)  s%global%plot_on = .false.

! And init everything.

if (.not. s%global%initialized) then
  call tao_init (err)
  if (err) then
    call out_io (s_fatal$, r_name, 'TAO INIT FILE NOT FOUND. STOPPING.')
    stop
  endif
  s%global%initialized = .true. 
endif

u => s%u(1)  ! Used for debugging

! MPI slave settings
!if (.not. s%mpi%master) then 
!  !Turn off screen output
!  s%com%print_to_terminal = .false.
!  call output_direct( do_print = .false.)
!endif

! Command loop
do
  err = .false.
  
  if (interactive) then
    call tao_get_user_input (cmd_line)
  else
    call tao_get_user_input (cmd_line, will_need_cmd_line_input = will_need_cmd_line_input)
  endif
  
 ! if (s%mpi%master) call tao_get_user_input (cmd_line)
  
  ! Broadcast command to slaves
  !if (s%mpi%on) call tao_broadcast_chars_mpi(cmd_line)
  
  if (s%com%single_mode) then
    ! single mode
    call tao_single_mode (cmd_line(1:1))
    ! Do the standard calculations and plotting after command execution.
    call tao_cmd_end_calc ()
  else
    ! command line mode
    call tao_hook_command (cmd_line, found)
    if (.not. found) call tao_command (cmd_line, err)
  endif
  if (.not. err) call tao_cmd_history_record (cmd_line)

  ! Non-interactive will exit
  if (will_need_cmd_line_input) exit
enddo

if (present(errcode) )then
  if (err) then
    errcode = 1
  else
    errcode = 0
  endif
endif

end subroutine

