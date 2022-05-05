!+
! Subroutine tao_top_level (command, errcode)
!
! Top level Tao routine.
!
! Basic idea: If the command argument is not present, this routine will go into an endless loop
! reading commands from the terminal. If the command argument is present, for example when 
! interfacing with Python, this routine handles a command and then returns.
! [When being called from Python, Python calls tao_c_command which in turn calls this routine.]
!
! Input:
!   command    -- character(*), optional: Tao command string. If present, getting user 
!                   input from the terminal is bypassed. This is used when interfacing to python.
!                                          
! Output:
!   errcode    -- integer, optional: Return error code: 0 => OK, Not 0 => Err.
!-

subroutine tao_top_level (command, errcode)

use tao_command_mod, dummy => tao_top_level
use tao_get_user_input_mod, only: tao_get_user_input
!use tao_mpi_mod

implicit none

type (tao_super_universe_struct), pointer :: s_ptr  ! For debug purposes
type (tao_universe_struct), pointer :: u

integer, optional :: errcode
integer n_lev

character(*), optional :: command
character(1000) cmd_out
character(16) :: r_name = 'tao_top_level'

logical found, err, need_input, err_is_fatal

! init

s_ptr => s       ! Used for debugging
if (present(errcode)) errcode = 1   ! Assume error

! MPI: Turn off plotting for slaves
! if (.not. s%mpi%master)  s%global%plot_on = .false.

! And init everything.

if (.not. s%initialized) then
  call tao_parse_command_args (err, command)
  if (err) return
  call tao_init (err)
  if (err) return
  s%initialized = .true.
  n_lev = s%com%cmd_file_level
  need_input = (s%com%saved_cmd_line == '' .and. (n_lev == 0 .or. s%com%cmd_file(n_lev)%paused) .and. &
                                                                                  .not. s%com%single_mode)
  if (present(command) .and. need_input) then
    if (present(errcode)) errcode = 0
    return
  endif
endif

u => s%u(1)  ! Used for debugging

! MPI: slave settings
! if (.not. s%mpi%master) then 
!   ! Turn off screen output
!   s%com%print_to_terminal = .false.
!   call output_direct( print_and_capture = .false.)
! endif

! Command loop

do
  err = .false.
  
  call tao_get_user_input (cmd_out, cmd_in = command)

  ! MPI: Broadcast command to slaves
  ! if (s%mpi%master) call tao_get_user_input (cmd_out)
  ! if (s%mpi%on) call tao_broadcast_chars_mpi(cmd_out)
  
  if (s%com%single_mode) then
    ! single mode
    call tao_single_mode (cmd_out(1:1))
    ! Do the standard calculations and plotting after command execution.
    call tao_cmd_end_calc ()
  else
    ! Command line mode
    call tao_hook_command (cmd_out, found)
    if (.not. found) call tao_command (cmd_out, err, err_is_fatal)
    if (err_is_fatal) exit
    if (err) call tao_abort_command_file()
    call tao_cmd_history_record (cmd_out)
  endif

  ! Exit if current command line parsing is finished (may not be if multiple commands were present) and 
  ! Tao is getting commands through the command argument
  n_lev = s%com%cmd_file_level
  need_input = (s%com%saved_cmd_line == '' .and. (n_lev == 0 .or. s%com%cmd_file(n_lev)%paused) .and. &
                                                                                  .not. s%com%single_mode)
  if (present(command) .and. need_input) exit
enddo

if (present(errcode)) then
  if (err) then
    errcode = 1
  else
    errcode = 0
  endif
endif

end subroutine

