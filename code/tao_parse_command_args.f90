!+
! Subroutine tao_parse_command_args (error, cmd_line)
!
! Subroutine to parse the command line arguments.
!
! Input:
!   cmd_words(:) -- Character(*), optional: If present then this is used
!                    in place of the command line.
! Output:
!   error -- Logical: Set True if there is an error. False otherwise.
!-

subroutine tao_parse_command_args (error, cmd_line)

use tao_interface, dummy => tao_parse_command_args
use tao_command_mod, only: tao_cmd_split

implicit none

character(*), optional :: cmd_line
character(200) :: cmd_words(12)
character(80) arg0, base, switch
character(24) :: r_name = 'tao_parse_command_args'

integer n_arg, i_arg, ix
logical error

! Get command line input

error = .false.

call tao_hook_parse_command_args()
if (.not. s%com%parse_cmd_args) return

if (present(cmd_line)) then
  call tao_cmd_split(cmd_line, 12, cmd_words, .false., error)
  if (error) return
  n_arg = size(cmd_words)
  if (cmd_words(1) == '') return
else
  n_arg = cesr_iargc()
  if (n_arg == 0) return
endif

! loop over all arguments

i_arg = 0

do 

  if (i_arg == n_arg) exit
  call get_next_arg (arg0)

  call match_word (arg0, [character(41):: '-?', '-init', '-noinit', '-beam_all', '-beam0', &
        '-noplot', '-lat', '-log_startup', '-beam', '-var', '-data', '-building_wall', '-plot', &
        '-startup', 'help', '-help', '?', '-geometry', '-rf_on', '-debug', '-disable_smooth_line_calc', &
        '-color_prompt', '-no_stopping', '-hook_init_file', '-beam_position0', '-silent_run', &
        '-beam_init_file_name'], &
              ix, .true., matched_name=switch)

  select case (switch)

  case ('-init')
    call get_next_arg (s%com%init_tao_file)
    s%com%init_tao_file_arg_set = .true.
    if (s%com%init_tao_file == '') then
      call out_io (s_fatal$, r_name, 'NO TAO INIT FILE NAME ON COMMAND LINE.')
      call err_exit
    endif
    ix = SplitFileName(s%com%init_tao_file, s%com%init_tao_file_path, base)

  case ('-beam')
    call get_next_arg (s%com%beam_file)

  case ('-beam_all')
    call get_next_arg (s%com%beam_all_file)

  case ('-beam_position0', '-beam0', '-beam_init_file_name')
    call get_next_arg (s%com%beam_init_file_name)
    if (switch /= '-beam_init_file_name') then
      call out_io (s_warn$, r_name, 'Note: Switch: ' // quote(switch) // &
                                ' has been replaced with "-beam_init_file_name"')
    endif

  case ('-building_wall')
    call get_next_arg (s%com%building_wall_file)

  case ('-color_prompt')
    s%global%prompt_color = 'BLUE'

  case ('-data')
    call get_next_arg (s%com%data_file)

  case ('-disable_smooth_line_calc')
    s%global%disable_smooth_line_calc = .true.

  case ('-debug')
    s%global%debug_on = .true.
    s%global%stop_on_error = .false.

  case ('-geometry')
    call get_next_arg (s%com%plot_geometry)

  case ('help', '-help', '?', '-?')
    call tao_print_command_line_info
    stop

  case ('-hook_init_file')
    call get_next_arg (s%com%hook_init_file)

  case ('-lat')
    call get_next_arg (s%com%lat_file)

  case ('-log_startup')
    s%com%log_startup = .true.

  case ('-no_stopping')
    s%global%stop_on_error = .false.

  case ('-noinit')
    s%com%init_tao_file = ''

  case ('-noplot')
    s%com%noplot_arg_set = .true.

  case ('-rf_on')
    s%global%rf_on = .true.

  case ('-plot')
    call get_next_arg (s%com%plot_file)

  case ('-silent_run')
    call tao_silent_run_set(.true.)

  case ('-startup')
    call get_next_arg (s%com%startup_file)

  case ('-var')
    call get_next_arg (s%com%var_file)

  case default
    call out_io (s_error$, r_name, 'BAD COMMAND LINE ARGUMENT: ' // arg0)
    call tao_print_command_line_info
    error = .true.
    if (s%global%stop_on_error) stop
    return
  end select

enddo

!-----------------------------
contains

subroutine get_next_arg(arg)

character(*) arg

!

if (i_arg == n_arg) then
  call out_io (s_error$, r_name, 'MISSING COMMAND LINE ARGUMENT FOR: ' // arg0)
  error = .true.
  return
endif

i_arg = i_arg + 1

if (present(cmd_line)) then
  arg = cmd_words(i_arg)
  if (cmd_words(i_arg+1) == '') n_arg = i_arg
else
  call cesr_getarg(i_arg, arg)
endif

end subroutine get_next_arg

end subroutine tao_parse_command_args

