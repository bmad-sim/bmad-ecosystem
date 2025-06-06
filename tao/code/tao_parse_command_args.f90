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
use tao_set_mod, only: tao_set_openmp_n_threads

implicit none

type (tao_common_struct) com_saved

character(*), optional :: cmd_line
character(200) :: cmd_words(12)
character(80) arg0, arg1, base, switch
character(*), parameter :: r_name = 'tao_parse_command_args'

integer n_arg, i_arg, ix, old_n_threads
logical error, negate

! Init global and common structs.
! If s%initialized = True then this is a reinit so save some param values in s%com.

if (s%initialized) then
  com_saved = s%com
  s%com = tao_common0
  s%com%saved_cmd_line = com_saved%saved_cmd_line
  s%com%ix_history     = com_saved%ix_history
  s%com%n_history      = com_saved%n_history
  s%com%cmd_file       = com_saved%cmd_file
  s%com%cmd_file_level = com_saved%cmd_file_level
  s%com%lev_loop       = com_saved%lev_loop
  call move_alloc(com_saved%do_loop, s%com%do_loop)
else
  s%com = tao_common0
endif
  
s%global  = tao_global_struct()
bmad_com  = bmad_common_struct()
space_charge_com = space_charge_common_struct()

! Get command line input

error = .false.
s%com%command_arg_has_been_executed = .false.

if (associated(tao_hook_parse_command_args_ptr)) call tao_hook_parse_command_args_ptr()
if (.not. s%init%parse_cmd_args) return

if (present(cmd_line)) then
  call tao_cmd_split(cmd_line, 12, cmd_words, .false., error)
  if (error) return
  n_arg = size(cmd_words)
  if (cmd_words(1) == '') return
else
  n_arg = command_argument_count()
  if (n_arg == 0) return
endif

! loop over all arguments

i_arg = 0

do 

  if (i_arg == n_arg) exit
  call get_next_arg ('', arg0, i_arg, n_arg)
  if (arg0(1:2) == '--') then
    arg1 = arg0(2:)
    negate = .true.
  else
    arg1 = arg0
    negate = .false.
  endif

  call match_word (arg1, [character(44):: '-?', '-init_file', '-noinit', '-nostartup', '-beam0', '-clear', &
        '-noplot', '-lattice_file', '-log_startup', '-beam_file', '-var_file', '-data_file', &
        '-startup_file', 'help', '-help', '?', '-geometry', '-rf_on', '-debug', '-disable_smooth_line_calc', &
        '-color_prompt', '-no_stopping', '-hook_init_file', '-beam_position0', '-silent_run', &
        '-beam_init_file_name', '-slice_lattice', '-start_branch_at', '-prompt_color', '-beam_init_position_file', &
        '-plot_file', '-external_plotting', '-quiet', '-no_rad_int', '-command', &
        '-symbol_import', '-building_wall_file', '-reverse'], ix, .true., matched_name=switch)

  if (negate) switch = '-' // switch

  !

  select case (switch)

  case ('-beam_file')
    call get_next_arg (arg0, s%init%beam_file_arg, i_arg, n_arg)

  case ('-beam_position0', '-beam0', '-beam_init_file_name')
    call get_next_arg (arg0, s%init%beam_init_position_file_arg, i_arg, n_arg)
    call out_io (s_warn$, r_name, 'Note: Switch: ' // quote(switch) // &
                                ' has been replaced with "-beam_init_position_file"')

  case ('-beam_init_position_file')
    call get_next_arg (arg0, s%init%beam_init_position_file_arg, i_arg, n_arg)

  case ('-building_wall_file')
    call get_next_arg (arg0, s%init%building_wall_file_arg, i_arg, n_arg)

  case ('-clear')
    s%init = tao_init_struct()
    s%com = tao_common0
    old_n_threads = s%global%n_threads
    s%global = tao_global_struct()
    s%global%n_threads = old_n_threads

  case ('-command')
    call get_next_arg (arg0, s%init%command_arg, i_arg, n_arg, .true.)

  case ('-data_file')
    call get_next_arg (arg0, s%init%data_file_arg, i_arg, n_arg)

  case ('-disable_smooth_line_calc')
    s%init%disable_smooth_line_calc_arg = present_str

  case ('-debug')
    s%init%debug_switch = .true.

  case ('-external_plotting')
    s%init%external_plotting_switch = .true.

  case ('-geometry')
    call get_next_arg (arg0, s%init%geometry_arg, i_arg, n_arg, .true.)

  case ('help', '-help', '?', '-?')
    call tao_print_command_line_info
    stop

  case ('-hook_init_file')
    call get_next_arg (arg0, s%init%hook_init_file_arg, i_arg, n_arg)

  case ('-init_file')
    call get_next_arg (arg0, s%init%init_file_arg, i_arg, n_arg)
    ix = SplitFileName(s%init%init_file_arg, s%init%init_file_arg_path, base)

  case ('-lattice_file')
    call get_next_arg (arg0, s%init%lattice_file_arg, i_arg, n_arg)
    s%init%noinit_arg = ''

  case ('-log_startup')
    s%init%log_startup_arg = present_str

  case ('-no_stopping')
    s%init%no_stopping_arg = present_str

  case ('-noinit')
    s%init%noinit_arg = present_str
    s%init%init_file_arg = ''

  case ('-noplot')
    s%init%noplot_arg = present_str

  case ('-nostartup')
    s%init%nostartup_arg = present_str

  case ('-no_rad_int')
    s%init%no_rad_int_arg = present_str

  case ('-plot_file')
    call get_next_arg (arg0, s%init%plot_file_arg, i_arg, n_arg)

  case ('-prompt_color', '-color_prompt')
    if (i_arg < n_arg) then
      call get_next_arg (arg0, s%init%prompt_color_arg, i_arg, n_arg)
      call upcase_string(s%init%prompt_color_arg)
      if (s%init%prompt_color_arg(1:1) == '-') then
        i_arg = i_arg - 1
        s%init%prompt_color_arg = 'BLUE'
      endif
    endif

  case ('-quiet', '-silent_run')       ! "-silent_run" is old syntax
    if (i_arg < n_arg) then
      call get_next_arg (arg0, s%init%quiet_arg, i_arg, n_arg)
      if (s%init%quiet_arg(1:1) == '-') then
        s%init%quiet_arg = 'all'
        i_arg = i_arg - 1
      endif
    else
      s%init%quiet_arg = 'all'
    endif

  case ('-reverse')
    s%init%reverse_arg = present_str

  case ('-rf_on')
    s%init%rf_on_arg = present_str

  case ('-slice_lattice')
    call get_next_arg (arg0, s%init%slice_lattice_arg, i_arg, n_arg, .true.)

  case ('-start_branch_at')
    call get_next_arg (arg0, s%init%start_branch_at_arg, i_arg, n_arg)

  case ('-startup_file')
    call get_next_arg (arg0, s%init%startup_file_arg, i_arg, n_arg)

  case ('-symbol_import')
    s%init%symbol_import_arg = present_str

  case ('-var_file')
    call get_next_arg (arg0, s%init%var_file_arg, i_arg, n_arg)

  case default
    if (ix == 0) then
      call out_io (s_error$, r_name, 'BAD COMMAND LINE ARGUMENT: ' // arg0)
    else
      call out_io (s_error$, r_name, 'THIS MATCHES TO MULTIPLE COMMAND LINE ARGUMENTS: ' // arg0)
    endif
    call tao_print_command_line_info
    error = .true.
    return

  ! Negate cases

  case ('--beam_file');                           s%init%beam_file_arg = ''
  case ('--beam_position0', '--beam0', '--beam_init_file_name'); 
                                                  s%init%beam_init_position_file_arg = ''
  case ('--beam_init_position_file');             s%init%beam_init_position_file_arg = ''
  case ('--building_wall_file');                  s%init%building_wall_file_arg = ''
  case ('--command');                             s%init%command_arg = ''
  case ('--data_file');                           s%init%data_file_arg = ''
  case ('--disable_smooth_line_calc');            s%init%disable_smooth_line_calc_arg = negated_str
  case ('--debug');                               s%init%debug_switch = .false.
  case ('--external_plotting');                   s%init%external_plotting_switch = .false.
  case ('--geometry');                            s%init%geometry_arg = ''
  case ('--hook_init_file');                      s%init%hook_init_file_arg = ''
  case ('--init_file');                           s%init%init_file_arg = ''; s%init%init_file_arg_path = ''
  case ('--lattice_file');                        s%init%lattice_file_arg = ''
  case ('--log_startup');                         s%init%log_startup_arg = negated_str
  case ('--no_stopping');                         s%init%no_stopping_arg = negated_str
  case ('--noinit');                              s%init%noinit_arg = ''
  case ('--noplot');                              s%init%noplot_arg = negated_str
  case ('--nostartup');                           s%init%nostartup_arg = ''
  case ('--no_rad_int');                          s%init%no_rad_int_arg = negated_str
  case ('--plot_file');                           s%init%plot_file_arg = ''
  case ('--prompt_color', '--color_prompt');      s%init%prompt_color_arg = 'DEFAULT'  ! read_a_line recognizes this.
  case ('--reverse');                             s%init%reverse_arg = negated_str
  case ('--rf_on');                               s%init%rf_on_arg = negated_str
  case ('--quiet', '--silent_run');               s%init%quiet_arg = 'off'
  case ('--slice_lattice');                       s%init%slice_lattice_arg = ''
  case ('--start_branch_at');                     s%init%start_branch_at_arg = ''
  case ('--startup_file');                        s%init%startup_file_arg = ''
  case ('--symbol_import');                       s%init%symbol_import_arg = ''
  case ('--var_file');                            s%init%var_file_arg = ''
  end select
enddo

!

if (s%init%debug_switch) then
  s%global%debug_on = .true.
  s%global%stop_on_error = .false.
endif

s%global%external_plotting = s%init%external_plotting_switch

!-----------------------------
contains

subroutine get_next_arg(arg0, arg_next, i_arg, n_arg, may_have_blanks)

integer i_arg, n_arg
character(*) arg0, arg_next
character(40) sub
logical, optional :: may_have_blanks

!

if (i_arg == n_arg) then
  call out_io (s_error$, r_name, 'MISSING COMMAND LINE ARGUMENT FOR: ' // arg0)
  error = .true.
  return
endif

i_arg = i_arg + 1

if (present(cmd_line)) then
  arg_next = cmd_words(i_arg)
  if (cmd_words(i_arg+1) == '') n_arg = i_arg
else
  call get_command_argument(i_arg, arg_next)
endif

! may_have_blanks = T means that the arg_next string may contain blank characters.

if (logic_option(.false., may_have_blanks)) then
  do
    if (i_arg == n_arg) return

    if (present(cmd_line)) then
      sub = cmd_words(i_arg+1)
    else
      call get_command_argument(i_arg+1, sub)
    endif

    if (sub(1:1) == '-') return
    arg_next = trim(arg_next) // sub
    i_arg = i_arg + 1

    if (present(cmd_line)) then
      if (cmd_words(i_arg+1) == '') n_arg = i_arg
    endif
  enddo
endif

end subroutine get_next_arg

end subroutine tao_parse_command_args

