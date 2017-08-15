!+
! Subroutine tao_command (command_line, err)
!
! Interface to all standard (non hook) tao commands. 
! This routine esentially breaks the command line into words
! and then calls the appropriate routine.
! Commands are case sensitive.
!
! Input:
!   command_line -- character(*): command line
!
!  Output:
!-

subroutine tao_command (command_line, err)

use tao_change_mod
use tao_command_mod
use tao_data_and_eval_mod
use tao_dmerit_mod
use tao_misalign_mod
use tao_cut_ring_mod
use tao_plot_window_mod
use tao_scale_mod
use tao_set_mod
use tao_show_mod
use tao_wave_mod
use tao_x_scale_mod

! MPI use tao_mpi_mod

implicit none

integer i, j, iu, ios, n_word, n_eq
integer ix, ix_line, ix_cmd, which
integer int1, int2, uni, wrt, n_level

real(rp) value1, value2, time

character(*) :: command_line
character(300) cmd_line
character(20) :: r_name = 'tao_command'
character(300) :: cmd_word(12)
character(40) gang_str, switch, word
character(16) cmd_name, set_word, axis_name

character(16) :: cmd_names(37) = [ &
    'quit         ', 'exit         ', 'show         ', 'plot         ', 'place        ', &
    'clip         ', 'scale        ', 'veto         ', 'use          ', 'restore      ', &
    'run_optimizer', 'flatten      ', 'change       ', 'set          ', 'cut_ring     ', &
    'call         ', 'ptc          ', 'alias        ', 'help         ', 'single_mode  ', &
    're_execute   ', 'reinitialize ', 'x_scale      ', 'x_axis       ', 'derivative   ', &
    'spawn        ', 'xy_scale     ', 'read         ', 'misalign     ', 'end_file     ', &
    'pause        ', 'continue     ', 'wave         ', 'timer        ', 'write        ', &
    'python       ', 'quiet        ']

character(16) :: cmd_names_old(6) = [&
    'x-scale      ', 'xy-scale     ', 'single-mode  ', 'x-axis       ', 'end-file     ', &
    'output       ']

character(16) :: set_names(21) = [&
    'data         ', 'var          ', 'lattice      ', 'global       ', 'plot_page    ', &
    'universe     ', 'curve        ', 'graph        ', 'beam_init    ', 'wave         ', &
    'plot         ', 'bmad_com     ', 'element      ', 'opti_de_param', 'ran_state    ', &
    'csr_param    ', 'floor_plan   ', 'lat_layout   ', 'geodesic_lm  ', 'default      ', &
    'key          ']

logical quit_tao, err, silent, gang, abort

! pause if single stepping

if (s%global%single_step) then
  n_level = s%com%cmd_file_level
  s%com%cmd_file(n_level)%paused = .true.
endif

! blank line => nothing to do

call string_trim (command_line, cmd_line, ix_line)
if (ix_line == 0 .or. cmd_line(1:1) == '!') return

! '/' denotes an option so put a space before it so it does not look like part of the command.

ix = index(cmd_line(1:ix_line), '/')
if (ix /= 0) then
  cmd_line = cmd_line(1:ix-1) // ' ' // trim(cmd_line(ix:))
  ix_line = ix - 1
endif

! strip the command line of comments

ix = index(cmd_line, '!')
if (ix /= 0) cmd_line = cmd_line(:ix-1)        ! strip off comments

! match first word to a command name

call match_word (cmd_line, cmd_names, ix_cmd, .true., matched_name = cmd_name)

if (ix_cmd == 0) then  ! Accept old-style names with "-" instead of "_".
  call match_word (cmd_line, cmd_names_old, ix_cmd, .true., matched_name = cmd_name)
  ix = index(cmd_name, '-')
  if (ix /= 0) cmd_name(ix:ix) = '_'
  if (cmd_name == 'output') cmd_name = 'write'
endif

if (ix_cmd == 0) then
  call out_io (s_error$, r_name, 'UNRECOGNIZED COMMAND: ' // cmd_line)
  call tao_abort_command_file()
  return
elseif (ix_cmd < 0) then
  call out_io (s_error$, r_name, 'AMBIGUOUS COMMAND')
  call tao_abort_command_file()
  return
endif

! Strip off command name from cmd_line and select the appropriate command.

call string_trim (cmd_line(ix_line+1:), cmd_line, ix_line)

select case (cmd_name)

!--------------------------------
! ALIAS

case ('alias')

  call tao_cmd_split(cmd_line, 2, cmd_word, .false., err); if (err) return
  call tao_alias_cmd (cmd_word(1), cmd_word(2))
  return

!--------------------------------
! CALL

case ('call')

  call tao_cmd_split(cmd_line, 10, cmd_word, .true., err); if (err) return
  call tao_call_cmd (cmd_word(1), cmd_word(2:10))
  return

!--------------------------------
! CHANGE

case ('change')

  call tao_cmd_split (cmd_line, 2, cmd_word, .false., err); if (err) return

  silent = .false.
  if (index('-silent', trim(cmd_word(1))) == 1) then
    silent = .true.
    call tao_cmd_split (cmd_word(2), 2, cmd_word, .false., err); if (err) return
  endif

  if (index ('variable', trim(cmd_word(1))) == 1) then
    call tao_cmd_split (cmd_word(2), 2, cmd_word, .false., err); if (err) return
    call tao_change_var (cmd_word(1), cmd_word(2), silent)
  elseif (index('element', trim(cmd_word(1))) == 1) then
    call tao_cmd_split (cmd_word(2), 3, cmd_word, .false., err); if (err) return
    call tao_change_ele (cmd_word(1), cmd_word(2), cmd_word(3))
  elseif (index(trim(cmd_word(1)), 'beam_start') /= 0) then     ! Could be "2@beam_start"
    word = cmd_word(1)
    call tao_cmd_split (cmd_word(2), 2, cmd_word, .false., err); if (err) return
    call tao_change_ele (word, cmd_word(1), cmd_word(2))
  else
    call out_io (s_error$, r_name, 'Change who? (should be: "element", "beam_start", or "variable")')
  endif


!--------------------------------
! CLIP

case ('clip')

  call tao_cmd_split (cmd_line, 4, cmd_word, .true., err); if (err) return

  gang = .false.
  if (index('-gang', trim(cmd_word(1))) == 1 .and. len_trim(cmd_word(1)) > 1) then
    gang = .true.
    cmd_word(1:3) = cmd_word(2:4)
  endif

  if (cmd_word(2) == ' ') then
    call tao_clip_cmd (gang, cmd_word(1), 0.0_rp, 0.0_rp) 
  else
    call tao_to_real (cmd_word(2), value1, err);  if (err) return
    if (cmd_word(3) /= ' ') then
      call tao_to_real (cmd_word(3), value2, err);  if (err) return
    else
      value2 = value1
      value1 = -value1
    endif
    call tao_clip_cmd (gang, cmd_word(1), value1, value2)
  endif

!--------------------------------
! CONTINUE

case ('continue')

  n_level = s%com%cmd_file_level
  if (s%com%cmd_file(n_level)%paused) then
    s%com%cmd_file(n_level)%paused = .false.
  else
    call out_io (s_error$, r_name, 'NO PAUSED COMMAND FILE HERE.')
  endif

  return

!--------------------------------
! CUT_RING

case ('cut_ring')

  call tao_cut_ring ()
  return

!--------------------------------
! DERIVATIVE

case ('derivative')

  call tao_dmodel_dvar_calc(.true.)
  call out_io (s_blank$, r_name, 'Derivative calculated')

  return

!--------------------------------
! END_FILE

case ('end_file')

  n_level = s%com%cmd_file_level
  if (n_level == 0) then
    call out_io (s_error$, r_name, 'END_FILE COMMAND ONLY ALLOWED IN A COMMAND FILE!')
    return
  endif

  close (s%com%cmd_file(n_level)%ix_unit)
  s%com%cmd_file(n_level)%ix_unit = 0 
  s%com%cmd_file_level = n_level - 1 ! signal that the file has been closed

  if (s%com%cmd_file(n_level-1)%paused) then
    call out_io (s_info$, r_name, 'To continue the paused command file type "continue".')
  endif

  return

!--------------------------------
! EXIT/QUIT

case ('exit', 'quit')

  call string_trim (command_line, cmd_line, ix)
  if (ix < 3) then
    call out_io (s_error$, r_name, &
            'SAFETY FEATURE: YOU NEED TO TYPE AT LEAST THREE CHARACTERS TO QUIT.')
    return
  endif

  if (s%global%plot_on) call tao_destroy_plot_window
  call out_io (s_dinfo$, r_name, "Stopping.")
  !MPI !Finalize MPI if it is on
  !MPI if (s%mpi%on) call tao_mpi_finalize()
  stop
 
!--------------------------------
! HELP

case ('help')

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
  call tao_help (cmd_word(1), cmd_word(2))
  return

!--------------------------------
! MISALIGN

case ('misalign')

  call tao_cmd_split (cmd_line, 5, cmd_word, .true., err); if (err) return
  call tao_misalign (cmd_word(1), cmd_word(2), cmd_word(3), cmd_word(4), &
                     cmd_word(5))


!--------------------------------
! PAUSE

case ('pause')

  time = 0
  call tao_cmd_split (cmd_line, 1, cmd_word, .true., err); if (err) return
  if (cmd_word(1) /= '') then
    read (cmd_word(1), *, iostat = ios) time
    if (ios /= 0) then
      call out_io (s_error$, r_name, 'TIME IS NOT A NUMBER.')
      return
    endif
  endif

  call tao_pause_cmd (time)
  return

!--------------------------------
! PLACE

case ('place')

  call tao_cmd_split (cmd_line, 3, cmd_word, .true., err); if (err) return
  if (cmd_word(3) /= ' ') then
    call out_io (s_error$, r_name, 'NOT RECOGNIZED: ' // cmd_word(3))
    return
  endif

  call tao_place_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! PLOT

case ('plot')

  call tao_cmd_split (cmd_line, 2, cmd_word, .false., err)
  call tao_plot_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! PTC

case ('ptc')

  call tao_cmd_split (cmd_line, 2, cmd_word, .false., err); if (err) return

  call tao_ptc_cmd (cmd_word(1), cmd_word(2))
  return

!--------------------------------
! PYTHON

case ('python')

  call tao_python_cmd (cmd_line)
  return

!--------------------------------
! QUIET PYTHON

case ('quiet')

if (s%com%cmd_file_level == 0) then 
  call out_io (s_error$, r_name, 'The "quiet" command may only be used in command files.')
else
  s%global%quiet = .true.
endif

!--------------------------------
! RE_EXECUTE

case ('re_execute')

  call tao_re_execute (cmd_line, err)
  return

!--------------------------------
! READ

case ('read')

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
  call tao_read_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! RESTORE, USE, VETO

case ('restore', 'use', 'veto')

  call tao_cmd_split(cmd_line, 2, cmd_word, .true., err)
  if (err) return
  
  call match_word (cmd_word(1), ["data    ", "variable"], which, .true., matched_name = switch)
  
  if (switch == 'data') then
    call tao_use_data (cmd_name, cmd_word(2))
  elseif (switch == 'variable') then
    call tao_use_var (cmd_name, cmd_word(2))
  else
    call out_io (s_error$, r_name, "Use/veto/restore what? data or variable?")
    return
  endif 

!--------------------------------
! REINITIALIZE

case ('reinitialize')

  call tao_cmd_split(cmd_line, 10, cmd_word, .false., err)
  if (err) return

  call match_word (cmd_word(1), ['data', 'tao ', 'beam'], ix, .true., matched_name=word)

  select case (word)

  case ('beam') 
    do i = lbound(s%u, 1), ubound(s%u, 1)
      s%u(i)%beam%init_beam0 = .true.
      s%u(i)%calc%lattice = .true.
    enddo

  case ('data') 
    s%u(:)%calc%lattice = .true.

  case ('tao') 

    call tao_parse_command_args (err, cmd_word(2:))
    if (err) return

    ! quit the plot window so it will be recreated    

    if (s%global%plot_on) then
      call tao_destroy_plot_window
      s%global%init_plot_needed = .true.
    endif
      
    ! quit the plot window so it will be recreated    
    call tao_destroy_plot_window
    s%global%init_plot_needed = .true.
    
    if (s%com%init_tao_file /= '') call out_io (s_info$, r_name, 'Reinitializing with: ' // s%com%init_tao_file)
    call tao_init (err)
    return

  case default
    call out_io (s_error$, r_name, 'Reinit what? Choices are: "beam", "data", or "tao".')
    return
    
  end select

!--------------------------------
! RUN, FLATTEN

case ('run_optimizer', 'flatten')

  call tao_cmd_split (cmd_line, 1, cmd_word, .true., err); if (err) return
  call tao_run_cmd (cmd_word(1), abort)

!--------------------------------
! SCALE

case ('scale')

  call tao_cmd_split (cmd_line, 5, cmd_word, .true., err); if (err) return

  axis_name = ''
  gang_str = ''

  do 
    word = cmd_word(1)
    if (word(1:1) /= '-') exit
    call match_word (word, ['-y     ', '-y2    ', '-nogang', '-gang  '], ix, .true., matched_name=switch)

    select case (switch)
    case ('-y', '-y2') 
     axis_name = switch(2:)
    case ('-gang', '-nogang')
      gang_str = switch(2:)
    case default
      call out_io (s_error$, r_name, 'BAD SWITCH: ' // switch)
      return
    end select
    cmd_word(1:4) = cmd_word(2:5)
  enddo

  if (cmd_word(2) == ' ') then
    call tao_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, axis_name, gang_str) 
  else
    call tao_to_real (cmd_word(2), value1, err);  if (err) return
    if (cmd_word(3) /= ' ') then
      call tao_to_real (cmd_word(3), value2, err);  if (err) return
    else
      value2 = value1
      value1 = -value1
    endif
    call tao_scale_cmd (cmd_word(1), value1, value2, axis_name, gang_str)
  endif

!--------------------------------
! SET

case ('set')

  call tao_cmd_split (cmd_line, 2, cmd_word, .false., err, '=')

  call match_word (cmd_word(1), set_names, ix, .true., matched_name=set_word)
  if (ix < 1) then
    call out_io (s_error$, r_name, 'NOT RECOGNIZED OR AMBIGUOUS: ' // cmd_word(1))
    return
  endif

  cmd_line = cmd_word(2)
  select case (set_word)
  case ('ran_state'); n_word = 2; n_eq = 1
  case ('beam_init', 'bmad_com', 'csr_param', 'data', 'global', 'lattice', 'default', &
        'opti_de_param', 'var', 'wave', 'floor_plan', 'lat_layout', 'geodesic_lm', 'key'); n_word = 3; n_eq = 2
  case ('universe'); n_word = 3; n_eq = 10
  case ('plot_page'); n_word = 4; n_eq = 2
  case ('curve', 'element', 'graph', 'plot'); n_word = 4; n_eq = 3
  end select

  call tao_cmd_split (cmd_line, n_word, cmd_word, .false., err, '=')

  if (set_word /= 'universe' .and. cmd_word(n_eq) /= '=') then
    call out_io (s_error$, r_name, 'SYNTAX PROBLEM. "=" NOT IN CORRECT PLACE.')
    return
  endif

  select case (set_word)
  case ('beam_init')
    call tao_set_beam_init_cmd (cmd_word(1), cmd_word(3))
  case ('bmad_com')
    call tao_set_bmad_com_cmd (cmd_word(1), cmd_word(3))
  case ('csr_param')
    call tao_set_csr_param_cmd (cmd_word(1), cmd_word(3))
  case ('curve')
    call tao_set_curve_cmd (cmd_word(1), cmd_word(2), cmd_word(4)) 
  case ('data')
    call tao_set_data_cmd (cmd_word(1), cmd_word(3))
  case ('default')
    call tao_set_default_cmd (cmd_word(1), cmd_word(3))
  case ('element')
    call tao_set_elements_cmd (cmd_word(1), cmd_word(2), cmd_word(4))
  case ('geodesic_lm')
    call tao_set_geodesic_lm_cmd (cmd_word(1), cmd_word(3))
  case ('global')
    call tao_set_global_cmd (cmd_word(1), cmd_word(3))
  case ('graph')
    call tao_set_graph_cmd (cmd_word(1), cmd_word(2), cmd_word(4))
  case ('key')
    call tao_set_key_cmd (cmd_word(1), cmd_word(3))    
  case ('lattice')
    call tao_set_lattice_cmd (cmd_word(1), cmd_word(3))
  case ('opti_de_param')
    call tao_set_opti_de_param_cmd (cmd_word(1), cmd_word(3))
  case ('plot ')
    call tao_set_plot_cmd (cmd_word(1), cmd_word(2), cmd_word(4))
  case ('plot_page')
    call tao_set_plot_page_cmd (cmd_word(1), cmd_word(3), cmd_word(4))
  case ('ran_state')
    call tao_set_ran_state_cmd (cmd_word(2))
  case ('floor_plan')
    call tao_set_drawing_cmd (s%plot_page%floor_plan, cmd_word(1), cmd_word(3))
  case ('lat_layout')
    call tao_set_drawing_cmd (s%plot_page%lat_layout, cmd_word(1), cmd_word(3))
  case ('universe')    
    call tao_set_universe_cmd (cmd_word(1), cmd_word(2), cmd_word(3))
  case ('var')
    call tao_set_var_cmd (cmd_word(1), cmd_word(3))
  case ('wave')
    call tao_set_wave_cmd (cmd_word(1), cmd_word(3), err)
    if (err) return
    call tao_cmd_end_calc
    call tao_show_cmd ('wave', '')
  end select

!--------------------------------
! SHOW

case ('show')

  call tao_cmd_split (cmd_line, 2, cmd_word, .false., err); if (err) return
  if (cmd_word(1) == '-file') then
    cmd_line = cmd_word(2)
    call tao_cmd_split (cmd_line, 3, cmd_word, .false., err); if (err) return
    iu = lunget()
    open (iu, file = cmd_word(1))
    call output_direct (iu)
    call tao_show_cmd (cmd_word(2), cmd_word(3))
    call output_direct(-1)
    close (iu)
    call out_io (s_blank$, r_name, 'Output written to file: ' // cmd_word(1))
    return
  endif

  call tao_show_cmd (cmd_word(1), cmd_word(2))
  return

!--------------------------------
! SINGLE-MODE

case ('single_mode')

  s%com%single_mode = .true.
  call out_io (s_blank$, r_name, 'Entering Single Mode...')
  return

!--------------------------------
! SPAWN

case ('spawn')

  call system_command (cmd_line)
  return

!--------------------------------
! timer

case ('timer')

  call tao_timer (cmd_line)
  return

!--------------------------------
! WAVE

case ('wave')

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
  call tao_wave_cmd (cmd_word(1), cmd_word(2), err); if (err) return
  call tao_cmd_end_calc
  call tao_show_cmd ('wave', '')
  return

!--------------------------------
! write

case ('write')

  call tao_write_cmd (cmd_line)
  return

!--------------------------------
! X_AXIS

case ('x_axis')

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
  call tao_x_axis_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! X_SCALE

case ('x_scale')

  call tao_cmd_split (cmd_line, 5, cmd_word, .true., err); if (err) return

  gang_str = ''
  do 
    if (cmd_word(1)(1:1) /= '-') exit
    call match_word (cmd_word(1), ['-nogang', '-gang  '], ix, .true., matched_name=switch)

    select case (switch)
    case ('-gang', '-nogang')
      gang_str = switch(2:)
    case default
      call out_io (s_error$, r_name, 'BAD SWITCH: ' // switch)
      return
    end select
    cmd_word(1:4) = cmd_word(2:5)
  enddo

  if (cmd_word(2) == ' ') then
    call tao_x_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, err, gang_str)
  else
    call tao_to_real (cmd_word(2), value1, err); if (err) return
    call tao_to_real (cmd_word(3), value2, err); if (err) return
    call tao_x_scale_cmd (cmd_word(1), value1, value2, err, gang_str)
  endif

!--------------------------------
! XY_SCALE

case ('xy_scale')

  call tao_cmd_split (cmd_line, 3, cmd_word, .true., err); if (err) return

  if (cmd_word(2) == ' ') then
    call tao_x_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, err)
    call tao_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp) 
  else
    call tao_to_real (cmd_word(2), value1, err);  if (err) return
    if (cmd_word(3) /= ' ') then
      call tao_to_real (cmd_word(3), value2, err);  if (err) return
    else
      value2 = value1
      value1 = -value1
    endif
    call tao_x_scale_cmd (cmd_word(1), value1, value2, err)
    call tao_scale_cmd (cmd_word(1), value1, value2)
  endif

!--------------------------------
! DEFAULT

case default

  call out_io (s_error$, r_name, 'INTERNAL COMMAND PARSING ERROR!')
  call err_exit

end select

!------------------------------------------------------------------------
! Do the standard calculations and plotting after command
! Note: wave command bypasses this.

call tao_cmd_end_calc

end subroutine tao_command




