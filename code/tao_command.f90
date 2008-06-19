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

#include "CESR_platform.inc"

subroutine tao_command (command_line, err)

use tao_mod
use quick_plot
use tao_command_mod
use tao_dmerit_mod
use tao_scale_mod
use tao_x_scale_mod
use tao_set_mod
use tao_plot_window_mod
use tao_show_mod
use tao_change_mod
use tao_misalign_mod

implicit none

integer i, j, iu
integer ix, ix_line, ix_cmd, which
integer int1, int2, uni, wrt, n_level

real(rp) value1, value2

character(*) :: command_line
character(140) cmd_line
character(20) :: r_name = 'tao_command'
character(80) :: cmd_word(12)
character(16) cmd_name, set_word, axis_name

character(16) :: cmd_names(30) = (/  &
    'quit        ', 'exit        ', 'show        ', 'plot        ', 'place       ', &
    'clip        ', 'scale       ', 'veto        ', 'use         ', 'restore     ', &
    'run         ', 'flatten     ', 'output      ', 'change      ', 'set         ', &
    'call        ', 'view        ', 'alias       ', 'help        ', 'history     ', &
    'single-mode ', 'reinitialize', 'x-scale     ', 'x-axis      ', 'derivative  ', &
    'spawn       ', 'xy-scale    ', 'read        ', 'misalign    ', 'end_file    ' /)

character(16) :: set_names(8) = (/ &
    'data        ', 'var         ', 'lattice     ', 'global      ', 'plot_page   ', &
    'universe    ', 'curve       ', 'graph       ' /)



logical quit_tao, err, silent, gang

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

call match_word (cmd_line, cmd_names, ix_cmd, .true.)
if (ix_cmd == 0) then
  call out_io (s_error$, r_name, 'UNRECOGNIZED COMMAND: ' // cmd_line)
  return
elseif (ix_cmd < 0) then
  call out_io (s_error$, r_name, 'AMBIGUOUS COMMAND')
  return
endif

cmd_name = cmd_names(ix_cmd)

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

!--------------------------------
! CHANGE

case ('change')

  call tao_cmd_split (cmd_line, 2, cmd_word, .false., err); if (err) return

  silent = .false.
  if (index('-silent', trim(cmd_word(1))) == 1) then
    silent = .true.
    call tao_cmd_split (cmd_word(2), 2, cmd_word, .false., err); if (err) return
  endif

  if (cmd_word(1) == 'var') then
    call tao_cmd_split (cmd_word(2), 2, cmd_word, .false., err); if (err) return
    call tao_change_var (cmd_word(1), cmd_word(2), silent)
  elseif (cmd_word(1) == 'ele') then
    call tao_cmd_split (cmd_word(2), 3, cmd_word, .false., err); if (err) return
    call tao_change_ele (cmd_word(1), cmd_word(2), cmd_word(3))
  elseif (cmd_word(1) == 'beam_start') then
    call tao_cmd_split (cmd_word(2), 2, cmd_word, .false., err); if (err) return
    call tao_change_ele ('beam_start', cmd_word(1), cmd_word(2))
  else
    call out_io (s_error$, r_name, &
             'Error: Change who? (should be: "ele", "bunch_start", or "var")')
  endif


!--------------------------------
! CLIP

case ('clip')

  call tao_cmd_split (cmd_line, 4, cmd_word, .true., err); if (err) return

  gang = .false.
  if (index('-gang', cmd_word(1)) == 1 .and. len_trim(cmd_word(1)) > 0) then
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
! DERIVATIVE

case ('derivative')

  call tao_dmodel_dvar_calc(.true.)
  call out_io (s_blank$, r_name, 'Derivative calculated')

!--------------------------------
! END_FILE

case ('end_file')

  n_level = tao_com%cmd_file_level
  if (n_level == 0) then
    call out_io (s_error$, r_name, 'END_FILE COMMAND ONLY ALLOWED IN A COMMAND FILE!')
    return
  endif
  close (tao_com%cmd_file(n_level)%ix_unit)
  tao_com%cmd_file(n_level)%ix_unit = 0 
  tao_com%cmd_file_level = n_level - 1 ! signal that the file has been closed

!--------------------------------
! EXIT/QUIT

case ('exit', 'quit')

  quit_tao = .false.
  call tao_query_logical ('y', 'n', 'Quit?', quit_tao)
  if (.not. quit_tao) return
  if (s%global%plot_on) call tao_destroy_plot_window
  call out_io (s_dinfo$, r_name, "Stopping.")
  stop
 
!--------------------------------
! HELP

case ('help')

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
  call tao_help (cmd_word(1), cmd_word(2))
  return

!--------------------------------
! HISTORY

case ('history')

  call tao_cmd_split (cmd_line, 1, cmd_word, .true., err); if (err) return
  call tao_history_cmd (cmd_word(1), err)
  return

!--------------------------------
! MISALIGN

case ('misalign')

  call tao_cmd_split (cmd_line, 5, cmd_word, .true., err); if (err) return
  call tao_misalign (cmd_word(1), cmd_word(2), cmd_word(3), cmd_word(4), &
                     cmd_word(5))


!--------------------------------
! OUTPUT

case ('output')

  call tao_output_cmd (cmd_line)
  return

!--------------------------------
! PLACE

case ('place')

  if (.not. s%global%plot_on) then
    call out_io (s_error$, r_name, "PLOTTING TURNED OFF!")
    return
  endif

  call tao_cmd_split (cmd_line, 3, cmd_word, .true., err); if (err) return
  if (cmd_word(3) /= ' ') then
    call out_io (s_error$, r_name, 'NOT RECOGNIZED: ' // cmd_word(3))
    return
  endif

  call tao_place_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! PLOT

case ('plot')

  if (.not. s%global%plot_on) then
    call out_io (s_error$, r_name, "PLOTTING TURNED OFF!")
    return
  endif

  call tao_cmd_split (cmd_line, 9, cmd_word, .false., err, '+-')

  i = 1; j = 1
  do 
    j = j + 1
    if (cmd_word(j) == ' ') exit
    i = i + 1
    if (cmd_word(j) == '+' .or. cmd_word(j) == '-') then
      cmd_word(i) = cmd_word(j)(1:1) // cmd_word(j+1)
      j = j + 1
    elseif (i == 2) then
      cmd_word(i) = '+' // cmd_word(j)
    else
      cmd_word(i) = cmd_word(j)
    endif
  enddo

  call tao_plot_cmd (cmd_word(1), cmd_word(2:i))

!--------------------------------
! VETO, RESTORE, USE

case ('use', 'veto', 'restore')

  call tao_cmd_split(cmd_line, 3, cmd_word, .false., err)
  
  call match_word (cmd_word(1), name$%data_or_var, which, .true.)
  
  if (which .eq. data$) then
    call tao_use_data (cmd_name, cmd_word(2))
  elseif (which .eq. variable$) then
    call tao_use_var (cmd_name, cmd_word(2))
  else
    call out_io (s_error$, r_name, "Use/veto/restore what? data or variable?")
    return
  endif 

!--------------------------------
! READ

case ('read')

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
  call tao_read_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! REINITIALIZE

case ('reinitialize')

  call tao_cmd_split(cmd_line, 10, cmd_word, .false., err)
  if (err) return

  select case (cmd_word(1))

  case ('beam') 
    tao_com%init_beam0 = .true.
    tao_com%lattice_recalc = .true.

  case ('data') 
    tao_com%lattice_recalc = .true.

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
    
    call out_io (s_info$, r_name, 'Reinitializing with ' // tao_com%init_tao_file)
    call tao_init ()
    return

  case default
    call out_io (s_error$, r_name, 'REINIT WHAT? CHOICES ARE: "beam" OR "tao".')
    return
    
  end select

!--------------------------------
! RUN, FLATTEN

case ('run', 'flatten')

  call tao_cmd_split (cmd_line, 1, cmd_word, .true., err); if (err) return
  call tao_run_cmd (cmd_word(1))

!--------------------------------
! SET

case ('set')

  call tao_cmd_split (cmd_line, 5, cmd_word, .false., err, '=')

  call match_word (cmd_word(1), set_names, ix, .true.)
  if (ix == 0) then
    call out_io (s_error$, r_name, 'NOT RECOGNIZED: ' // cmd_word(1))
    return
  endif

  set_word = set_names(ix)

  if ( (set_word == 'curve'     .and. cmd_word(4) /= '=') .or. &
       (set_word == 'data'      .and. cmd_word(3) /= '=') .or. &
       (set_word == 'var'       .and. cmd_word(3) /= '=') .or. &
       (set_word == 'global'    .and. cmd_word(3) /= '=') .or. &
       (set_word == 'plot_page' .and. cmd_word(3) /= '=') .or. &
       (set_word == 'graph'     .and. cmd_word(4) /= '=') .or. &
     (set_word == 'lattice'   .and. cmd_word(3) /= '=')) then
    call out_io (s_error$, r_name, 'SYNTAX PROBLEM. "=" NOT IN CORRECT PLACE.')
    return
  endif

  select case (set_word)
  case ('data')
    cmd_word(4) = trim(cmd_word(4)) // cmd_word(5)
    call tao_set_data_cmd (cmd_word(2), cmd_word(4))
  case ('var')
    cmd_word(4) = trim(cmd_word(4)) // cmd_word(5)
    call tao_set_var_cmd (cmd_word(2), cmd_word(4))
  case ('lattice')
    cmd_word(4) = trim(cmd_word(4)) // cmd_word(5)
    call tao_set_lattice_cmd (cmd_word(2), cmd_word(4)) 
  case ('curve')
    call tao_set_curve_cmd (cmd_word(2), cmd_word(3), cmd_word(5)) 
  case ('global')
    call tao_set_global_cmd (cmd_word(2), cmd_word(4))
  case ('plot_page')
    call tao_set_plot_page_cmd (cmd_word(2), cmd_word(4), cmd_word(5))
  case ('graph')
    call tao_set_graph_cmd (cmd_word(2), cmd_word(3), cmd_word(5))
  case ('universe')    
    call tao_set_uni_cmd (cmd_word(2), cmd_word(3))
  end select

!--------------------------------
! SCALE

case ('scale')

  if (.not. s%global%plot_on) then
    call out_io (s_error$, r_name, "PLOTTING TURNED OFF!")
    return
  endif

  call tao_cmd_split (cmd_line, 5, cmd_word, .true., err); if (err) return

  axis_name = ''
  if (cmd_word(1) == '-y' .or. cmd_word(1) == '-y2') then
   axis_name = cmd_word(1)
   axis_name = axis_name(2:)
   cmd_word(1:3) = cmd_word(2:4)
  endif 

  if (cmd_word(2) == ' ') then
    call tao_scale_cmd (cmd_word(1), axis_name, 0.0_rp, 0.0_rp) 
  else
    call tao_to_real (cmd_word(2), value1, err);  if (err) return
    if (cmd_word(3) /= ' ') then
      call tao_to_real (cmd_word(3), value2, err);  if (err) return
    else
      value2 = value1
      value1 = -value1
    endif
    call tao_scale_cmd (cmd_word(1), axis_name, value1, value2)
  endif

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
    call output_direct(0)
    close (iu)
    call out_io (s_blank$, r_name, 'Output written to file: ' // cmd_word(1))
    return
  endif

  call tao_show_cmd (cmd_word(1), cmd_word(2))
  return

!--------------------------------
! SINGLE-MODE

case ('single-mode')

  tao_com%single_mode = .true.
  call out_io (s_blank$, r_name, 'Entering Single Mode...')
  return

!--------------------------------
! SPAWN

case ('spawn')

  call system_command (cmd_line)
  return

!--------------------------------
! VIEW

case ('view')

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
  call tao_to_int (cmd_word(1), int1, err); if (err) return
  call tao_view_cmd (int1)

!--------------------------------
! X-AXIS

case ('x-axis')

  if (.not. s%global%plot_on) then
    call out_io (s_error$, r_name, "PLOTTING TURNED OFF!")
    return
  endif

  call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
  call tao_x_axis_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! X-SCALE

case ('x-scale')

  if (.not. s%global%plot_on) then
    call out_io (s_error$, r_name, "PLOTTING TURNED OFF!")
    return
  endif

  call tao_cmd_split (cmd_line, 4, cmd_word, .true., err); if (err) return
  if (cmd_word(2) == ' ') then
    call tao_x_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, err)
  else
    call tao_to_real (cmd_word(2), value1, err); if (err) return
    call tao_to_real (cmd_word(3), value2, err); if (err) return
    call tao_x_scale_cmd (cmd_word(1), value1, value2, err)
  endif

!--------------------------------
! XY-SCALE

case ('xy-scale')

  if (.not. s%global%plot_on) then
    call out_io (s_error$, r_name, "PLOTTING TURNED OFF!")
    return
  endif

  call tao_cmd_split (cmd_line, 3, cmd_word, .true., err); if (err) return
  if (cmd_word(2) == ' ') then
    call tao_x_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, err)
    call tao_scale_cmd (cmd_word(1), '', 0.0_rp, 0.0_rp) 
  else
    call tao_to_real (cmd_word(2), value1, err);  if (err) return
    if (cmd_word(3) /= ' ') then
      call tao_to_real (cmd_word(3), value2, err);  if (err) return
    else
      value2 = value1
      value1 = -value1
    endif
    call tao_x_scale_cmd (cmd_word(1), value1, value2, err)
    call tao_scale_cmd (cmd_word(1), '', value1, value2)
  endif

!--------------------------------
! DEFAULT

case default

  call out_io (s_error$, r_name, 'INTERNAL COMMAND PARSING ERROR!')
  call err_exit

end select

!------------------------------------------------------------------------
! Do the standard calculations and plotting after command

call tao_cmd_end_calc

end subroutine tao_command




