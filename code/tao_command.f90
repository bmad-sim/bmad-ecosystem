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

  implicit none

  integer i, j
  integer ix, ix_line, ix_cmd, which
  integer int1, int2, uni

  real(rp) value1, value2

  character(*) :: command_line
  character(140) cmd_line
  character(20) :: r_name = 'tao_command'
  character(80) :: cmd_word(12)
  character(16) cmd_name, set_word

  character(16) :: cmd_names(26) = (/  &
    'quit        ', 'exit        ', 'show        ', 'plot        ', 'place       ', &
    'clip        ', 'scale       ', 'veto        ', 'use         ', 'restore     ', &
    'run         ', 'flatten     ', 'output      ', 'change      ', 'set         ', &
    'call        ', 'view        ', 'alias       ', 'help        ', 'history     ', &
    'single-mode ', 'reinitialize', 'x-scale     ', 'x-axis      ', 'derivative  ', &
    'spawn       '/)

  character(16) :: set_names(8) = (/ &
    'data        ', 'var         ', 'lattice     ', 'global      ', 'plot_page   ', &
    'universe    ', 'curve       ', 'graph       ' /)



  logical quit_tao, err

! blank line => nothing to do

  call string_trim (command_line, cmd_line, ix_line)
  if (ix_line == 0 .or. cmd_line(1:1) == '!') return

! strip the command line of comments

  ix = index(cmd_line, '!')
  if (ix /= 0) cmd_line = cmd_line(:ix-1)        ! strip off comments

! match first word to a command name

  call match_word (cmd_line, cmd_names, ix_cmd)
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

    call tao_cmd_split (cmd_line, 4, cmd_word, .false., err)

    if (cmd_word(1) == 'var') then
      call tao_change_var (cmd_word(2), cmd_word(3))
    elseif (cmd_word(1) == 'ele') then
      call tao_change_ele (cmd_word(2), cmd_word(3), cmd_word(4))
    else
      call out_io (s_error$, r_name, &
             'ERROR: CHANGE WHO? (SHOULD BE "ele" OR "var")')
    endif


!--------------------------------
! CLIP

  case ('clip')

    call tao_cmd_split (cmd_line, 3, cmd_word, .true., err); if (err) return
    if (cmd_word(2) == ' ') then
      call tao_clip_cmd (cmd_word(1), 0.0_rp, 0.0_rp) 
    else
      call tao_to_real (cmd_word(2), value1, err);  if (err) return
      if (cmd_word(3) /= ' ') then
        call tao_to_real (cmd_word(3), value2, err);  if (err) return
      else
        value2 = value1
        value1 = -value1
      endif
      call tao_clip_cmd (cmd_word(1), value1, value2)
    endif

!--------------------------------
! DERIVATIVE

  case ('derivative')

    call tao_dmodel_dvar_calc(.true.)
    call out_io (s_blank$, r_name, 'Derivative calculated')

!--------------------------------
! EXIT/QUIT

  case ('exit', 'quit')

    quit_tao = .false.
    call tao_query_logical ('y', 'n', 'Quit?', quit_tao)
    if (.not. quit_tao) return
    if (s%global%plot_on) call tao_destroy_plot_window
    stop
 
!--------------------------------
! RUN, FLATTEN

  case ('run', 'flatten')

    call tao_cmd_split (cmd_line, 1, cmd_word, .true., err); if (err) return
    call tao_run_cmd (cmd_word(1))

!--------------------------------
! HELP

  case ('help')

    call tao_cmd_split (cmd_line, 1, cmd_word, .true., err); if (err) return
    call tao_help (cmd_word(1))
    return

!--------------------------------
! HISTORY

  case ('history')

    call tao_cmd_split (cmd_line, 1, cmd_word, .true., err); if (err) return
    call tao_history_cmd (cmd_word(1), err)
    return

!--------------------------------
! OUTPUT

  case ('output')

    call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
    call tao_output_cmd (cmd_word(1), cmd_word(2))
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
    
    call match_word (cmd_word(1), name$%data_or_var, which)
    
    if (which .eq. data$) then
      call tao_use_data (cmd_name, cmd_word(2))
    elseif (which .eq. variable$) then
      call tao_use_var (cmd_name, cmd_word(2))
    else
      call out_io (s_error$, r_name, "Use/veto/restore what? data or variable?")
      return
    endif 

!--------------------------------
! REINITIALIZE

  case ('reinitialize')

    call tao_cmd_split(cmd_line, 2, cmd_word, .false., err)

    call out_io (s_warn$, r_name, &
         "Use this command with a little caution. There is a small memory leak somewhere!")
    
    ! quit the plot window so it will be recreated    
    call tao_destroy_plot_window
    s%global%init_plot_needed = .true.
    
    if (cmd_word(1) .eq. ' ') then
      call out_io (s_info$, r_name, 'Reinitializing with ' // &
                                         s%global%current_init_file)
      call tao_init (s%global%current_init_file)
    elseif (index(cmd_word(1), 'default') .ne. 0) then
      call out_io (s_info$, r_name, 'Reinitializing with ' // &
                                         s%global%default_init_file)
      call tao_init (s%global%default_init_file)
      s%global%current_init_file = s%global%default_init_file
    else
      call out_io (s_info$, r_name, 'Reinitializing with ' // &
                                         cmd_word(1))
      call tao_init (cmd_word(1))
      s%global%current_init_file = cmd_word(1)
    endif
      
!--------------------------------
! SET

  case ('set')

    call tao_cmd_split (cmd_line, 5, cmd_word, .false., err, '=')

    call match_word (cmd_word(1), set_names, ix)
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
         (set_word == 'graph'     .and. cmd_word(3) /= '=') .or. &
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
      call tao_set_graph_cmd (cmd_word(2), cmd_word(4), cmd_word(5))
    case ('universe')
      
      if (cmd_word(4) .eq. "recalc") then
        call tao_set_uni_cmd (cmd_word(2), cmd_word(3), .true.)
      else
        call tao_set_uni_cmd (cmd_word(2), cmd_word(3), .false.)
      endif
    end select

!--------------------------------
! SCALE

  case ('scale')

    if (.not. s%global%plot_on) then
      call out_io (s_error$, r_name, "PLOTTING TURNED OFF!")
      return
    endif

    call tao_cmd_split (cmd_line, 3, cmd_word, .true., err); if (err) return
    if (cmd_word(2) == ' ') then
      call tao_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp) 
    else
      call tao_to_real (cmd_word(2), value1, err);  if (err) return
      if (cmd_word(3) /= ' ') then
        call tao_to_real (cmd_word(3), value2, err);  if (err) return
      else
        value2 = value1
        value1 = -value1
      endif
      call tao_scale_cmd (cmd_word(1), value1, value2)
    endif

!--------------------------------
! SHOW

  case ('show')

    call tao_cmd_split (cmd_line, 2, cmd_word, .false., err); if (err) return
    call tao_show_cmd (cmd_word(1), cmd_word(2))
    return

!--------------------------------
! SINGLE-MODE

  case ('single-mode')

    s%global%single_mode = .true.
    call out_io (s_blank$, r_name, 'Entering Single Mode...')
    return

!--------------------------------
! SPAWN

  case ('spawn')
#ifdef CESR_LINUX
      call system (cmd_line)
      return
#endif
#ifdef CESR_UNIX
      call system (cmd_line)
      return
#endif

  call out_io (s_error$, r_name, "Spawn only works on linux and unix")
 

!--------------------------------
! VIEW

  case ('view')

    call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
    call tao_to_int (cmd_word(1), int1, err); if (err) return
    call tao_view_cmd (int1)

!--------------------------------
! X_AXIS

  case ('x-axis')

    if (.not. s%global%plot_on) then
      call out_io (s_error$, r_name, "PLOTTING TURNED OFF!")
      return
    endif

    call tao_cmd_split (cmd_line, 2, cmd_word, .true., err); if (err) return
    call tao_x_axis_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! X_SCALE

  case ('x-scale')

    if (.not. s%global%plot_on) then
      call out_io (s_error$, r_name, "PLOTTING TURNED OFF!")
      return
    endif

    call tao_cmd_split (cmd_line, 3, cmd_word, .true., err); if (err) return
    if (cmd_word(2) == ' ') then
      call tao_x_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, err)
    else
      call tao_to_real (cmd_word(2), value1, err); if (err) return
      call tao_to_real (cmd_word(3), value2, err); if (err) return
      call tao_x_scale_cmd (cmd_word(1), value1, value2, err)
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




