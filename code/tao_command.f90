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

  use tao_mod
  use quick_plot
  use tao_cmd_history_mod
  use tao_dmerit_mod
  use tao_scale_mod
  use tao_x_scale_mod
  use tao_set_mod

  implicit none


  integer i, j
  integer ix, ix_line, ix_cmd, which
  integer int1, int2

  real(rp) value1, value2

  character(*) :: command_line
  character(140) cmd_line
  character(20) :: r_name = 'tao_command'
  character(40) :: cmd_word(12)
 
  character(16) cmd_name
  character(16) :: cmd_names(26) = (/  &
        'quit        ', 'exit        ', 'show        ', 'plot        ', 'place       ', &
        'clip        ', 'scale       ', 'veto        ', 'use         ', 'restore     ', &
        'run         ', 'flatten     ', 'output      ', 'change      ', 'set         ', &
        'call        ', 'view        ', 'alias       ', 'help        ', 'history     ', &
        'single-mode ', '            ', 'x-scale     ', 'x-axis      ', 'derivative  ', &
        'reinitialize' /)

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

    call cmd_split(2, .false., err); if (err) return
    call tao_alias_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! CALL

  case ('call')

    call cmd_split(10, .true., err); if (err) return
    call tao_call_cmd (cmd_word(1), cmd_word(2:10))

!--------------------------------
! CHANGE

  case ('change')

    call cmd_split (4, .false., err)
    call tao_change_cmd (cmd_word(1), cmd_word(2), cmd_word(3), cmd_word(4))

!--------------------------------
! CLIP

  case ('clip')

    call cmd_split (3, .true., err); if (err) return
    if (cmd_word(2) == ' ') then
      call tao_clip_cmd (cmd_word(1), 0.0_rp, 0.0_rp) 
    else
      call to_real (cmd_word(2), value1, err);  if (err) return
      if (cmd_word(3) /= ' ') then
        call to_real (cmd_word(3), value2, err);  if (err) return
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
    call qp_close_page
    stop
 
!--------------------------------
! RUN, FLATTEN

  case ('run', 'flatten')

    call cmd_split (1, .true., err); if (err) return
    call tao_run_cmd (cmd_word(1))

!--------------------------------
! HELP

  case ('help')

    call cmd_split (1, .true., err); if (err) return
    call tao_help (cmd_word(1))

!--------------------------------
! HISTORY

  case ('history')

    call cmd_split (1, .true., err); if (err) return
    call tao_history_cmd (cmd_word(1), err)

!--------------------------------
! OUTPUT

  case ('output')

    call cmd_split (1, .true., err); if (err) return
    call tao_output_cmd (cmd_word(1))

!--------------------------------
! PLACE

  case ('place')

    call cmd_split (3, .true., err); if (err) return
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

    call cmd_split (9, .false., err, '+-')

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

    call cmd_split(3, .false., err)
    
    call match_word (cmd_word(1), name$%data_or_var, which)
    
    if (which .eq. data$) then
      call tao_use_data (cmd_name, cmd_word(2), cmd_word(3))
    elseif (which .eq. variable$) then
      call tao_use_var (cmd_name, cmd_word(2), cmd_word(3))
    else
      call out_io (s_error$, r_name, "Use/veto/restore what? data or variable?")
      return
    endif 

!--------------------------------
! REINITIALIZE

  case ('reinitialize')

    call cmd_split(2, .false., err)
! quit the plot window so it will be recreated    
    call qp_close_page
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

    call cmd_split (6, .false., err, '=')

    if (((cmd_word(1) == 'data' .or. cmd_word(1) == 'var') .and. &
                  cmd_word(4) /= '=') .or. &
                  (cmd_word(1) == 'global' .and. cmd_word(3) /= '=') .or. &
                  (cmd_word(1) == 'plot'   .and. cmd_word(3) /= '=') .or. &
		  (cmd_word(1)(1:7) == 'lattice' .and. cmd_word(3) /= '=')) then
      call out_io (s_error$, r_name, 'SYNTAX PROBLEM. "=" NOT IN CORRECT PLACE.')
      return
    endif


    select case (cmd_word(1))
    case ('data')
      call tao_set_data_cmd (cmd_word(2), &
                              cmd_word(3), cmd_word(5), cmd_word(6)) 
    case ('var')
      call tao_set_var_cmd (cmd_word(2), &
                              cmd_word(3), cmd_word(5), cmd_word(6)) 
    case ('lattice')
      call tao_set_lattice_cmd (cmd_word(1), &
                              cmd_word(2), cmd_word(4)) 
    case ('global')
      call tao_set_global_cmd (cmd_word(2), cmd_word(4))
    case ('plot')
      call tao_set_plot_cmd (cmd_word(2), cmd_word(4), cmd_word(5))
    case default
      call out_io (s_error$, r_name, 'NOT RECOGNIZED: ' // cmd_word(1))
    end select

!--------------------------------
! SCALE

  case ('scale')

    call cmd_split (3, .true., err); if (err) return
    if (cmd_word(2) == ' ') then
      call tao_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp) 
    else
      call to_real (cmd_word(2), value1, err);  if (err) return
      if (cmd_word(3) /= ' ') then
        call to_real (cmd_word(3), value2, err);  if (err) return
      else
        value2 = value1
        value1 = -value1
      endif
      call tao_scale_cmd (cmd_word(1), value1, value2)
    endif

!--------------------------------
! SHOW

  case ('show')

    call cmd_split (3, .false., err)
    if (cmd_word(1) == ' ') then
      call out_io (s_error$, r_name, 'SHOW WHAT?')
      return
    endif

    call tao_show_cmd (cmd_word(1), cmd_word(2), cmd_word(3))

!--------------------------------
! SINGLE-MODE

  case ('single-mode')

    s%global%single_mode = .true.
    call out_io (s_blank$, r_name, 'Entering Single Mode...')

!--------------------------------
! VIEW

  case ('view')

    call cmd_split (2, .true., err); if (err) return
    call to_int (cmd_word(1), int1, err); if (err) return
    call tao_view_cmd (int1)

!--------------------------------
! X_AXIS

  case ('x-axis')

    call cmd_split (2, .true., err); if (err) return
    call tao_x_axis_cmd (cmd_word(1), cmd_word(2))

!--------------------------------
! X_SCALE

  case ('x-scale')

    call cmd_split (3, .true., err); if (err) return
    if (cmd_word(2) == ' ') then
      call tao_x_scale_cmd (cmd_word(1), 0.0_rp, 0.0_rp, err)
    else
      call to_real (cmd_word(2), value1, err); if (err) return
      call to_real (cmd_word(3), value2, err); if (err) return
      call tao_x_scale_cmd (cmd_word(1), value1, value2, err)
    endif

!--------------------------------
! DEFAULT

  case default

    call out_io (s_error$, r_name, 'INTERNAL COMMAND PARSING ERROR!')
    call err_exit

  end select

! do the standard calculations and plotting after command
  call tao_cmd_end_calc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains


! This routine splits the command into words.
!
! Input: 
!  n_word         -- integer: number of words to split command line into
!  no_extra_words -- logical: are extra words allowed at the end?
!  err            -- logical: error in splitting words
!  separator      -- character(*): a list of characters that, besides a blank space,
!                                  signify a word boundary. 
!
! Output:
!  cmd_word(n_word) -- character(40): the individual words
!
! For example: 
!   separator = '-+' 
!   cmd_line = 'model-design'
! Restult:
!   cmd_word(1) = 'model'
!   cmd_word(2) = '-'
!   cmd_word(3) = 'design'
!
! Anything between single or double quotes is treated as a single word.
!

subroutine cmd_split (n_word, no_extra_words, err, separator)

  integer i, n, n_word, ix_end_quote
  character(*), optional :: separator
  logical err
  logical no_extra_words

!

  cmd_word(:) = ' '

  if (present(separator)) then
    i = 0
    do 
      if (i == len(cmd_line)) exit
      i = i + 1
      if (index(separator, cmd_line(i:i)) /= 0) then
        cmd_line = cmd_line(:i-1) // ' ' // cmd_line(i:i) // ' ' // &
                                                            cmd_line(i+1:)
        i = i + 3
      endif
    enddo
  endif

  ix_line = 0
  ix_end_quote = 0
  
  do n = 1, n_word
    call string_trim (cmd_line(ix_line+1:), cmd_line, ix_line)
    if (cmd_line(1:1) .eq. '"') then
      ix_end_quote = index(cmd_line(2:), '"')
      if (ix_end_quote .eq. 0) ix_end_quote = len(cmd_line)
      cmd_word(n) = cmd_line(2:ix_end_quote)
    elseif (cmd_line(1:1) .eq. "'") then
      ix_end_quote = index(cmd_line(2:), "'")
      if (ix_end_quote .eq. 0) ix_end_quote = len(cmd_line)
      cmd_word(n) = cmd_line(2:ix_end_quote)
    else
      cmd_word(n) = cmd_line(:ix_line)
    endif 
    if (ix_line == 0) return
  enddo

  cmd_word(n_word) = cmd_line

  if (no_extra_words) then
    call string_trim (cmd_line(ix_line+1:), cmd_line, ix_line)
    if (ix_line /= 0) then
      call out_io (s_error$, r_name, 'EXTRA STUFF ON COMMAND LINE: ' // &
                                                                  cmd_line)
      err = .true.
    endif
  endif

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! contains
! This routine converts a string to a real number.

subroutine to_real (str, r_real, err)

character(*) str
real(rp) r_real
integer ios
logical err

!

err = .false.
read (str, *, iostat = ios) r_real

if (ios /= 0) then
  call out_io (s_error$, r_name, 'EXPECTING REAL NUMBER: ' // str)
  err = .true.
  return
endif

end subroutine

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! contains
! This routine converts a string to an integer

subroutine to_int (str, i_int, err)

character(*) str
integer ios, i_int
logical err

!

err = .false.
read (str, *, iostat = ios) i_int

if (ios /= 0) then
  call out_io (s_error$, r_name, 'EXPECTING INTEGER: ' // str)
  err = .true.
  return
endif

end subroutine


end subroutine tao_command




