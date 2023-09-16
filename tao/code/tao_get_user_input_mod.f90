module tao_get_user_input_mod

use tao_interface
use input_mod

character(5), parameter, private :: sub_str(9) = ['[[1]]', '[[2]]', '[[3]]', '[[4]]', '[[5]]', &
                            '[[6]]', '[[7]]', '[[8]]', '[[9]]']

character(5), parameter, private :: opt_sub_str(9) = ['[<1>]', '[<2>]', '[<3>]', '[<4>]', '[<5>]', &
                            '[<6>]', '[<7>]', '[<8>]', '[<9>]']

private tao_alias_translate

contains

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!+
! Subroutine tao_get_user_input (cmd_out, prompt_str, wait_flag, cmd_in)
!
! Subroutine to get the next Tao command. In order of precedence, input may come from:
!   1) s%com%cmd string (if s%com%use_cmd_here is set to True). 
!      Used for recalling commands from the history stack.
!   2) A saved command string. 
!   3) A command file.
!   4) The cmd_in argument (if present). Used, for example, when interfacing with Python.
!   5) The terminal.
!
! Note: A saved command string is present if a prior input string contained multiple commands. 
! For example, the following string is read from a command file or terminal or passed via cmd_in:
!         "show ele 1; set opti de; run"
! Then cmd_out would be "show ele 1" and "set opti de; run" would be saved for the next call to this routine.
!
! Note: In single character mode, the input precedence order is ignored and input is taken from the terminal.
! 
! Input:
!   prompt_str -- Character(*), optional: Primpt string to print at terminal. If not
!                   present then s%global%prompt_string will be used.
!   wait_flag  -- logical, optional: Used for single mode: Wait state for get_a_char call.
!   cmd_in     -- Character(*), optional: Command to be used in place getting user input. 
!
! Output:
!   cmd_out    -- Character(*): Command from the user.
!-

subroutine tao_get_user_input (cmd_out, prompt_str, wait_flag, cmd_in)

implicit none

integer i, j, ix, ix1, ix2
integer ios, n_level

character(*) :: cmd_out
character(*), optional :: prompt_str, cmd_in
character(80) prompt_string, prompt_string_with_color
character(40) name
character(*), parameter :: r_name = 'tao_get_user_input'

logical, optional :: wait_flag
logical err, wait, flush, boldit, using_saved_cmd
logical, save :: init_needed = .true.

! Init single char input

if (init_needed) then
#ifdef CESR_WINCVF
#else
  call init_tty_char
#endif
  init_needed = .false.
endif

! Init

cmd_out = ''

prompt_string = s%global%prompt_string
if (present(prompt_str)) prompt_string = prompt_str
prompt_string_with_color = prompt_string

select case (upcase(s%global%prompt_color))
case ('BLACK');   call add_color (prompt_string, prompt_string_with_color, black_color)
case ('RED');     call add_color (prompt_string, prompt_string_with_color, red_color)
case ('GREEN');   call add_color (prompt_string, prompt_string_with_color, green_color)
case ('YELLOW');  call add_color (prompt_string, prompt_string_with_color, yellow_color)
case ('BLUE');    call add_color (prompt_string, prompt_string_with_color, blue_color)
case ('MAGENTA'); call add_color (prompt_string, prompt_string_with_color, magenta_color)
case ('CYAN');    call add_color (prompt_string, prompt_string_with_color, cyan_color)
case ('GRAY');    call add_color (prompt_string, prompt_string_with_color, gray_color)
end select

! If single character input wanted then...

if (s%com%single_mode) then
  if (s%global%wait_for_CR_in_single_mode) then
    if (s%com%single_mode_buffer == '') then
      do
        read '(a)', s%com%single_mode_buffer
        if (s%com%single_mode_buffer /= '') exit
      enddo
    endif
    cmd_out(1:1) = s%com%single_mode_buffer(1:1)
    s%com%single_mode_buffer = s%com%single_mode_buffer(2:)
  else
    wait = logic_option(.true., wait_flag)
    call get_a_char (cmd_out(1:1), wait, [' '])  ! ignore blanks
    s%com%cmd_from_cmd_file = .false.
  endif
  return
endif

s%com%single_mode_buffer = '' ! Reset buffer when not in single mode

! If a command was recalled from the cmd history stack use it.

if (s%com%use_cmd_here) then
  cmd_out = s%com%cmd
  call out_io (s_blank$, r_name, '  ' // cmd_out)
  cmd_out = tao_alias_translate (cmd_out, err)
  s%com%use_cmd_here = .false.
  return
endif

! Check if we still have something from a line with multiple commands

using_saved_cmd = .false.
if (s%com%saved_cmd_line /= '') then
  call string_trim (s%com%saved_cmd_line, s%com%saved_cmd_line, ix)
  cmd_out = s%com%saved_cmd_line
  using_saved_cmd = .true.
endif

! If a command file is open then read a line from the file.

n_level = s%com%cmd_file_level
if (n_level == 0) call tao_quiet_set ('off')  ! verbose if not running from a command file

if (n_level /= 0 .and. .not. s%com%cmd_file(n_level)%paused) then

  if (s%global%single_step) then
    call read_a_line ('Single_step: Press <return> to continue...', name, &
                                      prompt_color = s%global%prompt_color, prompt_bold = boldit)

  endif

  call output_direct (print_and_capture = (s%global%quiet /= 'all'), min_level = s_blank$, max_level = s_dwarn$)

  if (cmd_out == '') then
    ix = 0
    do
      read (s%com%cmd_file(n_level)%ix_unit, '(a)', end = 8000, iostat = ios) cmd_out(ix+1:)
      if (ios /= 0) then
        call tao_quiet_set('off')
        call out_io (s_error$, r_name, 'CANNOT READ LINE FROM FILE: ' // s%com%cmd_file(n_level)%full_name)
        goto 8000
      endif

      s%com%cmd_file(n_level)%n_line = s%com%cmd_file(n_level)%n_line + 1
      s%com%cmd_from_cmd_file = .true.
      call string_trim (cmd_out, cmd_out, ix)
      ix = len_trim(cmd_out)
      if (ix == 0) cycle
      if (cmd_out(ix:ix) == '&') then
        cmd_out(ix:ix) = ' '
        ix = len_trim(cmd_out)
        cycle
      endif
      !!!! if (any(cmd_out(ix:ix) == [',', '(', '{', '[', '='])) cycle
      exit
    enddo

    ! Nothing more to do if an alias definition

    if (cmd_out(1:5) == 'alias') then
      call out_io (s_blank$, r_name, '', trim(prompt_string_with_color) // ': ' // trim(cmd_out))
      call tao_quiet_set(s%global%quiet)
      return
    endif

    ! Replace argument variables

    do i = 1, 9
      do j = 1, 10
        ix = max(index(cmd_out, sub_str(i)), index(cmd_out, opt_sub_str(i)))
        if (ix == 0) exit
        cmd_out = cmd_out(1:ix-1) // trim(s%com%cmd_file(n_level)%cmd_arg(i)) // cmd_out(ix+5:)
      enddo
    enddo

    loop1: do
      ix1 = index(cmd_out, '[[')
      ix2 = index(cmd_out, ']]')
      if (ix1 == 0) exit
      if (.not. allocated(s%com%do_loop) .or. ix2 < ix1) then
        call out_io (s_error$, r_name, 'MALFORMED LINE IN COMMAND FILE: ' // cmd_out)
        call tao_abort_command_file()
        cmd_out = ''
        return
      endif
      name = cmd_out(ix1+2:ix2-1)

      do i = 1, size(s%com%do_loop)
        if (name == s%com%do_loop(i)%name) then
          cmd_out = cmd_out(1:ix1-1) // int_str(s%com%do_loop(i)%value) // cmd_out(ix2+2:) 
          cycle loop1
        endif
      enddo

      call out_io (s_error$, r_name, 'CANNOT MATCH NAME IN [[...]] CONSTRUCT: ' // cmd_out)
      call tao_abort_command_file()
      cmd_out = ''
      return

    enddo loop1

    if (s%global%quiet /= 'all') then
      if (s%global%blank_line_between_commands) call out_io (s_blank$, r_name, '')
      call out_io (s_blank$, r_name, trim(prompt_string_with_color) // ': ' // trim(cmd_out))
    endif

    ! Check if in a do loop
    call parse_do_loop(s%com%lev_loop, s%com%do_loop, n_level, cmd_out)
  endif

  !

  cmd_out = tao_alias_translate (cmd_out, err)
  call check_for_multi_commands()

  if (using_saved_cmd .or. s%com%saved_cmd_line /= '') then
    call out_io (s_blank$, r_name, '', trim(prompt_string_with_color) // ': ' // trim(cmd_out))
  endif

  call tao_quiet_set(s%global%quiet)
  return
endif

!-------------------------------------------------------------------------
! Here if no command file is being read from...

! Get a command if cmd_out has not already been set due to stuff saved in s%com%saved_cmd_line.

if (cmd_out == '') then
  if (present(cmd_in)) then
    cmd_out = cmd_in
  elseif (s%init%command_arg /= '' .and. .not. s%com%command_arg_has_been_executed) then
    cmd_out = s%init%command_arg
    s%com%command_arg_has_been_executed = .true.
  else
    s%com%cmd_from_cmd_file = .false.
    boldit = (s%global%prompt_color /= '' .and. s%global%prompt_color /= 'DEFAULT')
    if (s%global%blank_line_between_commands) call out_io (s_blank$, r_name, '')
    call read_a_line (trim(prompt_string) // '> ', cmd_out, prompt_color = s%global%prompt_color, &
                            prompt_bold = boldit, history_file = s%global%history_file)
    if (cmd_out == achar(24)) cmd_out = 'exit'   ! Cntl-D pressed

  endif
endif

! Alias translate the command and do other bookkeeping.

cmd_out = tao_alias_translate (cmd_out, err)
call check_for_multi_commands ()

if (using_saved_cmd .or. s%com%saved_cmd_line /= '') then
  call out_io (s_blank$, r_name, '', trim(prompt_string_with_color) // ': ' // trim(cmd_out))
endif

return

!-------------------------------------------------------------------------
! Here on an end of file detected or a file read error.

8000 continue
call tao_close_command_file()
cmd_out = ''

! If there exists a still lower nested command file to complete then return.
if (s%com%cmd_file_level /= 0) then
  if (s%com%cmd_file(n_level-1)%paused) then
    call out_io (s_info$, r_name, 'To continue the paused command file type "continue".')
  else
    return 
  endif
endif

call tao_quiet_set('off')
return

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
contains

subroutine check_for_multi_commands ()

integer ix
character(1) quote

!

if (cmd_out(1:5) == 'alias') return

quote = ''   ! Not in quote string
do ix = 1, len(cmd_out)

  if (quote == '') then
    select case (cmd_out(ix:ix))
    case ('!')
      exit

    case (';')
      s%com%saved_cmd_line = cmd_out(ix+1:)
      cmd_out = cmd_out(:ix-1)
      return

    case ("'", '"')
     quote = cmd_out(ix:ix)
    end select

  else ! quote /= ''
    if (cmd_out(ix:ix) == quote .and. cmd_out(ix-1:ix-1) /= '\') quote = ''           ! '
  endif

enddo

s%com%saved_cmd_line = ' '

end subroutine check_for_multi_commands

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains

subroutine parse_do_loop (lev_loop, loop, n_level, cmd_out)

use tao_command_mod

type (do_loop_struct), allocatable :: loop(:)
integer lev_loop, n_level, ix

character(*) cmd_out
character(8) :: r_name = "do_loop"
character(20) cmd_word(9)

logical err

! Check if a "do" statement

call tao_cmd_split (cmd_out, 9, cmd_word, .false., err, '=,')
if (cmd_word(1) /= 'do' .and. cmd_word(1) /= 'enddo') return

!

if (cmd_word(1) == 'do') then

  if (.not. is_integer(cmd_word(4)) .or. .not. is_integer(cmd_word(6)) .or. &
                        (.not. is_integer(cmd_word(8)) .and. cmd_word(8) /= '')) then 
    call out_io (s_error$, r_name, 'DO LOOP VARIABLES MUST BE INTEGERS.')
    call tao_abort_command_file()
    return
  endif

  if (cmd_word(3) /= '=' .or. cmd_word(5) /= ',' .or. cmd_word(9) /= '' .or. &
      (cmd_word(7) == ',' .and. cmd_word(8) == '') .or. (cmd_word(7) /= ',' .and. cmd_word(7) /= '')) then
    call out_io (s_error$, r_name, 'MALFORMED DO STATEMENT.')
    call tao_abort_command_file()
    return
  endif

  call set_loop_level (lev_loop, lev_loop + 1, loop)
  loop(lev_loop)%name = cmd_word(2)
  read (cmd_word(4), *) loop(lev_loop)%start
  read (cmd_word(6), *) loop(lev_loop)%end
  loop(lev_loop)%step = 1
  if (cmd_word(7) /= '') read (cmd_word(8), *) loop(lev_loop)%step

  loop(lev_loop)%n_line_start = s%com%cmd_file(n_level)%n_line
  loop(lev_loop)%value = loop(lev_loop)%start
  call out_io (s_blank$, r_name, 'Loop: ' // trim(loop(lev_loop)%name) // ' = \i0\ ', &
                                                  i_array = (/ loop(lev_loop)%value /) )
  cmd_out = '' ! So tao_command will not try to parse this line.

! Check if hit 'enddo'.

elseif (cmd_word(1) == 'enddo') then
  cmd_out = ''  ! so tao_command will not try to process this.
  if (lev_loop == 0) then
    call out_io (s_error$, r_name, 'ENDDO FOUND WITHOUT CORRESPODING DO STATEMENT')
    call tao_abort_command_file()
    return
  endif

  loop(lev_loop)%value = loop(lev_loop)%value + loop(lev_loop)%step
  if ((loop(lev_loop)%value <= loop(lev_loop)%end .and. loop(lev_loop)%step > 0) .or. &
      (loop(lev_loop)%value >= loop(lev_loop)%end .and. loop(lev_loop)%step < 0)) then
    ! rewind
    do i = s%com%cmd_file(n_level)%n_line, loop(lev_loop)%n_line_start+1, -1
      backspace (s%com%cmd_file(s%com%cmd_file_level)%ix_unit)
      s%com%cmd_file(n_level)%n_line = s%com%cmd_file(n_level)%n_line - 1
    enddo
    call out_io (s_blank$, r_name, 'Loop: ' // trim(loop(lev_loop)%name) // ' = \i0\ ', &
                                                  i_array = (/ loop(lev_loop)%value /) )
  else
    ! looped correct number of times, now exit loop
    call set_loop_level (lev_loop, lev_loop-1, loop)
  endif
  ! read next line
endif
    
end subroutine parse_do_loop

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains

subroutine set_loop_level (lev_loop, lev_new, loop)

implicit none

type (do_loop_struct), allocatable :: loop(:)
type (do_loop_struct), allocatable :: temp(:)
integer lev_loop, lev_new

! If going down then only need to set lev_loop

if (lev_new < lev_loop) then
  lev_loop = lev_new
  return
endif

!

lev_loop = lev_new

if (allocated(loop)) then
  if (size(loop) >= lev_loop) return
  call move_alloc(loop, temp)
endif

allocate (loop(lev_new))
if (allocated(temp)) then
  loop(1:size(temp)) = temp
  deallocate(temp)
endif

end subroutine set_loop_level

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains

subroutine add_color (str, c_str, color)
character(*) str, c_str, color
c_str = trim(color) // trim(str) // trim(reset_color)
end subroutine add_color

end subroutine tao_get_user_input

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

recursive function tao_alias_translate (cmd_str_in, err, depth) result(cmd_str_out)

character(*) cmd_str_in
character(:), allocatable :: cmd_str_out, tail, alias_cmd

integer, optional :: depth
integer ic, i, j, ix, ix2, this_depth
logical err, found_a_dummy_arg, found_this_dummy_arg
character(*), parameter :: r_name = 'tao_alias_translate'

!

if (present(depth)) then
  this_depth = depth + 1
  if (this_depth > 100) then
    call out_io (s_error$, r_name, 'INFINITE RECURSION LOOP TRANSLATING ALIASES.')
    cmd_str_out = ''
    err = .true.
    return
  endif
else
  this_depth = 1
  err = .false.
endif

! Look multiple commands

ix = index(cmd_str_in, ';')
if (ix /= 0) then
  cmd_str_out = tao_alias_translate(cmd_str_in(1:ix-1), err, this_depth)
  tail = cmd_str_in(ix+1:)

  do 
    ix = index(cmd_str_out, ';')
    if (ix == 0) exit
    cmd_str_out = trim(cmd_str_out) // ';' // tao_alias_translate(tail(1:ix-1), err, depth)
    tail = tail(ix+1:)
  enddo

  cmd_str_out = trim(cmd_str_out) // ';' // tao_alias_translate(tail, err, depth)
  return
endif

! Here if there is just a single command to translate

tail = cmd_str_in
call string_trim (tail, tail, ic)

do i = 1, s%com%n_alias

  if (tail(1:ic) /= s%com%alias(i)%name) cycle

  ! We have a match...
  ! Now get the actual arguments and replace dummy args with actual args.

  alias_cmd = trim(s%com%alias(i)%expanded_str)
  found_a_dummy_arg = .false.

  outer_loop: do j = 1, 9
    call string_trim (tail(ic+1:), tail, ic)
    found_this_dummy_arg = .false.
    ! There can be multiple dummy args for a given actual arg so need to loop here.
    do
      ix = index (alias_cmd, sub_str(j))
      ix2 = index(alias_cmd, opt_sub_str(j))
      if (ix == 0 .and. ix2 == 0 .and. .not. found_this_dummy_arg) exit outer_loop
      if (ix == 0 .and. ix2 == 0) exit
      found_a_dummy_arg = .true.
      found_this_dummy_arg = .true.
      if (tail(1:ic) == '' .and. ix /= 0) then
        call out_io (s_error$, r_name, 'ALIAS COMMAND DEMANDS MORE ARGUMENTS! ' // cmd_str_in)
        err = .true.
        cmd_str_out = ''
        return
      endif
      if (ix /= 0) then
        alias_cmd = alias_cmd(1:ix-1) // trim(tail(1:ic)) // trim(alias_cmd(ix+5:))
      else
        alias_cmd = alias_cmd(1:ix2-1) // trim(tail(1:ic)) // trim(alias_cmd(ix2+5:))
      endif
    enddo
  enddo outer_loop

  if (tail /= '') then
    if (found_a_dummy_arg) then
      call out_io (s_error$, r_name, 'EXTRA ARGUMENTS PRESENT IN COMMAND! ' // cmd_str_in)
      err = .true.
      cmd_str_out = ''
      return
    ! If there are no dummy args present then just append the tail to the command
    else
      alias_cmd = trim(alias_cmd) // ' ' // trim(tail)
    endif
  endif

  ! The translation may need to be translated.

  cmd_str_out = tao_alias_translate (alias_cmd, err, this_depth) ! Translation is an alias?
  return
enddo

! No translation needed

cmd_str_out = trim(cmd_str_in)

end function tao_alias_translate

end module

