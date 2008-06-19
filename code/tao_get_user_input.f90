!+
! Subroutine tao_get_user_input (cmd_line, prompt_str)
!
! Subroutine to get input from the terminal.
!
! Input:
!   prompt_str -- Character(*), optional: Primpt string to print at terminal. If not
!                   present then s%global%prompt_string will be used.
!
! Output:
!   cmd_line -- Character(*): Command line from the user.
!-

subroutine tao_get_user_input (cmd_line, prompt_str)

use tao_mod
use input_mod

implicit none


integer i, ix, ix2
integer, allocatable, save :: indx(:), indx_start(:), indx_end(:), indx_step(:) ! for do loops
integer, allocatable, save :: loop_line_count(:) ! lines in each nested loop
integer, save :: in_loop = 0 ! in loop nest level
integer ios, n_level
character(10), save, allocatable :: indx_name(:) ! do loop index name
character(*) :: cmd_line
character(*), optional :: prompt_str
character(80) prompt_string

character(5) :: sub_str(9) = (/ '[[1]]', '[[2]]', '[[3]]', '[[4]]', '[[5]]', &
                            '[[6]]', '[[7]]', '[[8]]', '[[9]]' /)
character(40) tag
character(200), save :: saved_line
character(40) :: r_name = 'tao_get_user_input'

logical err, wait, flush
logical, save :: init_needed = .true.

! Init single char input

prompt_string = s%global%prompt_string
if (present(prompt_str)) prompt_string = prompt_str

if (init_needed) then
#ifndef CESR_WINCVF
  call init_tty_char
#endif
  init_needed = .false.
endif

! If single character input wanted then...

if (tao_com%single_mode) then
  call get_a_char (cmd_line(1:1), .true., (/ ' ' /))  ! ignore blanks
  tao_com%cmd_from_cmd_file = .false.
  return
endif

! check if we still have something from a line with multiple commands

if (tao_com%multi_commands_here) then
  call string_trim (saved_line, saved_line, ix)
  if (ix == 0) then
    tao_com%multi_commands_here = .false.
  else
    cmd_line = saved_line
  endif
endif

! If recalling a command from the cmd history stack...

if (tao_com%use_cmd_here) then
  cmd_line = tao_com%cmd
  call alias_translate (cmd_line, err)
  tao_com%use_cmd_here = .false.
  return
endif

! If a command file is open then read a line from the file.

if (tao_com%cmd_file_level /= 0) then

  call output_direct (do_print = s%global%command_file_print_on)

  n_level = tao_com%cmd_file_level
  if (.not. tao_com%multi_commands_here) then
    do
      read (tao_com%cmd_file(n_level)%ix_unit, '(a)', end = 8000) cmd_line
      tao_com%cmd_from_cmd_file = .true.
      call string_trim (cmd_line, cmd_line, ix)
      if (ix /= 0) exit
    enddo

    ! nothing more to do if an alias definition

    if (cmd_line(1:5) == 'alias') then
      call out_io (s_blank$, r_name, trim(prompt_string) // ': ' // trim(cmd_line))
      return
    endif

    ! replace argument variables

    do i = 1, 9
      ix = index (cmd_line, sub_str(i))
      if (ix /= 0) cmd_line = cmd_line(1:ix-1) // &
                          trim(tao_com%cmd_file(n_level)%cmd_arg(i)) // cmd_line(ix+5:)
    enddo
    
    call out_io (s_blank$, r_name, trim(prompt_string) // ': ' // trim(cmd_line))
    
    ! Check if in a do loop
    call do_loop()
    
  endif

  call alias_translate (cmd_line, err)
  call check_for_multi_commands

  return

  8000 continue
  close (tao_com%cmd_file(n_level)%ix_unit)
  tao_com%cmd_file(n_level)%ix_unit = 0 
  tao_com%cmd_file_level = n_level - 1 ! signal that the file has been closed
  if (tao_com%cmd_file_level /= 0) return ! still lower nested command file to complete
endif

! Here if no command file is being used.

if (.not. tao_com%multi_commands_here) then
  cmd_line = ' '
  tag = trim(prompt_string) // '> '
  tao_com%cmd_from_cmd_file = .false.
  call read_a_line (tag, cmd_line)
endif

call alias_translate (cmd_line, err)
call check_for_multi_commands

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
contains

subroutine alias_translate (cmd_line, err)

character(*) cmd_line
character(100) old_cmd_line

logical err, translated

!

old_cmd_line = cmd_line ! Save old command line for the command history.
translated = .false.    ! No translation done yet
call alias_translate2 (cmd_line, err, translated)

if (translated) then
  write (*, '(2a)') 'Alias: ', trim (cmd_line)
  cmd_line = trim(cmd_line) // "  ! " // trim(old_cmd_line)  
endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains

recursive subroutine alias_translate2 (cmd_line, err, translated)

character(*) cmd_line
character(100) string

integer ic, i, j
logical err, translated

! Look for a translation for the first word

call string_trim (cmd_line, cmd_line, ic)

do i = 1, tao_com%n_alias

  if (cmd_line(1:ic) /= tao_com%alias(i)%name) cycle

  ! We have a match...
  ! Now get the actual arguments and replace dummy args with actual args.

  string = cmd_line
  cmd_line = tao_com%alias(i)%string

  do j = 1, 9
    ix = index (cmd_line, sub_str(j))
    if (ix == 0) exit
    call string_trim (string(ic+1:), string, ic)
    cmd_line = cmd_line(1:ix-1) // &
                          trim(string(1:ic)) // cmd_line(ix+5:)
  enddo

  ! Append rest of string

  call string_trim (string(ic+1:), string, ic)
  cmd_line = trim(cmd_line) // ' ' // string
  call alias_translate2 (cmd_line, err, translated) ! Translation is an alias?
  translated = .true.

  return

enddo

translated = .false.

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains

subroutine check_for_multi_commands

  integer ix

  if (cmd_line(1:5) == 'alias') return

  ix = index (cmd_line, ';')
  if (ix /= 0) then
    tao_com%multi_commands_here = .true.
    saved_line = cmd_line(ix+1:)
    cmd_line = cmd_line(:ix-1)
  else
    saved_line = ' '
  endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains
!
! This routine has gotten ugly and should be re-written

subroutine do_loop

integer enddo_count ! for getting the line number count correct with nested loops
integer :: inner_loop_count =0

character(6) do_word ! 'do' or 'enddo'
character(15) indx_char
character(8) :: r_name = "do_loop"

1001 continue

  do_word = ' '
  call string_trim (cmd_line, cmd_line, ix)
  if (ix .le. len(do_word)) &
    call str_upcase(do_word(1:ix), cmd_line(1:ix))
  if (ix .eq. 2 .and. do_word(1:3) .eq. "DO ") then
    call set_loop_cmd_file_level (in_loop + 1)
    ! next word is loop index
    indx_name(in_loop) = ' '
    call string_trim (cmd_line(ix+1:), cmd_line, ix)
    indx_name(in_loop) = cmd_line(1:ix)
    ! now index start
    call string_trim (cmd_line(ix+1:), cmd_line, ix)
    read (cmd_line(1:ix), '(I)') indx_start(in_loop)
    ! now index end
    call string_trim (cmd_line(ix+1:), cmd_line, ix)
    read (cmd_line(1:ix), '(I)') indx_end(in_loop)
    indx(in_loop) = indx_start(in_loop) - 1 ! add one before first loop below
    ! finally index stepsize
    call string_trim (cmd_line(ix+1:), cmd_line, ix)
    if (ix /= 0) then ! specify index stepsize
      read (cmd_line(1:ix), '(I)') indx_step(in_loop)
    else
      indx_step(in_loop) = 1
    endif

    ! count loop statements so I know how many records to backspace on 'ENDDO"
    loop_line_count(in_loop) = 0
    enddo_count = 0
    inner_loop_count = 0
    do 
      read (tao_com%cmd_file(tao_com%cmd_file_level)%ix_unit, '(a)', end = 9000) cmd_line
      write (*, '(3a)') trim(prompt_string), ': ', trim(cmd_line)
      call string_trim (cmd_line, cmd_line, ix)
      do_word = ' '
      if (ix .le. len(do_word)) &
        call str_upcase(do_word(1:ix), cmd_line(1:ix))
        if (ix .eq. 5 .and. do_word(1:6) .eq. "ENDDO ") then
          if (enddo_count .eq. 0) then
            exit
          else
            enddo_count = enddo_count - 1
          endif
      endif
      if (ix .eq. 2 .and. do_word(1:3) .eq. "DO ") then
        inner_loop_count = inner_loop_count + 1
        enddo_count = enddo_count + 1
      endif
      loop_line_count(in_loop) = loop_line_count(in_loop) + 1
    enddo

    if (inner_loop_count /= 0) then
      ! There's an inner loop so rewind to beginning of first inner loop
      do i = 1, loop_line_count(in_loop) + 1
        backspace (tao_com%cmd_file(tao_com%cmd_file_level)%ix_unit)
        backspace (tao_com%cmd_file(tao_com%cmd_file_level)%ix_unit)
        read (tao_com%cmd_file(tao_com%cmd_file_level)%ix_unit, '(a)', end = 9000) cmd_line
        call string_trim (cmd_line, cmd_line, ix)
        do_word = ' '
        if (ix .le. len(do_word)) &
          call str_upcase(do_word(1:ix), cmd_line(1:ix))
        if (ix .eq. 2 .and. do_word(1:3) .eq. "DO ") then
          inner_loop_count = inner_loop_count - 1
          if (inner_loop_count .eq. 0) then
            indx(in_loop) = indx(in_loop) + indx_step(in_loop)
            goto 1001
          endif
        endif
        if (i .eq. loop_line_count(in_loop) + 1) then
          call out_io (s_error$, r_name, "Internal Error in routine!")
        endif
      enddo
    endif
  endif

  ! check if hit 'ENDDO'
  call string_trim (cmd_line, cmd_line, ix)
  do_word = ' '
  if (ix .le. len(do_word)) &
    call str_upcase(do_word(1:ix), cmd_line(1:ix))
  if (ix .eq. 5 .and. do_word(1:6) .eq. "ENDDO ") then
    if (in_loop .eq. 0) then
      call out_io (s_error$, r_name, &
                   "ENDDO found without correspoding DO statement")
      return
    endif
    indx(in_loop) = indx(in_loop) + indx_step(in_loop)
    if (indx(in_loop) .le. indx_end(in_loop) .and. indx_step(in_loop) .gt. 0) then
      ! rewind
      do i = 1, loop_line_count(in_loop) + 1
        backspace (tao_com%cmd_file(tao_com%cmd_file_level)%ix_unit)
      enddo
    elseif (indx(in_loop) .ge. indx_end(in_loop) .and. indx_step(in_loop) .lt. 0) then
      call out_io (s_error$, r_name, &
        "negative step size in loops not yet implemented: will ignore loop")
    else
      ! looped correct number of times, now exit loop
      call set_loop_cmd_file_level (in_loop - 1)
      if (in_loop .ge. 1) then
        read (tao_com%cmd_file(tao_com%cmd_file_level)%ix_unit, '(a)', end = 9000) cmd_line
        goto 1001
      else
      endif
    endif
    ! read next line
    read (tao_com%cmd_file(tao_com%cmd_file_level)%ix_unit, '(a)', end = 9100) cmd_line
    goto 1001
  endif
  
  ! insert index name variables for each loop
  if (in_loop /= 0) then
    do i = 1, in_loop
      do 
        ix = index (cmd_line, '[' // trim(indx_name(i)) // ']')
        if (ix /= 0) then
          write (indx_char, '(I)') indx(i)
          call string_trim(indx_char, indx_char, ix2)
          ix2 = len(trim(indx_name(i)))+2
          write (cmd_line, *) cmd_line(1:ix-1), trim(indx_char), &
                                    trim(cmd_line(ix+ix2:))
        else
          exit
        endif
      enddo
    enddo
  endif

  return

  ! No 'ENDDO' statement
  9000 continue
  call out_io (s_error$, r_name, &
       "No corresponding 'enddo' statment found, loop will be ignored")
  9100 continue
  close (tao_com%cmd_file(tao_com%cmd_file_level)%ix_unit)
  tao_com%cmd_file(tao_com%cmd_file_level)%ix_unit = 0 
  tao_com%cmd_file_level = tao_com%cmd_file_level - 1 ! signal that the file has been closed
  in_loop = 0

  ! don't send last 'ENDDO' to tao_command.f90
  call string_trim (cmd_line, cmd_line, ix)
  do_word = ' '
  if (ix .le. len(do_word)) &
    call str_upcase(do_word(1:ix), cmd_line(1:ix))
  if (ix .eq. 5 .and. do_word(1:6) .eq. "ENDDO ") cmd_line = " "
    
end subroutine do_loop

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains
!

subroutine set_loop_cmd_file_level (level)

implicit none

integer level
  
  in_loop = level
  call re_allocate (loop_line_count, in_loop)
  call re_allocate (indx, in_loop)
  call re_allocate (indx_start, in_loop)
  call re_allocate (indx_end, in_loop)
  call re_allocate (indx_step, in_loop)
  call re_allocate (indx_name, 10, in_loop)

end subroutine set_loop_cmd_file_level

end subroutine tao_get_user_input

