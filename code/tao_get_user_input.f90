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
use single_char_input_mod

implicit none


integer i, ix, ix2
integer, save, allocatable :: indx(:), indx_start(:), indx_end(:), indx_step(:) ! for do loops
integer, save, allocatable :: loop_line_count(:) ! lines in each nested loop
integer, save :: in_loop = 0 ! in loop nest level

character(10), save, allocatable :: indx_name(:) ! do loop index name
character(*) :: cmd_line
character(*), optional :: prompt_str
character(80) prompt_string

character(3) :: str(9) = (/ '[1]', '[2]', '[3]', '[4]', '[5]', &
                            '[6]', '[7]', '[8]', '[9]' /)
character(40) tag
character(200), save :: saved_line

logical err, wait, flush
logical, save :: init_needed = .true.
logical, save :: multi_commands_here = .false.

! Init single char input

prompt_string = s%global%prompt_string
if (present(prompt_str)) prompt_string = prompt_str

if (init_needed) then
  call init_tty_char
  init_needed = .false.
endif

! If single character input wanted then...

if (s%global%single_mode) then
  call get_a_char (cmd_line(1:1), .true., (/ ' ' /))  ! ignore blanks
  return
endif

! check if we still have something from a line with multiple commands

if (multi_commands_here) then
  call string_trim (saved_line, saved_line, ix)
  if (ix == 0) then
    multi_commands_here = .false.
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

if (tao_com%nest_level /= 0) then
  if (.not. multi_commands_here) then
    read (tao_com%lun_command_file(tao_com%nest_level), '(a)', end = 8000) cmd_line
    call string_trim (cmd_line, cmd_line, ix)

    ! replace argument variables
    if (cmd_line(1:5) == 'alias') return
    do i = 1, 9
      ix = index (cmd_line, str(i))
      if (ix /= 0) cmd_line = cmd_line(1:ix-1) // trim(tao_com%cmd_arg(i)) // &
                              cmd_line(ix+3:)
    enddo
    
    write (*, '(3a)') trim(prompt_string), ': ', trim(cmd_line)
    
    ! Check if in a do loop
    call do_loop()
    
  endif
  call alias_translate (cmd_line, err)
  call check_for_multi_commands

  return

  8000 continue
  close (tao_com%lun_command_file(tao_com%nest_level))
  tao_com%lun_command_file(tao_com%nest_level) = 0 
  tao_com%nest_level = tao_com%nest_level - 1 ! signal that the file has been closed
  if (tao_com%nest_level .ne. 0) return ! still lower nested command file to complete
endif

! Here if no command file is being used.

if (.not. multi_commands_here) then
  cmd_line = ' '
  tag = trim(prompt_string) // '> ' // achar(0)
  call read_line (trim(tag), cmd_line)
endif
call alias_translate (cmd_line, err)
call check_for_multi_commands

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
contains

subroutine alias_translate (cmd_line, err)

character(*) cmd_line
character(100) string

integer ic, i, j
logical err

!

call string_trim (cmd_line, cmd_line, ic)

do i = 1, tao_com%n_alias

  if (cmd_line(1:ic) /= tao_com%alias(i)%name) cycle

  ! get actual arguments and replace dummy args with actual args

  string = cmd_line
  cmd_line = tao_com%alias(i)%string

  do j = 1, 9
    ix = index (cmd_line, str(j))
    if (ix == 0) exit
    call string_trim (string(ic+1:), string, ic)
    cmd_line = cmd_line(1:ix-1) // &
                          trim(string(1:ic)) // cmd_line(ix+3:)
  enddo

  ! append rest of string

  call string_trim (string(ic+1:), string, ic)
  cmd_line = trim(cmd_line) // ' ' // string

  write (*, '(2a)') 'Alias: ', trim (cmd_line)
  return

enddo

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! contains

subroutine check_for_multi_commands

  integer ix

  if (cmd_line(1:5) == 'alias') return

  ix = index (cmd_line, ';')
  if (ix /= 0) then
    multi_commands_here = .true.
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
    call set_loop_nest_level (in_loop + 1)
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
    if (ix .ne. 0) then ! specify index stepsize
      read (cmd_line(1:ix), '(I)') indx_step(in_loop)
    else
      indx_step(in_loop) = 1
    endif

    ! count loop statements so I know how many records to backspace on 'ENDDO"
    loop_line_count(in_loop) = 0
    enddo_count = 0
    inner_loop_count = 0
    do 
      read (tao_com%lun_command_file(tao_com%nest_level), '(a)', end = 9000) cmd_line
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

    if (inner_loop_count .ne. 0) then
      ! There's an inner loop so rewind to beginning of first inner loop
      do i = 1, loop_line_count(in_loop) + 1
        backspace (tao_com%lun_command_file(tao_com%nest_level))
        backspace (tao_com%lun_command_file(tao_com%nest_level))
        read (tao_com%lun_command_file(tao_com%nest_level), '(a)', end = 9000) cmd_line
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
        backspace (tao_com%lun_command_file(tao_com%nest_level))
      enddo
    elseif (indx(in_loop) .ge. indx_end(in_loop) .and. indx_step(in_loop) .lt. 0) then
      call out_io (s_error$, r_name, &
        "negative step size in loops not yet implemented: will ignore loop")
    else
      ! looped correct number of times, now exit loop
      call set_loop_nest_level (in_loop - 1)
      if (in_loop .ge. 1) then
        read (tao_com%lun_command_file(tao_com%nest_level), '(a)', end = 9000) cmd_line
        goto 1001
      else
      endif
    endif
    ! read next line
    read (tao_com%lun_command_file(tao_com%nest_level), '(a)', end = 9100) cmd_line
    goto 1001
  endif
  
  ! insert index name variables for each loop
  if (in_loop .ne. 0) then
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
  close (tao_com%lun_command_file(tao_com%nest_level))
  tao_com%lun_command_file(tao_com%nest_level) = 0 
  tao_com%nest_level = tao_com%nest_level - 1 ! signal that the file has been closed
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

subroutine set_loop_nest_level (level)

implicit none

integer level
  
  in_loop = level
  call reallocate_integer (loop_line_count, in_loop)
  call reallocate_integer (indx, in_loop)
  call reallocate_integer (indx_start, in_loop)
  call reallocate_integer (indx_end, in_loop)
  call reallocate_integer (indx_step, in_loop)
  call reallocate_string (indx_name, 10, in_loop)

end subroutine set_loop_nest_level

end subroutine tao_get_user_input

