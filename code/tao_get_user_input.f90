!+
! Subroutine tao_get_user_input (s, cmd_line)
!
! Subroutine to get input from the terminal.
!
! Input:
!   s%global%prompt_string -- Prompt string.
!
! Output:
!   cmd_line -- Character(*): Command line from the user.
!-

subroutine tao_get_user_input (s, cmd_line)

use tao_mod
use tao_common
use tao_single_mod
use single_char_input_mod

implicit none

type (tao_super_universe_struct) :: s

integer i, ix

character(*) :: cmd_line
character(3) :: str(9) = (/ '[1]', '[2]', '[3]', '[4]', '[5]', &
                            '[6]', '[7]', '[8]', '[9]' /)
character(40) tag

logical err, wait, flush
logical, save :: init_needed = .true.

! Init single char input

if (init_needed) then
  call init_tty_char
  init_needed = .false.
endif

! If single character input wanted then...

if (s%global%single_mode) then
  call get_a_char (cmd_line(1:1), .true., (/ ' ' /))  ! ignore blanks
  return
endif

! If recalling a command from the cmd history stack...

if (tao_com%use_cmd_here) then
  cmd_line = tao_com%cmd
  call alias_translate (cmd_line, err)
  tao_com%use_cmd_here = .false.
  return
endif

! If a command file is open then read a line from the file.

if (s%global%lun_command_file /= 0) then
  read (s%global%lun_command_file, '(a)', end = 8000) cmd_line
  do i = 1, 9
    ix = index (cmd_line, str(i))
    if (ix /= 0) cmd_line = cmd_line(1:ix-1) // trim(tao_com%cmd_arg(i)) // &
                            cmd_line(ix+3:)
  enddo
  write (*, '(3a)') trim(s%global%prompt_string), ': ', trim(cmd_line)
  call alias_translate (cmd_line, err)
  return

  8000 continue
  close (s%global%lun_command_file)
  s%global%lun_command_file = 0 ! signal that the file has been closed
endif

! Here if no command file is being used.

!! print '(1x, 2a, $)', trim(s%global%prompt_string), '> '
!! read (*, '(a)') cmd_line

cmd_line = ' '
tag = trim(s%global%prompt_string) // '> ' // achar(0)
call read_line (trim(tag), cmd_line)
call alias_translate (cmd_line, err)

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

end subroutine tao_get_user_input

