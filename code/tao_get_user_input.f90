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

implicit none

type (tao_super_universe_struct) :: s

integer i, ix

character(*) :: cmd_line
character(3) :: str(9) = (/ '[1]', '[2]', '[3]', '[4]', '[5]', &
                            '[6]', '[7]', '[8]', '[9]' /)

logical err

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

write (*, '(2a)', advance = 'NO') trim(s%global%prompt_string), '> '
read (*, '(a)') cmd_line
call alias_translate (cmd_line, err)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
contains

subroutine alias_translate (cmd_line, err)

character(*) cmd_line
character(100) string
character(40) alias_arg(9)

integer ic, i, j
logical err

!

call string_trim (cmd_line, cmd_line, ic)

do i = 1, tao_com%n_alias
  if (cmd_line(1:ic) == tao_com%alias(i)%name) then

    alias_arg = ' '

    do j = 1, 9
      call string_trim (cmd_line(ic+1:), string, ic)
      if (ic == 0) exit
      alias_arg(j) = string(1:ic)
    enddo

    cmd_line = tao_com%alias(i)%string

    do j = 1, 9
      ix = index (cmd_line, str(j))
      if (ix /= 0) cmd_line = cmd_line(1:ix-1) // &
                          trim(alias_arg(i)) // cmd_line(ix+3:)

    enddo

    write (*, '(2a)') 'Alias: ', trim (cmd_line)
    return

  endif
enddo

end subroutine

end subroutine tao_get_user_input

