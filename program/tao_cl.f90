!+
! Program tao_cl
!
! Command Line version of tao.
!-

program tao_cl

  use tao_mod
  use tao_cmd_history_mod

  implicit none

  logical end_flag, found, err
  character(100) cmd_line

! init

  call tao_init ('tao.init')

! loop over commands.

  end_flag = .false.

  do
    err = .false.
    call tao_get_user_input (cmd_line)
    call tao_hook_command (cmd_line, found)
    if (.not. found) call tao_command (cmd_line, err)
    if (.not. err) call tao_cmd_history_record (cmd_line)
  enddo

end program
