!+
! Program tao_cl
!
! Command Line version of tao.
!-

program tao_cl

  use tao_mod
  use tao_cmd_history_mod

  implicit none

  type (tao_super_universe_struct) s
  logical end_flag, found, err
  character(100) cmd_line

! init

  call tao_init (s, 'tao.init')

! loop over commands.

  end_flag = .false.

  do
    err = .false.
    call tao_get_user_input (s, cmd_line)
    call tao_hook_command (s, cmd_line, found)
    if (.not. found) call tao_command (s, cmd_line, err)
    if (.not. err) call tao_cmd_history_record (cmd_line)
  enddo

end program
