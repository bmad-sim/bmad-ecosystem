!+
! Program tao_cl
!
! Command Line version of tao.
!-

program tao_cl

  use tao_mod
  use tao_command_mod

  implicit none

  logical end_flag, found, err
  character(100) cmd_line

! init

  call tao_init (s%global%default_init_file)

! loop over commands.

  end_flag = .false.

  do
    err = .false.
    call tao_get_user_input (cmd_line)
    if (s%global%single_mode) then
      ! single mode
      call tao_single_mode (cmd_line(1:1))
      ! Do the standard calculations and plotting after command execution.
      call tao_cmd_end_calc
    else
      ! command line mode
      call tao_hook_command (cmd_line, found)
      if (.not. found) call tao_command (cmd_line, err)
    endif
    if (.not. err) call tao_cmd_history_record (cmd_line)
  enddo

end program
