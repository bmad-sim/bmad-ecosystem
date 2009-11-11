!+
! Subroutine tao_cl (prompt_string)
!
! Command Line version of tao.
!
! Modules needed:
!   use tao_mod
!
! Input:
!   prompt_string -- Character(*), optional: Prompt string. Default is 'Tao'.
!-

subroutine tao_cl (prompt_string)

  use tao_mod
  use tao_command_mod

  implicit none

  type (tao_super_universe_struct), pointer :: s_ptr  ! For debug purposes
  type (tao_common_struct), pointer :: t_ptr          ! For debug purposes

  character(200) cmd_line
  character(16) :: r_name = 'tao_cl'
  character(*), optional :: prompt_string

  logical end_flag, found, err

! init

  s_ptr => s
  t_ptr => tao_com
  tao_com%init_tao_file  = tao_com%default_init_tao_file

  call tao_parse_command_args (err)
  if (err) stop

  call tao_init (err)
  if (err) then
    call out_io (s_fatal$, r_name, 'TAO INIT FILE NOT FOUND. STOPPING.')
    stop
  endif

  if (present(prompt_string)) s%global%prompt_string = prompt_string

! loop over commands.

  end_flag = .false.

  do
    err = .false.
    call tao_get_user_input (cmd_line)
    if (tao_com%single_mode) then
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

end subroutine
