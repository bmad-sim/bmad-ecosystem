!+
! Subroutine tao_hook_parse_command_args (is_set) 
!
! Hook routine to read the command arguments
!
! Output:
!   is_set -- Logical: Set True if the command arguments are parsed by this routine.
!                False otherwise. This tells the calling routine whether to skip
!                the standard command arg parsing.
!-

subroutine tao_hook_parse_command_args (is_set)

implicit none

logical is_set

!

is_set = .false.

end subroutine
