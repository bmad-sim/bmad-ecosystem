!+
! Subroutine tao_hook_command (s, cmd_line, found)
!
!
! Dummy subroutine that needs to be over written in order to implement
! custom commands.
!
! Input:
!   s          -- tao_super_universe_struct
!   cmd_line   -- Character(*): command line
!   found      -- Logical: Set True if the command not handled by this routine.
!
!  Output:
!   s          -- tao_super_universe_struct
!-

subroutine tao_hook_command (s, cmd_line, found)

  use tao_mod

  implicit none

  type (tao_super_universe_struct), target :: s
  character(*) cmd_line
  logical found

!

  found = .false.

  if (cmd_line(1:4) == 'load') call load

!----------------------------------------------
contains

subroutine load

  type (tao_d2_data_struct), pointer :: d2_ptr
  type (tao_data_struct), pointer :: d

  integer i, j
  logical err

  found = .true.
  call tao_find_data (err, s%u(1), 'phase', d2_ptr)
  if (err) call err_exit

  do i = 1, size(d2_ptr%d1)
    do j = lbound(d2_ptr%d1(i)%d, 1), ubound(d2_ptr%d1(i)%d, 1)
      d => d2_ptr%d1(i)%d(j)
      d%good_data = .true.
    enddo
  enddo

end subroutine

end subroutine 
