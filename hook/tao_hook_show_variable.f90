!+
! Subroutine tao_hook_show_variable (var, lines, nl)
! 
! Dummy routine to add information to the "show variable" command.
!
! Input:
!   var     -- Tao_var_struct: Variable.
!
! Output:
!   lines(:) -- Character(*): Output lines
!   nl       -- Integer: Current line number
!-

subroutine tao_hook_show_variable (var, lines, nl)

use tao_mod

implicit none

type (tao_var_struct) var

character(*) lines(:)
character(20) :: r_name = 'tao_hook_show_variable'

integer nl

!

end subroutine
