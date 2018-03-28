!+
! Subroutine tao_hook_show_cmd (what, result_id, lines, nl)
! 
! Dummy routine to modify the results of the "show" command.
!
! Input:
!   what      -- Character(*): What to show.
!   result_id -- Character(*): Set by tao_show_this routine
!   lines(:)  -- Character(*), allocatable: Output lines
!   nl        -- Integer: Current line number
!
!
! Output:
!   lines(:) -- Character(*): Output lines
!   nl       -- Integer: Current line number
!-

subroutine tao_hook_show_cmd (what, result_id, lines, nl)

use tao_interface, dummy => tao_hook_show_cmd

implicit none

type (tao_var_struct) var

character(*) what, result_id
character(*), allocatable :: lines(:)
character(*), parameter :: r_name = 'tao_hook_show_cmd'

integer nl

!

end subroutine
