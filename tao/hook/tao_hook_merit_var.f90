!+
! Subroutine tao_hook_merit_var (i_uni, j_var, var)
! 
! Dummy routine that needs to be over written in order to implement a
! custom merit calculation for variables.
!
! Input:
!   i_uni   -- Integer: Universe index for the variable in the s%u(:) array.
!   j_var   -- Integer: index of the variable in the u%var(:) array.
!   var     -- Tao_var_struct: Variable whose contribution to the merit function
!               is to be calculated.
!
! Output:
!   var      -- Tao_var_struct: Variable whose contribution to the merit function
!                is to be calculated.
!     %merit -- Real(rp): Contribution to the merit function.
!-

subroutine tao_hook_merit_var (i_uni, j_var, var)

use tao_interface

implicit none

type (tao_var_struct) var

integer, intent(in) :: i_uni, j_var
character(*), parameter :: r_name = 'tao_hook_merit_var'

!

call out_io (s_abort$, r_name, 'THIS ROUTINE SHOULD NOT HAVE BEEN CALLED')
call out_io (s_abort$, r_name, 'MERIT_TYPE NOT RECOGNIZED FOR VAR: ' &
        // tao_var1_name(var), 'MERIT_TYPE: ' // var%merit_type)

call err_exit

end subroutine
