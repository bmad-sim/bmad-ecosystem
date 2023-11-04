!+
! Subroutine tao_hook_optimizer (abort)
!
! Dummy subroutine that needs to be over written in order to implement a custom
! optimizer.
!
! Output:
!   abort -- Logical: Set True if an user stop signal detected.
!-

subroutine tao_hook_optimizer (abort)

use tao_interface

implicit none

character(*), parameter :: r_name = 'tao_hook_optimizer'
logical abort

!

call out_io (s_error$, r_name, 'THIS DUMMY ROUTINE NEEDS TO BE OVER WRITTEN!')

end subroutine
