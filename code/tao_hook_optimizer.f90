!+
! Subroutine tao_hook_optimizer (s)
!
! Dummy subroutine that needs to be over written in order to implement a custom
! optimizer.
!
! Input:
!   s  -- tao_super_universe_struct: 
!
! Output:
!   s  -- tao_super_universe_struct:
!-

subroutine tao_hook_optimizer (s)

use tao_mod
implicit none

type (tao_super_universe_struct), target :: s
character(20) :: r_name = 'tao_hook_optimizer'
!

  call out_io (s_error$, r_name, &
                      'THIS DUMMY ROUTINE NEEDS TO BE OVER WRITTEN!')

end subroutine
