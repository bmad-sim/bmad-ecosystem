!+
! Subroutine tao_hook_init_data (is_set) 
!
! Hook routine to initialize data.
!
! Output:
!   is_set -- Logical: Set True if data has been initialized by this routine.
!                False otherwise. This tells the calling routine whether to skip
!                the standard data initialization.
!-

subroutine tao_hook_init_data (is_set) 

use tao_mod

implicit none

logical is_set

!

is_set = .false.

end subroutine
