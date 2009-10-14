!+
! Subroutine tao_hook_init_var (is_set) 
!
! Hook routine to initialize Tao variables.
!
! Output:
!   is_set -- Logical: Set True if variables have been initialized by this routine.
!                False otherwise. This tells the calling routine whether to skip
!                the standard variable initialization.
!-

subroutine tao_hook_init_var (is_set) 

use tao_mod

implicit none

logical is_set

!

is_set = .false.

end subroutine
