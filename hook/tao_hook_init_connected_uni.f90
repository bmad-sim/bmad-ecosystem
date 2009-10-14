!+
! Subroutine tao_hook_init_connected_uni (is_set)
!
! Hook routine to initialize connected universe information.
!
! Output:
!   is_set -- Logical: Set True if connected universe info has been initialized 
!                by this routine. False otherwise. This tells the calling routine 
!                whether to skip the standard connected universe initialization.
!-

subroutine tao_hook_init_connected_uni (is_set)

implicit none

logical is_set

!

is_set = .false.

end subroutine
