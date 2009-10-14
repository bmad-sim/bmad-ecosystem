!+
! Subroutine tao_hook_init_beam (is_set) 
!
! Hook routine to initialize the beam for beam tracking
!
! Output:
!   is_set -- Logical: Set True if the beam has been initialized by this routine.
!                False otherwise. This tells the calling routine whether to skip
!                the standard beam initialization.
!-

subroutine tao_hook_init_beam (is_set)

implicit none

logical is_set

!

is_set = .false.

end subroutine
