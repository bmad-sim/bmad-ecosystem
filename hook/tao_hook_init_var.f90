!+
! Subroutine tao_hook_init_var (do_standard_setup) 
!
! Hook routine to initialize Tao variables.
!
! Output:
!   do_standard_setup -- logical: Set False to prevent the standard var init code from running.
!-

subroutine tao_hook_init_var (do_standard_setup) 

use tao_init_variables_mod, dummy => tao_hook_init_var

implicit none

logical do_standard_setup

!

do_standard_setup = .true. ! Change to .false. when doing a custom setup.

end subroutine
