!+
! Subroutine tao_hook_init_data (do_standard_setup) 
!
! Hook routine to initialize data.
!
! Output:
!   do_standard_setup -- logical: Set False to prevent the standard data init code from running.
!-

subroutine tao_hook_init_data (do_standard_setup) 

use tao_init_data_mod, dummy => tao_hook_init_data

implicit none

logical do_standard_setup

!

do_standard_setup = .true. ! Change to .false. when doing a custom setup.

end subroutine
