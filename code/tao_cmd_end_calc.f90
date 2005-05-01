!+
! Subroutine tao_cmd_end_calc ()
! 
! After every command this will do the standard lattice calculations and
! regenerate the plotting window
!
! Input:
!
! Output:
! s        -- tao_super_universe_struct: lattice calculations and plotting
!                                        update
!
!-

subroutine tao_cmd_end_calc

use tao_mod
use tao_plot_mod

implicit none

real(rp) this_merit !not really used here

! Note: tao_merit calls tao_lattice_calc.

this_merit =  tao_merit()         
call tao_plot_data_setup()       ! transfer data to the plotting structures
call tao_hook_plot_data_setup()
call tao_plot_out()              ! Update the plotting window

end subroutine tao_cmd_end_calc


