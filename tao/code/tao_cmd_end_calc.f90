!+
! Subroutine tao_cmd_end_calc ()
! 
! After every command this will do the standard lattice calculations and
! regenerate the plotting window
!
! Input:
!
! Output:
! s        -- tao_super_universe_struct: lattice calculations and plotting update.
!-

subroutine tao_cmd_end_calc

use tao_struct
use tao_plot_mod, only: tao_draw_plots
use tao_scale_mod
use tao_x_scale_mod

implicit none

type (tao_plot_region_struct), pointer :: r

real(rp) this_merit !not really used here
integer i
logical err

! Note: tao_merit calls tao_lattice_calc.

this_merit =  tao_merit()         

! update variable values to reflect lattice values

call tao_plot_setup()       ! transfer data to the plotting structures
call tao_hook_plot_setup()

do i = 1, size(s%plot_page%region)
  r => s%plot_page%region(i)
  if (.not. r%visible) cycle
  if (r%plot%autoscale_x)   call tao_x_scale_plot (r%plot, 0.0_rp, 0.0_rp)
  if (r%plot%autoscale_y)   call tao_scale_plot (r%plot, 0.0_rp, 0.0_rp)
enddo

call tao_draw_plots()              ! Update the plotting window

end subroutine tao_cmd_end_calc


