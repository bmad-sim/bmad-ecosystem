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

real(rp) tim, this_merit
integer i
logical err
character(*), parameter :: r_name = 'tao_cmd_end_calc'

! Note: tao_merit calls tao_lattice_calc.

this_merit = tao_merit()         

! Update variable values to reflect lattice values

call run_timer('START')
call tao_plot_setup()       ! transfer data to the plotting structures
if (associated(tao_hook_plot_setup_ptr)) call tao_hook_plot_setup_ptr()

do i = 1, size(s%plot_page%region)
  r => s%plot_page%region(i)
  if (.not. r%visible) cycle
  if (r%plot%autoscale_x)   call tao_x_scale_plot (r%plot, 0.0_rp, 0.0_rp)
  if (r%plot%autoscale_y)   call tao_scale_plot (r%plot, 0.0_rp, 0.0_rp)
enddo

call tao_draw_plots()              ! Update the plotting window

call run_timer('READ', tim)
if (s%global%max_plot_time > 0 .and. tim > s%global%max_plot_time) then
  call out_io(s_blank$, r_name, 'Time to plot is: ' // real_str(tim, n_decimal = 1) // ' seconds.', &
                                'To reduce plotting time try "set plot_page n_curve_pts = N" where N is less than the current ' // int_str(s%plot_page%n_curve_pts), &
                                'To turn off this message, "set global max_plot_time" to something negative or something greater than the plot time.')  
endif

end subroutine tao_cmd_end_calc


