!+
! Subroutine tao_x_axis_cmd (where, what)
!
! Routine to set the type of axis for a plot.
! 
! Input:
!   where -- Character(*): Region to axis. Eg: "top"
!   what  -- Character(*): "s" or "index"
!-

subroutine tao_x_axis_cmd (where, what)

use tao_mod
use quick_plot
use tao_graph_setup_mod
use tao_x_scale_mod

implicit none

type (tao_plot_array_struct), allocatable, save :: plot(:)

integer i, j

character(*) where, what
character(16) :: r_name = 'tao_x_axis_cmd'

logical err

! check what

if (what /= "s" .and. what /= "index" .and. what /= "ele_index") then
  call out_io (s_error$, r_name, 'I DO NOT UNDERSTAND THIS: ' // what)
  return
endif

! If the where argument is blank then axis all graphs

if (len_trim(where) == 0 .or. where == '*' .or. where == 'all') then
  do j = 1, size(s%plot_region)
    if (.not. s%plot_region(j)%visible) cycle
    s%plot_region(j)%plot%x_axis_type = what
    call tao_x_scale_plot (s%plot_region(j)%plot, 0.0_rp, 0.0_rp)
  enddo
  return
endif

! locate the plot by the region name given by the where argument.

call tao_find_plots (err, where, 'REGION', plot)
if (err) return
do i = 1, size(plot)
  plot(i)%p%x_axis_type = what
  call tao_x_scale_plot (plot(i)%p, 0.0_rp, 0.0_rp)
enddo


end subroutine 






