!+
! Subroutine tao_scale_cmd (where, y_min, y_max)
!
! Routine to scale a plot. If y_min = y_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   where -- Character(*): Region to scale. Eg: "top:x"
!   y_min -- Real(rp): Plot y-axis min value.
!   y_max -- Real(rp): Plot y-axis max value.
!
!  Output:
!-

subroutine tao_scale_cmd (where, y_min, y_max)

use tao_mod
use quick_plot

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

real(rp) y_min, y_max

integer i, j, ix, places

character(*) where

real(rp) this_min, this_max

logical err

! If the where argument is blank then scale all graphs

if (len_trim(where) == 0 .or. where(1:3) == 'all') then
  do j = 1, size(s%plot_page%plot)
    plot => s%plot_page%plot(j)
    if (.not. plot%visible) cycle
    call scale_plot (plot)
  enddo
  return
endif

! locate the plot by the region name given by the where argument.
! If where has a ':' then we are dealing with just one graph of the plot.
! Otherwise we scale all the graphs of the plot.

call tao_find_plot (err, s%plot_page%plot, 'BY_REGION', where, plot, graph)
if (err) return

ix = index(where, ':')
if (ix == 0) then                ! If all the graphs of a plot...
  call scale_plot (plot)
else                          ! else just the one graph...
  call scale_graph (graph)
endif

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
contains

subroutine scale_plot (plot)

type (tao_plot_struct) plot

! If we scale a whole plot with auto scale then at the end all graphs
! are adjusted to have the same scale such that all the data fits on
! all the graphs.

do i = 1, size(plot%graph)
  call scale_graph (plot%graph(i))
enddo

if (y_min == y_max) then  ! if auto scale was done...
  do i = 1, size(plot%graph)
    graph => plot%graph(i)
    call qp_calc_axis_scale (this_min, this_max, graph%y%major_div, &
                graph%y%bounds, graph%y%places, graph%y%min, graph%y%max)
    if (graph%y%places < 0) graph%y%places = 0
  enddo
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! contains

subroutine scale_graph (graph)

type (tao_graph_struct) graph
integer i

! If y_min = y_max then autoscale: That is we need to find the 
! min/max so all the data points are within bounds.

if (y_min == y_max) then

  do i = 1, size(graph%curve)
    if (associated(graph%curve(i)%y_symb)) then
      this_max = max(this_max, maxval(graph%curve(i)%y_symb))
      this_min = min(this_min, minval(graph%curve(i)%y_symb))
    endif
  enddo

!  if (this_max-this_min < s%global%y_axis_plot_dmin) then
!    this_min = (this_min + this_max - s%global%y_axis_plot_dmin) / 2 
!    this_max = this_min + s%global%y_axis_plot_dmin
!  endif 

  call qp_calc_axis_scale (this_min, this_max, graph%y%major_div, &
                graph%y%bounds, graph%y%places, graph%y%min, graph%y%max)
  if (graph%y%places < 0) graph%y%places = 0
  return

endif

! If specific min/max values are given then life is easy.

graph%y%min = y_min
graph%y%max = y_max
call qp_calc_axis_places (y_min, y_max, graph%y%major_div, graph%y%places)
if (graph%y%places < 0) graph%y%places = 0

end subroutine

end subroutine 






