module tao_scale_mod

use tao_mod
use quick_plot

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
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

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph

real(rp) y_min, y_max

integer i, j, ix, places

character(*) where

logical err

! If the where argument is blank or 'all' then scale all plots.

if (len_trim(where) == 0 .or. where(1:3) == 'all') then
  do j = 1, size(s%plot_page%region)
    plot => s%plot_page%region(j)%plot
    if (.not. s%plot_page%region(j)%visible) cycle
    call tao_scale_plot (plot, y_min, y_max)
  enddo
  return
endif

! locate the plot by the region name given by the where argument.
! If where has a ':' then we are dealing with just one graph of the plot.
! Otherwise we scale all the graphs of the plot.

call tao_find_plot_by_region (err, where, plot, graph)
if (err) return

ix = index(where, ':')
if (ix == 0) then                ! If all the graphs of a plot...
  call tao_scale_plot (plot, y_min, y_max)
else                          ! else just the one graph...
  call tao_scale_graph (graph, y_min, y_max)
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_scale_plot (plot, y_min, y_max)

type (tao_plot_struct) plot
type (tao_graph_struct), pointer :: graph
real(rp) y_min, y_max, this_min, this_max
integer i

! If we scale a whole plot with auto scale then at the end all graphs
! are adjusted to have the same scale such that all the data fits on
! all the graphs.

if (.not. associated (plot%graph)) return

do i = 1, size(plot%graph)
  call tao_scale_graph (plot%graph(i), y_min, y_max)
enddo

if (y_min == y_max .and. .not. plot%independent_graphs) then  ! if auto scale was done...
  this_min = minval (plot%graph(:)%y%min)
  this_max = maxval (plot%graph(:)%y%max)
  do i = 1, size(plot%graph)
    graph => plot%graph(i)
    call qp_calc_axis_scale (this_min, this_max, graph%y%major_div, &
                graph%y%bounds, graph%y%places, graph%y%min, graph%y%max)
  enddo
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_scale_graph (graph, y_min, y_max)

type (tao_graph_struct) graph
real(rp) y_min, y_max, this_min, this_max
integer i

! If y_min = y_max then autoscale: That is we need to find the 
! min/max so all the data points are within bounds.

if (.not. associated (graph%curve)) return

if (y_min == y_max) then

  this_min =  1e30
  this_max = -1e30

  do i = 1, size(graph%curve)
    if (associated(graph%curve(i)%y_symb)) then
      this_min = min(this_min, minval(graph%curve(i)%y_symb))
      this_max = max(this_max, maxval(graph%curve(i)%y_symb))
    endif
  enddo

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

end module
