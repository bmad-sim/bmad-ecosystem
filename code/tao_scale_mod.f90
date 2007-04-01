module tao_scale_mod

use tao_mod
use quick_plot

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_scale_cmd (where, axis, y_min, y_max)
!
! Routine to scale a plot. If y_min = y_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   where -- Character(*): Region to scale. Eg: "top:x"
!   axis  -- Character(*): 'y', 'y2', or '' (both).
!   y_min -- Real(rp): Plot y-axis min value.
!   y_max -- Real(rp): Plot y-axis max value.
!
!  Output:
!-

subroutine tao_scale_cmd (where, axis, y_min, y_max)

implicit none

type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)

real(rp) y_min, y_max

integer i, j, ix, places

character(*) where, axis
character(20) :: r_name = 'tao_scale_cmd'

logical err

! Error check

if (axis /= 'y' .and. axis /= 'y2' .and. axis /= '') then
  call out_io (s_error$, r_name, 'BAD AXIS NAME: ' // axis)
  return
endif

! If the where argument is blank or 'all' then scale all plots.

if (len_trim(where) == 0 .or. where(1:3) == 'all') then
  do j = 1, size(s%plot_page%region)
    if (.not. s%plot_page%region(j)%visible) cycle
    call tao_scale_plot (s%plot_page%region(j)%plot, axis, y_min, y_max)
  enddo
  return
endif

! locate the plot by the region name given by the where argument.
!If no graph is specified then we scale all the graphs of the plot.

call tao_find_plots (err, where, 'REGION', plot, graph)
if (err) return

if (allocated(graph)) then                ! If all the graphs of a plot...
  do j = 1, size(graph)
    call tao_scale_graph (graph(j)%g, axis, y_min, y_max)
  enddo
else                          ! else just the one graph...
  do i = 1, size(plot)
    call tao_scale_plot (plot(i)%p, axis, y_min, y_max)
  enddo
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_scale_plot (plot, axis, y_min, y_max)

type (tao_plot_struct), target :: plot
type (tao_graph_struct), pointer :: graph
real(rp) y_min, y_max, this_min, this_max
character(*) axis
integer i

! If we scale a whole plot with auto scale then at the end all graphs
! are adjusted to have the same scale such that all the data fits on
! all the graphs.

if (.not. allocated (plot%graph)) return

do i = 1, size(plot%graph)
  call tao_scale_graph (plot%graph(i), axis, y_min, y_max)
enddo

if (y_min == y_max .and. .not. plot%independent_graphs) then  ! if auto scale was done...

  if (axis == '' .or. axis == 'y') then
    this_min = minval (plot%graph(:)%y%min)
    this_max = maxval (plot%graph(:)%y%max)
    do i = 1, size(plot%graph)
      graph => plot%graph(i)
      call qp_calc_axis_scale (this_min, this_max, graph%y)
    enddo
  endif

  if (axis == '' .or. axis == 'y2') then
    this_min = minval (plot%graph(:)%y2%min)
    this_max = maxval (plot%graph(:)%y2%max)
    do i = 1, size(plot%graph)
      graph => plot%graph(i)
      call qp_calc_axis_scale (this_min, this_max, graph%y2)
    enddo
  endif

endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_scale_graph (graph, axis, y_min, y_max)

type (tao_graph_struct) graph
real(rp) y_min, y_max, this_min, this_max, this_min2, this_max2
integer i
character(*) axis
logical found_data, found_data2

! If specific min/max values are given then life is easy.

if (y_min /= y_max) then

  if (axis == '' .or. axis == 'y') then
    graph%y%min = y_min
    graph%y%max = y_max
    call qp_calc_axis_places (graph%y)
  endif

  if (axis == '' .or. axis == 'y2') then
    graph%y2%min = y_min
    graph%y2%max = y_max
    call qp_calc_axis_places (graph%y2)
  endif

  return

endif

! Since y_min = y_max then autoscale: That is we need to find the 
! min/max so all the data points are within bounds.

! For a floor plan the min and max have been stored in graph%y_min, graph%y_max

if (graph%type == 'floor_plan') then
  call qp_calc_axis_scale (graph%y_min, graph%y_max, graph%y)
  return
endif

!

if (.not. allocated (graph%curve)) return

this_min =  1e30
this_max = -1e30
this_min2 =  1e30
this_max2 = -1e30
found_data = .false.
found_data2 = .false.

do i = 1, size(graph%curve)

  if (allocated(graph%curve(i)%y_symb)) then
    if (size(graph%curve(i)%y_symb) > 0) then
      if (graph%curve(i)%use_y2) then
        this_min2 = min(this_min2, minval(graph%curve(i)%y_symb))
        this_max2 = max(this_max2, maxval(graph%curve(i)%y_symb))
        found_data2 = .true.
      else
        this_min = min(this_min, minval(graph%curve(i)%y_symb))
        this_max = max(this_max, maxval(graph%curve(i)%y_symb))
        found_data = .true.
      endif
    endif
  endif

  if (allocated(graph%curve(i)%y_line)) then
    if (size(graph%curve(i)%y_line) > 0) then
      if (graph%curve(i)%use_y2) then
        this_min2 = min(this_min2, minval(graph%curve(i)%y_line))
        this_max2 = max(this_max2, maxval(graph%curve(i)%y_line))
        found_data2 = .true.
      else
        this_min = min(this_min, minval(graph%curve(i)%y_line))
        this_max = max(this_max, maxval(graph%curve(i)%y_line))
        found_data = .true.
      endif
    endif
  endif

enddo

if (.not. found_data) then
  this_max = 10
  this_min = -10
endif

if (this_max >  1d252) this_max =  1d252
if (this_min < -1d252) this_min = -1d252

if (.not. found_data2) then
  this_max2 = 10
  this_min2 = -10
endif

if (this_max2 >  1d252) this_max2 =  1d252
if (this_min2 < -1d252) this_min2 = -1d252

if (axis == '' .or. axis == 'y') call qp_calc_axis_scale (this_min, this_max, graph%y)
if (axis == '' .or. axis == 'y2') call qp_calc_axis_scale (this_min2, this_max2, graph%y2)

end subroutine

end module
