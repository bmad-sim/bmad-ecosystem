module tao_x_scale_mod

use tao_mod
use quick_plot

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_cmd (where, x_min, x_max, err, gang)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   where -- Character(*): Region to scale. Eg: "top"
!   x_min -- Real(rp): Plot x-axis min value.
!   x_max -- Real(rp): Plot x-axis max value.
!   gang  -- Character(*), optional: 'gang', 'nogang', ''. Default = ''.
!
!  Output:
!   err -- Logical: Set to True if the plot cannot be found. False otherwise.
!-

subroutine tao_x_scale_cmd (where, x_min, x_max, err, gang)

implicit none

type (tao_plot_struct), pointer :: p, p2
type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)

real(rp) x_min, x_max

integer i, j, n, ix, places, im

character(*) where
character(*), optional :: gang

logical err

! find plots to scale

if (allocated(plot)) deallocate(plot)
if (allocated(graph)) deallocate(graph)

if (len_trim(where) == 0 .or. where == 'all' .or. where == 's') then
  n = 0
  do j = 1, size(s%plot_region)
    if (.not. s%plot_region(j)%visible) cycle
    if (where == 's' .and. s%plot_region(j)%plot%x_axis_type /= 's') cycle
    n = n + 1
  enddo
  allocate (plot(n))
  n = 0
  do j = 1, size(s%plot_region)
    if (.not. s%plot_region(j)%visible) cycle
    if (where == 's' .and. s%plot_region(j)%plot%x_axis_type /= 's') cycle
    n = n + 1
    plot(n)%p => s%plot_region(j)%plot
  enddo
else
  call tao_find_plots (err, where, 'REGION', plot, graph)
  if (err) return
endif

! If autoscaling then figure out the limits of the data by calling tao_plot_setup.

if (x_min == x_max) then
  if (allocated(graph)) then
    do i = 1, size(graph)
      graph(i)%g%x%min = -1e30  ! So no cliping of points
      graph(i)%g%x%max = 1e30   ! So no cliping of points
    enddo
  else
    do i = 1, size(plot)
      do j = 1, size(plot(i)%p%graph)
        plot(i)%p%graph(j)%x%min = -1e30  ! So no cliping of points
        plot(i)%p%graph(j)%x%max = 1e30   ! So no cliping of points
      enddo
    enddo
  endif
  call tao_plot_setup()
endif

! Set max and min

if (allocated(graph)) then
  do i = 1, size(graph)
    call tao_x_scale_graph (graph(i)%g, x_min, x_max)
  enddo
else
  do i = 1, size(plot)
    call tao_x_scale_plot (plot(i)%p, x_min, x_max, gang)
  enddo
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_plot (plot, x_min, x_max, gang)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   plot  -- Tao_plot_struct: Plot to scale. Eg: "top"
!   x_min -- Real(rp): Plot x-axis min value.
!   x_max -- Real(rp): Plot x-axis max value.
!   gang  -- Character(*), optional: 'gang', 'nogang', ''. Default = ''.
!-

subroutine tao_x_scale_plot (plot, x_min, x_max, gang)

type (tao_plot_struct), target :: plot

real(rp) x_min, x_max, this_min, this_max
integer i, j, p1, p2
character(*), optional :: gang
character(16) :: r_name = 'tao_x_scale_plot'
logical do_gang

! Check if the thing exists

if (.not. allocated (plot%graph)) return

!

do j = 1, size(plot%graph)
  call tao_x_scale_graph (plot%graph(j), x_min, x_max)
enddo

! if auto scale was done...

do_gang = plot%autoscale_gang_x
if (present(gang)) then
  if (gang == 'gang') do_gang = .true.
  if (gang == 'nogang') do_gang = .false.
  if (gang /= '') then
    call out_io (s_error$, r_name, 'BAD GANG SWITCH: ' // gang)
    call err_exit
  endif
endif

if (x_min == x_max .and. do_gang) then
  this_min = minval (plot%graph(:)%x%min)
  this_max = maxval (plot%graph(:)%x%max)
  p1 = nint(0.7 * plot%x%major_div_nominal)  
  p2 = nint(1.3 * plot%x%major_div_nominal)
  call qp_calc_and_set_axis ('X', this_min, this_max, p1, p2, 'GENERAL', plot%x%type)
  call qp_get_axis ('X', plot%x%min, plot%x%max, plot%x%major_div, plot%x%places)
  do i = 1, size(plot%graph)
    plot%graph(i)%x = plot%x
  enddo
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_x_scale_graph (graph, x_min, x_max)

type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (floor_position_struct) end
type (lat_struct), pointer :: lat

integer i, j, k, n, p1, p2, iu, ix
real(rp) x_min, x_max
real(rp) this_min, this_max
logical curve_here

! If specific min/max values are given then life is easy.

if (graph%type == 'key_table') return

if (x_max /= x_min) then

  p1 = nint(0.6 * graph%x%major_div_nominal)  
  p2 = nint(1.4 * graph%x%major_div_nominal)  
  graph%x%min = x_min
  graph%x%max = x_max
  call qp_calc_axis_divisions (x_min, x_max, p1, p2, graph%x%major_div)
  call qp_calc_axis_places (graph%x)
  call qp_set_axis ('X', graph%x%min, graph%x%max, graph%x%major_div, graph%x%places)

  return

endif

! Auto scale

this_min =  1e20
this_max = -1e20
curve_here = .false.

if (graph%type == 'floor_plan') then
  ix = tao_universe_number(graph%ix_universe)
  lat => s%u(ix)%model%lat
  this_min = 1e30
  this_max = -1e30
  do i = 0, lat%n_ele_track
    call floor_to_screen_coords (lat%ele(i)%floor, end)
    this_min = min(this_min, end%x)
    this_max = max(this_max, end%x)
  enddo
  curve_here = .true.
else if (graph%p%x_axis_type == 's') then
  iu = tao_universe_number(graph%ix_universe)
  this_min = min (this_min, s%u(iu)%model%lat%ele(0)%s)
  ix = s%u(iu)%model%lat%n_ele_track
  this_max = max (this_max, s%u(iu)%model%lat%ele(ix)%s)
  curve_here = .true.
else
  if (.not. allocated(graph%curve)) return
  do k = 1, size(graph%curve)
    curve => graph%curve(k)
    if (allocated (curve%x_symb)) then
      curve_here = .true.
      this_min = min (this_min, minval(curve%x_symb(:)))
      this_max = max (this_max, maxval(curve%x_symb(:)))
    endif
    if (allocated (curve%x_line)) then
      curve_here = .true.
      this_min = min (this_min, minval(curve%x_line(:)))
      this_max = max (this_max, maxval(curve%x_line(:)))
    endif
  enddo
endif

if (.not. curve_here) then
  this_max = 100
  this_min = 0
endif

p1 = nint(0.7 * graph%x%major_div_nominal)  
p2 = nint(1.3 * graph%x%major_div_nominal)
call qp_calc_and_set_axis ('X', this_min, this_max, p1, p2, 'GENERAL', graph%x%type)
call qp_get_axis ('X', graph%x%min, graph%x%max, graph%x%major_div, graph%x%places)
  
end subroutine 

end module
