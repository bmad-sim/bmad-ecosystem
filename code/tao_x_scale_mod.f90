module tao_x_scale_mod

use tao_mod
use quick_plot
use tao_graph_setup_mod

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_cmd (where, x_min_in, x_max_in, err, gang, turn_autoscale_off)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   where    -- Character(*): Region to scale. Eg: "top"
!   x_min_in -- Real(rp): Plot x-axis min value.
!   x_max_in -- Real(rp): Plot x-axis max value.
!   gang     -- Character(*), optional: 'gang', 'nogang', ''. Default = ''.
!   turn_autoscale_off
!            -- Logical, optional: If True (default) then turn off plot%autoscale_x logical
!               for all plots that are scaled.
!
!  Output:
!   err -- Logical: Set to True if the plot cannot be found. False otherwise.
!-

subroutine tao_x_scale_cmd (where, x_min_in, x_max_in, err, gang, turn_autoscale_off)

implicit none

type (tao_plot_struct), pointer :: p, p2
type (tao_plot_array_struct), allocatable, save :: plot(:)
type (tao_graph_array_struct), allocatable, save :: graph(:)

real(rp) x_min_in, x_max_in, x_min, x_max

integer i, j, n, ix, places, im

character(*) where
character(*), optional :: gang
character(20) :: r_name = 'tao_x_scale_cmd'

logical, optional :: turn_autoscale_off
logical err, all_same

! Use local vars for x_min and x_max in case the actual args are something 
! like graph%x%min, etc.

x_min = x_min_in
x_max = x_max_in

! find plots to scale

if (allocated(plot)) deallocate(plot)
if (allocated(graph)) deallocate(graph)

if (len_trim(where) == 0 .or. where == '*' .or. where == 'all' .or. where == 's') then
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

! Set max and min.

all_same = .true.

if (allocated(graph)) then
  do i = 1, size(graph)
    call tao_x_scale_graph (graph(i)%g, x_min, x_max)
    if (graph(i)%g%x%max /= graph(1)%g%x%max) all_same = .false.
    if (graph(i)%g%x%min /= graph(1)%g%x%min) all_same = .false.
    if (logic_option(.true., turn_autoscale_off)) graph(i)%g%p%autoscale_x = .false.
  enddo
else
  do i = 1, size(plot)
    call tao_x_scale_plot (plot(i)%p, x_min, x_max, gang)
    if (plot(i)%p%x%max /= plot(1)%p%x%max) all_same = .false.
    if (plot(i)%p%x%min /= plot(1)%p%x%min) all_same = .false.
    if (logic_option(.true., turn_autoscale_off)) plot(i)%p%autoscale_x = .false.
  enddo
endif

! Issue warning if not all plots have the same scale

if (.not. all_same) call out_io (s_warn$, r_name, &
      'Note: Not all plots have the same min/max due to different ', &
      '      x-axis major_div or major_div_nominal values.')

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_plot (plot, x_min_in, x_max_in, gang)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   plot     -- Tao_plot_struct: Plot to scale. Eg: "top"
!   x_min_in -- Real(rp): Plot x-axis min value.
!   x_max_in -- Real(rp): Plot x-axis max value.
!   gang     -- Character(*), optional: 'gang', 'nogang', ''. Default = ''.
!-

subroutine tao_x_scale_plot (plot, x_min_in, x_max_in, gang)

type (tao_plot_struct), target :: plot
type (tao_graph_struct), pointer :: graph

real(rp) x_min_in, x_max_in, x_min, x_max, this_min, this_max, major_div_nominal
integer i, j, p1, p2
character(*), optional :: gang
character(16) :: r_name = 'tao_x_scale_plot'
logical do_gang

! Check if the thing exists

if (.not. allocated (plot%graph)) return
if (size(plot%graph) == 0) return

! Use local vars in case the actual args are something like graph%x%min, etc.

x_min = x_min_in
x_max = x_max_in

!

do j = 1, size(plot%graph)
  call tao_x_scale_graph (plot%graph(j), x_min, x_max)
enddo

! if auto scale is needed...

do_gang = plot%autoscale_gang_x
if (present(gang)) then
  if (gang == 'gang') do_gang = .true.
  if (gang == 'nogang') do_gang = .false.
  if (gang /= '') then
    call out_io (s_error$, r_name, 'BAD GANG SWITCH: ' // gang)
    call err_exit
  endif
endif

if (do_gang) then

  if (x_min == x_max) then
    this_min = minval (plot%graph(:)%x%min)
    this_max = maxval (plot%graph(:)%x%max)
    major_div_nominal = real(sum(plot%graph(:)%x%major_div_nominal)) / size(plot%graph)
    if (major_div_nominal > 0) then
      p1 = nint(0.7 * major_div_nominal)  
      p2 = nint(1.3 * major_div_nominal)
    else
      p1 = real(sum(plot%graph(:)%x%major_div)) / size(plot%graph)
      p2 = p1
    endif
    do i = 1, size(plot%graph)
      graph => plot%graph(i)
      if (.not. graph%visible) cycle
      call qp_calc_and_set_axis ('X', this_min, this_max, p1, p2, 'GENERAL', graph%x%type)
      call qp_get_axis_attrib ('X', graph%x%min, graph%x%max, graph%x%major_div, graph%x%places)
    enddo
  endif

  plot%x = plot%graph(1)%x

endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_x_scale_graph (graph, x_min, x_max)

implicit none

type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (floor_position_struct) end
type (lat_struct), pointer :: lat

integer i, j, k, n, p1, p2, iu, ix, ib
real(rp) x_min, x_max
real(rp) this_min, this_max
logical curve_here

character(24) :: r_name = 'tao_x_scale_graph'

! If specific min/max values are given then life is easy.

if (graph%type == 'key_table') return
if (.not. graph%visible) return

if (x_max /= x_min) then

  if (graph%x%major_div_nominal> 0) then
    p1 = nint(0.7 * graph%x%major_div_nominal)  
    p2 = nint(1.3 * graph%x%major_div_nominal)  
  else
    p1 = graph%x%major_div
    p2 = p1
  endif
  graph%x%min = x_min
  graph%x%max = x_max
  call qp_calc_axis_divisions (x_min, x_max, p1, p2, graph%x%major_div)
  call qp_calc_axis_places (graph%x)
  call qp_set_axis ('X', graph%x%min, graph%x%max, graph%x%major_div, graph%x%places)

  return

endif

! Auto scale

this_min =  1e30
this_max = -1e30
curve_here = .false.

if (graph%type == 'floor_plan') then
  ix = tao_universe_number(graph%ix_universe)
  lat => s%u(ix)%model%lat
  do ib = 0, ubound(lat%branch, 1)
    do i = 0, lat%branch(ib)%n_ele_track
      call floor_to_screen_coords (lat%branch(ib)%ele(i)%floor, end)
      this_min = min(this_min, end%x)
      this_max = max(this_max, end%x)
    enddo
  enddo
  curve_here = .true.

else if (graph%p%x_axis_type == 's') then
  if (allocated(graph%curve)) then
    do i = 1, size(graph%curve)
      iu = tao_universe_number(graph%curve(i)%ix_universe)
      ib = graph%curve(i)%ix_branch
      this_min = min (this_min, s%u(iu)%model%lat%branch(ib)%ele(0)%s)
      ix = s%u(iu)%model%lat%branch(ib)%n_ele_track
      this_max = max (this_max, s%u(iu)%model%lat%branch(ib)%ele(ix)%s)
    enddo
  else
    ib = graph%ix_branch
    iu = tao_universe_number(graph%ix_universe)
    this_min = min (this_min, s%u(iu)%model%lat%branch(ib)%ele(0)%s)
    ix = s%u(iu)%model%lat%branch(ib)%n_ele_track
    this_max = max (this_max, s%u(iu)%model%lat%branch(ib)%ele(ix)%s)
  endif
  curve_here = .true.

else if (graph%p%x_axis_type == 'var' .or. graph%p%x_axis_type == 'lat') then
  call out_io (s_error$, r_name, 'CANNOT AUTO X-SCALE A PLOT WITH X_AXIS_TYPE OF: ' // graph%p%x_axis_type)
  return

else
  if (.not. allocated(graph%curve)) return
  if (x_min == x_max) then
    graph%x%min = -1e30  ! So no cliping of points
    graph%x%max = 1e30   ! So no cliping of points
    call tao_graph_setup(graph%p, graph)
  endif
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

if (graph%x%major_div_nominal > 0) then
  p1 = nint(0.7 * graph%x%major_div_nominal)  
  p2 = nint(1.3 * graph%x%major_div_nominal)
else
  p1 = graph%x%major_div
  p2 = p1
endif

call qp_calc_and_set_axis ('X', this_min, this_max, p1, p2, 'GENERAL', graph%x%type)
call qp_get_axis_attrib ('X', graph%x%min, graph%x%max, graph%x%major_div, graph%x%places)

end subroutine 

end module
