module tao_x_scale_mod

use tao_interface
use quick_plot
use tao_graph_setup_mod

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_cmd (where, x_min_in, x_max_in, err, include_wall, gang, exact, turn_autoscale_off)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   where              -- Character(*): Region to scale. Eg: "top"
!   x_min_in           -- Real(rp): Plot x-axis min value.
!   x_max_in           -- Real(rp): Plot x-axis max value.
!   include_wall       -- logical, optional: Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in.
!                           If present and True include the building wall position will be included in determining the the scale.
!   gang               -- Character(*), optional: 'gang', 'nogang', ''. Default = ''.
!   exact              -- logical, optional: Exact plot y_max, y_min to correspond to y_min_in, y_max_in?
!                           Default is False. Only relavent when y_min_in /= y_max_in.
!   turn_autoscale_off -- Logical, optional: If True (default) then turn off plot%autoscale_x logical for all plots that are scaled.
!
!  Output:
!   err                -- Logical: Set to True if the plot cannot be found. False otherwise.
!-

subroutine tao_x_scale_cmd (where, x_min_in, x_max_in, err, include_wall, gang, exact, turn_autoscale_off)

implicit none

type (tao_plot_struct), pointer :: p, p2
type (tao_plot_array_struct), allocatable :: plot(:)
type (tao_graph_array_struct), allocatable :: graph(:)

real(rp) x_min_in, x_max_in, x_min, x_max, x0, x1

integer i, j, n, ix, places, im

character(*) where
character(*), optional :: gang
character(20) :: r_name = 'tao_x_scale_cmd'

logical, optional :: include_wall, exact, turn_autoscale_off
logical err, all_same, have_scaled

! Use local vars for x_min and x_max in case the actual args are something 
! like graph%x%min, etc.

x_min = x_min_in
x_max = x_max_in

! find plots to scale

if (len_trim(where) == 0 .or. where == '*' .or. where == 'all' .or. where == 's') then
  n = 0
  do j = 1, size(s%plot_page%region)
    if (.not. s%plot_page%region(j)%visible) cycle
    if (where == 's' .and. s%plot_page%region(j)%plot%x_axis_type /= 's') cycle
    n = n + 1
  enddo
  allocate (plot(n))
  n = 0
  do j = 1, size(s%plot_page%region)
    if (.not. s%plot_page%region(j)%visible) cycle
    if (where == 's' .and. s%plot_page%region(j)%plot%x_axis_type /= 's') cycle
    n = n + 1
    plot(n)%p => s%plot_page%region(j)%plot
  enddo
  allocate(graph(0))

else
  call tao_find_plots (err, where, 'REGION', plot, graph, only_visible = .false.)
  if (err) return
endif

! Set max and min.

all_same = .true.

if (size(graph) > 0) then
  x0 = real_garbage$
  do i = 1, size(graph)
    call set_this_exact (graph(i)%g, x_min_in, x_max_in, exact)
    call tao_x_scale_graph (graph(i)%g, x_min, x_max, include_wall, have_scaled)
    if (.not. have_scaled) cycle
    if (graph(i)%g%p%type == 'floor_plan') cycle

    if (x0 == real_garbage$) then
      x0 = graph(i)%g%x%min
      x1 = graph(i)%g%x%max
    else
      if (graph(i)%g%x%min /= x0 .or. graph(i)%g%x%max /= x1) all_same = .false.
    endif

    if (logic_option(.true., turn_autoscale_off)) graph(i)%g%p%autoscale_x = .false.
  enddo

else
  x0 = real_garbage$
  do i = 1, size(plot)
    if (.not. allocated(plot(i)%p%graph)) cycle
    do j = 1, size(plot(i)%p%graph)
      call set_this_exact (plot(i)%p%graph(j), x_min_in, x_max_in, exact)
      call tao_x_scale_graph (plot(i)%p%graph(j), x_min, x_max, include_wall, have_scaled)
    enddo
    call tao_x_scale_plot (plot(i)%p, x_min, x_max, include_wall, gang, have_scaled)
    if (.not. have_scaled) cycle
    if (plot(i)%p%type == 'floor_plan') cycle

    do j = 1, size(plot(i)%p%graph)
      if (x0 == real_garbage$) then
        x0 = plot(i)%p%graph(j)%x%min
        x1 = plot(i)%p%graph(j)%x%max
      else
        if (plot(i)%p%graph(j)%x%min /= x0 .or. plot(i)%p%graph(j)%x%max /= x1) all_same = .false.
      endif
    enddo

    if (logic_option(.true., turn_autoscale_off)) plot(i)%p%autoscale_x = .false.
  enddo
endif

! Issue warning if not all plots have the same scale

if (.not. all_same) call out_io (s_warn$, r_name, &
      'Note: Not all plots have the same min/max due to different ', &
      '      x-axis major_div or major_div_nominal values.')

!---------------------------------------------------------------------------
contains

subroutine set_this_exact (graph, x_min_in, x_max_in, exact)

type (tao_graph_struct) graph
real(rp) x_min_in, x_max_in
logical, optional :: exact

!

if (x_min_in == x_max_in) then
  if (graph%x%bounds == 'EXACT') graph%x%bounds = 'GENERAL'
else
  if (logic_option(.false., exact)) graph%x%bounds = 'EXACT'
endif

end subroutine set_this_exact

end subroutine tao_x_scale_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_x_scale_plot (plot, x_min_in, x_max_in, include_wall, gang, have_scaled)
!
! Routine to scale a plot. If x_min = x_max
! Then the scales will be chosen to show all the data.
! 
! Input:
!   plot          -- Tao_plot_struct: Plot to scale. Eg: "top"
!   x_min_in      -- Real(rp): Plot x-axis min value.
!   x_max_in      -- Real(rp): Plot x-axis max value.
!   include_wall  -- logical, optional: Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in.
!                      If present and True include the building wall position will be included in determining the the scale.
!   gang          -- Character(*), optional: 'gang', 'nogang', ''. Default = ''.
!
! Output:
!   have_scaled   -- Logical, optional: Has a graph been scaled?
!-

subroutine tao_x_scale_plot (plot, x_min_in, x_max_in, include_wall, gang, have_scaled)

type (tao_plot_struct), target :: plot
type (tao_graph_struct), pointer :: graph

real(rp) x_min_in, x_max_in, x_min, x_max, this_min, this_max, major_div_nominal
integer i, j, p1, p2, n
character(*), optional :: gang
character(16) :: r_name = 'tao_x_scale_plot'
logical, optional :: include_wall, have_scaled
logical do_gang, scaled, valid

! Check if the thing exists

if (present(have_scaled)) have_scaled = .false.
if (.not. allocated (plot%graph)) return
if (size(plot%graph) == 0) return

! Use local vars in case the actual args are something like graph%x%min, etc.

x_min = x_min_in
x_max = x_max_in

!

valid = .true.
do j = 1, size(plot%graph)
  graph => plot%graph(j)
  call tao_x_scale_graph (graph, x_min, x_max, include_wall, scaled)
  if (present(have_scaled)) have_scaled = (have_scaled .or. scaled)
  if (.not. graph%is_valid) valid = .false.
  if (allocated(graph%curve)) then
    if (.not. any(graph%curve%valid)) valid = .false.
  endif
enddo

! if ganging is needed...

do_gang = plot%autoscale_gang_x
if (present(gang)) then
  select case (gang)
  case ('gang');   do_gang = .true.
  case ('nogang'); do_gang = .false.
  case ('')
  case default
    call out_io (s_error$, r_name, 'BAD GANG SWITCH: ' // gang)
    return
  end select
endif

if (do_gang .and. valid) then
  if (x_min == x_max) then
    this_min = 1e30; this_max = -1e30
    n = 0; major_div_nominal = 0
    do i = 1, size(plot%graph)
      graph => plot%graph(i)
      if (graph%type == 'key_table') cycle
      n = n + 1
      this_min = min (this_min, graph%x%min)
      this_max = max (this_max, graph%x%max)
      major_div_nominal = major_div_nominal + graph%x%major_div_nominal
    enddo
    if (n == 0) return
    major_div_nominal = major_div_nominal / n

    if (major_div_nominal > 0) then
      p1 = nint(0.7 * major_div_nominal)  
      p2 = nint(1.3 * major_div_nominal)
    else
      p1 = real(sum(plot%graph(:)%x%major_div)) / size(plot%graph)
      p2 = p1
    endif

    do i = 1, size(plot%graph)
      graph => plot%graph(i)
      if (graph%type == 'key_table') cycle
      call qp_calc_axis_params (this_min, this_max, p1, p2, graph%x)
      graph%x%eval_min = this_min
      graph%x%eval_max = this_max
    enddo
  endif
endif

end subroutine tao_x_scale_plot

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_x_scale_graph (graph, x_min, x_max, include_wall, have_scaled)

implicit none

type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (tao_ele_shape_struct), pointer :: shape
type (tao_universe_struct), pointer :: u
type (floor_position_struct) floor, end
type (lat_struct), pointer :: lat
type (tao_building_wall_point_struct) pt

integer i, j, k, n, p1, p2, ix, ib
real(rp) x_min, x_max
real(rp) this_min, this_max, del
logical, optional :: include_wall, have_scaled
logical curve_here

character(*), parameter :: r_name = 'tao_x_scale_graph'

! If specific min/max values are given then life is easy.

if (present(have_scaled)) have_scaled = .false.
if (graph%type == 'key_table') return
if (present(have_scaled)) have_scaled = .true.

if (x_max /= x_min) then
  if (graph%x%major_div_nominal > 0) then
    p1 = nint(0.7 * graph%x%major_div_nominal)  
    p2 = nint(1.3 * graph%x%major_div_nominal)
  else
    p1 = graph%x%major_div
    p2 = p1
  endif
  graph%x%min = x_min
  graph%x%max = x_max
  graph%x%eval_min = x_min
  graph%x%eval_max = x_max
  call qp_calc_axis_params (x_min, x_max, p1, p2, graph%x)
  return
endif

! Auto scale.

this_min =  1e30
this_max = -1e30
curve_here = .false.

if (graph%type == 'floor_plan') then
  u => tao_pointer_to_universe(graph%ix_universe, .true.)
  lat => u%model%lat

  do ib = 0, ubound(lat%branch, 1)
    do i = 0, lat%branch(ib)%n_ele_track
      call tao_floor_to_screen_coords (graph, lat%branch(ib)%ele(i)%floor, end)
      this_min = min(this_min, end%r(1))
      this_max = max(this_max, end%r(1))
    enddo
  enddo

  if (logic_option(.false., include_wall) .and. allocated(s%building_wall%section)) then
    do i = 1, size(s%plot_page%floor_plan%ele_shape)
      shape => s%plot_page%floor_plan%ele_shape(i)
      if (shape%ele_id(1:15) /= 'building_wall::') cycle
      if (.not. shape%draw) cycle
      do j = 1, size(s%building_wall%section)
        do k = 1, size(s%building_wall%section(j)%point)
          pt = tao_oreint_building_wall_pt(s%building_wall%section(j)%point(k))
          floor%r(1) = pt%x
          floor%r(2) = 0
          floor%r(3) = pt%z
          floor%theta = 0
          call tao_floor_to_screen_coords (graph, floor, end)
          this_min = min(this_min, end%r(1))
          this_max = max(this_max, end%r(1))
        enddo
      enddo
    enddo
  endif

  curve_here = .true.

elseif (graph%type == 'histogram') then
  curve_here = .false.
  do i = 1, size(graph%curve)
    curve => graph%curve(i)
    if (.not. allocated(curve%x_line)) cycle
    curve_here = .true.
    this_min = min(this_min, curve%hist%minimum)
    this_max = max(this_max, curve%hist%maximum)
  enddo
  if (.not. curve_here) return

else if (graph%p%x_axis_type == 's') then
  if (allocated(graph%curve)) then
    do i = 1, size(graph%curve)
      u => tao_pointer_to_universe(tao_curve_ix_uni(graph%curve(i)), .true.)
      if (.not. associated(u)) cycle
      ib = tao_branch_index(graph%curve(i)%ix_branch)
      this_min = min (this_min, u%model%lat%branch(ib)%ele(0)%s)
      ix = u%model%lat%branch(ib)%n_ele_track
      this_max = max (this_max, u%model%lat%branch(ib)%ele(ix)%s)
    enddo
  else
    ib = tao_branch_index(graph%ix_branch)
    u => tao_pointer_to_universe(graph%ix_universe, .true.)
    if (.not. associated(u)) then
      graph%is_valid = .false.
      graph%why_invalid = 'Bad universe index.'
      return
    endif
    this_min = min (this_min, u%model%lat%branch(ib)%ele(0)%s)
    ix = u%model%lat%branch(ib)%n_ele_track
    this_max = max (this_max, u%model%lat%branch(ib)%ele(ix)%s)
  endif

  curve_here = .true.

else if (graph%p%x_axis_type == 'var' .or. graph%p%x_axis_type == 'lat') then
  call out_io (s_error$, r_name, 'CANNOT AUTO X-SCALE A PLOT WITH X_AXIS_TYPE OF: ' // graph%p%x_axis_type, &
                                 'YOU, THE USER, MUST SET THIS.')
  graph%is_valid = .false.
  graph%why_invalid = 'CANNOT AUTO X-SCALE A PLOT WITH X_AXIS_TYPE OF: ' // graph%p%x_axis_type
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
else
  del = this_max - this_min
  this_min = this_min - graph%scale_margin%x1 * del
  this_max = this_max + graph%scale_margin%x2 * del
endif

if (graph%x%major_div_nominal > 0) then
  p1 = nint(0.7 * graph%x%major_div_nominal)  
  p2 = nint(1.3 * graph%x%major_div_nominal)
else
  p1 = graph%x%major_div
  p2 = p1
endif

graph%x%eval_min = this_min
graph%x%eval_max = this_max

call qp_calc_axis_params (this_min, this_max, p1, p2, graph%x)

end subroutine tao_x_scale_graph

end module
