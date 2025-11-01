module tao_scale_mod

use tao_interface
use quick_plot

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_scale_cmd (where, y_min_in, y_max_in, axis, include_wall, gang, exact, turn_autoscale_off)
!
! Routine to scale a plot. 
! If y_min = y_max, the scales will be chosen to show all the data.
! 
! Input:
!   where              -- character(*): Region to scale. Eg: "top:x"
!   y_min_in           -- real(rp): Plot y-axis min value.
!   y_max_in           -- real(rp): Plot y-axis max value.
!   axis               -- character(*), optional: 'y', 'y2', or '' (both). Default = ''.
!   include_wall       -- logical, optional: Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in.
!                           If present and True include the building wall position will be included in determining the the scale.
!   gang               -- Character(*), optional: 'gang', 'nogang', ''. Default = ''.
!   exact              -- logical, optional: Exact plot y_max, y_min to correspond to y_min_in, y_max_in?
!                           Default is False. Only relavent when y_min_in /= y_max_in.
!   turn_autoscale_off -- logical, optional: If True (default) then turn off plot%autoscale_y logical
!                           for all plots that are scaled.
!-

subroutine tao_scale_cmd (where, y_min_in, y_max_in, axis, include_wall, gang, exact, turn_autoscale_off)

implicit none

type (tao_plot_array_struct), allocatable, target :: plot(:)
type (tao_plot_struct), pointer :: plt
type (tao_graph_array_struct), allocatable, target :: graph(:)
type (tao_graph_struct), pointer :: gph

real(rp) y_min_in, y_max_in, y_min, y_max

integer i, j, ix, places, p1, p2

character(*) where
character(*), optional :: axis, gang
character(*), parameter :: r_name = 'tao_scale_cmd'

logical, optional :: turn_autoscale_off, include_wall, exact
logical err

! Error check

if (present(axis)) then
  if (axis /= 'y' .and. axis /= 'y2' .and. axis /= '') then
    call out_io (s_error$, r_name, 'BAD AXIS NAME: ' // axis)
    return
  endif
endif

! Use local vars for y_min and y_max in case the actual args are something 
! like graph%y%min, etc.

y_min = y_min_in
y_max = y_max_in

! If the where argument is blank or '*', then scale all plots.
! Exception: lat_layout

if (len_trim(where) == 0 .or. where == '*' .or. where == 'all') then
  do i = 1, size(s%plot_page%region)
    plt => s%plot_page%region(i)%plot
    if (.not. s%plot_page%region(i)%visible)cycle
    if (.not. allocated (plt%graph)) cycle
    do j = 1, size(plt%graph)
      call set_this_exact(plt%graph(j), y_min_in, y_max_in, exact, axis)
    enddo
    call tao_scale_plot (plt, y_min, y_max, axis, include_wall, gang, .true.)
  enddo
  return
endif

! Locate the plot by the region name given by the where argument.
! If no graph is specified then we scale all the graphs of the plot.

call tao_find_plots (err, where, 'REGION', plot, graph)
if (err) return

if (size(graph) > 0) then                ! If all the graphs of a plot...
  do j = 1, size(graph)
    gph => graph(j)%g
    call set_this_exact (gph, y_min_in, y_max_in, exact, axis)
    call tao_scale_graph (gph, y_min, y_max, axis, include_wall)
    if (logic_option(.true., turn_autoscale_off)) gph%p%autoscale_y = .false.
  enddo
else                          ! else just the one graph...
  do i = 1, size(plot)
    plt => plot(i)%p
    if (.not. allocated (plt%graph)) cycle
    do j = 1, size(plt%graph)
      call set_this_exact(plt%graph(j), y_min_in, y_max_in, exact, axis)
    enddo
    call tao_scale_plot (plt, y_min, y_max, axis, include_wall, gang)
    if (logic_option(.true., turn_autoscale_off)) plt%autoscale_y = .false.
  enddo
endif

!---------------------------------------------------------
contains

subroutine set_this_exact (graph, y_min_in, y_max_in, exact, axis)

type (tao_graph_struct) graph
real(rp) y_min_in, y_max_in
logical, optional :: exact
character(*), optional :: axis

!

if (.not. present(axis)) then
  call set_this_bounds (graph%y, y_min_in, y_max_in, exact)
  call set_this_bounds (graph%y2, y_min_in, y_max_in, exact)
  return
endif

select case (axis)
case ('')
  call set_this_bounds (graph%y, y_min_in, y_max_in, exact)
  call set_this_bounds (graph%y2, y_min_in, y_max_in, exact)
case ('y');     call set_this_bounds (graph%y, y_min_in, y_max_in, exact)
case ('y2');    call set_this_bounds (graph%y2, y_min_in, y_max_in, exact)
end select

end subroutine set_this_exact

!---------------------------------------------------------
! contains

subroutine set_this_bounds (axis, y_min_in, y_max_in, exact)

type (qp_axis_struct) axis
real(rp) y_min_in, y_max_in
logical, optional :: exact

!

if (y_min_in == y_max_in) then
  if (axis%bounds == 'EXACT') axis%bounds = 'GENERAL'
else
  if (logic_option(.false., exact)) axis%bounds = 'EXACT'
endif

end subroutine set_this_bounds

end subroutine tao_scale_cmd

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_scale_plot (plot, y_min_in, y_max_in, axis, include_wall, gang, skip_lat_layout)
!
! Routine to scale the y-axis and/or y2-axis of the graphs of the plot.
! If y_min_in = y_max_in then autoscaling will be done and the particular value
! of y_min_in and y_max_in is ignored.
!
! Input:
!   plot            -- tao_plot_struct: Plot with graphs to be scaled.
!   y_min_in        -- real(rp): Axis [min, max] must cover [y_min_in, y_max_in] if not autoscaling.
!   y_max_in        -- real(rp): Axis [min, max] must cover [y_min_in, y_max_in] if not autoscaling.
!   axis            -- character(*), optional: Axis to scale.
!                         ''   -> scale y and y2 (default).
!                         'y'  -> scale y-axis.
!                         'y2' -> scale y2-axis
!   include_wall    -- logical, optional: Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in.
!                        If present and True include the building wall position will be included in determining the the scale.
!   gang            -- character(*), optional: If autoscale then make all graph y-axes the same and/or
!                         make all y2-axes the same? 
!                         ''        -> (default) Use setting of plot%autoscale_gang_y
!                         'gang'    -> Gang graphs.
!                         'nogang'  -> Do not gang graphs.
!   skip_lat_layout -- logical, optional: If True, skip scaling any lat_layout graphs. Default is false.
!
! Output:
!   plot            -- tao_plot_struct: Plot with scaled graphs.
!-

subroutine tao_scale_plot (plot, y_min_in, y_max_in, axis, include_wall, gang, skip_lat_layout)

type (tao_plot_struct), target :: plot
type (tao_graph_struct), pointer :: graph
type (qp_axis_struct) axis_save

real(rp) y_min_in, y_max_in, y_min, y_max
real(rp) y_range(2), y2_range(2)

character(*), optional :: axis, gang
character(16) this_axis
character(*), parameter :: r_name = 'tao_scale_plot'
integer i, p1, p2

logical, optional :: include_wall, skip_lat_layout
logical do_gang, found_one

! Use local vars in case the actual args are something like graph%y%min, etc.

y_min = y_min_in
y_max = y_max_in

! If we scale a whole plot with auto scale then at the end all graphs
! are adjusted to have the same scale such that all the data fits on
! all the graphs.

if (.not. allocated (plot%graph)) return

y_range = [1d30, -1d30]
y2_range = [1d30, -1d30]

found_one = .false.
do i = 1, size(plot%graph)
  if (logic_option(.false., skip_lat_layout) .and. plot%graph(i)%type == 'lat_layout') cycle
  found_one = .true.
  call tao_scale_graph (plot%graph(i), y_min, y_max, axis, include_wall, y_range, y2_range)
enddo

if (.not. found_one) return

! If auto scale was done...

this_axis = string_option ('', axis)

do_gang = plot%autoscale_gang_y
if (present(gang)) then
  if (gang == 'gang') then
    do_gang = .true.
  elseif (gang == 'nogang') then
    do_gang = .false.
  elseif (gang /= '') then
    call out_io (s_error$, r_name, 'BAD GANG SWITCH: ' // gang)
    return
  endif
endif

if (y_min == y_max .and. do_gang) then

  if (this_axis == '' .or. this_axis == 'y') then
    do i = 1, size(plot%graph)
      graph => plot%graph(i)
      if (graph%y%major_div_nominal > 0) then
        p1 = nint(0.7 * graph%y%major_div_nominal)  
        p2 = nint(1.3 * graph%y%major_div_nominal)
        call qp_calc_axis_params (y_range(1), y_range(2), p1, p2, graph%y)
      else
        call qp_calc_axis_scale (y_range(1), y_range(2), graph%y)
      endif
    enddo
  endif

  if (this_axis == '' .or. this_axis == 'y2') then
    do i = 1, size(plot%graph)
      graph => plot%graph(i)
      if (graph%y2_mirrors_y) then
        axis_save = graph%y2
        graph%y2 = graph%y
        graph%y2%label        = axis_save%label
        graph%y2%draw_label   = axis_save%draw_label
        graph%y2%draw_numbers = axis_save%draw_numbers
      else
        call qp_calc_axis_scale (y2_range(1), y2_range(2), graph%y2)
        if (graph%y2%major_div_nominal > 0) then
          p1 = nint(0.7 * graph%y2%major_div_nominal)  
          p2 = nint(1.3 * graph%y2%major_div_nominal)
          call qp_calc_axis_params (y2_range(1), y2_range(2), p1, p2, graph%y2)
        else
          call qp_calc_axis_scale (y2_range(1), y2_range(2), graph%y2)
        endif
      endif
    enddo
  endif

endif

end subroutine tao_scale_plot

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine tao_scale_graph (graph, y_min, y_max, axis, include_wall, y_range, y2_range)
!
! Routine to scale the y-axis and/or y2-axis of a graph
! If y_min = y_max then autoscaling will be done and the particular value of y_min and y_max is ignored.
! Note: y_min/y_max is ignored if scaling y2-axis and graph%y2_mirrors_y = T.
!
! Input:
!   graph           -- tao_graph_struct: Graph with axis/axes to be scaled.
!   y_min           -- real(rp): Axis [min, max] must cover [y_min, y_max] if not autoscaling.
!   y_max           -- real(rp): Axis [min, max] must cover [y_min, y_max] if not autoscaling.
!   axis            -- character(*), optional: Axis to scale.
!                        ''   -> scale y and y2 (default).
!                        'y'  -> scale y-axis.
!                        'y2' -> scale y2-axis
!   include_wall    -- logical, optional: Used for floor_plan plots where a building wall is drawn and y_min_in = y_max_in.
!                        If present and True include the building wall position will be included in determining the the scale.
!
! Output:
!   graph           -- tao_graph_struct: Graph with scaled axis/axes.
!   y_range(2)      -- real(rp), optional: Only used by tao_scale_plot when ganging graphs.
!   y2_range(2)     -- real(rp), optional: Only used by tao_scale_plot when ganging graphs.
!-

subroutine tao_scale_graph (graph, y_min, y_max, axis, include_wall, y_range, y2_range)

type (tao_graph_struct), target :: graph
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele
type (tao_ele_shape_struct), pointer :: shape
type (tao_curve_struct), pointer :: curve
type (floor_position_struct) floor, end
type (qp_axis_struct) axis_save
type (tao_building_wall_point_struct) pt

real(rp), optional :: y_range(2), y2_range(2)
real(rp) y_min, y_max, this_min, this_max, this_min2, this_max2, del, y1, y2

integer i, j, k, ix, ib, p1, p2, iu
logical, optional :: include_wall
logical found_data, found_data2

character(*), optional :: axis
character(4) this_axis
character(40) label_name
character(*), parameter :: r_name = 'tao_scale_graph'

! If specific min/max values are given then life is easy.

this_axis = string_option ('', axis)

if (y_min /= y_max) then

  if (this_axis == '' .or. this_axis == 'y') then
    if (graph%y%major_div_nominal> 0) then
      p1 = nint(0.7 * graph%y%major_div_nominal)  
      p2 = nint(1.3 * graph%y%major_div_nominal)  
    else
      p1 = graph%y%major_div
      p2 = p1
    endif
    call qp_calc_axis_params (y_min, y_max, p1, p2, graph%y)
  endif

  ! y2

  if (graph%y2_mirrors_y) then
    axis_save = graph%y2
    graph%y2 = graph%y
    graph%y2%label        = axis_save%label
    graph%y2%draw_label   = axis_save%draw_label
    graph%y2%draw_numbers = axis_save%draw_numbers

  elseif (this_axis == '' .or. this_axis == 'y2') then
    if (graph%y2%major_div_nominal> 0) then
      p1 = nint(0.7 * graph%y2%major_div_nominal)  
      p2 = nint(1.3 * graph%y2%major_div_nominal)  
    else
      p1 = graph%y2%major_div
      p2 = p1
    endif
    call qp_calc_axis_params (y_min, y_max, p1, p2, graph%y2)
  endif

  return
endif

! Here if y_min = y_max. Need to autoscale.
! That is we need to find the min/max so all the data points are within bounds.

! For a floor plan 

if (graph%type == 'floor_plan') then
  ix = tao_universe_index(graph%ix_universe)
  this_min = 1e30
  this_max = -1e30
  found_data = .false.

  do iu = 1, ubound(s%u, 1)
    if (ix /= -2 .and. ix /= iu) cycle
    lat => s%u(iu)%model%lat
    do ib = 0, ubound(lat%branch, 1)
      do i = 0, lat%branch(ib)%n_ele_track
        ele => lat%branch(ib)%ele(i)
        call tao_floor_to_screen_coords (graph, ele%floor, end)
        if (end%r(1) > graph%x%max .or. end%r(1) < graph%x%min) cycle
        this_min = min(this_min, end%r(2))
        this_max = max(this_max, end%r(2))
        found_data = .true.
        ! For a linac, the shape can extend past the plot. So take into account the shape size.
        call tao_ele_shape_info (iu, ele, s%plot_page%floor_plan%ele_shape, shape, label_name, y1, y2)
        if (associated(shape)) then
          y1 = y1 * s%plot_page%floor_plan_shape_scale
          y2 = y2 * s%plot_page%floor_plan_shape_scale
          this_min = min(this_min, end%r(2) - y1, end%r(2) - y2)
          this_max = max(this_min, end%r(2) + y1, end%r(2) + y2)
        endif
      enddo
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
          if (end%r(1) > graph%x%max .or. end%r(1) < graph%x%min) cycle
          this_min = min(this_min, end%r(2))
          this_max = max(this_max, end%r(2))
          found_data = .true.
        enddo
      enddo
    enddo
  endif

  if (.not. found_data) then
    call out_io (s_error$, r_name, 'CANNOT Y-SCALE FLOOR_PLAN SINCE NO PART OF THE LATTICE/BUILDING_WALL', &
                                   'IS WITHIN THE GRAPH X-AXIS RANGE. CONSIDER USING xy_scale.')
    if (graph%y%min == graph%y%max) then
      this_min = -10  
      this_max = 10
    else
      if (present(y_range)) y_range = [graph%y%min, graph%y%max]
      return
    endif
  endif

! For a lat_layout: Default is [-1, 1]

elseif (graph%type == 'lat_layout') then
    this_min = -1
    this_max = 1

! Else not a floor_plan nor lat_layout.

else
  if (.not. allocated (graph%curve)) return
  if (.not. any(graph%curve%valid) .or. .not. graph%is_valid) return  ! Do not scale until there is valid data

  this_min =  1e30
  this_max = -1e30
  this_min2 =  1e30
  this_max2 = -1e30
  found_data = .false.
  found_data2 = .false.

  do i = 1, size(graph%curve)
    curve => graph%curve(i)

    if (allocated(curve%y_symb) .and. curve%draw_symbols) then
      if (size(curve%y_symb) > 0) then
        if (curve%use_y2) then
          this_min2 = min(this_min2, minval(curve%y_symb))
          this_max2 = max(this_max2, maxval(curve%y_symb))
          found_data2 = .true.
        else
          this_min = min(this_min, minval(curve%y_symb))
          this_max = max(this_max, maxval(curve%y_symb))
          found_data = .true.
        endif
      endif
    endif

    if (allocated(curve%y_line) .and. curve%draw_line) then
      if (size(curve%y_line) > 0) then
        if (curve%use_y2) then
          this_min2 = min(this_min2, minval(curve%y_line))
          this_max2 = max(this_max2, maxval(curve%y_line))
          found_data2 = .true.
        else
          this_min = min(this_min, minval(curve%y_line))
          this_max = max(this_max, maxval(curve%y_line))
          found_data = .true.
        endif
      endif
    endif
  enddo

  if (.not. found_data) then
    call out_io (s_error$, r_name, 'CANNOT SCALE GRAPH ' // trim(tao_graph_name(graph)) // &
                                                       ' SINCE NO DATA IS WITHIN THE GRAPH X-AXIS RANGE.', &
                                   'USE x_scale TO FIRST SCALE THE X-AXIS.')
    if (graph%y%min /= graph%y%max) return
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

endif

if (graph%y%major_div_nominal > 0) then
  p1 = nint(0.7 * graph%y%major_div_nominal)  
  p2 = nint(1.3 * graph%y%major_div_nominal)
else
  p1 = graph%y%major_div
  p2 = p1
endif

del = this_max - this_min
this_min = this_min - graph%scale_margin%y1 * del
this_max = this_max + graph%scale_margin%y2 * del

if (this_axis == '' .or. this_axis == 'y') then
  call qp_calc_axis_params (this_min, this_max, p1, p2, graph%y)
  if (present(y_range)) y_range = [min(y_range(1), this_min), max(y_range(2), this_max)]
endif

! y2

if (graph%y2_mirrors_y) then
  axis_save = graph%y2
  graph%y2 = graph%y
  graph%y2%label        = axis_save%label
  graph%y2%draw_label   = axis_save%draw_label
  graph%y2%draw_numbers = axis_save%draw_numbers

elseif (this_axis == '' .or. this_axis == 'y2') then
  call qp_calc_axis_params (this_min, this_max, p1, p2, graph%y2)
  if (present(y2_range)) y2_range = [min(y2_range(1), this_min), max(y2_range(2), this_max)]
endif

end subroutine tao_scale_graph

end module
