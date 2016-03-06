module tao_plot_mod

use tao_mod
use quick_plot
use tao_plot_window_mod

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_plots (do_clear)
!
! Subroutine to draw the plots on the plot window.
!
! Input:
!   do_clear -- Logical, optional: If present and False then call qp_clear_page.
!                 This argument is used when drawing PS or GIF.
!-

subroutine tao_draw_plots (do_clear)

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (qp_rect_struct) border1, border2
type (tao_data_array_struct), allocatable, save :: d_array(:)

real(rp) location(4), dx, dy, h

integer i, j, k, ic, id

character(80) text
character(16) :: r_name = 'tao_draw_plots'
character(3) view_str

logical, optional :: do_clear
logical found, err, beam_source

! inits

if (.not. s%global%plot_on) return
call tao_create_plot_window () ! This routine knows not to create multiple windows.
if (logic_option(.true., do_clear)) call qp_clear_page

h = s%plot_page%text_height
call qp_set_text_attrib ('TEXT', height = h)
call qp_set_text_attrib ('MAIN_TITLE', height = h * s%plot_page%main_title_text_scale)
call qp_set_text_attrib ('GRAPH_TITLE', height = h * s%plot_page%graph_title_text_scale)
call qp_set_text_attrib ('LEGEND', height = h * s%plot_page%legend_text_scale)
call qp_set_text_attrib ('AXIS_NUMBERS', height = h * s%plot_page%axis_number_text_scale)
call qp_set_text_attrib ('AXIS_LABEL', height = h * s%plot_page%axis_label_text_scale)

! print the title 

h = s%plot_page%text_height * s%plot_page%main_title_text_scale
do i = 1, size(s%plot_page%title)
  if (s%plot_page%title(i)%draw_it)                                         &
    call qp_draw_text (s%plot_page%title(i)%string, s%plot_page%title(i)%x, &
                     s%plot_page%title(i)%y, s%plot_page%title(i)%units,    &
                     s%plot_page%title(i)%justify, height = h)
enddo

! Draw view universe

if (size(s%u) > 1) then
  write (view_str, '(i3)') s%com%default_universe
  call qp_draw_text ('View Universe:' // view_str, -2.0_rp, -2.0_rp, 'POINTS/PAGE/RT', 'RT')
endif

! loop over all plots

do i = 1, size(s%plot_page%region)

  if (.not. s%plot_page%region(i)%visible) cycle
  plot => s%plot_page%region(i)%plot

  ! set the s%plot_page border for this particular region

  location = s%plot_page%region(i)%location
  border1%units = '%PAGE'
  call qp_convert_rectangle_rel (s%plot_page%border, border1)
  dx = 1 - (border1%x2 - border1%x1)
  dy = 1 - (border1%y2 - border1%y1)
  border2%x1 = border1%x1 + dx * location(1)
  border2%x2 = border1%x2 + dx * (1 - location(2))
  border2%y1 = border1%y1 + dy * location(3)
  border2%y2 = border1%y2 + dy * (1 - location(4))
  border2%units = '%PAGE'
  call qp_set_layout (page_border = border2)

  ! loop over all the graphs of the plot and draw them.

  g_loop: do j = 1, size(plot%graph)

    graph => plot%graph(j)
    if (.not. graph%visible) cycle

    ! For a non-valid graph just print a message

    if (.not. graph%valid) then
      call qp_set_layout (box = graph%box)
      text = 'Error In The Plot Calculation'
      if (graph%why_invalid /= '') text = graph%why_invalid
      call qp_draw_text (text, 0.5_rp, 0.5_rp, '%BOX', color = red$, justify = 'CC')
    endif

    ! Now we can draw the graph

    call tao_hook_draw_graph (plot, graph, found)
    if (found) cycle

    select case (graph%type)
    case ('data', 'phase_space', 'dynamic_aperture')
      call tao_plot_data (plot, graph)
    case ('wave.0', 'wave.a', 'wave.b')
      call tao_plot_wave (plot, graph)
     case ('lat_layout')
      call tao_draw_lat_layout (plot, graph)
    case ('key_table')
      call tao_plot_key_table (plot, graph)
    case ('floor_plan')
      call tao_draw_floor_plan (plot, graph)
    case ('histogram')
      call tao_plot_histogram (plot, graph)
    case default
      call out_io (s_fatal$, r_name, 'UNKNOWN GRAPH TYPE: ' // graph%type)
    end select

    ! Draw a rectangle so the box and graph boundries can be seen.

    if (s%global%box_plots) then
      call qp_draw_rectangle (0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp, '%GRAPH/LB')
      call qp_draw_rectangle (0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp, '%BOX/LB')
    endif

  enddo g_loop

enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_plot_histogram (plot, graph)
!
! Routine to draw one graph for the histogram analysis plot.
!
! Input:
!   plot  -- Tao_plot_struct: Plot containing the graph.
!   graph -- Tao_graph_struct: Graph to plot.
!-

subroutine tao_plot_histogram (plot, graph)

type (tao_plot_struct) :: plot
type (tao_graph_struct), target :: graph

integer k
logical have_data

! Draw the graph outline.

call tao_draw_data_graph (plot, graph)
if (.not. graph%valid) return

! loop over all the curves of the graph and draw them

have_data = .false.

do k = 1, size(graph%curve)
  call tao_draw_histogram_data (plot, graph, graph%curve(k), have_data)
enddo

if (.not. have_data) call qp_draw_text ('**No Plottable Data**', &
                            0.18_rp, -0.15_rp, '%/GRAPH/LT', color = red$) 

end subroutine tao_plot_histogram

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_plot_wave (plot, graph)
!
! Routine to draw one graph for the wave analysis plot.
!
! Input:
!   plot  -- Tao_plot_struct: Plot containing the graph.
!   graph -- Tao_graph_struct: Graph to plot.
!-

subroutine tao_plot_wave (plot, graph)

type (tao_plot_struct) :: plot
type (tao_graph_struct), target :: graph
type (qp_axis_struct), pointer :: y

real(rp) y0, y1

! Draw the data

call tao_plot_data (plot, graph)

! Now draw the rectangles of the fit regions.

y => graph%y
y0 = y%min + 0.1 * (y%max - y%min)
y1 = y%max - 0.1 * (y%max - y%min)

if (graph%type == 'wave.0' .or. graph%type == 'wave.a') then
  call qp_draw_rectangle (1.0_rp * s%wave%ix_a1, 1.0_rp * s%wave%ix_a2, y0, y1, color = blue$, width = 2)
endif

if (graph%type == 'wave.0' .or. graph%type == 'wave.b') then
  call qp_draw_rectangle (1.0_rp * s%wave%ix_b1, 1.0_rp * s%wave%ix_b2, y0, y1, color = blue$, width = 2)
endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_plot_key_table (plot, graph)
!
! Routine to draw a key table graph.
!
! Input:
!   plot  -- Tao_plot_struct: Plot containing the graph.
!   graph -- Tao_graph_struct: Graph to plot.
!-


subroutine tao_plot_key_table (plot, graph)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (tao_var_struct), pointer :: var

integer i, j, k, ix_var, i_off
real(rp) :: y_here, x1, x2, y1, y2
real(rp) :: height

character(120) str, header
character(7) prefix
character(24) :: r_name = 'tao_plot_key_table'

!

call qp_set_layout (box = graph%box, margin = graph%margin)

call qp_get_layout_attrib ('GRAPH', x1, x2, y1, y2, 'POINTS/GRAPH')
y_here = y2  ! start from the top of the graph
height = s%plot_page%text_height * s%plot_page%key_table_text_scale


i_off = s%com%ix_key_bank
call tao_key_info_to_str (1, i_off+1, i_off+10, str, header)
call qp_draw_text ('   Ix  ' // header, 0.0_rp, y_here, 'POINTS/GRAPH', &
                             height = height, uniform_spacing = .true.)
  

do i = 1, 10

  k = i + i_off
  if (k > ubound(s%key, 1)) return

  prefix = ''
  j = mod(i, 10)
  if (i == 1) then
    write (prefix, '(i2, a, i2)') s%com%ix_key_bank, ':', j
  else
    write (prefix(4:), '(i2)') j
  endif

  ix_var = s%key(k)

  call tao_key_info_to_str (i+i_off, i_off+1, i_off+10, str, header)

  y_here = y_here - 1.1 * height
  call qp_draw_text (prefix // str, 0.0_rp, y_here, 'POINTS/GRAPH', &
                          height = height, uniform_spacing = .true.)

enddo

end subroutine tao_plot_key_table

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_floor_plan (plot, graph)
!
! Routine to draw a floor plan graph.
!
! Input:
!   plot  -- Tao_plot_struct: Plot containing the graph.
!   graph -- Tao_graph_struct: Graph to plot.
!-

subroutine tao_draw_floor_plan (plot, graph)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (qp_axis_struct) x_ax, y_ax
type (lat_struct), pointer :: lat
type (tao_ele_shape_struct), pointer :: ele_shape
type (tao_lattice_branch_struct), pointer :: lat_branch
type (floor_position_struct) end1, end2, floor
type (tao_building_wall_point_struct), pointer :: pt(:)
type (ele_struct), pointer :: ele, slave
type (branch_struct), pointer :: branch

real(rp) theta, v_vec(3), theta1, dtheta, dat_var_value
real(rp) x_bend(0:1000), y_bend(0:1000)

integer i, j, k, n, n_bend, isu, ic, ib, icol, ix_shape, ix_shape_min

character(40) dat_var_name
character(20) :: r_name = 'tao_draw_floor_plan'

logical err

! Each graph is a separate floor plan plot (presumably for different universes). 
! setup the placement of the graph on the plot page.

call qp_save_state(.false.)

call set_this_axis(graph%x, x_ax, 'X')
call set_this_axis(graph%y, y_ax, 'Y')

call qp_set_layout (x_axis = x_ax, y_axis = y_ax, y2_axis = graph%y2, x2_mirrors_x = .true., &
                    y2_mirrors_y = .true., box = graph%box, margin = graph%margin)

call qp_set_axis ('X2', draw_numbers = .true.)

if (graph%correct_xy_distortion) call qp_eliminate_xy_distortion

!

if (graph%draw_axes) then
  if (graph%title == '') then
    call qp_set_graph (title = '')
  else
    call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
  endif
  call qp_draw_axes (draw_grid = graph%draw_grid)
endif

isu = tao_universe_number(graph%ix_universe)
lat => s%u(isu)%model%lat

if (.not. graph%valid) return

! loop over all elements in the lattice. 

do n = 0, ubound(lat%branch, 1)
  branch => lat%branch(n)
  branch%ele%logic = .false.  ! Used to mark as drawn.
  do i = 0, branch%n_ele_max
    ele => branch%ele(i)

    ix_shape_min = 1
    do
      ele_shape => tao_pointer_to_ele_shape (isu, ele, s%plot_page%floor_plan%ele_shape(ix_shape_min:), &
                                                                             dat_var_name, dat_var_value, ix_shape)
      if (ele%ix_ele > branch%n_ele_track .and. .not. associated(ele_shape)) exit   ! Nothing to draw
      if (ele%lord_status == multipass_lord$) then
        do j = 1, ele%n_slave
          slave => pointer_to_slave(ele, j)
          call tao_draw_ele_for_floor_plan (plot, graph, isu, lat, slave, dat_var_name, dat_var_value, ele_shape)
        enddo
      else
        call tao_draw_ele_for_floor_plan (plot, graph, isu, lat, ele, dat_var_name, dat_var_value, ele_shape)
      endif
      if (.not. associated(ele_shape)) exit
      if (.not. ele_shape%multi) exit
      ix_shape_min = ix_shape_min + ix_shape
    enddo

  enddo
enddo

! Draw the building wall

if (allocated(s%building_wall%section)) then
  do i = 1, size(s%plot_page%floor_plan%ele_shape)
    ele_shape => s%plot_page%floor_plan%ele_shape(i)
    if (ele_shape%ele_id /= 'wall::building') cycle
    if (.not. ele_shape%draw) cycle
    icol =qp_translate_to_color_index (ele_shape%color)

    do ib = 1, size(s%building_wall%section)
      pt => s%building_wall%section(ib)%point

      do j = 2, size(pt)
        if (pt(j)%radius == 0) then   ! line
          call floor_to_screen (graph, [pt(j-1)%x, 0.0_rp, pt(j-1)%z], end1%r(1), end1%r(2))
          call floor_to_screen (graph, [pt(j)%x, 0.0_rp, pt(j)%z], end2%r(1), end2%r(2))
          call qp_draw_line(end1%r(1), end2%r(1), end1%r(2), end2%r(2), color = icol)

        else                    ! arc
          theta1 = atan2(pt(j-1)%x - pt(j)%x_center, pt(j-1)%z - pt(j)%z_center)
          dtheta = atan2(pt(j)%x - pt(j)%x_center, pt(j)%z - pt(j)%z_center) - theta1
          if (abs(dtheta) > pi) dtheta = modulo2(dtheta, pi)
          n_bend = abs(50 * dtheta) + 1
          do k = 0, n_bend
            theta = theta1 + k * dtheta / n_bend
            v_vec(1) = pt(j)%x_center + abs(pt(j)%radius) * sin(theta)
            v_vec(2) = 0
            v_vec(3) = pt(j)%z_center + abs(pt(j)%radius) * cos(theta)
            call floor_to_screen (graph, v_vec, x_bend(k), y_bend(k))
          enddo
          call qp_draw_polyline(x_bend(:n_bend), y_bend(:n_bend), color = icol)
        endif
      enddo

    enddo
    exit
  enddo
end if

! Draw any data curves and beam chamber wall curve

do i = 1, size(graph%curve)
  ! ... needs to be filled in ..
enddo

! Hook routine for more plotting if desired...

call tao_hook_draw_floor_plan (plot, graph)

call qp_restore_state

!-------------------------------------------------------------------------
contains

subroutine set_this_axis (axis_in, axis_out, which)

type (qp_axis_struct) axis_in, axis_out
real(rp) f
integer irot
character(*) which
character(1) x_str, y_str
character(2) label(0:3)

!

axis_out = axis_in
if (axis_out%label /= 'SMART LABEL') return

f = modulo(4*graph%floor_plan_rotation, 1.0_rp)
f = f - fraction(f)

! If rotation is not multiple of 90 degrees then must use blank label.

if (abs(f) > 0.01) then
  axis_out%label = ''
  return
endif

! Normal case

x_str = upcase(graph%floor_plan_view(1:1))
y_str = upcase(graph%floor_plan_view(2:2))
label = [x_str // ' ', '-' // y_str, '-' // x_str, y_str // ' ']

irot = modulo(nint(4*graph%floor_plan_rotation), 4)
if (which == 'Y') irot = modulo(irot-1, 4)
axis_out%label = label(irot)

end subroutine set_this_axis

end subroutine tao_draw_floor_plan 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_ele_for_floor_plan (plot, graph, ix_uni, lat, ele, dat_var_name, dat_var_value, ele_shape)
!
! Routine to draw one lattice element or one datum location for the floor plan graph. 
!
! Input:
!   plot            -- tao_plot_struct: Plot containing the graph.
!   graph           -- tao_graph_struct: Graph to plot.
!   ix_uni          -- integer: Universe index.
!   lat             -- lat_struct: Lattice containing the element.
!   ele             -- ele_struct: Element to draw.
!   dat_var_name    -- Character(*): If not blank then name to print beside the element.
!   dat_var_value   -- real(rp): Use for vvar_box and asym_vvar_box.
!   ele_shape       -- tao_ele_shape_struct: Shape to draw from s%plot_page%floor_plan%ele_shape(:) array.
!                       Will be NULL if no associated shape for this element.
!-

recursive subroutine tao_draw_ele_for_floor_plan (plot, graph, ix_uni, lat, ele, dat_var_name, dat_var_value, ele_shape)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (lat_struct) :: lat
type (ele_struct) :: ele
type (ele_struct) :: drift
type (ele_struct), pointer :: ele1, ele2, lord
type (floor_position_struct) end1, end2, floor, x_ray
type (tao_building_wall_point_struct), pointer :: pt(:)
type (tao_ele_shape_struct), pointer :: ele_shape, branch_shape

integer ix_uni, i, j, k, icol, n_bend, n, ix, ic, n_mid

real(rp) off, off1, off2, angle, rho, dx1, dy1, dx2, dy2, ang, length
real(rp) dt_x, dt_y, x_center, y_center, dx, dy, theta, e_edge, dat_var_value
real(rp) x_bend(0:1000), y_bend(0:1000), dx_bend(0:1000), dy_bend(0:1000)
real(rp) v_old(3), w_old(3,3), r_vec(3), dr_vec(3), v_vec(3), dv_vec(3)
real(rp) cos_t, sin_t, cos_p, sin_p, cos_a, sin_a, height
real(rp) x_inch, y_inch, x0, y0, x1, x2, y1, y2, e1_factor, e2_factor
real(rp) r0_plus(2), r0_minus(2), dr2_p(2), dr2_m(2), dr_p(2), dr_m(2)
real(rp) x_min, x_max, y_min, y_max

character(*) dat_var_name
character(80) str
character(40) name
character(40) :: r_name = 'tao_draw_ele_for_floor_plan'
character(8) :: draw_units
character(2) justify

logical is_data_or_var, is_there
logical shape_has_box, is_bend

!

call find_element_ends (ele, ele1, ele2)
if (.not. associated(ele1)) return

is_data_or_var = .false.
if (associated(ele_shape)) is_data_or_var = (ele_shape%ele_id(1:6) == 'data::' .or. ele_shape%ele_id(1:5) == 'var::')

if (is_data_or_var) then  ! pretend this is zero length element
  ele1 => ele2
  is_bend = .false.
else
  is_bend = (ele%key == sbend$)
endif

call floor_to_screen_coords (graph, ele1%floor, end1)
call floor_to_screen_coords (graph, ele2%floor, end2)

! Only draw those element that have at least one point in bounds.

! If qp_eliminate_distortion has been called then the min/max
! values of the actual plot are different from graph%x%min, etc.
! In this case, use the actual min/max values

call qp_get_axis_attrib ('X', x_min, x_max)
call qp_get_axis_attrib ('Y', y_min, y_max)

if ((end1%r(1) < x_min .or. x_max < end1%r(1) .or. end1%r(2) < y_min .or. y_max < end1%r(2)) .and. &
    (end2%r(1) < x_min .or. x_max < end2%r(1) .or. end2%r(2) < y_min .or. y_max < end2%r(2))) return

! Bends can be tricky if they are not in the X-Z plane. 
! Bends are parameterized by a set of points (x_bend, y_bend) along their  
! centerline and a set of vectors (dx_bend, dy_bend) tangent to the centerline.

if (is_bend) then

  floor = ele1%floor
  v_old = floor%r
  call floor_angles_to_w_mat (floor%theta, floor%phi, 0.0_rp, w_old)

  n_bend = min(abs(int(100 * ele%value(angle$))) + 1, ubound(x_bend, 1))
  ang    = ele%value(angle$) * ele%orientation
  length = ele%value(l$)     * ele%orientation
  do j = 0, n_bend
    angle = j * ang / n_bend
    cos_t = cos(ele%value(ref_tilt_tot$))
    sin_t = sin(ele%value(ref_tilt_tot$))
    cos_a = cos(angle)
    sin_a = sin(angle)
    if (ele%value(g$) == 0) then
      r_vec = length * j * [0, 0, 1]
    else
      r_vec = ele%value(rho$) * [cos_t * (cos_a - 1), sin_t * (cos_a - 1), sin_a]
    endif
    dr_vec = [-cos_t * sin_a, -sin_t * sin_a, cos_a]
    ! This keeps dr_vec pointing to the inside (important for the labels).
    if (cos_t < 0) dr_vec = -dr_vec
    v_vec = matmul (w_old, r_vec) + v_old
    dv_vec = matmul (w_old, dr_vec) 
    call floor_to_screen (graph, v_vec, x_bend(j), y_bend(j))
    call floor_to_screen (graph, dv_vec, dx_bend(j), dy_bend(j))

    ! Correct for e1 and e2 face angles which are a rotation of the faces about
    ! the local y-axis.

    if (j == 0) then
      e_edge = ele%value(e1$)
      if (ele%orientation == -1) e_edge = -ele%value(e2$)
      dr_vec = tan(e_edge) * [cos_t * cos_a, sin_t * cos_a, sin_a]
      dv_vec = matmul (w_old, dr_vec) 
      call floor_to_screen (graph, dv_vec, dx1, dy1)
      dx_bend(j) = dx_bend(j) - dx1
      dy_bend(j) = dy_bend(j) - dy1
      e1_factor = sqrt(dx_bend(j)**2 + dy_bend(j)**2)
    endif

    if (j == n_bend) then
      e_edge = ele%value(e2$)
      if (ele%orientation == -1) e_edge = -ele%value(e1$)
      dr_vec = tan(e_edge) * [cos_t * cos_a, sin_t * cos_a, sin_a]
      dv_vec = matmul (w_old, dr_vec) 
      call floor_to_screen (graph, dv_vec, dx1, dy1)
      dx_bend(j) = dx_bend(j) + dx1
      dy_bend(j) = dy_bend(j) + dy1
      e2_factor = sqrt(dx_bend(j)**2 + dy_bend(j)**2)
    endif
  enddo

endif

! Only those elements with an associated ele_shape are to be drawn in full.
! All others are drawn with a simple line or arc

is_there = .false.
if (associated(ele_shape)) is_there = ele_shape%draw
if (.not. is_there) then
  if (is_bend) then
    call qp_draw_polyline(x_bend(:n_bend), y_bend(:n_bend))
  else
    call qp_draw_line(end1%r(1), end2%r(1), end1%r(2), end2%r(2))
  endif
  return
endif

! Here if element is to be drawn...

icol = qp_translate_to_color_index (ele_shape%color)

off = ele_shape%size * s%plot_page%floor_plan_shape_scale 
off1 = off
off2 = off

select case (ele_shape%shape)
case ('VAR_BOX', 'ASYM_VAR_BOX')
  select case (ele%key)
  case (sbend$)
    off1 = off * ele%value(g$)
  case (quadrupole$)
    off1 = off * ele%value(k1$)
  case (sextupole$)
    off1 = off * ele%value(k2$)
  case (octupole$)
    off1 = off * ele%value(k3$)
  case (solenoid$)
    off1 = off * ele%value(ks$)
  end select
  off2 = off1
  if (ele_shape%shape == 'ASYM_VAR_BOX') off1 = 0
case ('VVAR_BOX')
  off1 = dat_var_value
  off2 = dat_var_value
case ('ASYM_VVAR_BOX')
  off1 = 0
  off2 = dat_var_value
endselect

if (s%plot_page%floor_plan_size_is_absolute) then
  draw_units = 'DATA'
else
  draw_units = 'POINTS'
endif

! x-ray line parameters if present

if (attribute_index(ele, 'X_RAY_LINE_LEN') > 0 .and. ele%value(x_ray_line_len$) > 0) then
  call init_ele(drift)
  drift%key = drift$
  drift%value(l$) = ele%value(x_ray_line_len$)
  call ele_geometry (ele2%floor, drift, drift%floor) 
  call floor_to_screen_coords (graph, drift%floor, x_ray)
  call qp_convert_point_abs (x_ray%r(1), x_ray%r(2), 'DATA', x_ray%r(1), x_ray%r(2), draw_units)
endif

! Draw the shape. Since the conversion from floor coords to screen coords can
! be different along x and y, we convert to screen coords to make sure that rectangles
! remain rectangular.

call qp_convert_point_abs (end1%r(1), end1%r(2), 'DATA', end1%r(1), end1%r(2), draw_units)
call qp_convert_point_abs (end2%r(1), end2%r(2), 'DATA', end2%r(1), end2%r(2), draw_units)

! dx1, etc. are offsets perpendicular to the refernece orbit

call qp_convert_point_rel (cos(end1%theta), sin(end1%theta), 'DATA', dt_x, dt_y, draw_units)
dx1 =  off1 * dt_y / sqrt(dt_x**2 + dt_y**2)
dy1 = -off1 * dt_x / sqrt(dt_x**2 + dt_y**2)

call qp_convert_point_rel (cos(end2%theta), sin(end2%theta), 'DATA', dt_x, dt_y, draw_units)
dx2 =  off2 * dt_y / sqrt(dt_x**2 + dt_y**2)
dy2 = -off2 * dt_x / sqrt(dt_x**2 + dt_y**2)

if (is_bend) then
  do j = 0, n_bend
    call qp_convert_point_abs (x_bend(j), y_bend(j), 'DATA', x_bend(j), y_bend(j), draw_units)
    call qp_convert_point_rel (dx_bend(j), dy_bend(j), 'DATA', dt_x, dt_y, draw_units)
    dx_bend(j) =  off * dt_y / sqrt(dt_x**2 + dt_y**2)
    dy_bend(j) = -off * dt_x / sqrt(dt_x**2 + dt_y**2)
  enddo
  dx_bend(0) = dx_bend(0) * e1_factor
  dy_bend(0) = dy_bend(0) * e1_factor
  dx_bend(n_bend) = dx_bend(n_bend) * e2_factor
  dy_bend(n_bend) = dy_bend(n_bend) * e2_factor

  ! Finite e1 or e2 may mean some points extend beyound the outline of the bend.
  ! Throw out these points.

  ! First look at the first half of the bend for points that are beyound due to e1.

  r0_plus  = [x_bend(0) + dx_bend(0), y_bend(0) + dy_bend(0)]
  r0_minus = [x_bend(0) - dx_bend(0), y_bend(0) - dy_bend(0)]
  n_mid = n_bend/2
  dr2_p = [x_bend(n_mid) + dx_bend(n_mid), y_bend(n_mid) + dy_bend(n_mid)] - r0_plus
  dr2_m = [x_bend(n_mid) - dx_bend(n_mid), y_bend(n_mid) - dy_bend(n_mid)] - r0_minus

  j = n_bend/2 
  do 
    if (j == 0) exit
    dr_p  = [x_bend(j) + dx_bend(j), y_bend(j) + dy_bend(j)] - r0_plus
    dr_m  = [x_bend(j) - dx_bend(j), y_bend(j) - dy_bend(j)] - r0_minus
    ! If one of the points is outside then exit
    if (dot_product(dr_p, dr2_p) < 0 .or.dot_product(dr_m, dr2_m) < 0) exit
    j = j - 1
  enddo

  ! If there are points outside, delete them

  if (j > 0) then
    x_bend(1:n_bend-j)  = x_bend(j+1:n_bend)
    y_bend(1:n_bend-j)  = y_bend(j+1:n_bend)
    dx_bend(1:n_bend-j) = dx_bend(j+1:n_bend)
    dy_bend(1:n_bend-j) = dy_bend(j+1:n_bend)
    n_bend = n_bend - j
  endif

  ! Now look at the last half of the bend for points that are beyound due to e2.

  r0_plus  = [x_bend(n_bend) + dx_bend(n_bend), y_bend(n_bend) + dy_bend(n_bend)]
  r0_minus = [x_bend(n_bend) - dx_bend(n_bend), y_bend(n_bend) - dy_bend(n_bend)]
  n_mid = n_bend/2
  dr2_p = [x_bend(n_mid) + dx_bend(n_mid), y_bend(n_mid) + dy_bend(n_mid)] - r0_plus
  dr2_m = [x_bend(n_mid) - dx_bend(n_mid), y_bend(n_mid) - dy_bend(n_mid)] - r0_minus

  j = n_bend/2 
  do 
    if (j == n_bend) exit
    dr_p  = [x_bend(j) + dx_bend(j), y_bend(j) + dy_bend(j)] - r0_plus
    dr_m  = [x_bend(j) - dx_bend(j), y_bend(j) - dy_bend(j)] - r0_minus
    ! If one of the points is outside then exit
    if (dot_product(dr_p, dr2_p) < 0 .or.dot_product(dr_m, dr2_m) < 0) exit
    j = j + 1
  enddo

  ! If there are points outside, delete them

  if (j < n_bend) then
    x_bend(j)  = x_bend(n_bend)
    y_bend(j)  = y_bend(n_bend)
    dx_bend(j) = dx_bend(n_bend)
    dy_bend(j) = dy_bend(n_bend)
    n_bend = j
  endif

endif

! Draw the element...

! Draw x-ray line

if (attribute_index(ele, 'X_RAY_LINE_LEN') > 0 .and. ele%value(x_ray_line_len$) > 0) then
  drift%key = photon_fork$
  drift%name = ele%name
  branch_shape => tao_pointer_to_ele_shape (ix_uni, drift, s%plot_page%floor_plan%ele_shape)
  if (associated(branch_shape)) then
    if (branch_shape%draw) then
      call qp_draw_line (x_ray%r(1), end2%r(1), x_ray%r(2), end2%r(2), units = draw_units, &
                                        color = qp_translate_to_color_index (branch_shape%color))
    endif
  endif
endif

shape_has_box = (index(ele_shape%shape, 'BOX') /= 0)

! Draw diamond

if (ele_shape%shape == 'DIAMOND') then
  if (is_bend) then
    n = n_bend / 2
    x1 = (x_bend(n) + dx_bend(n)) / 2
    x2 = (x_bend(n) - dx_bend(n)) / 2
    y1 = (y_bend(n) + dy_bend(n)) / 2
    y2 = (y_bend(n) - dy_bend(n)) / 2
  else
    x1 = ((end1%r(1) + end2%r(1)) + (dx1 + dx2)) / 2
    x2 = ((end1%r(1) + end2%r(1)) - (dx1 + dx2)) / 2
    y1 = ((end1%r(2) + end2%r(2)) + (dy1 + dy2)) / 2
    y2 = ((end1%r(2) + end2%r(2)) - (dy1 + dy2)) / 2
  endif
  call qp_draw_line (end1%r(1), x1, end1%r(2), y1, units = draw_units, color = icol)
  call qp_draw_line (end1%r(1), x2, end1%r(2), y2, units = draw_units, color = icol)
  call qp_draw_line (end2%r(1), x1, end2%r(2), y1, units = draw_units, color = icol)
  call qp_draw_line (end2%r(1), x2, end2%r(2), y2, units = draw_units, color = icol)
endif

! Draw a circle.

if (ele_shape%shape == 'CIRCLE') then
  call qp_draw_circle ((end1%r(1)+end2%r(1))/2, (end1%r(2)+end2%r(2))/2, off, &
                                                  units = draw_units, color = icol)
endif

! Draw an X.

if (ele_shape%shape == 'X') then
  if (is_bend) then
    n = n_bend / 2
    x0  = x_bend(n)
    y0  = y_bend(n)
    dx1 = dx_bend(n)
    dy1 = dy_bend(n)
  else
    x0 = (end1%r(1) + end2%r(1)) / 2
    y0 = (end1%r(2) + end2%r(2)) / 2
  endif
  call qp_draw_line (x0 - dx1, x0 + dx1, y0 - dy1, y0 + dy1, units = draw_units, color = icol) 
  call qp_draw_line (x0 - dx1, x0 + dx1, y0 + dy1, y0 - dy1, units = draw_units, color = icol) 
endif

! Draw top and bottom of boxes and bow_tiw

if (ele_shape%shape == 'BOW_TIE' .or. shape_has_box) then
  if (is_bend) then
    call qp_draw_polyline(x_bend(:n_bend) + dx_bend(:n_bend), &
                          y_bend(:n_bend) + dy_bend(:n_bend), units = draw_units, color = icol)
    call qp_draw_polyline(x_bend(:n_bend) - dx_bend(:n_bend), &
                          y_bend(:n_bend) - dy_bend(:n_bend), units = draw_units, color = icol)

  else
    call qp_draw_line (end1%r(1)+dx1, end2%r(1)+dx1, end1%r(2)+dy1, end2%r(2)+dy1, &
                                                    units = draw_units, color = icol)
    call qp_draw_line (end1%r(1)-dx2, end2%r(1)-dx2, end1%r(2)-dy2, end2%r(2)-dy2, &
                                                    units = draw_units, color = icol)
  endif
endif

! Draw sides of boxes

if (shape_has_box) then
  if (is_bend) then
    call qp_draw_line (x_bend(0)-dx_bend(0), x_bend(0)+dx_bend(0), &
                       y_bend(0)-dy_bend(0), y_bend(0)+dy_bend(0), units = draw_units, color = icol)
    n = n_bend
    call qp_draw_line (x_bend(n)-dx_bend(n), x_bend(n)+dx_bend(n), &
                       y_bend(n)-dy_bend(n), y_bend(n)+dy_bend(n), units = draw_units, color = icol)
  else
    call qp_draw_line (end1%r(1)+dx1, end1%r(1)-dx2, end1%r(2)+dy1, end1%r(2)-dy2, &
                                                  units = draw_units, color = icol)
    call qp_draw_line (end2%r(1)+dx1, end2%r(1)-dx2, end2%r(2)+dy1, end2%r(2)-dy2, &
                                                  units = draw_units, color = icol)
  endif
endif

! Draw X for xbox or bow_tie

if (ele_shape%shape == 'XBOX' .or. ele_shape%shape == 'BOW_TIE') then
  call qp_draw_line (end1%r(1)+dx1, end2%r(1)-dx2, end1%r(2)+dy1, end2%r(2)-dy2, &
                                                  units = draw_units, color = icol)
  call qp_draw_line (end1%r(1)-dx2, end2%r(1)+dx1, end1%r(2)-dy1, end2%r(2)+dy2, &
                                                  units = draw_units, color = icol)
endif

! Draw the label.
! Since multipass slaves are on top of one another, just draw the multipass lord's name.
! Also place a bend's label to the outside of the bend.

if (ele_shape%label == 'name') then
  if (dat_var_name /= '') then
    name = dat_var_name
  elseif (ele%slave_status == multipass_slave$) then
    lord => pointer_to_lord(ele, 1)
    name = lord%name
  else
    name = ele%name
  endif
elseif (ele_shape%label == 's') then
  write (name, '(f16.2)') ele%s - ele%value(l$) / 2
  call string_trim (name, name, ix)
elseif (ele_shape%label /= 'none') then
  call out_io (s_error$, r_name, 'BAD ELEMENT LABEL: ' // ele_shape%label)
  call err_exit
endif 

if (ele_shape%label /= 'none') then
  if (ele%key /= sbend$ .or. ele%value(g$) == 0) then
    x_center = (end1%r(1) + end2%r(1)) / 2 
    y_center = (end1%r(2) + end2%r(2)) / 2 
    dx = -2 * dt_y / sqrt(dt_x**2 + dt_y**2)
    dy =  2 * dt_x / sqrt(dt_x**2 + dt_y**2)
  else
    n = n_bend / 2
    x_center = x_bend(n) 
    y_center = y_bend(n) 
    dx = -2 * dx_bend(n) / sqrt(dx_bend(n)**2 + dy_bend(n)**2)
    dy = -2 * dy_bend(n) / sqrt(dx_bend(n)**2 + dy_bend(n)**2)
  endif
  ! The extra factors of 2 are to shift the branch cut away from +/- 90.
  ! This is done since many lattices have elements with theta at +/- 90.
  theta = modulo2 (2 + atan2(dy, dx) * 180 / pi, 90.0_rp) - 2
  if (dx > 0) then
    justify = 'LC'
  else
    justify = 'RC'
  endif
  height = s%plot_page%text_height * s%plot_page%legend_text_scale
  call qp_draw_text (name, x_center+dx*off2, y_center+dy*off2, units = draw_units, &
                               height = height, justify = justify, ANGLE = theta)    
endif

end subroutine tao_draw_ele_for_floor_plan

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_lat_layout (plot, graph)
!
! Routine to draw a lattice layout graph.
!
! Input:
!   plot  -- Tao_plot_struct: Plot containing the graph.
!   graph -- Tao_graph_struct: Graph to plot.
!-

subroutine tao_draw_lat_layout (plot, graph)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (tao_lattice_branch_struct), pointer :: lat_branch
type (tao_ele_shape_struct), pointer :: ele_shape
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele, ele1, ele2
type (branch_struct), pointer :: branch, branch2
type (tao_data_array_struct), allocatable, target :: d_array(:)
type (tao_var_array_struct), allocatable, target :: v_array(:)
type (tao_logical_array_struct), allocatable :: logic_array(:)
type (tao_data_struct), pointer :: datum
type (tao_var_struct), pointer :: var

real(rp) x1, x2, y1, y2, y, s_pos, x0, y0
real(rp) lat_len, height, dx, dy, key_number_height, dummy, l2

integer i, j, k, n, kk, ix, ix1, isu
integer ix_var, ixv

logical shape_has_box, err, have_data

character(80) str
character(40) name
character(20) :: r_name = 'tao_draw_lat_layout'
character(20) shape_name

! Init

if (.not. graph%valid) return

select case (plot%x_axis_type)
case ('index')
  call out_io (s_error$, r_name, '"index" x-axis type not valid with lat_layout')
  graph%valid = .false.
  return

case ('ele_index')
  call out_io (s_error$, r_name, '"ele_index" x-axis type not valid with lat_layout')
  graph%valid = .false.
  return

case ('s')

case default
  call out_io (s_warn$, r_name, "Unknown x_axis_type")
  graph%valid = .false.
  return
end select

isu = tao_universe_number(graph%ix_universe)
lat => s%u(isu)%model%lat
branch => lat%branch(graph%ix_branch)
lat_branch => s%u(isu)%model%lat_branch(graph%ix_branch)

lat_len = branch%param%total_length
  
! Setup the placement of the graph on the plot page.

call qp_set_layout (x_axis = graph%x, y_axis = graph%y, x2_mirrors_x = .false., &
                    box = graph%box, margin = graph%margin)

call qp_draw_line (graph%x%min, graph%x%max, 0.0_rp, 0.0_rp)

! loop over all elements in the branch. Only draw those element that
! are within bounds.

do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  if (ele%slave_status == super_slave$) cycle
  call draw_ele_for_lat_layout (ele, graph)
enddo

! Loop over all control elements.

do i = lat%n_ele_track+1, lat%n_ele_max
  ele => branch%ele(i)
  if (ele%lord_status == multipass_lord$) cycle
  branch2 => pointer_to_branch(ele)
  if (branch2%ix_branch /= branch%ix_branch) cycle
  call draw_ele_for_lat_layout (ele, graph)
enddo

! Draw x-axis min max

if (graph%x%draw_numbers) then
  call qp_to_axis_number_text (graph%x, 0, str)
  call qp_convert_point_abs (graph%x%min, -10.0_rp, 'DATA', x1, y1, 'POINTS')
  call qp_draw_text (trim(str) // '-|', x1, y1, 'POINTS', justify = 'RT')
  call qp_to_axis_number_text (graph%x, graph%x%major_div, str)
  call qp_convert_point_abs (graph%x%max, -10.0_rp, 'DATA', x1, y1, 'POINTS')
  call qp_draw_text ('|-' // trim(str), x1, y1, 'POINTS', justify = 'LT')
endif

! This is for drawing the key numbers under the appropriate elements

key_number_height = 10

if (s%global%label_keys) then
  do kk = 1, 10
    k = kk + 10*s%com%ix_key_bank
    if (k > ubound(s%key, 1)) cycle
    ix_var = s%key(k)
    if (ix_var < 1) cycle
    write (str, '(i1)') mod(kk, 10)
    var => s%var(ix_var)
    do ixv = 1, size(var%this)
      if (var%this(ixv)%ix_uni /= isu) cycle
      ele => pointer_to_ele(lat, var%this(ixv)%ix_ele, var%this(ixv)%ix_branch)
      if (ele%n_slave /= 0 .and. ele%lord_status /= super_lord$) then
        do j = 1, ele%n_slave
          ele1 => pointer_to_slave(ele, j)
          l2 = ele1%value(l$) / 2
          s_pos = ele1%s - l2
          if (s_pos > graph%x%max .and. s_pos-lat_len > graph%x%min) s_pos = s_pos - lat_len
          if (s_pos + l2 < graph%x%min .or. s_pos - l2 > graph%x%max) cycle
          call qp_draw_text (trim(str), s_pos, graph%y%max, justify = 'CT', height = key_number_height)  
        enddo
      else
        l2 = ele%value(l$) / 2
        s_pos = ele%s - l2
        if (s_pos > graph%x%max .and. s_pos-lat_len > graph%x%min) s_pos = s_pos - lat_len
        if (s_pos + l2 < graph%x%min .or. s_pos - l2 > graph%x%max) cycle
        call qp_draw_text (trim(str), s_pos, graph%y%max, justify = 'CT', height = key_number_height)  
      endif
    enddo
  enddo
endif

! Draw data and beam_chamber curves

if (allocated(graph%curve)) then
  do i = 1, size(graph%curve)
    call tao_draw_curve_data (plot, graph, graph%curve(i), have_data)
  enddo
endif

!--------------------------------------------------------------------------------------------------
contains 

subroutine draw_ele_for_lat_layout (ele, graph)

type (ele_struct) ele
type (tao_graph_struct) graph
type (ele_struct), pointer :: ele1, ele2
type (tao_ele_shape_struct), pointer :: ele_shape

real(rp) dat_var_value
integer section_id, icol, ix_shape_min, ix_shape

character(40) dat_var_name, this_name

! Draw element shape...

ix_shape_min = 1
do
  ele_shape => tao_pointer_to_ele_shape (isu, ele, s%plot_page%lat_layout%ele_shape(ix_shape_min:), &
                                                                      dat_var_name, dat_var_value, ix_shape)
  if (.not. associated(ele_shape)) return
  if (.not. ele_shape%draw) return

  shape_name = ele_shape%shape

  call find_element_ends (ele, ele1, ele2)
  if (.not. associated(ele1)) return
  if (ele1%ix_branch /= graph%ix_branch) return
  x1 = ele1%s
  x2 = ele2%s
  ! If out of range then try a negative position
  if (branch%param%geometry == closed$ .and. x1 > graph%x%max) then
    x1 = x1 - lat_len
    x2 = x2 - lat_len
  endif
    
  if (x1 > graph%x%max) return
  if (x2 < graph%x%min) return

  ! Here if element is to be drawn...
  ! r1 and r2 are the scale factors for the lines below and above the center line.

  y = ele_shape%size * s%plot_page%lat_layout_shape_scale 
  y1 = -y
  y2 =  y

  select case (ele_shape%shape)
  case ('VAR_BOX', 'ASYM_VAR_BOX')
    select case (ele%key)
    case (sbend$)
      y2 = y * ele%value(g$)
    case (quadrupole$)
      y2 = y * ele%value(k1$)
    case (sextupole$)
      y2 = y * ele%value(k2$)
    case (octupole$)
      y2 = y * ele%value(k3$)
    case (solenoid$)
      y2 = y * ele%value(ks$)
    end select
    y1 = -y2
    if (shape_name == 'ASYM_VAR_BOX') y1 = 0
  case ('VVAR_BOX')
    y1 = -dat_var_value
    y2 = dat_var_value
  case ('ASYM_VVAR_BOX')
    y1 = 0
    y2 = dat_var_value
  endselect

  y1 = max(graph%y%min, min(y1, graph%y%max))
  y2 = max(graph%y%min, min(y2, graph%y%max))

  this_name = dat_var_name
  if (this_name == '') this_name = ele%name
  call draw_shape_for_lat_layout (this_name, ele%s - ele%value(l$) / 2, ele_shape)

  if (.not. ele_shape%multi) return
  ix_shape_min = ix_shape_min + ix_shape
enddo

end subroutine draw_ele_for_lat_layout

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! contains

subroutine draw_shape_for_lat_layout (name_in, s_pos, ele_shape)

type (tao_ele_shape_struct) ele_shape
real(rp) s_pos, r_dum, y_off
integer icol
character(*) name_in
character(20) shape_name

!

shape_name = ele_shape%shape
shape_has_box = (index(shape_name, 'BOX') /= 0)
icol = qp_translate_to_color_index (ele_shape%color)

! Draw the shape

if (shape_name == 'DIAMOND') then
  call qp_draw_line (x1, (x1+x2)/2, 0.0_rp, y1, color = icol)
  call qp_draw_line (x1, (x1+x2)/2, 0.0_rp, y2, color = icol)
  call qp_draw_line (x2, (x1+x2)/2, 0.0_rp, y1, color = icol)
  call qp_draw_line (x2, (x1+x2)/2, 0.0_rp, y2, color = icol)
endif

if (shape_name == 'CIRCLE') then
  call qp_convert_point_abs ((x1+x2)/2, (y1+y2)/2, 'DATA', x0, y0, 'POINTS')
  call qp_convert_point_rel (r_dum, y1, 'DATA', r_dum, y1, 'POINTS')
  call qp_draw_circle (x0, y0, abs(y1), units = 'POINTS', color = icol)
endif

if (shape_name == 'X') then
  call qp_convert_point_abs ((x1+x2)/2, (y1+y2)/2, 'DATA', x0, y0, 'POINTS')
  call qp_convert_point_rel (x1, y1, 'DATA', x1, y1, 'POINTS')
  call qp_convert_point_rel (x2, y2, 'DATA', x2, y2, 'POINTS')
  call qp_draw_line (x0-y1, x0+y1, y0-y1, y0+y1, units = 'POINTS', color = icol)
  call qp_draw_line (x0-y1, x0+y1, y0+y1, y0-y1, units = 'POINTS', color = icol)
endif

if (shape_name == 'BOW_TIE') then
  call qp_draw_line (x1, x2, y1, y1, color = icol)
  call qp_draw_line (x1, x2, y2, y2, color = icol)
endif

if (shape_has_box) then
  call qp_draw_rectangle (x1, x2, y1, y2, color = icol)
endif

! Draw X for XBOX or BOW_TIE

if (shape_name == 'XBOX' .or. shape_name == 'BOW_TIE') then
  call qp_draw_line (x1, x2, y2, y1, color = icol)
  call qp_draw_line (x1, x2, y1, y2, color = icol)
endif

! Put on a label

if (s%global%label_lattice_elements .and. ele_shape%label /= 'none') then

  call qp_from_inch_rel (0.0_rp, graph%y%label_offset, r_dum, y_off, 'DATA')

  if (ele_shape%label == 'name') then
    name = name_in
  elseif (ele_shape%label == 's') then
    write (name, '(f16.2)') s_pos
    call string_trim (name, name, ix)
  else
    call out_io (s_error$, r_name, 'BAD ELEMENT LABEL: ' // ele_shape%label)
    call err_exit
  endif 

  if (s_pos > graph%x%max .and. s_pos-lat_len > graph%x%min) s_pos = s_pos - lat_len
  height = s%plot_page%text_height * s%plot_page%legend_text_scale
  call qp_draw_text (name, s_pos, graph%y%min-y_off, height = height, justify = 'LC', ANGLE = 90.0_rp)

endif

end subroutine draw_shape_for_lat_layout

end subroutine tao_draw_lat_layout

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_beam_chamber_wall (plot, graph)
!
! Routine to draw the beam chamber wall.
!
! Input:
!   plot  -- Tao_plot_struct: Plot containing the graph.
!   graph -- Tao_graph_struct: Graph to plot.
!-

subroutine tao_draw_beam_chamber_wall (plot, graph)

implicit none

type (tao_plot_struct) plot
type (ele_struct), pointer :: ele
type (tao_graph_struct), target :: graph
type (lat_struct), pointer :: lat
type (tao_lattice_branch_struct), pointer :: lat_branch
type (branch_struct), pointer :: branch, branch2

real(rp) lat_len, dummy

integer i, isu

!

isu = tao_universe_number(graph%ix_universe)
lat => s%u(isu)%model%lat
branch => lat%branch(graph%ix_branch)
lat_branch => s%u(isu)%model%lat_branch(graph%ix_branch)

lat_len = branch%param%total_length
  
! Setup the placement of the graph on the plot page.

call qp_set_layout (x_axis = graph%x, y_axis = graph%y, x2_mirrors_x = .false., &
                    box = graph%box, margin = graph%margin)
call qp_draw_line (graph%x%min, graph%x%max, 0.0_rp, 0.0_rp)

! Loop over all elements

do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  if (ele%slave_status == super_slave$) cycle
  call draw_ele_beam_chamber (ele)
enddo

! Loop over all control elements.

do i = lat%n_ele_track+1, lat%n_ele_max
  ele => branch%ele(i)
  if (ele%lord_status == multipass_lord$) cycle
  branch2 => pointer_to_branch(ele)
  if (branch2%ix_branch /= branch%ix_branch) cycle
  call draw_ele_beam_chamber (ele)
enddo

!------------------------------------------------------------------------
contains

subroutine draw_ele_beam_chamber (ele)

type (ele_struct) ele

real(rp) y1_plus, y1_minus, y2_plus, y2_minus, x1, x2, y1, y2
integer section_id, icol

! Draw beam chamber wall. 

icol = black$
if (allocated (graph%curve)) icol = graph%curve(1)%line%color

if (associated(ele%wall3d)) then
  call calc_wall_radius (ele%wall3d(1)%section(1)%v,  1.0_rp, 0.0_rp,  y1_plus, dummy)
  call calc_wall_radius (ele%wall3d(1)%section(1)%v, -1.0_rp, 0.0_rp,  y1_minus, dummy)
  x1 = ele%s - ele%value(l$) + ele%wall3d(1)%section(1)%s

  ! Skip points so close to the last point that the points have negligible spacing

  do section_id = 2, size(ele%wall3d(1)%section)
    x2 = ele%s - ele%value(l$) + ele%wall3d(1)%section(section_id)%s
    if (section_id /= size(ele%wall3d(1)%section) .and. &
            (x2 - x1) < (graph%x%max - graph%x%min) / s%plot_page%n_curve_pts) cycle
    call calc_wall_radius (ele%wall3d(1)%section(section_id)%v,  1.0_rp, 0.0_rp,  y2_plus, dummy)
    call calc_wall_radius (ele%wall3d(1)%section(section_id)%v, -1.0_rp, 0.0_rp,  y2_minus, dummy)
    !scale wall
    call qp_draw_line (x1, x2, y1_plus, y2_plus, color = icol)
    call qp_draw_line (x1, x2, -y1_minus, -y2_minus, color = icol)
    x1       = x2
    y1_plus  = y2_plus
    y1_minus = y2_minus 
  end do
endif

end subroutine draw_ele_beam_chamber

end subroutine tao_draw_beam_chamber_wall

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_plot_data (plot, graph)
!
! Routine to draw a graph with data and/or variable curves. 
!
! Input:
!   plot  -- Tao_plot_struct: Plot containing the graph.
!   graph -- Tao_graph_struct: Graph to plot.
!-

subroutine tao_plot_data (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph

integer k
logical have_data

! Draw the graph outline.

call tao_draw_data_graph (plot, graph)
if (.not. graph%valid) return

! loop over all the curves of the graph and draw them

have_data = .false.

do k = 1, size(graph%curve)
  call tao_draw_curve_data (plot, graph, graph%curve(k), have_data)
enddo

if (.not. have_data) call qp_draw_text ('**No Plottable Data**', &
                            0.18_rp, -0.15_rp, '%/GRAPH/LT', color = red$) 

end subroutine tao_plot_data

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_curve_data (plot, graph, curve, have_data)
!
! Routine to draw a graph with data and/or variable curves. 
!
! Input:
!   plot      -- Tao_plot_struct: Plot containing the graph.
!   graph     -- Tao_graph_struct: Graph containing the curve.
!   curve     -- Tao_curve_struct: Curve to draw.
!   have_data -- Logical: Intitial state.
! Output:
!    have_data -- Logical: Is there any data to plot? Set True if so.
!                   But never reset to False. 
!-

subroutine tao_draw_curve_data (plot, graph, curve, have_data)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct) :: curve

integer i
logical have_data
character(16) num_str

!

if (curve%use_y2) call qp_use_axis (y = 'Y2')
call qp_set_symbol (curve%symbol)

if (curve%draw_symbols .and. allocated(curve%x_symb)) then
  if (size(curve%x_symb) > 0) have_data = .true.
  if (graph%symbol_size_scale > 0) then
    do i = 1, size(curve%x_symb), max(1, curve%symbol_every)
      call qp_draw_symbol (curve%x_symb(i), curve%y_symb(i), height = curve%symb_size(i), clip = graph%clip)
    enddo
  else
    call qp_draw_symbols (curve%x_symb, curve%y_symb, symbol_every = curve%symbol_every, clip = graph%clip)
  endif
endif

if (curve%draw_symbol_index .and. allocated(curve%ix_symb)) then
  if (size(curve%ix_symb) > 0) have_data = .true.
  do i = 1, size(curve%ix_symb)
    if (graph%clip) then
      if (curve%x_symb(i) < graph%x%min .or. curve%x_symb(i) > graph%x%max)  cycle
      if (curve%y_symb(i) < graph%y%min .or. curve%y_symb(i) > graph%y%max) cycle
    endif
    write (num_str, '(i0)') curve%ix_symb(i)
    call qp_draw_text (num_str, curve%x_symb(i), curve%y_symb(i))
  enddo
endif

if (curve%draw_line .and. allocated(curve%x_line)) then
  if (size(curve%x_line) > 0) have_data = .true.
  call qp_set_line ('PLOT', curve%line) 
  call qp_draw_polyline (curve%x_line, curve%y_line, clip = graph%clip, style = 'PLOT')
endif

call qp_use_axis (y = 'Y')  ! reset

end subroutine tao_draw_curve_data

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_histogram_data (plot, graph, curve, have_data)
!
! Routine to draw a graph with data and/or variable histograms. 
!
! Input:
!   plot      -- Tao_plot_struct: Plot containing the graph.
!   graph     -- Tao_graph_struct: Graph containing the histogram.
!   curve     -- Tao_curve_struct: Histogram to draw.
!   have_data -- Logical: Intitial state.
!
! Output:
!    have_data -- Logical: Is there any data to plot? Set True if so.
!                   But never reset to False. 
!-

subroutine tao_draw_histogram_data (plot, graph, curve, have_data)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct) :: curve

integer i
logical have_data
character(16) num_str


! Draw
call qp_draw_histogram (curve%x_line, curve%y_line, line_color = curve%line%color, &
	                    fill_color = curve%line%color, fill_pattern = curve%line%pattern) !, fill_color, fill_pattern, line_color, clip)
have_data = .true.

end subroutine tao_draw_histogram_data

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_data_graph (plot, graph)
!
! Routine to draw a just the graph part of a data graph.
! The calling routine takes care of drawing any curves.
!
! Input:
!   plot  -- Tao_plot_struct: Plot containing the graph.
!   graph -- Tao_graph_struct: Graph to plot.
!-

subroutine tao_draw_data_graph (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (qp_line_struct), allocatable :: line(:)
type (qp_symbol_struct), allocatable :: symbol(:)

integer i, j, k, n
real(rp) x, y, x1
character(100), allocatable :: text(:)

! Set scales, margens, etc

call qp_set_layout (box = graph%box, margin = graph%margin)
call qp_set_layout (x_axis = graph%x, x2_mirrors_x = .false.)
call qp_set_layout (y_axis = graph%y, y2_axis = graph%y2, y2_mirrors_y = graph%y2_mirrors_y)
if (graph%title == '') then
  call qp_set_graph (title = '')
else
  call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
endif
call qp_draw_axes (draw_grid = graph%draw_grid)

! Draw the default x-axis label if there is none. 

if (graph%x%draw_label .and. graph%x%label == '') then
  select case (plot%x_axis_type) 
  case ('index', 'ele_index', 's')
    call qp_to_inch_rel (1.0_rp, 0.0_rp, x1, y, '%GRAPH')
    x = x1 * (graph%x%major_div - 0.5) / graph%x%major_div
    y = -graph%x%number_offset
    call qp_draw_text (plot%x_axis_type, x, y, 'INCH', justify = 'CT')
  end select
endif

!

if (.not. graph%valid) return

if (graph%limited .and. graph%clip .and. s%global%draw_curve_off_scale_warn) &
  call qp_draw_text ('**Curve Off Scale**', -0.30_rp, -0.15_rp, '%/GRAPH/RT', color = red$) 


! Draw the text legend if there is one

if (any(graph%text_legend /= ' ')) call qp_draw_text_legend (graph%text_legend, &
       graph%text_legend_origin%x, graph%text_legend_origin%y, graph%text_legend_origin%units)

! Draw the curve legend if needed

n = size(graph%curve)
allocate (text(n), symbol(n), line(n))

do i = 1, n
  curve => graph%curve(i)
  text(i) = curve%legend_text
  if (text(i) == '') text(i) = curve%data_type
  symbol(i) = curve%symbol
  if (size(curve%x_symb) == 0) symbol(i)%type = -1 ! Do not draw
  if (.not. curve%draw_symbols) symbol(i)%type = -1
  line(i) = curve%line
  if (size(curve%x_line) == 0) line(i)%width = -1 ! Do not draw
  if (.not. curve%draw_line) line(i)%width = -1
enddo

if (graph%draw_curve_legend .and. n > 1) then
  call qp_draw_curve_legend (graph%curve_legend_origin%x, graph%curve_legend_origin%y, &
            graph%curve_legend_origin%units, line, s%plot_page%curve_legend_line_len, &
            symbol, text, s%plot_page%curve_legend_text_offset)
endif

! Draw any curve info messages

j = 0
do i = 1, n
  curve => graph%curve(i)
  if (curve%message_text == '') cycle
  j = j + 1
  text(j) = curve%message_text
enddo

if (j > 1) then
  call qp_draw_text_legend (text(1:j), 0.50_rp, 0.95_rp, '%GRAPH/LB')
endif

!

deallocate (text, symbol, line)

end subroutine tao_draw_data_graph

end module
