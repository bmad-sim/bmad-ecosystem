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
  write (view_str, '(i3)') s%global%u_view
  call qp_draw_text ('View Universe:' // view_str, -2.0_rp, -2.0_rp, 'POINTS/PAGE/RT', 'RT')
endif

! loop over all plots

do i = 1, size(s%plot_region)

  if (.not. s%plot_region(i)%visible) cycle
  plot => s%plot_region(i)%plot

  ! set the s%plot_page border for this particular region

  location = s%plot_region(i)%location
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
    case ('data', 'phase_space', 'data_slice')
      call tao_plot_data (plot, graph)
    case ('wave.0', 'wave.a', 'wave.b')
      call tao_plot_wave (plot, graph)
     case ('lat_layout')
      call tao_plot_lat_layout (plot, graph)
    case ('key_table')
      call tao_plot_key_table (plot, graph)
    case ('floor_plan')
      call tao_plot_floor_plan (plot, graph)
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

subroutine tao_plot_wave (plot, graph)

type (tao_plot_struct) :: plot
type (tao_graph_struct), target :: graph
type (qp_axis_struct), pointer :: y

real(rp) y0, y1

!

call tao_plot_data (plot, graph)

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


i_off = tao_com%ix_key_bank
call tao_key_info_to_str (1, i_off+1, i_off+10, str, header)
call qp_draw_text ('   Ix  ' // header, 0.0_rp, y_here, 'POINTS/GRAPH', &
                             height = height, uniform_spacing = .true.)
  

do i = 1, 10

  k = i + i_off
  if (k > ubound(s%key, 1)) return

  prefix = ''
  j = mod(i, 10)
  if (i == 1) then
    write (prefix, '(i2, a, i2)') tao_com%ix_key_bank, ':', j
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

subroutine tao_plot_floor_plan (plot, graph)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (lat_struct), pointer :: lat
type (floor_position_struct) end1, end2, floor
type (tao_wall_point_struct), pointer :: pt(:)
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (tao_ele_shape_struct), pointer :: ele_shape
type (tao_data_struct), pointer :: datum
type (tao_data_array_struct), allocatable, target :: d_array(:)
type (tao_logical_array_struct), allocatable :: logic_array(:)

real(rp) theta, v_vec(3)
real(rp) x_bend(0:1000), y_bend(0:1000)

integer i, j, k, n, n_bend, isu, ic, ix_shape

character(20) :: r_name = 'tao_plot_floor_plan'

logical err

! Each graph is a separate floor plan plot (presumably for different universes). 
! setup the placement of the graph on the plot page.

call qp_set_layout (x_axis = graph%x, y_axis = graph%y, y2_axis = graph%y2, &
                                        box = graph%box, margin = graph%margin)

! Adjust the margins if there is to be no distortion of the drawing

if (graph%correct_xy_distortion) call qp_eliminate_xy_distortion

!

if (graph%draw_axes) then
  if (graph%title == '') then
    call qp_set_graph (title = '')
  else
    call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
  endif
  call qp_draw_axes
endif

isu = tao_universe_number(graph%ix_universe)
lat => s%u(isu)%model%lat

if (.not. graph%valid) return

! loop over all elements in the lattice. 

do n = 0, ubound(lat%branch, 1)
  branch => lat%branch(n)
  branch%ele%logic = .false.  ! Used to mark as drawn.
  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    call tao_find_ele_shape (ele, tao_com%ele_shape_floor_plan, ix_shape)
    if (ele%ix_ele > lat%n_ele_track .and. ix_shape == 0) cycle   ! Nothing to draw
    if (ele%lord_status == multipass_lord$) then
      do j = ele%ix1_slave, ele%ix2_slave
        ic = lat%control(j)%ix_slave
        call tao_draw_ele_for_floor_plan (plot, graph, lat, branch%ele(ic), ix_shape, .false.)
      enddo
    else
      call tao_draw_ele_for_floor_plan (plot, graph, lat, ele, ix_shape, .false.)
    endif
  enddo
enddo

! Draw data

do i = 1, size(tao_com%ele_shape_floor_plan)
  ele_shape => tao_com%ele_shape_floor_plan(i)
  if (ele_shape%ele_name(1:5) /= 'dat::') cycle
  call tao_find_data (err, ele_shape%ele_name, d_array = d_array, log_array = logic_array)
  if (err) cycle
  do j = 1, size(d_array)
    datum => d_array(j)%d
    if (datum%ix_branch /= graph%ix_branch) cycle
    if (size(logic_array) /= 0) then
      if (.not. logic_array(j)%l) cycle
    endif
    ele => pointer_to_ele (lat, datum%ix_branch, datum%ix_ele)
    call tao_draw_ele_for_floor_plan (plot, graph, lat, ele, i, .true.)
  enddo
enddo

! Draw the tunnel wall

if (allocated(s%wall)) then
  do i = 1, size(s%wall)
    pt => s%wall(i)%point

    do j = 1, size(pt)

      select case (pt(j)%type)
      case (point$)
        if (j == 1) cycle
        call floor_to_screen (pt(j-1)%x, 0.0_rp, pt(j-1)%z, end1%x, end1%y)
        call floor_to_screen (pt(j)%x, 0.0_rp, pt(j)%z, end2%x, end2%y)
        call qp_draw_line(end1%x, end2%x, end1%y, end2%y)

      case (arc$)
        n_bend = abs(int(100 * (pt(j)%theta2 - pt(j)%theta1))) + 1
        do k = 0, n_bend
          theta = pt(j)%theta1 + k * (pt(j)%theta2 - pt(j)%theta1) / n_bend
          v_vec(1) = pt(j)%x + pt(j)%r * sin(theta)
          v_vec(2) = 0
          v_vec(3) = pt(j)%z + pt(j)%r * cos(theta)
          call floor_to_screen (v_vec(1), v_vec(2), v_vec(3), x_bend(j), y_bend(j))
        enddo
        call qp_draw_polyline(x_bend(:n_bend), y_bend(:n_bend))
      end select

    enddo

  end do
end if

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Input:
!   plot
!   graph
!   lat
!   ele
!   ix_shape
!-

recursive subroutine tao_draw_ele_for_floor_plan (plot, graph, lat, ele, ix_shape, is_data)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (lat_struct) :: lat
type (ele_struct) :: ele
type (ele_struct) :: drift
type (ele_struct), pointer :: ele1, ele2
type (floor_position_struct) end1, end2, floor, x_ray
type (tao_wall_point_struct), pointer :: pt(:)
type (tao_ele_shape_struct), pointer :: ele_shape

integer i, j, k, ix_shape, icol, isu, n_bend, n, ix, ixs, ic

real(rp) off, off1, off2, angle, rho, dx1, dy1, dx2, dy2
real(rp) dt_x, dt_y, x_center, y_center, dx, dy, theta
real(rp) x_bend(0:1000), y_bend(0:1000), dx_bend(0:1000), dy_bend(0:1000)
real(rp) v_old(3), w_old(3,3), r_vec(3), dr_vec(3), v_vec(3), dv_vec(3)
real(rp) cos_t, sin_t, cos_p, sin_p, cos_a, sin_a, height
real(rp) x_inch, y_inch, x0, y0, x1, x2, y1, y2

character(80) str
character(40) name
character(40) :: r_name = 'tao_draw_ele_for_floor_plan'
character(16) shape
character(2) justify

logical is_data
logical shape_has_box, is_bend

!

call find_element_ends (lat, ele, ele1, ele2)
if (.not. associated(ele1)) return

if (is_data) then  ! pretend this is zero length element
  ele1 => ele2
  is_bend = .false.
else
  is_bend = (ele%key == sbend$)
endif

call floor_to_screen_coords (ele1%floor, end1)
call floor_to_screen_coords (ele2%floor, end2)

! Only draw those element that have at least one point in bounds.
  
if ((end1%x < graph%x%min .or. graph%x%max < end1%x .or. &
    end1%y < graph%y%min .or. graph%y%max < end1%y) .and. &
    (end2%x < graph%x%min .or. graph%x%max < end2%x .or. &
    end2%y < graph%y%min .or. graph%y%max < end2%y)) return

! Bends can be tricky if they are not in the X-Z plane. 
! Bends are parameterized by a set of points (x_bend, y_bend) along their  
! centerline and a set of vectors (dx_bend, dy_bend) perpendicular to the centerline.

if (is_bend) then

  if (ele%value(g$) == 0) then
    n_bend = 1
    x_bend(0) = end1%x; y_bend(0) = end1%y
    x_bend(1) = end2%x; y_bend(1) = end2%y

  else
    floor = ele1%floor
    v_old = (/ floor%x, floor%y, floor%z /)
    cos_t = cos(floor%theta)
    sin_t = sin(floor%theta)
    cos_p = cos(floor%phi)
    sin_p = sin(floor%phi)
    w_old(1, 1:3) = (/  cos_t,  -sin_t * sin_p, sin_t * cos_p /)
    w_old(2, 1:3) = (/ 0.0_rp,   cos_p,         sin_p /)
    w_old(3, 1:3) = (/ -sin_t,  -cos_t * sin_p, cos_t * cos_p /)

    rho = ele%value(rho$)

    n_bend = min(abs(int(100 * ele%value(angle$))) + 1, ubound(x_bend, 1))
    do j = 0, n_bend
      angle = j * ele%value(angle$) / n_bend
      cos_t = cos(ele%value(tilt$))
      sin_t = sin(ele%value(tilt$))
      cos_a = cos(angle)
      sin_a = sin(angle)
      r_vec = rho * (/ cos_t * (cos_a - 1), sin_t * (cos_a - 1), sin_a /)
      dr_vec = rho * (/ -cos_t * sin_a, -sin_t * sin_a, cos_a /)
      ! This keeps dr_vec pointing to the inside (important for the labels).
      if (cos_t < 0) dr_vec = -dr_vec
      v_vec = matmul (w_old, r_vec) + v_old
      dv_vec = matmul (w_old, dr_vec) 
      call floor_to_screen (v_vec(1), v_vec(2), v_vec(3), x_bend(j), y_bend(j))
      call floor_to_screen (dv_vec(1), dv_vec(2), dv_vec(3), dx_bend(j), dy_bend(j))
    enddo
  endif

endif

! Only those elements with ix_shape > 0 are to be drawn in full.
! All others are drawn with a simple line or arc

if (ix_shape < 1) then
  if (is_bend) then
    call qp_draw_polyline(x_bend(:n_bend), y_bend(:n_bend))
  else
    call qp_draw_line(end1%x, end2%x, end1%y, end2%y)
  endif
  return
endif

! Here if element is to be drawn...

ele_shape => tao_com%ele_shape_floor_plan(ix_shape)
shape = ele_shape%shape

call qp_translate_to_color_index (ele_shape%color, icol)

off = ele_shape%dy_pix
off1 = off
off2 = off
if (shape == 'VAR_BOX' .or. shape == 'ASYM_VAR_BOX') then
  select case (ele%key)
  case (quadrupole$)
    off1 = off * ele%value(k1$)
  case (sextupole$)
    off1 = off * ele%value(k2$)
  case (octupole$)
    off1 = off * ele%value(k3$)
  case (solenoid$)
    off1 = off * ele%value(ks$)
  end select
  off1 = max(-s%plot_page%shape_height_max, min(off1, s%plot_page%shape_height_max))
  off2 = off1
  if (shape == 'ASYM_VAR_BOX') off1 = 0
endif

! x-ray line parameters if present

if (attribute_index(ele, 'X_RAY_LINE_LEN') > 0 .and. ele%value(x_ray_line_len$) > 0) then
  call init_ele(drift)
  drift%key = drift$
  drift%value(l$) = ele%value(x_ray_line_len$)
  call ele_geometry (ele2%floor, drift, drift%floor) 
  call floor_to_screen_coords (drift%floor, x_ray)
  call qp_convert_point_abs (x_ray%x, x_ray%y, 'DATA', x_ray%x, x_ray%y, 'POINTS')
endif

! Draw the shape. Since the conversion from floor coords to screen coords can
! be different along x and y, we convert to screen coords to make sure that rectangles
! remain rectangular.

call qp_convert_point_abs (end1%x, end1%y, 'DATA', end1%x, end1%y, 'POINTS')
call qp_convert_point_abs (end2%x, end2%y, 'DATA', end2%x, end2%y, 'POINTS')

! dx1, etc. are offsets perpendicular to the refernece orbit

call qp_convert_point_rel (cos(end1%theta), sin(end1%theta), 'DATA', dt_x, dt_y, 'POINTS')
dx1 =  off1 * dt_y / sqrt(dt_x**2 + dt_y**2)
dy1 = -off1 * dt_x / sqrt(dt_x**2 + dt_y**2)

call qp_convert_point_rel (cos(end2%theta), sin(end2%theta), 'DATA', dt_x, dt_y, 'POINTS')
dx2 =  off2 * dt_y / sqrt(dt_x**2 + dt_y**2)
dy2 = -off2 * dt_x / sqrt(dt_x**2 + dt_y**2)

if (is_bend) then
  do j = 0, n_bend
    call qp_convert_point_abs (x_bend(j), y_bend(j), 'DATA', x_bend(j), y_bend(j), 'POINTS')
    call qp_convert_point_rel (dx_bend(j), dy_bend(j), 'DATA', dt_x, dt_y, 'POINTS')
    dx_bend(j) =  off * dt_y / sqrt(dt_x**2 + dt_y**2)
    dy_bend(j) = -off * dt_x / sqrt(dt_x**2 + dt_y**2)
  enddo
endif

! Draw the element...

! Draw x-ray line

if (attribute_index(ele, 'X_RAY_LINE_LEN') > 0 .and. ele%value(x_ray_line_len$) > 0) then
  drift%key = photon_branch$
  drift%name = ele%name
  call tao_find_ele_shape (drift, tao_com%ele_shape_floor_plan, ixs)
  if (ixs > 0) then
    call qp_translate_to_color_index (tao_com%ele_shape_floor_plan(ixs)%color, ic)
    call qp_draw_line (x_ray%x, end2%x, x_ray%y, end2%y, units = 'POINTS', color = ic)
  endif
endif

shape_has_box = (index(shape, 'BOX') /= 0)

! Draw diamond

if (shape == 'DIAMOND') then
  if (is_bend) then
    n = n_bend / 2
    x1 = (x_bend(n) + dx_bend(n)) / 2
    x2 = (x_bend(n) - dx_bend(n)) / 2
    y1 = (y_bend(n) + dy_bend(n)) / 2
    y2 = (y_bend(n) - dy_bend(n)) / 2
  else
    x1 = ((end1%x + end2%x) + (dx1 + dx2)) / 2
    x2 = ((end1%x + end2%x) - (dx1 + dx2)) / 2
    y1 = ((end1%y + end2%y) + (dy1 + dy2)) / 2
    y2 = ((end1%y + end2%y) - (dy1 + dy2)) / 2
  endif
  call qp_draw_line (end1%x, x1, end1%y, y1, units = 'POINTS', color = icol)
  call qp_draw_line (end1%x, x2, end1%y, y2, units = 'POINTS', color = icol)
  call qp_draw_line (end2%x, x1, end2%y, y1, units = 'POINTS', color = icol)
  call qp_draw_line (end2%x, x2, end2%y, y2, units = 'POINTS', color = icol)
endif

! Draw a circle.

if (shape == 'CIRCLE') then
  call qp_draw_circle ((end1%x+end2%x)/2, (end1%y+end2%y)/2, off, &
                                                  units = 'POINTS', color = icol)
endif

! Draw an X.

if (shape == 'X') then
  if (is_bend) then
    n = n_bend / 2
    x0  = x_bend(n)
    y0  = y_bend(n)
    dx1 = dx_bend(n)
    dy1 = dy_bend(n)
  else
    x0 = (end1%x + end2%x) / 2
    y0 = (end1%y + end2%y) / 2
  endif
  call qp_draw_line (x0 - dx1, x0 + dx1, y0 - dy1, y0 + dy1, units = 'POINTS', color = icol) 
  call qp_draw_line (x0 - dx1, x0 + dx1, y0 + dy1, y0 - dy1, units = 'POINTS', color = icol) 
endif

! Draw top and bottom of boxes and bow_tiw

if (shape == 'BOW_TIE' .or. shape_has_box) then
  if (is_bend) then
    call qp_draw_polyline(x_bend(:n_bend) + dx_bend(:n_bend), &
                          y_bend(:n_bend) + dy_bend(:n_bend), units = 'POINTS', color = icol)
    call qp_draw_polyline(x_bend(:n_bend) - dx_bend(:n_bend), &
                          y_bend(:n_bend) - dy_bend(:n_bend), units = 'POINTS', color = icol)

  else
    call qp_draw_line (end1%x+dx1, end2%x+dx1, end1%y+dy1, end2%y+dy1, &
                                                    units = 'POINTS', color = icol)
    call qp_draw_line (end1%x-dx2, end2%x-dx2, end1%y-dy2, end2%y-dy2, &
                                                    units = 'POINTS', color = icol)
  endif
endif

! Draw sides of boxes

if (shape_has_box) then
  if (is_bend) then
    call qp_draw_line (x_bend(0)-dx_bend(0), x_bend(0)+dx_bend(0), &
                       y_bend(0)-dy_bend(0), y_bend(0)+dy_bend(0), units = 'POINTS', color = icol)
    n = n_bend
    call qp_draw_line (x_bend(n)-dx_bend(n), x_bend(n)+dx_bend(n), &
                       y_bend(n)-dy_bend(n), y_bend(n)+dy_bend(n), units = 'POINTS', color = icol)
  else
    call qp_draw_line (end1%x+dx1, end1%x-dx2, end1%y+dy1, end1%y-dy2, &
                                                  units = 'POINTS', color = icol)
    call qp_draw_line (end2%x+dx1, end2%x-dx2, end2%y+dy1, end2%y-dy2, &
                                                  units = 'POINTS', color = icol)
  endif
endif

! Draw X for xbox or bow_tie

if (shape == 'XBOX' .or. shape == 'BOW_TIE') then
  call qp_draw_line (end1%x+dx1, end2%x-dx2, end1%y+dy1, end2%y-dy2, &
                                                  units = 'POINTS', color = icol)
  call qp_draw_line (end1%x-dx2, end2%x+dx1, end1%y-dy1, end2%y+dy2, &
                                                  units = 'POINTS', color = icol)
endif

! Draw the label.
! Since multipass slaves are on top of one another, just draw the multipass lord's name.
! Also place a bend's label to the outside of the bend.

if (ele_shape%label_type == 'name') then
  if (ele%slave_status == multipass_slave$) then
    ix = ele%ic1_lord
    ix = lat%control(lat%ic(ix))%ix_lord
    name = lat%ele(ix)%name
  else
    name = ele%name
  endif
elseif (ele_shape%label_type == 's') then
  write (name, '(f16.2)') ele%s - ele%value(l$) / 2
  call string_trim (name, name, ix)
elseif (ele_shape%label_type /= 'none') then
  call out_io (s_error$, r_name, 'BAD ELEMENT LABEL: ' // ele_shape%label_type)
  call err_exit
endif 

if (ele_shape%label_type /= 'none') then
  if (ele%key /= sbend$ .or. ele%value(g$) == 0) then
    x_center = (end1%x + end2%x) / 2 
    y_center = (end1%y + end2%y) / 2 
    dx = -2 * dt_y / sqrt(dt_x**2 + dt_y**2)
    dy =  2 * dt_x / sqrt(dt_x**2 + dt_y**2)
  else
    n = n_bend / 2
    x_center = x_bend(n) 
    y_center = y_bend(n) 
    dx = -2 * dx_bend(n) / sqrt(dx_bend(n)**2 + dy_bend(n)**2)
    dy = -2 * dy_bend(n) / sqrt(dx_bend(n)**2 + dy_bend(n)**2)
  endif
  theta = modulo2 (atan2(dy, dx) * 180 / pi, 90.0_rp)
  if (dx > 0) then
    justify = 'LC'
  else
    justify = 'RC'
  endif
  height = s%plot_page%text_height * s%plot_page%legend_text_scale
  call qp_draw_text (name, x_center+dx*off2, y_center+dy*off2, units = 'POINTS', &
                               height = height, justify = justify, ANGLE = theta)    
endif

end subroutine tao_draw_ele_for_floor_plan

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_plot_lat_layout (plot, graph)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele, ele1, ele2
type (tao_ele_shape_struct), pointer :: ele_shapes(:), ele_shape
type (branch_struct), pointer :: branch
type (tao_data_array_struct), allocatable, target :: d_array(:)
type (tao_var_array_struct), allocatable, target :: v_array(:)
type (tao_logical_array_struct), allocatable :: logic_array(:)
type (tao_data_struct), pointer :: datum
type (tao_var_struct), pointer :: var

real(rp) x1, x2, y1, y2, y, s_pos, y_off, y_bottom, y_top, x0, y0
real(rp) lat_len, height, dx, dy, key_number_height, dummy, l2

integer i, j, ix_shape, k, kk, ix, ix1, isu
integer ix_var, ixv

logical shape_has_box, err

character(80) str
character(40) name
character(20) :: r_name = 'tao_plot_lat_layout'
character(20) shape_name

! Init

if (.not. graph%valid) return

isu = tao_universe_number(graph%ix_universe)
lat => s%u(isu)%model%lat
branch => lat%branch(graph%ix_branch)

lat_len = branch%param%total_length
key_number_height = 10
  
! Each graph is a separate lattice layout (presumably for different universes). 
! Setup the placement of the graph on the plot page.
! Without labels: Vertically Center the graph.
! With labels: Shift the graph up to give room for the labels.
! Center line is at y = 0

call qp_set_layout (x_axis = graph%x, box = graph%box, margin = graph%margin)

call qp_get_layout_attrib ('GRAPH', x1, x2, y1, y2, 'POINTS/GRAPH')
dy = (y2 - y1) 
y_bottom = -dy / 2
y_top = dy / 2

if (s%global%label_lattice_elements) then
  y_off = s%plot_page%shape_height_max
  if (s%global%label_keys) y_off = y_off + key_number_height
  if (2*y_off < dy) then
    y_top = y_off
    y_bottom = y_top - dy
  endif
endif

call qp_set_axis ('Y', y_bottom, y_top, 1, 0)

! Figure out x axis
! If it's not 's' then just draw a vertical line at the proper index

if (plot%x_axis_type == 's') then
  ! continue
elseif (plot%x_axis_type == 'index') then
  ! cannot plot layout in this case!
  graph%valid = .false.
  return
elseif (plot%x_axis_type == 'ele_index') then
  ! plot vertical line at ele_index
  ! temporarily turn this off until I get scaling working
  graph%valid = .false.
  return
else
  call out_io (s_warn$, r_name, "Unknown x_axis_type")
  graph%valid = .false.
  return
endif
    
call qp_draw_line (graph%x%min, graph%x%max, 0.0_rp, 0.0_rp)

! loop over all elements in the lattice. Only draw those element that
! are within bounds.

ele_shapes => tao_com%ele_shape_lat_layout
height = s%plot_page%text_height * s%plot_page%legend_text_scale

do i = 1, branch%n_ele_max

  ele => branch%ele(i)
  call tao_find_ele_shape (ele, ele_shapes, ix_shape)

  if (ele%lord_status == multipass_lord$) cycle
  if (ele%slave_status == super_slave$) cycle
  if (i > branch%n_ele_track .and. ix_shape < 1) cycle

  if (plot%x_axis_type == 's') then
    call find_element_ends (lat, ele, ele1, ele2)
    if (.not. associated(ele1)) cycle
    if (ele1%ix_branch /= graph%ix_branch) cycle
    x1 = ele1%s
    x2 = ele2%s
    ! If out of range then try a negative position
    if (branch%param%lattice_type == circular_lattice$ .and. x1 > graph%x%max) then
      x1 = x1 - lat_len
      x2 = x2 - lat_len
    endif

  elseif (plot%x_axis_type == 'index') then
    ! shouldn't be here!
    call out_io (s_error$, r_name, "Shouldn't be here!")
    graph%valid = .false.
    return

  elseif (plot%x_axis_type == 'ele_index') then
    x1 = i
    x2 = i

  else
    call out_io (s_warn$, r_name, "Unknown x_axis_type")
    graph%valid = .false.
    return
  endif
    
  if (x1 > graph%x%max) cycle
  if (x2 < graph%x%min) cycle

  ! Only those elements with ix_shape > 0 are to be drawn.
  ! All others have the zero line drawn through them.


  if (ix_shape < 1) cycle
  ele_shape => ele_shapes(ix_shape)
  shape_name = ele_shape%shape

  ! Here if element is to be drawn...

  ! r1 and r2 are the scale factors for the lines below and above the center line.

  y = ele_shape%dy_pix
  y1 = -y
  y2 =  y
  if (shape_name == 'VAR_BOX' .or. shape_name == 'ASYM_VAR_BOX') then
    select case (ele%key)
    case (quadrupole$)
      y2 = y * ele%value(k1$)
    case (sextupole$)
      y2 = y * ele%value(k2$)
    case (octupole$)
      y2 = y * ele%value(k3$)
    case (solenoid$)
      y2 = y * ele%value(ks$)
    end select
    y2 = max(-s%plot_page%shape_height_max, min(y2, s%plot_page%shape_height_max))
    y1 = -y2
    if (shape_name == 'ASYM_VAR_BOX') y1 = 0
  end if

  y1 = max(y_bottom, min(y1, y_top))
  y2 = max(y_bottom, min(y2, y_top))

  call draw_this_shape (ele%name, ele%s - ele%value(l$) / 2, ele_shape)

enddo

! Draw data

do i = 1, size(ele_shapes)
  if (plot%x_axis_type /= 's') exit
  ele_shape => ele_shapes(i)
  if (ele_shape%ele_name(1:5) /= 'dat::') cycle
  call tao_find_data (err, ele_shape%ele_name, d_array = d_array, log_array = logic_array)
  if (err) cycle
  do j = 1, size(d_array)
    datum => d_array(j)%d
    if (datum%ix_branch /= graph%ix_branch) cycle
    if (size(logic_array) /= 0) then
      if (.not. logic_array(j)%l) cycle
    endif
    x0 = datum%s 
    if (x0 > graph%x%max) cycle
    if (x0 < graph%x%min) cycle
    y1 = ele_shape%dy_pix
    y1 = max(y_bottom, min(y1, y_top))
    y2 = -y1
    call qp_convert_point_rel (dummy, y1, 'DATA', dummy, y, 'INCH') 
    call qp_convert_point_rel (y, dummy, 'INCH', dx, dummy, 'DATA')
    x1 = x0 - dx
    x2 = x0 + dx
    call draw_this_shape (tao_datum_name(datum), datum%s, ele_shape)
  enddo
enddo

! Draw variables

do i = 1, size(ele_shapes)
  if (plot%x_axis_type /= 's') exit
  ele_shape => ele_shapes(i)
  if (ele_shape%ele_name(1:5) /= 'var::') cycle
  call tao_find_var (err, ele_shape%ele_name, v_array = v_array, log_array = logic_array)

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

if (s%global%label_keys) then
  do kk = 1, 10
    k = kk + 10*tao_com%ix_key_bank
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
          ele1 => pointer_to_slave (lat, ele, j)
          l2 = ele1%value(l$) / 2
          s_pos = ele1%s - l2
          if (s_pos > graph%x%max .and. s_pos-lat_len > graph%x%min) s_pos = s_pos - lat_len
          if (s_pos + l2 < graph%x%min .or. s_pos - l2 > graph%x%max) cycle
          call qp_draw_text (trim(str), s_pos, y_top, justify = 'CT', height = key_number_height)  
        enddo
      else
        s_pos = ele%s - ele%value(l$)/2
        if (s_pos > graph%x%max .and. s_pos-lat_len > graph%x%min) s_pos = s_pos - lat_len
        if (s_pos + l2 < graph%x%min .or. s_pos - l2 > graph%x%max) cycle
        call qp_draw_text (trim(str), s_pos, y_top, justify = 'CT', height = key_number_height)  
      endif
    enddo
  enddo
endif

!-----------------------------------------------------------
contains 
subroutine draw_this_shape (name_in, s_pos, ele_shape)

type (tao_ele_shape_struct) ele_shape
real(rp) s_pos
integer icol
character(*) name_in
character(20) shape_name

!

shape_name = ele_shape%shape
shape_has_box = (index(shape_name, 'BOX') /= 0)
call qp_translate_to_color_index (ele_shape%color, icol)


! Draw the shape

if (shape_name == 'DIAMOND') then
  call qp_draw_line (x1, (x1+x2)/2, 0.0_rp, y1, color = icol)
  call qp_draw_line (x1, (x1+x2)/2, 0.0_rp, y2, color = icol)
  call qp_draw_line (x2, (x1+x2)/2, 0.0_rp, y1, color = icol)
  call qp_draw_line (x2, (x1+x2)/2, 0.0_rp, y2, color = icol)
endif

if (shape_name == 'CIRCLE') then
  call qp_convert_point_abs ((x1+x2)/2, (y1+y2)/2, 'DATA', x0, y0, 'POINTS')
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

if (s%global%label_lattice_elements .and. ele_shape%label_type /= 'none') then
  
  if (ele_shape%label_type == 'name') then
    name = name_in
  elseif (ele_shape%label_type == 's') then
    write (name, '(f16.2)') s_pos
    call string_trim (name, name, ix)
  else
    call out_io (s_error$, r_name, 'BAD ELEMENT LABEL: ' // ele_shape%label_type)
    call err_exit
  endif 

  y_off = y_bottom   
  if (s_pos > graph%x%max .and. s_pos-lat_len > graph%x%min) s_pos = s_pos - lat_len
  call qp_draw_text (name, s_pos, y_off, height = height, justify = 'LC', ANGLE = 90.0_rp)

endif

end subroutine draw_this_shape

end subroutine tao_plot_lat_layout

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+

subroutine tao_plot_data (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (qp_line_struct), allocatable :: line(:)
type (qp_symbol_struct), allocatable :: symbol(:)

integer i, j, k, n
logical have_data
real(rp) x, y, x1

character(16) num_str
character(100), allocatable :: text(:)

! Set scales, margens, etc

call qp_set_layout (box = graph%box, margin = graph%margin)
call qp_set_layout (x_axis = graph%x, x2_mirrors_x = .true.)
call qp_set_layout (y_axis = graph%y, y2_axis = graph%y2, &
                                                y2_mirrors_y = graph%y2_mirrors_y)
if (graph%title == '') then
  call qp_set_graph (title = '')
else
  call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
endif
call qp_draw_axes

! Draw the default x-axis label if there is none. 

if (graph%x%draw_label .and. graph%x%label == '') then
  call qp_to_inch_rel (1.0_rp, 0.0_rp, x1, y, '%GRAPH')
  x = x1 * (graph%x%major_div - 0.5) / graph%x%major_div
  y = -graph%x%number_offset
  call qp_draw_text (plot%x_axis_type, x, y, 'INCH', justify = 'CT')
endif

!

if (.not. graph%valid) return

if (graph%limited .and. graph%clip .and. s%global%draw_curve_off_scale_warn) &
  call qp_draw_text ('**Curve Off Scale**', -0.30_rp, -0.15_rp, '%/GRAPH/RT', color = red$) 

! loop over all the curves of the graph and draw them

have_data = .false.

do k = 1, size(graph%curve)
  curve => graph%curve(k)
  if (curve%use_y2) call qp_use_axis (y = 'Y2')
  call qp_set_symbol (curve%symbol)
  call qp_set_line ('PLOT', curve%line) 

  if (curve%draw_symbols .and. allocated(curve%x_symb)) then
    if (size(curve%x_symb) > 0) have_data = .true.
    call qp_draw_data (curve%x_symb, curve%y_symb, .false., curve%symbol_every, graph%clip)
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
    call qp_draw_data (curve%x_line, curve%y_line, curve%draw_line, 0, graph%clip)
  endif

  call qp_use_axis (y = 'Y')  ! reset

enddo

if (.not. have_data) call qp_draw_text ('**No Plottable Data**', &
                            0.18_rp, -0.15_rp, '%/GRAPH/LT', color = red$) 
! Draw the text legend if there is one

if (any(graph%text_legend /= ' ')) call qp_draw_text_legend (graph%text_legend, &
       graph%text_legend_origin%x, graph%text_legend_origin%y, graph%text_legend_origin%units)

! Draw the curve legend if needed

n = size(graph%curve)
allocate (text(n), symbol(n), line(n))

do i = 1, n
  text(i) = graph%curve(i)%legend_text
  if (text(i) == '') text(i) = graph%curve(i)%data_type
  symbol(i) = graph%curve(i)%symbol
  if (size(graph%curve(i)%x_symb) == 0) symbol(i)%type = -1 ! Do not draw
  if (.not. graph%curve(i)%draw_symbols) symbol(i)%type = -1
  line(i) = graph%curve(i)%line
  if (size(graph%curve(i)%x_line) == 0) line(i)%width = -1 ! Do not draw
  if (.not. graph%curve(i)%draw_line) line(i)%width = -1
enddo

if (graph%draw_curve_legend .and. n > 1) then
  call qp_draw_curve_legend (graph%curve_legend_origin, &
            s%plot_page%curve_legend_text_offset, s%plot_page%curve_legend_line_len, &
            line, symbol, text)
endif

deallocate (text, symbol, line)

end subroutine tao_plot_data

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_find_ele_shape (ele, ele_shapes, ix_shape)

implicit none

type (ele_struct) ele
type (tao_ele_shape_struct) :: ele_shapes(:)

integer k, ix_shape, n_ele_track, ix_class

character(20) :: r_name = 'tao_find_ele_shape'
character(40) ele_name

logical err

!

ix_shape = 0

if (ele%lord_status == group_lord$) return
if (ele%lord_status == overlay_lord$) return
if (ele%slave_status == super_slave$) return

do k = 1, size(ele_shapes)

  if (ele_shapes(k)%ele_name == '') cycle
  if (ele_shapes(k)%ele_name(1:5) == 'dat::') cycle

  call tao_string_to_element_id (ele_shapes(k)%ele_name, ix_class, ele_name, err, .false.)
  if (err) then
    call out_io (s_error$, r_name, 'BAD ELEMENT KEY IN SHAPE: ' // ele_shapes(k)%ele_name)
    cycle
  endif

  if (ix_class /= 0 .and. ix_class /= ele%key) cycle
  if (.not. match_wild(ele%name, ele_name)) cycle

  ix_shape = k
  return
enddo

end subroutine tao_find_ele_shape

end module
