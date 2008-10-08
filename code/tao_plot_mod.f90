module tao_plot_mod

use tao_mod
use quick_plot
use tao_plot_window_mod


contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_plots ()
!
! Subroutine to draw the plots on the plot window.
!
! Input:
!-

subroutine tao_draw_plots ()

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

logical found, err, beam_source

! inits

if (.not. s%global%plot_on) return
call tao_create_plot_window () ! This routine knows not to create multiple windows.
call qp_clear_page

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

  plot => s%plot_region(i)%plot
  if (.not. s%plot_region(i)%visible) cycle

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

    ! For a non-valid graph just print a message

    if (.not. graph%valid) then
      call qp_set_layout (box = graph%box)
      text = 'Error In The Plot Calculation'
      if (graph%why_invalid /= '') text = graph%why_invalid
      call qp_draw_text (text, 0.1_rp, 0.5_rp, '%BOX', color = red$)
      cycle
    endif

    ! Now we can draw the graph

    call tao_hook_draw_graph (plot, graph, found)
    if (found) cycle

    select case (graph%type)
    case ('data', 'phase_space', 'data_slice')
      call tao_plot_data (plot, graph)
    case ('lat_layout')
      call tao_plot_lat_layout (plot, graph)
    case ('key_table')
      call tao_plot_key_table (plot, graph)
    case ('floor_plan')
      call tao_plot_floor_plan (plot, graph)
    case default
      call out_io (s_fatal$, r_name, 'UNKNOWN GRAPH TYPE: ' // graph%type)
    end select

  enddo g_loop

enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_plot_key_table (plot, graph)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (tao_var_struct), pointer :: var

integer i, j, k, n, m, p, j_ele, j_att, ix_var
real(rp) :: y_here, norm, v, x1, x2, y1, y2
real(rp) :: height
character(120) str, str2
character(60) fmt, fmt2
character(12) model_str, val0_str, delta_str
character(4) exp_str
character(24) :: r_name = 'tao_plot_key_table'

!

j_ele = 4
j_att = 5

call qp_set_layout (box = graph%box, margin = graph%margin)

call qp_get_layout_attrib ('GRAPH', x1, x2, y1, y2, 'POINTS/GRAPH')
y_here = y2  ! start from the top of the graph

height = s%plot_page%text_height * s%plot_page%key_table_text_scale

do k = 10*tao_com%ix_key_bank+1, 10*tao_com%ix_key_bank+10
  if (k > ubound(s%key, 1)) cycle
  ix_var = s%key(k)
  if (ix_var == 0) cycle
  j_ele = max(j_ele, len_trim(s%var(ix_var)%ele_name))
  j_att = max(j_att, len_trim(s%var(ix_var)%attrib_name))
enddo

write (fmt, '(a, i5, a, i2, a)') &
                        '(4x, a, ', j_ele-2, 'x, a, ', j_att-1, 'x, a)'
write (str, fmt) 'ix  Name', 'Attrib', &
                         'Value      Value0       Delta   Uni  useit_opt'

call qp_draw_text (str, 5.0_rp, y_here, 'POINTS/GRAPH', &
                             height = height, uniform_spacing = .true.)
  
write (fmt, '(a, i2.2, a, i2.2, a)') &
        '(i2, 2x, a', j_ele, ', 2x, a', j_att, ', 3a12, 2x, a, 3x, l)'

write (str, '(i2, a)') tao_com%ix_key_bank, ':'
y_here = y_here - 1.1 * height
call qp_draw_text (str, 5.0_rp, y_here, 'POINTS/GRAPH', &
                          height = height, uniform_spacing = .true.)

do i = 1, 10

  k = i + 10*tao_com%ix_key_bank
  if (k > ubound(s%key, 1)) cycle
  ix_var = s%key(k)
  j = mod(i, 10)

  if (ix_var == 0) then
    write (str, '(i2)') j
  else
    var => s%var(ix_var)
    str2 = tao_var_uni_string(var)
    v = maxval(abs( (/ var%model_value, var%key_val0, var%key_delta /) ))
    if (v == 0) then
      n = 0
      m = 3
    else
      m = 1.001 * log10(v)
      n = 3 * floor(m/3.0)
      p = 3 - (m - n)
    endif

    if (m >= -1 .and. m <= 1) then
      fmt2 = '(f12.4, a0)'
      n = 0
    elseif (m == 2) then
      fmt2 = '(f12.2, a0)'
      n = 0
    else
      write (fmt2, '(a, i1, a)') '(f8.', p, ', a)'
      write (exp_str, '(a, i3.2)') 'E', n
      if (exp_str(2:2) == '0') exp_str(2:2) = '+'
    endif

    write (model_str, fmt2) var%model_value / 10.0**n, exp_str
    write (val0_str,  fmt2) var%key_val0 / 10.0**n, exp_str
    write (delta_str, fmt2) var%key_delta / 10.0**n, exp_str

    write (str, fmt) j, var%ele_name, var%attrib_name, model_str, &
                        val0_str, delta_str, trim(str2), var%useit_opt
  endif

  call qp_draw_text (str, 25.0_rp, y_here, 'POINTS/GRAPH', &
                          height = height, uniform_spacing = .true.)
  y_here = y_here - 1.1 * height
enddo

end subroutine

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

real(rp) theta, v_vec(3)
real(rp) x_bend(0:1000), y_bend(0:1000)

integer i, j, k, n, n_bend, isu, ic, ix_shape

character(20) :: r_name = 'tao_plot_floor_plan'

! Each graph is a separate lattice layout (presumably for different universes). 
! setup the placement of the graph on the plot page.

call qp_set_layout (x_axis = graph%x, y_axis = graph%y, y2_axis = graph%y2, &
                                        box = graph%box, margin = graph%margin)

! Adjust the margins if there is to be no distortion of the drawing

if (graph%correct_xy_distortion) call qp_eliminate_xy_distortion

!

if (graph%draw_axes) then
  call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
  call qp_draw_axes
endif

isu = tao_universe_number(graph%ix_universe)
lat => s%u(isu)%model%lat

! loop over all elements in the lattice. 

do n = 0, ubound(lat%branch, 1)
  lat%branch(n)%ele%logic = .false.  ! Used to mark as drawn.
  do i = 1, lat%branch(n)%n_ele_max
    ele => lat%branch(n)%ele(i)
    call tao_find_ele_shape (ele, tao_com%ele_shape_floor_plan, lat%n_ele_track, ix_shape)
    if (ele%ix_ele > lat%n_ele_track .and. ix_shape == 0) cycle   ! Nothing to draw
    if (ele%control_type == multipass_lord$) then
      do j = ele%ix1_slave, ele%ix2_slave
        ic = lat%control(j)%ix_slave
        call tao_draw_ele_for_floor_plan (plot, graph, lat, lat%branch(n)%ele(ic), ix_shape)
      enddo
    else
      call tao_draw_ele_for_floor_plan (plot, graph, lat, ele, ix_shape)
    endif
  enddo
enddo

! Draw the tunnel wall

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

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

recursive subroutine tao_draw_ele_for_floor_plan (plot, graph, lat, ele, ix_shape)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (lat_struct) :: lat
type (ele_struct) :: ele
type (ele_struct), pointer :: ele1, ele2
type (floor_position_struct) end1, end2, floor
type (tao_wall_point_struct), pointer :: pt(:)
type (tao_ele_shape_struct), pointer :: ele_shape

integer i, j, k, ix_shape, icol, isu, n_bend, n, ix

real(rp) off, off1, off2, angle, rho, x0, y0, dx1, dy1, dx2, dy2
real(rp) dt_x, dt_y, x_center, y_center, dx, dy, theta
real(rp) x_bend(0:1000), y_bend(0:1000), dx_bend(0:1000), dy_bend(0:1000)
real(rp) v_old(3), w_old(3,3), r_vec(3), dr_vec(3), v_vec(3), dv_vec(3)
real(rp) cos_t, sin_t, cos_p, sin_p, cos_a, sin_a, height
real(rp) x_inch, y_inch, x1, x2, y1, y2

character(80) str
character(40) name
character(40) :: r_name = 'tao_draw_ele_for_floor_plan'
character(16) this_shape
character(2) justify

!

call find_element_ends (lat, ele, ele1, ele2)
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

if (ele%key == sbend$) then

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

    n_bend = abs(int(100 * ele%value(angle$))) + 1
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
  if (ele%key == sbend$) then
    call qp_draw_polyline(x_bend(:n_bend), y_bend(:n_bend))
  else
    call qp_draw_line(end1%x, end2%x, end1%y, end2%y)
  endif
  return
endif

! Here if element is to be drawn...

ele_shape => tao_com%ele_shape_floor_plan(ix_shape)

select case (ele_shape%shape)
case ('BOX', 'VAR_BOX', 'ASYM_VAR_BOX', 'XBOX', 'DIAMOND', 'BOW_TIE')
case default
  print *, 'ERROR: UNKNOWN SHAPE: ', ele_shape%shape
  call err_exit
end select

call qp_translate_to_color_index (ele_shape%color, icol)

off = ele_shape%dy_pix
off1 = off
off2 = off
if (ele_shape%shape == 'VAR_BOX' .or. ele_shape%shape == 'ASYM_VAR_BOX') then
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
  if (ele_shape%shape == 'ASYM_VAR_BOX') off1 = 0
endif

! Draw the shape. Since the conversion from floor coords and screen pixels can
! be different along x and y we convert to pixels to make sure that rectangles
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

if (ele%key == sbend$) then
  do j = 0, n_bend
    call qp_convert_point_abs (x_bend(j), y_bend(j), 'DATA', x_bend(j), y_bend(j), 'POINTS')
    call qp_convert_point_rel (dx_bend(j), dy_bend(j), 'DATA', dt_x, dt_y, 'POINTS')
    dx_bend(j) =  off * dt_y / sqrt(dt_x**2 + dt_y**2)
    dy_bend(j) = -off * dt_x / sqrt(dt_x**2 + dt_y**2)
  enddo
endif

! Draw the element...
! Draw top and bottom

if (ele_shape%shape /= 'DIAMOND') then
  if (ele%key == sbend$) then
    call qp_draw_polyline(x_bend(:n_bend) + dx_bend(:n_bend), &
                          y_bend(:n_bend) + dy_bend(:n_bend), units = 'POINTS', color = icol)
    call qp_draw_polyline(x_bend(:n_bend) - dx_bend(:n_bend), &
                          y_bend(:n_bend) - dy_bend(:n_bend), units = 'POINTS', color = icol)

  else
    call qp_draw_line (end1%x+dx1, end2%x+dx2, end1%y+dy1, end2%y+dy2, &
                                                    units = 'POINTS', color = icol)
    call qp_draw_line (end1%x-dx1, end2%x-dx2, end1%y-dy1, end2%y-dy2, &
                                                    units = 'POINTS', color = icol)
  endif
endif

if (ele_shape%shape == 'DIAMOND') then
  if (ele%key == sbend$) then
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

elseif (ele_shape%shape /= 'BOW_TIE') then
  ! Draw sides
  call qp_draw_line (end1%x+dx1, end1%x-dx1, end1%y+dy1, end1%y-dy1, &
                                                  units = 'POINTS', color = icol)
  call qp_draw_line (end2%x+dx2, end2%x-dx2, end2%y+dy2, end2%y-dy2, &
                                                  units = 'POINTS', color = icol)
endif

! Draw X for xbox or bow_tie

if (ele_shape%shape == 'XBOX' .or. ele_shape%shape == 'BOW_TIE') then
  call qp_draw_line (end1%x+dx1, end2%x-dx2, end1%y+dy1, end2%y-dy2, &
                                                  units = 'POINTS', color = icol)
  call qp_draw_line (end1%x-dx1, end2%x+dx2, end1%y-dy1, end2%y+dy2, &
                                                  units = 'POINTS', color = icol)
endif

! Draw the label.
! Since multipass slaves are on top of one another, just draw the multipass lord's name.
! Also place a bend's label to the outside of the bend.

if (ele_shape%label_type == 'name') then
  if (ele%control_type == multipass_slave$) then
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
  theta = modulo2 (atan2d(dy, dx), 90.0_rp)
  if (dx > 0) then
    justify = 'LC'
  else
    justify = 'RC'
  endif
  height = s%plot_page%text_height * s%plot_page%legend_text_scale
  call qp_draw_text (name, x_center+dx*off2, y_center+dy*off2, units = 'POINTS', &
                               height = height, justify = justify, ANGLE = theta)    
endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_plot_lat_layout (plot, graph)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele, ele1, ele2
type (tao_ele_shape_struct), pointer :: ele_shape(:)

real(rp) x1, x2, y1, y2, y, s_pos, y_off, y_bottom, y_top
real(rp) lat_len, height, dy, key_number_height

integer i, j, ix_shape, k, kk, ix, ix1, isu
integer icol, ix_var, ixv

character(80) str
character(40) name
character(20) :: r_name = 'tao_plot_lat_layout'
character(20) this_shape

! Init

isu = tao_universe_number(graph%ix_universe)
lat => s%u(isu)%model%lat
lat_len = lat%param%total_length
key_number_height = 10
  
! Each graph is a separate lattice layout (presumably for different universes). 
! Setup the placement of the graph on the plot page.
! Without labels: Vertically Center the graph.
! With labels: Shift the graph up to give room for the labels.

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
    
! loop over all elements in the lattice. Only draw those element that
! are within bounds.

ele_shape => tao_com%ele_shape_lat_layout
height = s%plot_page%text_height * s%plot_page%legend_text_scale

do i = 1, lat%n_ele_max

  ele => lat%ele(i)
  call tao_find_ele_shape (ele, tao_com%ele_shape_lat_layout, lat%n_ele_track, ix_shape)

  if (ele%control_type == multipass_lord$) cycle
  if (ele%control_type == super_slave$) cycle
  if (i > lat%n_ele_track .and. ix_shape < 1) cycle

  if (plot%x_axis_type == 's') then
    call find_element_ends (lat, ele, ele1, ele2)
    x1 = ele1%s
    x2 = ele2%s
    ! If out of range then try a negative position
    if (lat%param%lattice_type == circular_lattice$ .and. x1 > graph%x%max) then
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

  if (ix_shape < 1) then
    call qp_draw_line (x1, x2, 0.0_rp, 0.0_rp)
    cycle
  endif

  ! Here if element is to be drawn...

  this_shape = ele_shape(ix_shape)%shape

  select case (this_shape)
  case ('BOX', 'VAR_BOX', 'ASYM_VAR_BOX', 'XBOX', 'DIAMOND', 'BOW_TIE')
  case default
    print *, 'ERROR: UNKNOWN SHAPE: ', this_shape
    call err_exit
  end select

  call qp_translate_to_color_index (ele_shape(ix_shape)%color, icol)

  ! r1 and r2 are the scale factors for the lines below and above the center line.

  y = ele_shape(ix_shape)%dy_pix
  y1 = -y
  y2 =  y
  if (this_shape == 'VAR_BOX' .or. this_shape == 'ASYM_VAR_BOX') then
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
    if (this_shape == 'ASYM_VAR_BOX') y1 = 0
  end if

  y1 = max(y_bottom, min(y1, y_top))
  y2 = max(y_bottom, min(y2, y_top))

  ! Draw the shape

  if (this_shape == 'DIAMOND') then
    call qp_draw_line (x1, (x1+x2)/2, 0.0_rp, y1, color = icol)
    call qp_draw_line (x1, (x1+x2)/2, 0.0_rp, y2, color = icol)
    call qp_draw_line (x2, (x1+x2)/2, 0.0_rp, y1, color = icol)
    call qp_draw_line (x2, (x1+x2)/2, 0.0_rp, y2, color = icol)
  elseif (this_shape == 'BOW_TIE') then
    call qp_draw_line (x1, x2, y1, y1, color = icol)
    call qp_draw_line (x1, x2, y2, y2, color = icol)
  else
    call qp_draw_rectangle (x1, x2, y1, y2, color = icol)
  endif

  if (this_shape == 'XBOX' .or. this_shape == 'BOW_TIE') then
    call qp_draw_line (x1, x2, y2, y1, color = icol)
    call qp_draw_line (x1, x2, y1, y2, color = icol)
  endif

  ! Put on a label
  
  if (s%global%label_lattice_elements) then

    
    if (ele_shape(ix_shape)%label_type == 'name') then
      name = ele%name
    elseif (ele_shape(ix_shape)%label_type == 's') then
      write (name, '(f16.2)') ele%s - ele%value(l$) / 2
      call string_trim (name, name, ix)
    elseif (ele_shape(ix_shape)%label_type /= 'none') then
      call out_io (s_error$, r_name, 'BAD ELEMENT LABEL: ' // ele_shape(ix_shape)%label_type)
      call err_exit
    endif 

    if (ele_shape(ix_shape)%label_type /= 'none') then
      y_off = y_bottom   
      s_pos = ele%s - ele%value(l$)/2
      if (s_pos > graph%x%max .and. s_pos-lat_len > graph%x%min) s_pos = s_pos - lat_len
      call qp_draw_text (name, s_pos, y_off, &
                                   height = height, justify = 'LC', ANGLE = 90.0_rp)
    endif

  endif

  call qp_draw_line (x1, x2, 0.0_rp, 0.0_rp)

enddo

! This is for drawing the key numbers under the appropriate elements

if (s%global%label_keys) then
  do kk = 1, 10
    k = kk + 10*tao_com%ix_key_bank
    if (k > ubound(s%key, 1)) cycle
    ix_var = s%key(k)
    if (ix_var < 1) cycle
    do ixv = 1, size(s%var(ix_var)%this)
      if (s%var(ix_var)%this(ixv)%ix_uni /= isu) cycle
      ix = s%var(ix_var)%this(ixv)%ix_ele
      write (str, '(i1)') mod(kk, 10)
      if (ix > lat%n_ele_track) then
        do j = lat%ele(ix)%ix1_slave, lat%ele(ix)%ix2_slave
          ix1 = lat%control(j)%ix_slave
          s_pos = ele1%s - ele1%value(l$)/2
          if (s_pos > graph%x%max .and. s_pos-lat_len > graph%x%min) s_pos = s_pos - lat_len
          call qp_draw_text (trim(str), s_pos, y_top, &
                              justify = 'CT', height = key_number_height)  
        enddo
      else
        s_pos = lat%ele(ix)%s - lat%ele(ix)%value(l$)/2
        call qp_draw_text (trim(str), s_pos, y_top, &
                                justify = 'CT', height = key_number_height)  
      endif
    enddo
  enddo
endif

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+

subroutine tao_plot_data (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve

integer i, k
logical have_data
character(16) num_str
!

call qp_set_layout (box = graph%box, margin = graph%margin)
call qp_set_layout (x_axis = graph%x, x2_mirrors_x = .true.)
call qp_set_layout (y_axis = graph%y, y2_axis = graph%y2, &
                                                y2_mirrors_y = graph%y2_mirrors_y)
call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
call qp_draw_axes

if (graph%limited .and. graph%clip) &
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

if (graph%draw_curve_legend .and. size(graph%curve) > 1) then
  call qp_draw_curve_legend (graph%curve_legend_origin, s%plot_page%curve_legend_text_offset, &
            s%plot_page%curve_legend_line_len, graph%curve%line, graph%curve%symbol, &
            graph%curve%data_type, graph%curve(1)%draw_line, graph%curve(1)%draw_symbols)
endif

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

function tao_var_uni_string (var) result (str)

implicit none

type (tao_var_struct) var
character(10) str
integer i, ct
logical uni(lbound(s%u, 1):ubound(s%u, 1))

!

uni = .false.

do i = 1, size (var%this)
  uni(var%this(i)%ix_uni) = .true.
enddo

ct = count(uni)

if (ct == 1) then
  write (str, '(i3)') var%this(1)%ix_uni
elseif (ct == size(s%u)) then
  str = 'All'
else
  str = '?'
endif

end function

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine tao_find_ele_shape (ele, ele_shapes, n_ele_track, ix_shape)

implicit none

type (ele_struct) ele
type (tao_ele_shape_struct) :: ele_shapes(:)

integer k, ix_shape, n_ele_track, ix_class

character(20) :: r_name = 'tao_find_ele_shape'
character(40) ele_name

logical err

!

ix_shape = 0

if (ele%control_type == group_lord$) return
if (ele%control_type == overlay_lord$) return
if (ele%control_type == super_slave$) return
if (ele%ix_ele > n_ele_track .and. ele%control_type == free$) return

do k = 1, size(ele_shapes)

  if (ele_shapes(k)%ele_name == '') cycle

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

end subroutine

end module
