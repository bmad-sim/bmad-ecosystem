module tao_plot_mod

use tao_mod
use quick_plot
use tao_plot_window_mod


contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_plot_out ()
!
! Subroutine to draw the plots on the plot window.
!
! Input:
!-

subroutine tao_plot_out ()

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (qp_rect_struct) border1, border2

real(rp) location(4), dx, dy, h

integer i, j, k

character(16) :: r_name = 'tao_plot_out'
character(3) view_str

logical found

! inits

if (.not. s%global%plot_on) return

h = s%plot_page%text_height
call qp_set_text_attrib ('TEXT', height = h)
call qp_set_text_attrib ('MAIN_TITLE', height = h * s%plot_page%main_title_text_scale)
call qp_set_text_attrib ('GRAPH_TITLE', height = h * s%plot_page%graph_title_text_scale)
call qp_set_text_attrib ('LEGEND', height = h * s%plot_page%legend_text_scale)
call qp_set_text_attrib ('AXIS_NUMBERS', height = h * s%plot_page%axis_number_text_scale)
call qp_set_text_attrib ('AXIS_LABEL', height = h * s%plot_page%axis_label_text_scale)

call tao_create_plot_window () ! This routine knows not to create multiple windows.

call qp_clear_page

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

  do j = 1, size(plot%graph)
    graph => plot%graph(j)

    ! For a non-valid graph just print a message

    if (.not. graph%valid) then
      call qp_set_layout (box = graph%box)
      call qp_draw_text ('Error In The Plot Calculation', &
                                                       0.1_rp, 0.5_rp, '%BOX')
      cycle
    endif

    call tao_hook_plot_graph (plot, graph, found)
    if (found) cycle

    select case (graph%type)
    case ('data', 'phase_space')
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
  enddo

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
type (tao_keyboard_struct), pointer :: key

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

do k = tao_com%ix_key_bank+1, tao_com%ix_key_bank+10
  if (k > ubound(s%key, 1)) cycle
  ix_var = s%key(k)%ix_var
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

write (str, '(i2, a)') tao_com%ix_key_bank/10, ':'
y_here = y_here - 1.1 * height
call qp_draw_text (str, 5.0_rp, y_here, 'POINTS/GRAPH', &
                          height = height, uniform_spacing = .true.)

do i = 1, 10

  k = i + tao_com%ix_key_bank
  if (k > ubound(s%key, 1)) cycle
  key => s%key(k)
  ix_var = key%ix_var
  j = mod(i, 10)

  if (key%ix_var == 0) then
    write (str, '(i2)') j
  else
    var => s%var(ix_var)
    str2 = tao_var_uni_string(var)
    v = maxval(abs( (/ var%model_value, key%val0, key%delta /) ))
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
    write (val0_str,  fmt2) key%val0 / 10.0**n, exp_str
    write (delta_str, fmt2) key%delta / 10.0**n, exp_str

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
type (ele_struct), pointer :: ele
type (floor_position_struct) end1, end2, floor

integer i, j, ix_ptr, icol, ix1, ix2, isu, n_bend, n, ix

real(rp) off, off1, off2, angle, rho, x0, y0, dx1, dy1, dx2, dy2
real(rp) dt_x, dt_y, x_center, y_center, dx, dy, theta
real(rp) x_bend(0:1000), y_bend(0:1000), dx_bend(0:1000), dy_bend(0:1000)
real(rp) v_old(3), w_old(3,3), r_vec(3), dr_vec(3), v_vec(3), dv_vec(3)
real(rp) cos_t, sin_t, cos_p, sin_p, cos_a, sin_a, height
real(rp) x_inch, y_inch

character(80) str
character(40) name
character(20) :: r_name = 'tao_plot_floor_plan'
character(16) shape
character(2) justify

! Each graph is a separate lattice layout (presumably for different universes). 
! setup the placement of the graph on the plot page.

call qp_set_layout (x_axis = plot%x, y_axis = graph%y, &
                    box = graph%box, margin = graph%margin)

! Adjust the margins if there is to be no distortion of the drawing

if (graph%correct_xy_distortion) call qp_eliminate_xy_distortion

!

if (graph%draw_axes) then
  call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
  call qp_draw_axes
endif

isu = graph%ix_universe
! if garph%ix_universe .eq. 0 then graph currently viewed universe
if (isu .eq. 0) then
  lat => s%u(s%global%u_view)%model%lat
else
  lat => s%u(isu)%model%lat
endif

graph%x_max = -1e20
graph%x_min =  1e20
graph%y_max = -1e20
graph%y_min =  1e20
  
! loop over all elements in the lattice. 

do i = 1, lat%n_ele_max

  ele => lat%ele(i)
  ix_ptr = ele%ix_pointer

  if (ele%control_type == super_slave$) cycle
  if (i > lat%n_ele_track .and. ix_ptr < 1) cycle

  call find_element_ends (lat, i, ix1, ix2)
  call floor_to_screen_coords (lat%ele(ix1)%floor, end1)
  call floor_to_screen_coords (lat%ele(ix2)%floor, end2)

  ! Record min and max

  graph%x_max = max(graph%x_max, end1%x, end2%x)
  graph%x_min = min(graph%x_min, end1%x, end2%x)
  graph%y_max = max(graph%y_max, end1%y, end2%y)
  graph%y_min = min(graph%y_min, end1%y, end2%y)

  ! Only draw those element that have at least one point in bounds.
  
  if ((end1%x < plot%x%min .or. plot%x%max < end1%x .or. &
      end1%y < graph%y%min .or. graph%y%max < end1%y) .and. &
      (end2%x < plot%x%min .or. plot%x%max < end2%x .or. &
      end2%y < graph%y%min .or. graph%y%max < end2%y)) cycle

  ! Bends can be tricky if they are not in the X-Z plane. 
  ! Bends are parameterized by a set of points (x_bend, y_bend) along their  
  ! centerline and a set of vectors (dx_bend, dy_bend) perpendicular to the centerline.

  if (ele%key == sbend$) then
    if (ele%value(g$) == 0) then
      n_bend = 1
      x_bend(0) = end1%x; y_bend(0) = end1%y
      x_bend(1) = end2%x; y_bend(1) = end2%y

    else
      floor = lat%ele(ix1)%floor
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

  ! Only those elements with ele%ix_pointer > 0 are to be drawn.
  ! All others are drawn with a line or arc

  if (ix_ptr < 1) then
    if (ele%key == sbend$) then
      call qp_draw_polyline(x_bend(:n_bend), y_bend(:n_bend))
    else
      call qp_draw_line(end1%x, end2%x, end1%y, end2%y)
    endif
    cycle
  endif

  ! Here if element is to be drawn...

  select case (s%plot_page%ele_shape(ix_ptr)%shape)
  case ('BOX', 'VAR_BOX', 'ASYM_VAR_BOX', 'XBOX')
  case default
    print *, 'ERROR: UNKNOWN SHAPE: ', s%plot_page%ele_shape(ix_ptr)%shape
    call err_exit
  end select

  call qp_translate_to_color_index (s%plot_page%ele_shape(ix_ptr)%color, icol)

  shape = s%plot_page%ele_shape(ix_ptr)%shape

  off = s%plot_page%ele_shape(ix_ptr)%dy_pix/2
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
    off2 = off1
    if (shape == 'ASYM_VAR_BOX') off1 = 0
  end if

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

  ! Draw the element.

  call qp_draw_line (end1%x+dx1, end1%x-dx1, end1%y+dy1, end1%y-dy1, &
                                                    units = 'POINTS', color = icol)
  call qp_draw_line (end2%x+dx2, end2%x-dx2, end2%y+dy2, end2%y-dy2, &
                                                    units = 'POINTS', color = icol)

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

  if (s%plot_page%ele_shape(ix_ptr)%shape == 'XBOX') then
    call qp_draw_line (end1%x+dx1, end2%x-dx2, end1%y+dy1, end2%y-dy2, &
                                                    units = 'POINTS', color = icol)
    call qp_draw_line (end1%x-dx1, end2%x+dx2, end1%y-dy1, end2%y+dy2, &
                                                    units = 'POINTS', color = icol)
  endif

  ! Draw the label.
  ! Since multipass slaves are on top of one another, just draw the multipass lord's name.
  ! Also place a bend's label to the outside of the bend.

  if (s%plot_page%ele_shape(ix_ptr)%draw_name) then
    if (ele%control_type == multipass_slave$) then
      ix = ele%ic1_lord
      ix = lat%control(lat%ic(ix))%ix_lord
      name = lat%ele(ix)%name
    else
      name = ele%name
    endif

    if (ele%key /= sbend$ .or. ele%value(g$) == 0) then
      x_center = (end1%x + end2%x) / 2 
      y_center = (end1%y + end2%y) / 2 
      dx = -2 * off2 * dt_y / sqrt(dt_x**2 + dt_y**2)
      dy =  2 * off2 * dt_x / sqrt(dt_x**2 + dt_y**2)
    else
      n = n_bend / 2
      x_center = x_bend(n) 
      y_center = y_bend(n) 
      dx = -2 * off2 * dx_bend(n) / sqrt(dx_bend(n)**2 + dy_bend(n)**2)
      dy = -2 * off2 * dy_bend(n) / sqrt(dx_bend(n)**2 + dy_bend(n)**2)
    endif
    theta = modulo2 (atan2d(dy, dx), 90.0_rp)
    if (dx > 0) then
      justify = 'LC'
    else
      justify = 'RC'
    endif
    height = s%plot_page%text_height * s%plot_page%legend_text_scale
    call qp_draw_text (name, x_center+dx, y_center+dy, units = 'POINTS', &
                                 height = height, justify = justify, ANGLE = theta)    
  endif

enddo


!--------------------------------------------------------------------------
contains

subroutine floor_to_screen_coords (floor, screen)

type (floor_position_struct) floor, screen

!

call floor_to_screen (floor%x, floor%y, floor%z, screen%x, screen%y)
screen%theta = pi + floor%theta

end subroutine

!--------------------------------------------------------------------------
! contains

subroutine floor_to_screen (x_floor, y_floor, z_floor, x_screen, y_screen)

real(rp) x_floor, y_floor, z_floor, x_screen, y_screen

! Mapping from floor coords to screen coords is:
!   Floor   Screen 
!    z   ->  -x
!    x   ->  -y

x_screen = -z_floor
y_screen = -x_floor

end subroutine

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_plot_lat_layout (plot, graph)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele

real(rp) x1, x2, y1, y2, y, s_pos, y_off, y_bottom, y_top
real(rp) lat_len, height

integer i, j, ix_ptr, k, kk, ix, ix1, ix2, isu
integer icol, ix_var, ixv, j_label

character(80) str
character(20) :: r_name = 'tao_plot_lat_layout'
character(16) shape

! Each graph is a separate lattice layout (presumably for different universes). 
! setup the placement of the graph on the plot page.

call qp_set_layout (x_axis = plot%x, box = graph%box, margin = graph%margin)
call qp_get_layout_attrib ('GRAPH', x1, x2, y1, y2, 'POINTS/GRAPH')
y_bottom = -(y2 + 12 * (s%global%n_lat_layout_label_rows - 1)) / 2
y_top = y_bottom + y2
call qp_set_axis ('Y', y_bottom, y_top, 1, 0)
height = s%plot_page%text_height * s%plot_page%legend_text_scale

isu = graph%ix_universe
! if garph%ix_universe .eq. 0 then graph currently viewed universe
if (isu == 0) then
  lat => s%u(s%global%u_view)%model%lat
else
  lat => s%u(isu)%model%lat
endif
lat_len = lat%param%total_length
  
! Figure out x axis
! If it's not 's' then just draw a vertical line at the proper index

  if (plot%x_axis_type .eq. 's') then
    ! continue
  elseif (plot%x_axis_type .eq. 'index') then
    ! cannot plot layout in this case!
    graph%valid = .false.
    return
  elseif (plot%x_axis_type .eq. 'ele_index') then
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

j_label = 0

do i = 1, lat%n_ele_max

  ele => lat%ele(i)
  ix_ptr = ele%ix_pointer

  if (ele%control_type == multipass_lord$) cycle
  if (ele%control_type == super_slave$) cycle
  if (i > lat%n_ele_track .and. ix_ptr < 1) cycle

  if (plot%x_axis_type .eq. 's') then
    call find_element_ends (lat, i, ix1, ix2)
    x1 = lat%ele(ix1)%s
    x2 = lat%ele(ix2)%s
    ! If out of range then try a negative position
    if (lat%param%lattice_type == circular_lattice$ .and. x1 > plot%x%max) then
      x1 = x1 - lat_len
      x2 = x2 - lat_len
    endif

  elseif (plot%x_axis_type .eq. 'index') then
    ! shouldn't be here!
    call out_io (s_error$, r_name, "Shouldn't be here!")
    graph%valid = .false.
    return

  elseif (plot%x_axis_type .eq. 'ele_index') then
    x1 = i
    x2 = i

  else
    call out_io (s_warn$, r_name, "Unknown x_axis_type")
    graph%valid = .false.
    return
  endif
    

  if (x1 > plot%x%max) cycle
  if (x2 < plot%x%min) cycle

  ! Only those elements with ele%ix_pointer > 0 are to be drawn.
  ! All others have the zero line drawn through them.

  if (ix_ptr < 1) then
    call qp_draw_line (x1, x2, 0.0_rp, 0.0_rp)
    cycle
  endif

  ! Here if element is to be drawn...

  select case (s%plot_page%ele_shape(ix_ptr)%shape)
  case ('BOX', 'VAR_BOX', 'ASYM_VAR_BOX', 'XBOX')
  case default
    print *, 'ERROR: UNKNOWN SHAPE: ', s%plot_page%ele_shape(ix_ptr)%shape
    call err_exit
  end select

  call qp_translate_to_color_index (s%plot_page%ele_shape(ix_ptr)%color, icol)

  shape = s%plot_page%ele_shape(ix_ptr)%shape

  ! r1 and r2 are the scale factors for the lines below and above the center line.

  y = s%plot_page%ele_shape(ix_ptr)%dy_pix/2
  y1 = y
  y2 = y
  if (shape == 'VAR_BOX' .or. shape == 'ASYM_VAR_BOX') then
    select case (ele%key)
    case (quadrupole$)
      y1 = y * ele%value(k1$)
    case (sextupole$)
      y1 = y * ele%value(k2$)
    case (octupole$)
      y1 = y * ele%value(k3$)
    case (solenoid$)
      y1 = y * ele%value(ks$)
    end select
    y2 = y1
    if (shape == 'ASYM_VAR_BOX') y1 = 0
  end if

  ! Draw the shape

  call qp_draw_rectangle (x1, x2,  -y1, y2, color = icol)

  if (s%plot_page%ele_shape(ix_ptr)%shape == 'XBOX') then
    call qp_draw_line (x1, x2,  y2, -y1, color = icol)
    call qp_draw_line (x1, x2, -y1,  y2, color = icol)
  endif

  ! Put on a label
  
  if (s%global%label_lattice_elements .and. s%plot_page%ele_shape(ix_ptr)%draw_name) then
    j_label = j_label + 1
    if (j_label == s%global%n_lat_layout_label_rows) j_label = 0
    y_off = y_bottom   ! + 12.0_rp * j_label 
    s_pos = ele%s - ele%value(l$)/2
    if (s_pos > plot%x%max .and. s_pos-lat_len > plot%x%min) s_pos = s_pos - lat_len
    call qp_draw_text (ele%name, s_pos, y_off, &
                                 height = height, justify = 'CC', ANGLE = 90.0_rp)
  endif

  call qp_draw_line (x1, x2, 0.0_rp, 0.0_rp)

enddo

! This is for drawing the key numbers under the appropriate elements

if (s%global%label_keys) then
  do k = tao_com%ix_key_bank+1, tao_com%ix_key_bank+10
    if (k > ubound(s%key, 1)) cycle
    ix_var = s%key(k)%ix_var
    if (ix_var < 1) cycle
    do ixv = 1, size(s%var(ix_var)%this)
      if (s%var(ix_var)%this(ixv)%ix_uni /= isu) cycle
      ix = s%var(ix_var)%this(ixv)%ix_ele
      kk = mod(k - tao_com%ix_key_bank, 10)
      write (str, '(i1)') kk
      if (ix > lat%n_ele_track) then
        do j = lat%ele(ix)%ix1_slave, lat%ele(ix)%ix2_slave
          ix1 = lat%control(j)%ix_slave
          s_pos = lat%ele(ix1)%s - lat%ele(ix1)%value(l$)/2
          if (s_pos > plot%x%max .and. s_pos-lat_len > plot%x%min) s_pos = s_pos - lat_len
          call qp_draw_text (trim(str), s_pos, y_top, &
                              justify = 'CT', height = 10.0_rp)  
        enddo
      else
        s_pos = lat%ele(ix)%s - lat%ele(ix)%value(l$)/2
        call qp_draw_text (trim(str), s_pos, y_top, &
                                justify = 'CT', height = 10.0_rp)  
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

integer k
logical have_data

!

call qp_set_layout (box = graph%box, margin = graph%margin)
call qp_set_layout (x_axis = plot%x, x2_mirrors_x = .true.)
call qp_set_layout (y_axis = graph%y, y2_axis = graph%y2, &
                                                y2_mirrors_y = graph%y2_mirrors_y)
call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
call qp_draw_axes

if (graph%limited .and. graph%clip) &
  call qp_draw_text ('**Limited**', -0.18_rp, -0.15_rp, '%/GRAPH/RT', color = red$) 

! loop over all the curves of the graph and draw them

have_data = .false.

do k = 1, size(graph%curve)
  curve => graph%curve(k)
  if (curve%use_y2) call qp_use_axis (y = 'Y2')
  call qp_set_symbol (curve%symbol)
  call qp_set_line ('PLOT', curve%line) 
  if (curve%draw_symbols .and. allocated(curve%x_symb)) then
    if (size(curve%x_symb) > 0) then
      call qp_draw_data (curve%x_symb, curve%y_symb, .false., curve%symbol_every, graph%clip)
      have_data = .true.
    endif
  endif
  if (curve%draw_line .and. allocated(curve%x_line)) then
    if (size(curve%x_line) > 0) then
      call qp_draw_data (curve%x_line, curve%y_line, curve%draw_line, 0, graph%clip)
      have_data = .true.
    endif
  endif
  call qp_use_axis (y = 'Y')  ! reset
enddo

if (.not. have_data) call qp_draw_text ('**No Plottable Data**', &
                            0.18_rp, -0.15_rp, '%/GRAPH/LT', color = red$) 
! draw the legend if there is one

if (any(graph%legend /= ' ')) call qp_draw_legend (graph%legend, &
          graph%legend_origin%x, graph%legend_origin%y, graph%legend_origin%units)

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

function tao_var_uni_string (var) result (str)

implicit none

type (tao_var_struct) var
character(10) str
integer i, iu, ct
logical uni(100)

!

iu = size(s%u)
uni = .false.

do i = 1, size (var%this)
  uni(var%this(i)%ix_uni) = .true.
enddo

ct = count(uni(1:iu))

if (ct == 1) then
  write (str, '(i2)') var%this(1)%ix_uni
elseif (ct == iu) then
  str = 'All'
else
  str = '?'
endif

end function

end module
