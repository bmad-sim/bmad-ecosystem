module tao_plot_mod

use quick_plot
use tao_plot_window_mod
use attribute_mod, only: attribute_index

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
type (tao_data_array_struct), allocatable :: d_array(:)

real(rp) location(4), dx, dy, h

integer i, j, k, ic, id, nb

character(80) text
character(*), parameter :: r_name = 'tao_draw_plots'
character(3) default_uni

logical, optional :: do_clear
logical found, err, beam_source

! inits

if (.not. s%global%plot_on) return
call tao_create_plot_window () ! This routine knows not to create multiple windows.
if (logic_option(.true., do_clear)) call qp_clear_page

! It can be helpful for debugging purposes not wait to flush all at the end

if (.not. s%global%debug_on) call qp_wait_to_flush(.true.)

h = s%plot_page%text_height
call qp_set_text_attrib ('TEXT', height = h)
call qp_set_text_attrib ('MAIN_TITLE', height = h * s%plot_page%main_title_text_scale)
call qp_set_text_attrib ('GRAPH_TITLE', height = h * s%plot_page%graph_title_text_scale)
call qp_set_text_attrib ('LEGEND', height = h * s%plot_page%legend_text_scale)
call qp_set_text_attrib ('AXIS_NUMBERS', height = h * s%plot_page%axis_number_text_scale)
call qp_set_text_attrib ('AXIS_LABEL', height = h * s%plot_page%axis_label_text_scale)

! print the title 

h = s%plot_page%text_height * s%plot_page%main_title_text_scale
if (s%plot_page%title%draw_it .and. s%plot_page%title%string /= '') then
  call qp_draw_text (s%plot_page%title%string, s%plot_page%title%x, s%plot_page%title%y, &
                     s%plot_page%title%units, s%plot_page%title%justify, height = h)
endif

if (s%plot_page%subtitle%draw_it .and. s%plot_page%subtitle%string /= '') then
  call qp_draw_text (s%plot_page%subtitle%string, s%plot_page%subtitle%x, s%plot_page%subtitle%y, &
                     s%plot_page%subtitle%units, s%plot_page%subtitle%justify, height = h)
endif

! Draw view universe

if (s%plot_page%draw_graph_title_suffix) then
  nb = 0
  do i = 1, size(s%u)
    nb = max(nb, size(s%u(1)%model%lat%branch))
  enddo
  text = ','
  if (size(s%u) > 1) text = ', Default Universe:' // int_str(s%global%default_universe)
  if (nb > 1) text = trim(text) // ', Default Branch: ' // int_str(s%global%default_branch)
  text = text(3:)
  if (text /= '') call qp_draw_text (text, -2.0_rp, -2.0_rp, 'POINTS/PAGE/RT', 'RT')
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

    if (associated(tao_hook_draw_graph_ptr)) then
      call tao_hook_draw_graph_ptr (plot, graph, found)
      if (found) cycle
    endif

    ! For a non-valid graph or curves print a message

    call qp_set_layout (box = graph%box)

    if (.not. graph%is_valid) then
      call tao_draw_graph_axes(plot, graph)
      call qp_draw_text ('Graph: ' // trim(plot%name) // '.' // trim(graph%name) // '  ' // graph%title, &
                                                       0.5_rp, 0.6_rp, '%BOX', color = 'red', justify = 'CC')
      call qp_draw_text (graph%why_invalid, 0.5_rp, 0.4_rp, '%BOX', color = 'red', justify = 'CC')
      cycle
    endif

    if (allocated(graph%curve)) then
      dy = -0.2
      do ic = 1, size(graph%curve)
        if (graph%curve(ic)%valid) cycle
        if (graph%curve(ic)%why_invalid == 'IGNORE') cycle
        call qp_draw_text (graph%curve(ic)%why_invalid, 0.2_rp, dy, '%BOX/LT', color = 'red', justify = 'LC') 
        dy = dy - 0.15
      enddo
    endif

    ! Draw a rectangle so the box and graph boundries can be seen.

    if (s%global%box_plots) then
      call qp_draw_rectangle (0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp, '%GRAPH/LB')
      call qp_draw_rectangle (0.0_rp, 1.0_rp, 0.0_rp, 1.0_rp, '%BOX/LB')
    endif

    ! If y%min = y%max then was not able to scale the graph due to some problem.
    ! In this case it is not possible to draw the data.

    if (graph%y%min == graph%y%max) then
      call out_io (s_error$, r_name, 'SOMETHING IS WRONG. GRAPH MIN = MAX FOR: ' // tao_graph_name(graph))
      cycle g_loop
    endif

    ! Now we can draw the graph

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
      call out_io (s_error$, r_name, 'UNKNOWN GRAPH TYPE: ' // graph%type)
    end select

  enddo g_loop

enddo

if (.not. s%global%debug_on) call qp_wait_to_flush(.false.)

end subroutine tao_draw_plots

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

call tao_draw_graph_axes (plot, graph)

! loop over all the curves of the graph and draw them

have_data = .false.

do k = 1, size(graph%curve)
  if (.not. graph%curve(k)%valid) cycle
  call tao_draw_histogram_data (plot, graph, graph%curve(k), have_data)
enddo

if (.not. have_data) call qp_draw_text ('**No Plottable Data**', &
                            0.18_rp, -0.15_rp, '%/GRAPH/LT', color = 'red') 

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

real(rp) y0, y1, x_a1, x_a2, x_b1, x_b2

! Draw the data

call tao_plot_data (plot, graph)

! Now draw the rectangles of the fit regions.

y => graph%y
y0 = y%min + 0.1 * (y%max - y%min)
y1 = y%max - 0.1 * (y%max - y%min)

x_a1 = x_val(s%wave%ix_a1)
x_a2 = x_val(s%wave%ix_a2)
x_b1 = x_val(s%wave%ix_b1)
x_b2 = x_val(s%wave%ix_b2)

if (graph%type == 'wave.0' .or. graph%type == 'wave.a') then
  call qp_draw_rectangle (x_a1, x_a2, y0, y1, color = 'blue', width = 2)
endif

if (graph%type == 'wave.0' .or. graph%type == 'wave.b') then
  call qp_draw_rectangle (x_b1, 1.0_rp * x_b2, y0, y1, color = 'blue', width = 2)
endif

!---------------------------------------------------
contains

function x_val(ix_dat) result (x)

type (tao_curve_struct), pointer :: c
real(rp) x
integer ix_dat, i, i0, i1

!

if (plot%x_axis_type == 'index') then
  x = ix_dat
  return
endif

!

c => graph%curve(1)

i0 = -1; i1 = -1
do i = 1, size(c%ix_symb)
  if (c%ix_symb(i) <= ix_dat) i0 = i
  if (c%ix_symb(i) >= ix_dat) then
    i1 = i
    exit
  endif
enddo

if (i0 == -1) then
  x = c%x_symb(i1)
elseif (i1 == -1) then
  x = c%x_symb(i0)
elseif (i0 == i1) then
  x = c%x_symb(i0)
else
  x = (c%x_symb(i0) * (i1 - ix_dat) + c%x_symb(i1) * (ix_dat - i0)) / (i1 - i0)
endif

end function x_val

end subroutine tao_plot_wave

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
character(*), parameter :: r_name = 'tao_plot_key_table'

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

integer isu

character(*), parameter :: r_name = 'tao_draw_floor_plan'

logical err, found

! Each graph is a separate floor plan plot (presumably for different universes). 
! setup the placement of the graph on the plot page.

call qp_save_state(.false.)

call tao_set_floor_plan_axis_label(graph, graph%x, x_ax, 'X')
call tao_set_floor_plan_axis_label(graph, graph%y, y_ax, 'Y')

call qp_set_layout (x_axis = x_ax, y_axis = y_ax, x2_axis = graph%x2, y2_axis = graph%y2, &
                    x2_mirrors_x = .true., y2_mirrors_y = .true., box = graph%box, margin = graph%margin)

if (graph%floor_plan%correct_distortion) call qp_eliminate_xy_distortion

!

if (graph%draw_axes) then
  if (graph%title == '' .or. .not. graph%draw_title) then
    call qp_set_graph (title = '')
  else
    if (s%plot_page%draw_graph_title_suffix) then
      call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
    else
      call qp_set_graph (title = graph%title)
    endif
  endif
  call qp_draw_axes (draw_grid = graph%draw_grid)
endif

! Draw for a particular universe

if (graph%ix_universe == -2) then
  do isu = 1, size(s%u)
    call draw_this_floor_plan(isu)
  enddo
else
  isu = tao_universe_index(graph%ix_universe)
  call draw_this_floor_plan(isu)
endif

! Hook routine for more plotting if desired...

if (associated(tao_hook_draw_floor_plan_ptr)) call tao_hook_draw_floor_plan_ptr (plot, graph)

call qp_restore_state

!-------------------------------------------------------------
contains

subroutine draw_this_floor_plan(isu)

type (tao_ele_shape_struct), pointer :: ele_shape, ele_shape2
type (tao_lattice_struct), pointer :: tao_lat
type (tao_building_wall_point_struct) pt0, pt1
type (floor_position_struct) end1, end2, floor
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, slave
type (lat_struct), pointer :: lat

real(rp) x_min, x_max, y_min, y_max, y1, y2
real(rp) theta, v_vec(3), theta1, dtheta, dat_var_value
real(rp) x_bend(0:400), y_bend(0:400)

integer i, j, k, n, is, ix, n_bend, isu, ic, ib, ix_shape_min
integer ix_pass, n_links, iwidth
logical err

character(40) label_name

!

select case(graph%floor_plan%orbit_lattice)
case ('model');   tao_lat => s%u(isu)%model
case ('design');  tao_lat => s%u(isu)%design
case ('base');    tao_lat => s%u(isu)%base
case default;
  call out_io (s_error$, r_name, 'Bad floor_plan%orbit_lattice: ' // graph%floor_plan%orbit_lattice, &
                                 'Should be one of: "model", "design", or "base"', &
                                 'Will default to "model"')
  tao_lat => s%u(isu)%model
end select

lat => tao_lat%lat

! loop over all elements in the lattice. 

do n = 0, ubound(lat%branch, 1)
  branch => lat%branch(n)
  branch%ele%logic = .false.  ! Used to mark as drawn.
  do i = 0, branch%n_ele_max
    ele => branch%ele(i)
    if (ele%slave_status == super_slave$) cycle

    ix_shape_min = 1
    do
      call tao_ele_shape_info (isu, ele, s%plot_page%floor_plan%ele_shape, ele_shape, label_name, y1, y2, ix_shape_min)
      if (.not. associated(ele_shape) .and. (ele%key == overlay$ .or. &
                                             ele%key == group$ .or. ele%key == girder$)) exit   ! Nothing to draw

      if (graph%floor_plan%draw_only_first_pass .and. ele%slave_status == multipass_slave$) then
        call multipass_chain (ele, ix_pass, n_links)
        if (ix_pass > 1) exit
      endif

      if (ele%lord_status == multipass_lord$) then
        do j = 1, ele%n_slave
          if (graph%floor_plan%draw_only_first_pass .and. j > 1) exit
          slave => pointer_to_slave(ele, j)
          ele_shape2 => tao_pointer_to_ele_shape (isu, slave, s%plot_page%floor_plan%ele_shape)
          if (associated(ele_shape2)) cycle ! Already drawn. Do not draw twice
          call tao_draw_ele_for_floor_plan (plot, graph, tao_lat, slave, ele_shape, label_name, y1, y2)
        enddo
      else
        call tao_draw_ele_for_floor_plan (plot, graph, tao_lat, ele, ele_shape, label_name, y1, y2)
      endif
      if (.not. associated(ele_shape)) exit
      if (.not. ele_shape%multi) exit
    enddo

  enddo
enddo

! Draw the building wall

if (allocated(s%building_wall%section) .and. graph%floor_plan%draw_building_wall) then
  found = .false.

  do ib = 1, size(s%building_wall%section)
    ele_shape => tao_pointer_to_building_wall_shape(s%building_wall%section(ib)%name)
    if (.not. associated(ele_shape)) cycle

    found = .true.
    iwidth = ele_shape%line_width

    do j = 2, size(s%building_wall%section(ib)%point)
      pt0 = tao_oreint_building_wall_pt(s%building_wall%section(ib)%point(j-1))
      pt1 = tao_oreint_building_wall_pt(s%building_wall%section(ib)%point(j))

      if (pt1%radius == 0) then   ! line
        call tao_floor_to_screen (graph, [pt0%x, 0.0_rp, pt0%z], end1%r(1), end1%r(2))
        call tao_floor_to_screen (graph, [pt1%x, 0.0_rp, pt1%z], end2%r(1), end2%r(2))

        ix = max(1, index(ele_shape%shape, '_LINE'))
        call qp_draw_line(end1%r(1), end2%r(1), end1%r(2), end2%r(2), &
                    line_pattern = ele_shape%shape(1:ix-1), width = iwidth, color = ele_shape%color, clip = .true.)

      else                    ! arc
        theta1 = atan2(pt0%x - pt1%x_center, pt0%z - pt1%z_center)
        dtheta = atan2(pt1%x - pt1%x_center, pt1%z - pt1%z_center) - theta1
        if (abs(dtheta) > pi) dtheta = modulo2(dtheta, pi)
        n_bend = abs(50 * dtheta) + 1
        do k = 0, n_bend
          theta = theta1 + k * dtheta / n_bend
          v_vec(1) = pt1%x_center + abs(pt1%radius) * sin(theta)
          v_vec(2) = 0
          v_vec(3) = pt1%z_center + abs(pt1%radius) * cos(theta)
          call tao_floor_to_screen (graph, v_vec, x_bend(k), y_bend(k))
        enddo
        call qp_draw_polyline(x_bend(:n_bend), y_bend(:n_bend), width = iwidth, color = ele_shape%color, clip = .true.)
      endif
    enddo
  enddo  ! wall section

  if (.not. found .and. size(s%building_wall%section) > 0) then
    call out_io (s_info$, r_name, 'Building wall defined but no building wall shape(s) for floor_plan present.')
  endif
endif

! Draw any data curves and beam chamber wall curve

do i = 1, size(graph%curve)
  ! ... needs to be filled in ..
enddo

end subroutine draw_this_floor_plan

end subroutine tao_draw_floor_plan 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine tao_set_floor_plan_axis_label (graph, axis_in, axis_out, which)

type (tao_graph_struct) graph
type (qp_axis_struct) axis_in, axis_out

real(rp) f
integer irot

character(*) which
character(1) x_str, y_str
character(2) label(0:3)

!

axis_out = axis_in
if (axis_out%label /= 'SMART LABEL') return

f = modulo(4*graph%floor_plan%rotation, 1.0_rp)
f = f - fraction(f)

! If rotation is not multiple of 90 degrees then must use blank label.

if (abs(f) > 0.01) then
  axis_out%label = ''
  return
endif

! Normal case

x_str = upcase(graph%floor_plan%view(1:1))
y_str = upcase(graph%floor_plan%view(2:2))
label = [x_str // ' ', '-' // y_str, '-' // x_str, y_str // ' ']

irot = modulo(nint(4*graph%floor_plan%rotation), 4)
if (which == 'Y') irot = modulo(irot-1, 4)
axis_out%label = label(irot)

end subroutine tao_set_floor_plan_axis_label

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_ele_for_floor_plan (plot, graph, tao_lat, ele, ele_shape, label_name, offset1, offset2)
!
! Routine to draw one lattice element or one datum location for the floor plan graph. 
!
! Input:
!   plot              -- tao_plot_struct: Plot containing the graph.
!   graph             -- tao_graph_struct: Graph to plot.
!   tao_lat           -- tao_lattice_struct: Lattice containing the element.
!   ele               -- ele_struct: Element to draw.
!   ele_shape         -- tao_ele_shape_struct: Shape to draw from s%plot_page%floor_plan%ele_shape(:) array.
!                         Will be NULL if no associated shape for this element.
!   label_name        -- character(*): Shape label.
!   offset1, offset2  -- real(rp): Transverse distances used to scale the drawing of the element shape.
!-

recursive subroutine tao_draw_ele_for_floor_plan (plot, graph, tao_lat, ele, ele_shape, label_name, offset1, offset2)

implicit none

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (tao_lattice_struct), target :: tao_lat
type (branch_struct), pointer :: branch
type (lat_struct), pointer :: lat
type (ele_struct) :: ele
type (ele_struct) :: drift
type (ele_struct), pointer :: ele0, ele1, ele2, lord
type (floor_position_struct) end1, end2, f_orb, floor, x_ray, floor1, floor2
type (tao_building_wall_point_struct), pointer :: pt(:)
type (tao_ele_shape_struct), pointer :: ele_shape, branch_shape
type (coord_struct), pointer :: orbit(:)
type (coord_struct) orb_here, orb_start, orb_end
type (tao_shape_pattern_struct), pointer :: pat

integer, parameter :: n_bend_extra = 40, l1 = -n_bend_extra, l2 = 200 + n_bend_extra
integer i, j, k, n_bend, n, ix, ic, n_mid, min1_bend, min2_bend, max1_bend, max2_bend
integer n1, n2

real(rp) offset1, offset2
real(rp) off, off1, off2, angle, rho, dx1, dy1, dx2, dy2, ang, length
real(rp) x, y, dt_x, dt_y, x_center, y_center, dx, dy, theta, e_edge
real(rp) x_bend(l1:l2), y_bend(l1:l2), dx_bend(l1:l2), dy_bend(l1:l2)
real(rp) dx_orbit(0:100), dy_orbit(0:100)
real(rp) dx_e1, dy_e1, dx_e2, dy_e2, x_tangent, y_tangent, scale, bend_scale
real(rp) v_old(3), w_old(3,3), r_vec(3), dr_vec(3), v_vec(3), dv_vec(3)
real(rp) cos_t, sin_t, cos_p, sin_p, cos_a, sin_a, height, x_inch, y_inch, x0, y0
real(rp) x_00, x_01, x_02, x_10, x_12, x_20, x_21, x_22
real(rp) y_00, y_01, y_02, y_10, y_12, y_20, y_21, y_22
real(rp) x_min, x_max, y_min, y_max, xb, yb, s_here, r0(2), r1(2), dr1(2), dr2(2)

character(*) label_name
character(80) str, shape
character(*), parameter :: r_name = 'tao_draw_ele_for_floor_plan'
character(16) prefix
character(8) :: draw_units
character(2) justify

logical is_data_or_var, is_there
logical is_bend, can_test

!

call find_element_ends (ele, ele1, ele2)
if (.not. associated(ele1)) return

orbit => tao_lat%tao_branch(ele1%ix_branch)%orbit

orb_start = orbit(ele1%ix_ele)
orb_end = orbit(ele2%ix_ele)

is_data_or_var = .false.
if (associated(ele_shape)) is_data_or_var = (ele_shape%ele_id(1:6) == 'data::' .or. ele_shape%ele_id(1:5) == 'var::')

is_bend = (ele%key == sbend$ .and. .not. is_data_or_var)

if (ele%key == overlay$ .or. ele%key == group$) then
  floor%r = [0.0_rp, 0.0_rp, ele1%value(l$)]
  floor1 = coords_local_curvilinear_to_floor (floor, ele, .false.)

  floor%r = [0.0_rp, 0.0_rp, ele2%value(l$)]
  floor2 = coords_local_curvilinear_to_floor (floor, ele, .false.)

else
  floor%r = [0.0_rp, 0.0_rp, 0.0_rp]
  floor1 = coords_local_curvilinear_to_floor (floor, ele, .true.)

  floor%r = [0.0_rp, 0.0_rp, ele%value(l$)]
  floor2 = coords_local_curvilinear_to_floor (floor, ele, .true.)
endif

if (is_data_or_var) floor1 = floor2 ! pretend this is zero length element

call tao_floor_to_screen_coords (graph, floor1, end1)
call tao_floor_to_screen_coords (graph, floor2, end2)

! Only draw those element that have at least one point in bounds.

! If qp_eliminate_distortion has been called then the min/max
! values of the actual plot are different from graph%x%min, etc.
! In this case, use the actual min/max values

call qp_get_axis_attrib ('X', x_min, x_max)
call qp_get_axis_attrib ('Y', y_min, y_max)

if ((end1%r(1) < x_min .or. x_max < end1%r(1) .or. end1%r(2) < y_min .or. y_max < end1%r(2)) .and. &
    (end2%r(1) < x_min .or. x_max < end2%r(1) .or. end2%r(2) < y_min .or. y_max < end2%r(2))) return

!

if (graph%floor_plan%size_is_absolute) then
  draw_units = 'DATA'
else
  draw_units = 'POINTS'
endif

! Bends can be tricky if they are not in the X-Z plane. 
! Bends are parameterized by a set of points (x_bend, y_bend) along their  
! centerline and a set of vectors (dx_bend, dy_bend) perpendicular to the centerline.

if (is_bend) then
  ! Start at entrance end (not upstream end)
  if (ele%orientation == 1) then
    floor = floor1
  else
    floor = floor2
  endif

  v_old = floor%r
  call floor_angles_to_w_mat (floor%theta, floor%phi, 0.0_rp, w_old)

  ang    = ele%value(angle$) * ele%orientation
  length = ele%value(l$)     * ele%orientation
  bend_scale = 0

  ! Extra points are calculated since finite e1 or e2 will lengthen or shorten the side curves.
  n_bend = int(100*abs(ele%value(angle$))/pi) + 1
  n_bend = min(n_bend, ubound(x_bend, 1)-n_bend_extra)
  do j = -n_bend_extra, n_bend+n_bend_extra
    angle = j * ang / n_bend
    cos_t = cos(ele%value(ref_tilt_tot$))
    sin_t = sin(ele%value(ref_tilt_tot$))
    cos_a = cos(angle)
    sin_a = sin(angle)
    if (ele%value(g$) == 0) then
      r_vec = length * j * [0, 0, 1]
    else
      r_vec = ele%value(rho$) * [cos_t * cos_one(angle), sin_t * cos_one(angle), sin_a]
    endif
    dr_vec = [-cos_t * sin_a, -sin_t * sin_a, cos_a]   ! Tangent vector
    ! This keeps dr_vec pointing to the inside (important for the labels).
    if (cos_t < 0) dr_vec = -dr_vec
    v_vec = matmul (w_old, r_vec) + v_old
    dv_vec = matmul (w_old, dr_vec) 
    call tao_floor_to_screen (graph, v_vec, x_bend(j), y_bend(j))
    call tao_floor_to_screen (graph, dv_vec, x_tangent, y_tangent)

    ! Construct normals to the e1 and e2 faces.
    ! Notice that unlike (dx_bend, dy_bend), these normals are not normalized to 1.

    if (j == 0) then
      e_edge = ele%value(e1$)
      if (ele%orientation == -1) e_edge = -ele%value(e2$)
      dr_vec = tan(e_edge) * [cos_t * cos_a, sin_t * cos_a, sin_a]  ! Radial vector
      dv_vec = matmul (w_old, dr_vec) 
      call tao_floor_to_screen (graph, dv_vec, dx1, dy1)
      x = x_tangent - dx1
      y = y_tangent - dy1
      call qp_convert_point_rel (x, y, 'DATA', dt_x, dt_y, draw_units)
      dx_e1 =  norm2([x,y]) * dt_y / norm2([dt_x, dt_y])
      dy_e1 = -norm2([x,y]) * dt_x / norm2([dt_x, dt_y])
    endif

    if (j == n_bend) then
      e_edge = ele%value(e2$)
      if (ele%orientation == -1) e_edge = -ele%value(e1$)
      dr_vec = tan(e_edge) * [cos_t * cos_a, sin_t * cos_a, sin_a]
      dv_vec = matmul (w_old, dr_vec) 
      call tao_floor_to_screen (graph, dv_vec, dx1, dy1)
      x = x_tangent + dx1
      y = y_tangent + dy1
      call qp_convert_point_rel (x, y, 'DATA', dt_x, dt_y, draw_units)
      dx_e2 =  norm2([x,y]) * dt_y / norm2([dt_x, dt_y])
      dy_e2 = -norm2([x,y]) * dt_x / norm2([dt_x, dt_y])
    endif

    ! Convert tangent vector to perpendicular

    call qp_convert_point_abs (x_bend(j), y_bend(j), 'DATA', x_bend(j), y_bend(j), draw_units)
    call qp_convert_point_rel (x_tangent, y_tangent, 'DATA', dt_x, dt_y, draw_units)
    scale = norm2([dt_x, dt_y])
    dx_bend(j) =  dt_y / scale 
    dy_bend(j) = -dt_x / scale
    bend_scale = max(bend_scale, scale)
  enddo
endif

! Draw orbit?

call qp_save_state(.false.)
call qp_set_line_attrib ('STD', graph%floor_plan%orbit_width, graph%floor_plan%orbit_color, graph%floor_plan%orbit_pattern)

if (graph%floor_plan%orbit_scale /= 0 .and. ele%value(l$) /= 0) then
  if (is_bend) then
    n = int(100*abs(ele%value(angle$))/pi) + int(100 * graph%floor_plan%orbit_scale * &
                  (abs(orb_end%vec(2) - orb_start%vec(2)) + abs(orb_end%vec(4) - orb_start%vec(4))))
    n = max(min(n, ubound(dx_orbit, 1)), 1)
    do j = 0, n
      s_here = j * ele%value(l$) / n
      call twiss_and_track_intra_ele (ele, ele%branch%param, 0.0_rp, s_here, &
                                                       .true., .true., orb_start, orb_here)
      f_orb%r(1:2) = graph%floor_plan%orbit_scale * orb_here%vec(1:3:2)
      f_orb%r(3) = s_here
      f_orb = coords_local_curvilinear_to_floor (f_orb, ele, .false., relative_to = upstream_end$)
      call tao_floor_to_screen (graph, f_orb%r, dx_orbit(j), dy_orbit(j))
    enddo
    call qp_draw_polyline(dx_orbit(0:n), dy_orbit(0:n))

  elseif (ele%key == patch$) then
    ele0 => pointer_to_next_ele (ele, -1)
    floor%r(1:2) = graph%floor_plan%orbit_scale * orb_start%vec(1:3:2)
    floor%r(3) = ele0%value(l$)
    floor1 = coords_local_curvilinear_to_floor (floor, ele0, .false.)
    call tao_floor_to_screen_coords (graph, floor1, f_orb)
    dx_orbit(0) = f_orb%r(1)
    dy_orbit(0) = f_orb%r(2)

    floor%r(1:2) = graph%floor_plan%orbit_scale * orb_end%vec(1:3:2)
    floor%r(3) = ele%value(l$)
    floor1 = coords_local_curvilinear_to_floor (floor, ele, .false., relative_to = downstream_end$)
    call tao_floor_to_screen_coords (graph, floor1, f_orb)
    dx_orbit(1) = f_orb%r(1)
    dy_orbit(1) = f_orb%r(2)

    call qp_draw_polyline(dx_orbit(0:1), dy_orbit(0:1))

  else
    n = int(100 * (abs(orb_end%vec(2) - orb_start%vec(2)) + abs(orb_end%vec(4) - orb_start%vec(4)))) + &
                      int(ele%value(num_steps$)) + 1
    n = min(n, ubound(dx_orbit, 1))
    n = nint(min(1.0*n, 1 + 0.1 * ele%value(l$) / bmad_com%significant_length))
    do ic = 0, n
      s_here = ic * ele%value(l$) / n
      call twiss_and_track_intra_ele (ele, ele%branch%param, 0.0_rp, s_here, &
                                                 .true., .true., orb_start, orb_here)
      floor%r(1:2) = graph%floor_plan%orbit_scale * orb_here%vec(1:3:2)
      floor%r(3) = s_here
      floor1 = coords_local_curvilinear_to_floor (floor, ele, .false., relative_to = upstream_end$)
      call tao_floor_to_screen_coords (graph, floor1, f_orb)
      dx_orbit(ic) = f_orb%r(1)
      dy_orbit(ic) = f_orb%r(2)
    enddo

    call qp_draw_polyline(dx_orbit(0:n), dy_orbit(0:n))
  endif
endif

! coords_local_curvilinear_to_floor does not handle patch elements 
! correctly (this will be fixed) so just ignore patch elements.

if (ele%key == patch$) then
  call qp_restore_state
  return
endif

! Only those elements with an associated ele_shape are to be drawn in full.
! All others are drawn with a simple line or arc.

call qp_set_line_attrib ('STD', 1, 'black', 'solid')

is_there = .false.
if (associated(ele_shape)) is_there = ele_shape%draw
if (.not. is_there) then
  if (is_bend) then
    call qp_draw_polyline(x_bend(1:n_bend), y_bend(1:n_bend))
  else
    call qp_draw_line(end1%r(1), end2%r(1), end1%r(2), end2%r(2))
  endif
  call qp_restore_state
  return
endif

! Here if element is to be drawn...

off = ele_shape%size * s%plot_page%floor_plan_shape_scale 
off1 = offset1 * s%plot_page%floor_plan_shape_scale
off2 = offset2 * s%plot_page%floor_plan_shape_scale

! Draw the shape. Since the conversion from floor coords to screen coords can
! be different along x and y, we convert to screen coords to make sure that rectangles
! remain rectangular.

call qp_convert_point_abs (end1%r(1), end1%r(2), 'DATA', end1%r(1), end1%r(2), draw_units)
call qp_convert_point_abs (end2%r(1), end2%r(2), 'DATA', end2%r(1), end2%r(2), draw_units)

if (is_bend) then
  if (ele%value(g$) > 0) then
    if (off1 * ele%value(g$) > scale) off1 = 0.9 * scale * ele%value(rho$)
  elseif (ele%value(g$) < 0) then
    if (-off2 * ele%value(g$) > scale) off2 = -0.9 * scale * ele%value(rho$)
  endif

  ! Find bounds for drawing bend top and bottom.
  ! First look at top curve.

  r0 = [x_bend(0) + off1 * dx_e1, y_bend(0) + off1 * dy_e1]
  dr2 = [x_bend(1) - x_bend(-1), y_bend(1) - y_bend(-1)]

  do j = -n_bend_extra, n_bend
    dr1  = [x_bend(j) + off1 * dx_bend(j), y_bend(j) + off1 * dy_bend(j)] - r0
    if (j == -n_bend_extra) can_test = (dot_product(dr1, dr2) <= 0)
    if (can_test .and. dot_product(dr1, dr2) > 0) exit
    can_test = (dot_product(dr1, dr2) <= 0)
  enddo
  min1_bend = j

  r0 = [x_bend(n_bend) + off1 * dx_e2, y_bend(n_bend) + off1 * dy_e2]
  dr2 = [x_bend(n_bend-1) - x_bend(n_bend+1), y_bend(n_bend-1) - y_bend(n_bend+1)]

  do j = n_bend+n_bend_extra, 0, -1
    dr1  = [x_bend(j) + off1 * dx_bend(j), y_bend(j) + off1 * dy_bend(j)] - r0
    if (j == n_bend+n_bend_extra) can_test = (dot_product(dr1, dr2) <= 0)
    if (can_test .and. dot_product(dr1, dr2) > 0) exit
    can_test = (dot_product(dr1, dr2) <= 0)
  enddo
  max1_bend = j

  ! Now look at the bottom curve

  r0 = [x_bend(0) - off2 * dx_e1, y_bend(0) - off2 * dy_e1]
  dr2 = [x_bend(1) - x_bend(-1), y_bend(1) - y_bend(-1)]

  do j = -n_bend_extra, n_bend
    dr1  = [x_bend(j) - off2 * dx_bend(j), y_bend(j) - off2 * dy_bend(j)] - r0
    if (j == -n_bend_extra) can_test = (dot_product(dr1, dr2) <= 0)
    if (can_test .and. dot_product(dr1, dr2) > 0) exit
    can_test = (dot_product(dr1, dr2) <= 0)
  enddo
  min2_bend = j

  r0 = [x_bend(n_bend) - off2 * dx_e2, y_bend(n_bend) - off2 * dy_e2]
  dr2 = [x_bend(n_bend-1) - x_bend(n_bend+1), y_bend(n_bend-1) - y_bend(n_bend+1)]

  do j = n_bend+n_bend_extra, 0, -1
    dr1  = [x_bend(j) - off2 * dx_bend(j), y_bend(j) - off2 * dy_bend(j)] - r0
    if (j == n_bend+n_bend_extra) can_test = (dot_product(dr1, dr2) <= 0)
    if (can_test .and. dot_product(dr1, dr2) > 0) exit
    can_test = (dot_product(dr1, dr2) <= 0)
  enddo
  max2_bend = j
endif

! dx1, etc. are offsets perpendicular to the refernece orbit

call qp_convert_point_rel (cos(end1%theta), sin(end1%theta), 'DATA', dt_x, dt_y, draw_units)
dx1 =  off1 * dt_y / norm2([dt_x, dt_y])
dy1 = -off1 * dt_x / norm2([dt_x, dt_y])

call qp_convert_point_rel (cos(end2%theta), sin(end2%theta), 'DATA', dt_x, dt_y, draw_units)
dx2 =  off2 * dt_y / norm2([dt_x, dt_y])
dy2 = -off2 * dt_x / norm2([dt_x, dt_y])

! Draw the element...

call qp_set_line_attrib ('STD', ele_shape%line_width, ele_shape%color)

ix = index(ele_shape%shape, ':')
if (ix == 0) then
  prefix = ''
  shape = ele_shape%shape
else
  prefix = ele_shape%shape(:ix-1)
  shape = ele_shape%shape(ix+1:)
endif

! Draw diamond, etc

if (shape == 'diamond' .or. shape(3:) == 'triangle') then
  if (is_bend) then
    x_00 = x_bend(1) - off * dx_bend(1);       y_00 = y_bend(1) - off * dy_bend(1)
    x_01 = x_bend(1);                          y_01 = y_bend(1)
    x_02 = x_bend(1) + off * dx_bend(1);       y_02 = y_bend(1) + off * dy_bend(1)
    n = n_bend / 2
    x_10 = (x_bend(n) - off * dx_bend(n)) / 2; y_10 = (y_bend(n) - off * dy_bend(n)) / 2
    x_12 = (x_bend(n) + off * dx_bend(n)) / 2; y_12 = (y_bend(n) + off * dy_bend(n)) / 2
    n = n_bend
    x_20 = x_bend(n) - off * dx_bend(n);       y_20 = y_bend(n) - off * dy_bend(n)
    x_21 = x_bend(n);                          y_21 = y_bend(n)
    x_22 = x_bend(n) + off * dx_bend(n);       y_22 = y_bend(n) + off * dy_bend(n)

  else
    x_00 = end1%r(1) - dx1;                                y_00 = end1%r(2) - dy1
    x_01 = end1%r(1);                                      y_01 = end1%r(2)
    x_02 = end1%r(1) + dx1;                                y_02 = end1%r(2) + dy1
    x_10 = ((end1%r(1) + end2%r(1)) - (dx1 + dx2)) / 2;    y_10 = ((end1%r(2) + end2%r(2)) - (dy1 + dy2)) / 2
    x_12 = ((end1%r(1) + end2%r(1)) + (dx1 + dx2)) / 2;    y_12 = ((end1%r(2) + end2%r(2)) + (dy1 + dy2)) / 2
    x_20 = end2%r(1) - dx1;                                y_20 = end2%r(2) - dy1
    x_21 = end2%r(1);                                      y_21 = end2%r(2)
    x_22 = end2%r(1) + dx1;                                y_22 = end2%r(2) + dy1
  endif

  select case (shape)
  case ('diamond')
    call qp_draw_line (x_01, x_12, y_01, y_12, units = draw_units)
    call qp_draw_line (x_01, x_10, y_01, y_10, units = draw_units)
    call qp_draw_line (x_21, x_12, y_21, y_12, units = draw_units)
    call qp_draw_line (x_21, x_10, y_21, y_10, units = draw_units)
  case ('r_triangle')
    call qp_draw_line (x_00, x_21, y_00, y_21, units = draw_units)
    call qp_draw_line (x_02, x_21, y_02, y_21, units = draw_units)
    call qp_draw_line (x_00, x_02, y_00, y_02, units = draw_units)
  case ('l_triangle')
    call qp_draw_line (x_20, x_01, y_20, y_01, units = draw_units)
    call qp_draw_line (x_22, x_01, y_22, y_01, units = draw_units)
    call qp_draw_line (x_20, x_22, y_20, y_22, units = draw_units)
  case ('u_triangle')
    call qp_draw_line (x_00, x_12, y_00, y_12, units = draw_units)
    call qp_draw_line (x_20, x_12, y_20, y_12, units = draw_units)
    call qp_draw_line (x_00, x_20, y_00, y_20, units = draw_units)
  case ('d_triangle')
    call qp_draw_line (x_02, x_10, y_02, y_10, units = draw_units)
    call qp_draw_line (x_22, x_10, y_22, y_10, units = draw_units)
    call qp_draw_line (x_02, x_22, y_02, y_22, units = draw_units)
  end select
endif

! Draw a circle.

if (shape == 'circle') then
  call qp_draw_circle ((end1%r(1)+end2%r(1))/2, (end1%r(2)+end2%r(2))/2, off, units = draw_units)
endif

! Draw an X.

if (shape == 'x') then
  if (is_bend) then
    n = n_bend / 2
    x0  = x_bend(n)
    y0  = y_bend(n)
    dx1 = off * dx_bend(n)
    dy1 = off * dy_bend(n)
  else
    x0 = (end1%r(1) + end2%r(1)) / 2
    y0 = (end1%r(2) + end2%r(2)) / 2
  endif
  call qp_draw_line (x0 - dx1, x0 + dx1, y0 - dy1, y0 + dy1, units = draw_units)
  call qp_draw_line (x0 - dx1, x0 + dx1, y0 + dy1, y0 - dy1, units = draw_units)
endif

! Draw top and bottom of boxes and rbow_tie

if (shape == 'rbow_tie' .or. shape == 'box' .or. shape == 'xbox') then
  if (is_bend) then
    n = n_bend; n1 = min1_bend; n2 = max1_bend
    call qp_draw_polyline([x_bend(0)+off1*dx_e1, x_bend(n1:n2)+off1*dx_bend(n1:n2), x_bend(n)+off1*dx_e2], &
                          [y_bend(0)+off1*dy_e1, y_bend(n1:n2)+off1*dy_bend(n1:n2), y_bend(n)+off1*dy_e2], units = draw_units)
    n = n_bend; n1 = min2_bend; n2 = max2_bend
    call qp_draw_polyline([x_bend(0)-off2*dx_e1, x_bend(n1:n2)-off2*dx_bend(n1:n2), x_bend(n)-off2*dx_e2], &
                          [y_bend(0)-off2*dy_e1, y_bend(n1:n2)-off2*dy_bend(n1:n2), y_bend(n)-off2*dy_e2], units = draw_units)

  else
    call qp_draw_line (end1%r(1)+dx1, end2%r(1)+dx1, end1%r(2)+dy1, end2%r(2)+dy1, units = draw_units)
    call qp_draw_line (end1%r(1)-dx2, end2%r(1)-dx2, end1%r(2)-dy2, end2%r(2)-dy2, units = draw_units)
  endif
endif

! Draw sides of boxes

if (shape == 'bow_tie' .or. shape == 'box' .or. shape == 'xbox') then
  if (is_bend) then
    call qp_draw_line (x_bend(0)-off2*dx_e1, x_bend(0)+off1*dx_e1, &
                       y_bend(0)-off2*dy_e1, y_bend(0)+off1*dy_e1, units = draw_units)
    n = n_bend
    call qp_draw_line (x_bend(n)-off2*dx_e2, x_bend(n)+off1*dx_e2, &
                       y_bend(n)-off2*dy_e2, y_bend(n)+off1*dy_e2, units = draw_units)
  else
    call qp_draw_line (end1%r(1)+dx1, end1%r(1)-dx2, end1%r(2)+dy1, end1%r(2)-dy2, units = draw_units)
    call qp_draw_line (end2%r(1)+dx1, end2%r(1)-dx2, end2%r(2)+dy1, end2%r(2)-dy2, units = draw_units)
  endif
endif

! Draw X for xbox or bow_tie

if (shape == 'xbox' .or. shape == 'bow_tie' .or. shape == 'rbow_tie' .or. shape == 'x') then
  call qp_draw_line (end1%r(1)+dx1, end2%r(1)-dx2, end1%r(2)+dy1, end2%r(2)-dy2, units = draw_units)
  call qp_draw_line (end1%r(1)-dx2, end2%r(1)+dx1, end1%r(2)-dy1, end2%r(2)+dy2, units = draw_units)
endif

! Custom pattern

if (prefix == 'pattern') then
  do i = 1, size(s%plot_page%pattern)
    if (shape /= s%plot_page%pattern(i)%name) cycle
    pat => s%plot_page%pattern(i)
    do j = 1, size(pat%pt)
      r1 = end1%r(1:2) + pat%pt(j)%s * (end2%r(1:2) - end1%r(1:2)) - pat%pt(j)%y * [dx1, dy1]
      if (j > 1) call qp_draw_line (r0(1), r1(1), r0(2), r1(2), units = draw_units)
      r0 = r1
    enddo
  enddo
endif

! Draw the label.
! Place a bend's label to the outside of the bend.

if (label_name /= '') then
  if (ele%key /= sbend$ .or. ele%value(g$) == 0) then
    x_center = (end1%r(1) + end2%r(1)) / 2 
    y_center = (end1%r(2) + end2%r(2)) / 2 
    dx = -1.5 * dt_y / norm2([dt_x, dt_y])
    dy =  1.5 * dt_x / norm2([dt_x, dt_y])
  else
    n = n_bend / 2
    x_center = x_bend(n) 
    y_center = y_bend(n) 
    dx = -1.5 * dx_bend(n) / norm2([dx_bend(n), dy_bend(n)])
    dy = -1.5 * dy_bend(n) / norm2([dx_bend(n), dy_bend(n)])
  endif

  if (graph%floor_plan%flip_label_side) then
    dx = -dx
    dy = -dy
  endif

  ! The extra factors of 2 are to shift the branch cut away from +/- 90.
  ! This is done since many lattices have elements with theta at +/- 90.
  theta = modulo2 (2 + atan2(dy, dx) * 180 / pi, 90.0_rp) - 2
  if (dx*cos(pi*theta/180)+dy*sin(pi*theta/180) > 0) then
    justify = 'LC'
  else
    justify = 'RC'
  endif
  height = s%plot_page%text_height * s%plot_page%legend_text_scale * s%plot_page%floor_plan_text_scale
  call qp_draw_text (label_name, x_center+dx*abs(off2), y_center+dy*abs(off2), units = draw_units, &
                               height = height, justify = justify, ANGLE = theta)    
endif

call qp_restore_state

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
type (tao_lattice_branch_struct), pointer :: tao_branch
type (tao_ele_shape_struct), pointer :: ele_shape
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele, ele1, ele2
type (branch_struct), pointer :: branch, branch2
type (tao_data_array_struct), allocatable, target :: d_array(:)
type (tao_var_array_struct), allocatable, target :: v_array(:)
type (tao_logical_array_struct), allocatable :: logic_array(:)
type (tao_data_struct), pointer :: datum
type (tao_var_struct), pointer :: var

real(rp) x1, x2, y1, y2, y, s_pos, x0, y0, s_lat_min, s_lat_max
real(rp) lat_len, height, dx, dy, key_number_height, dummy, l2

integer i, j, k, n, kk, ix, ix1, isu, ixe, ib
integer ix_var, ixv

logical err, have_data

character(80) str
character(40) name
character(*), parameter :: r_name = 'tao_draw_lat_layout'

! Init

select case (plot%x_axis_type)
case ('index')
  call out_io (s_error$, r_name, '"index" x-axis type not valid with lat_layout. Must be "s"')
  call qp_draw_text ('"index" x-axis type not valid with lat_layout', 0.1_rp, 0.5_rp, '%BOX', color = 'red', justify = 'LC')
  return

case ('ele_index')
  call qp_draw_text ('"ele_index" x-axis type not valid with lat_layout', 0.1_rp, 0.5_rp, '%BOX', color = 'red', justify = 'LC')
  call out_io (s_error$, r_name, '"ele_index" x-axis type not valid with lat_layout. Must be "s"')
  return

case ('s')

case default
  call qp_draw_text ('Unknown x-axis type not valid with lat_layout', 0.1_rp, 0.5_rp, '%BOX', color = 'red', justify = 'LC')
  call out_io (s_warn$, r_name, 'Unknown x_axis_type. Must be "s".')
  return
end select

isu = tao_universe_index(graph%ix_universe, .true.)
lat => s%u(isu)%model%lat
ib = tao_branch_index(graph%ix_branch)
branch => lat%branch(ib)
tao_branch => s%u(isu)%model%tao_branch(ib)

lat_len = branch%param%total_length

! Setup the placement of the graph on the plot page.

call qp_set_layout (x_axis = graph%x, y_axis = graph%y, x2_mirrors_x = .false., &
                    box = graph%box, margin = graph%margin)

s_lat_min = branch%ele(0)%s
s_lat_max = branch%ele(branch%n_ele_track)%s

if (branch%param%geometry == open$) then
  x1 = max(s_lat_min, graph%x%min)
  x2 = min(s_lat_max, graph%x%max)
else
  x1 = graph%x%min
  x2 = min(graph%x%max, graph%x%min + lat_len)
endif

if (x2 <= x1) then
  call qp_draw_text ('Graph: ' // trim(plot%name) // '.' // trim(graph%name) // '  ' // graph%title, &
                                                       0.5_rp, 0.6_rp, '%BOX', color = 'red', justify = 'CC')
  call qp_draw_text ('Horzontal range does not overlap lattice', 0.5_rp, 0.4_rp, '%BOX', color = 'red', justify = 'CC')
  return
endif

call qp_draw_line (x1, x2, 0.0_rp, 0.0_rp)

! loop over all elements in the branch. Only draw those element that
! are within bounds.

do i = 1, branch%n_ele_track
  ele => branch%ele(i)
  if (ele%slave_status == super_slave$) cycle
  call draw_ele_for_lat_layout (ele, graph)
enddo

! Loop over all control elements.

do i = lat%n_ele_track+1, lat%n_ele_max
  ele => lat%ele(i)
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
    do ixv = 1, size(var%slave)
      if (var%slave(ixv)%ix_uni /= isu .and. graph%ix_universe /= -2) cycle
      ixe = var%slave(ixv)%ix_ele
      if (var%ele_name == 'PARTICLE_START') ixe = 0
      ele => pointer_to_ele(lat, ixe, var%slave(ixv)%ix_branch)
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

real(rp) y1, y2, x_lab
integer section_id, ix_shape_min, ix_shape

character(40) label_name

! Draw element shape...

ix_shape_min = 1
do
  call tao_ele_shape_info (graph%ix_universe, ele, s%plot_page%lat_layout%ele_shape, ele_shape, label_name, y1, y2, ix_shape_min)

  if (.not. associated(ele_shape)) return
  if (.not. ele_shape%draw) return

  call find_element_ends (ele, ele1, ele2)
  if (.not. associated(ele1)) return
  if (ele1%ix_branch /= tao_branch_index(graph%ix_branch)) return
  x1 = ele1%s
  x2 = ele2%s
  ! If out of range then try to shift by lat_len to get in range
  if (branch%param%geometry == closed$ .and. x1 > graph%x%max .and. x2 > graph%x%max) then
    x1 = x1 - lat_len
    x2 = x2 - lat_len
  endif
  if (branch%param%geometry == closed$ .and. x1 < graph%x%min .and. x2 < graph%x%min) then
    x1 = x1 + lat_len
    x2 = x2 + lat_len
  endif
    
  if (x1 > graph%x%max .and. x2 > graph%x%max) return
  if (x1 < graph%x%min .and. x2 < graph%x%min) return

  ! Here if element is to be drawn...
  ! y1 and y2 are the offsets for the lines below and above the center line.

  y1 = max(graph%y%min, min(y1 * s%plot_page%lat_layout_shape_scale, graph%y%max))
  y2 = max(graph%y%min, min(-y2 * s%plot_page%lat_layout_shape_scale, graph%y%max))

  ! Does this element wrap around?

  if (x2 < x1 .and. ele%value(l$) > 0) then
    if (x1 > graph%x%min .and. x1 < graph%x%max) then
      call draw_shape_for_lat_layout (label_name, x1, x1 + ele%value(l$), y1, y2, min(graph%x%max, x1+ele%value(l$)/2), ele_shape)
    endif
    if (x2 > graph%x%min .and. x2 < graph%x%max) then
      call draw_shape_for_lat_layout (label_name, x2-ele%value(l$), x2, y1, y2, max(graph%x%min, x2-ele%value(l$)/2), ele_shape)
    endif

  else
    x_lab = min(max((x1+x2)/2, graph%x%min), graph%x%max)
    call draw_shape_for_lat_layout (label_name, x1, x2, y1, y2, x_lab, ele_shape)
  endif

  if (.not. ele_shape%multi) return
enddo

end subroutine draw_ele_for_lat_layout

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! contains

subroutine draw_shape_for_lat_layout (label_name, x1, x2, y1, y2, s_pos, ele_shape)

type (tao_ele_shape_struct) ele_shape
type (tao_shape_pattern_struct), pointer :: pat

real(rp) :: s_pos, y_off, r_dum = 0, x1, x2, y1, y2, r0(2), r1(2)
integer ix, iwidth

character(*) label_name
character(16) color, prefix
character(40) shape

!

ix = index(ele_shape%shape, ':')
if (ix == 0) then
  prefix = ''
  shape = ele_shape%shape
else
  prefix = ele_shape%shape(:ix-1)
  shape = ele_shape%shape(ix+1:)
endif

color = ele_shape%color
iwidth = ele_shape%line_width

! Draw the shape

if (shape == 'diamond') then
  call qp_draw_line (x1, (x1+x2)/2, 0.0_rp, y1, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x1, (x1+x2)/2, 0.0_rp, y2, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x2, (x1+x2)/2, 0.0_rp, y1, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x2, (x1+x2)/2, 0.0_rp, y2, width = iwidth, color = color, clip = .true.)
endif

if (shape == 'circle') then
  call qp_convert_point_abs ((x1+x2)/2, (y1+y2)/2, 'DATA', x0, y0, 'POINTS')
  call qp_convert_point_rel (r_dum, y1, 'DATA', r_dum, y1, 'POINTS')
  call qp_draw_circle (x0, y0, abs(y1), units = 'POINTS', width = iwidth, color = color)
endif

if (shape == 'x') then
  call qp_convert_point_abs ((x1+x2)/2, (y1+y2)/2, 'DATA', x0, y0, 'POINTS')
  call qp_convert_point_rel (x1, y1, 'DATA', x1, y1, 'POINTS')
  call qp_convert_point_rel (x2, y2, 'DATA', x2, y2, 'POINTS')
  call qp_draw_line (x0-y1, x0+y1, y0-y1, y0+y1, units = 'POINTS', width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x0-y1, x0+y1, y0+y1, y0-y1, units = 'POINTS', width = iwidth, color = color, clip = .true.)
endif

if (shape == 'bow_tie') then
  call qp_draw_line (x1, x1, y1, y2, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x2, x2, y1, y2, width = iwidth, color = color, clip = .true.)
endif

if (shape == 'rbow_tie') then
  call qp_draw_line (x1, x2, y1, y1, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x1, x2, y2, y2, width = iwidth, color = color, clip = .true.)
endif

if (shape == 'box' .or. shape == 'xbox') then
  call qp_draw_rectangle (x1, x2, y1, y2, width = iwidth, color = color, clip = .true.)
endif

if (shape == 'r_triangle') then
  call qp_draw_line (x1, x2, y1, 0.0_rp, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x1, x2, y2, 0.0_rp, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x1, x1, y1, y2, width = iwidth, color = color, clip = .true.)
endif

if (shape == 'l_triangle') then
  call qp_draw_line (x1, x2, 0.0_rp, y1, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x1, x2, 0.0_rp, y2, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x2, x2, y1, y2, width = iwidth, color = color, clip = .true.)
endif

if (shape == 'u_triangle') then
  call qp_draw_line (x1, x2, y2, y2, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x1, (x1+x2)/2, y2, y1, width = iwidth, color = color, clip = .true.)
  call qp_draw_line ((x1+x2)/2, x2, y1, y2, width = iwidth, color = color, clip = .true.)
endif

if (shape == 'd_triangle') then
  call qp_draw_line (x1, x2, y1, y1, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x1, (x1+x2)/2, y1, y2, width = iwidth, color = color, clip = .true.)
  call qp_draw_line ((x1+x2)/2, x2, y2, y1, width = iwidth, color = color, clip = .true.)
endif

! Draw X for XBOX or BOW_TIE

if (shape == 'xbox' .or. shape == 'bow_tie' .or. shape == 'rbow_tie') then
  call qp_draw_line (x1, x2, y2, y1, width = iwidth, color = color, clip = .true.)
  call qp_draw_line (x1, x2, y1, y2, width = iwidth, color = color, clip = .true.)
endif

! Custom pattern

if (prefix == 'pattern') then
  do i = 1, size(s%plot_page%pattern)
    if (shape /= s%plot_page%pattern(i)%name) cycle
    pat => s%plot_page%pattern(i)
    do j = 1, size(pat%pt)
      r1 = [x1, 0.0_rp] + [pat%pt(j)%s, pat%pt(j)%y] * [x2-x1, y1]
      if (j > 1) call qp_draw_line (r0(1), r1(1), r0(2), r1(2), width = iwidth, color = color, clip = .true.)
      r0 = r1
    enddo
  enddo
endif

! Put on a label

if (s%global%label_lattice_elements .and. label_name /= '') then

  call qp_from_inch_rel (0.0_rp, graph%y%label_offset, r_dum, y_off, 'DATA')

  if (s_pos > graph%x%max .and. s_pos-lat_len > graph%x%min) s_pos = s_pos - lat_len
  height = s%plot_page%text_height * s%plot_page%legend_text_scale * s%plot_page%lat_layout_text_scale
  call qp_draw_text (label_name, s_pos, graph%y%min-y_off, height = height, justify = 'LC', ANGLE = 90.0_rp)

endif

end subroutine draw_shape_for_lat_layout

end subroutine tao_draw_lat_layout

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

! NOTE: THIS ROUTINE IS NOT CURRENTLY ACITVE (NOT CALLED BY ANY OTHER ROUTINE).

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

use wall3d_mod, only: calc_wall_radius

implicit none

type (tao_plot_struct) plot
type (ele_struct), pointer :: ele
type (tao_graph_struct), target :: graph
type (lat_struct), pointer :: lat
type (tao_lattice_branch_struct), pointer :: tao_branch
type (branch_struct), pointer :: branch, branch2

real(rp) lat_len, dummy

integer i, ib, isu

!

isu = tao_universe_index(graph%ix_universe)
lat => s%u(isu)%model%lat
ib = graph%ix_branch
branch => lat%branch(ib)
tao_branch => s%u(isu)%model%tao_branch(ib)

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
integer section_id, n_curve_pts
character(16) color

! Draw beam chamber wall. 

color = 'black'
if (allocated (graph%curve)) color = graph%curve(1)%line%color

n_curve_pts = s%plot_page%n_curve_pts
if (plot%n_curve_pts > 0) n_curve_pts = plot%n_curve_pts

if (associated(ele%wall3d)) then
  call calc_wall_radius (ele%wall3d(1)%section(1)%v,  1.0_rp, 0.0_rp,  y1_plus, dummy)
  call calc_wall_radius (ele%wall3d(1)%section(1)%v, -1.0_rp, 0.0_rp,  y1_minus, dummy)
  x1 = ele%s_start + ele%wall3d(1)%section(1)%s

  ! Skip points so close to the last point that the points have negligible spacing

  do section_id = 2, size(ele%wall3d(1)%section)
    x2 = ele%s_start + ele%wall3d(1)%section(section_id)%s
    if (section_id /= size(ele%wall3d(1)%section) .and. &
            (x2 - x1) < (graph%x%max - graph%x%min) / n_curve_pts) cycle
    call calc_wall_radius (ele%wall3d(1)%section(section_id)%v,  1.0_rp, 0.0_rp,  y2_plus, dummy)
    call calc_wall_radius (ele%wall3d(1)%section(section_id)%v, -1.0_rp, 0.0_rp,  y2_minus, dummy)
    !scale wall
    call qp_draw_line (x1, x2, y1_plus, y2_plus, color = color)
    call qp_draw_line (x1, x2, -y1_minus, -y2_minus, color = color)
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

call tao_draw_graph_axes (plot, graph)

! loop over all the curves of the graph and draw them

do k = 1, size(graph%curve)
  if (.not. graph%curve(k)%valid) cycle
  call tao_draw_curve_data (plot, graph, graph%curve(k), have_data)
enddo

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
!   have_data -- Logical: Is there any data to plot? Set True if so.
!                   But never reset to False. 
!-

subroutine tao_draw_curve_data (plot, graph, curve, have_data)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct) :: curve

real(rp) :: dz, dx, dummy
integer i, color
logical have_data
character(16) num_str

!

if (.not. curve%valid) return
if (curve%use_y2) call qp_use_axis (y = 'Y2')

call qp_set_symbol (curve%symbol)

if (curve%draw_symbols .and. allocated(curve%x_symb)) then
  if (size(curve%x_symb) > 0) have_data = .true.
  
  if (graph%symbol_size_scale > 0) then  
    do i = 1, size(curve%x_symb), max(1, curve%symbol_every)
      call qp_draw_symbol (curve%x_symb(i), curve%y_symb(i), height = curve%symb_size(i), clip = graph%clip)
    enddo 
    
  ! Color by z  
  elseif (curve%z_color%is_on) then  
    if (curve%z_color%autoscale) then
      curve%z_color%min = minval(curve%z_symb(:))
      curve%z_color%max = maxval(curve%z_symb(:))
    endif

    dz = (curve%z_color%max - curve%z_color%min)

    do i = 1, size(curve%x_symb), max(1, curve%symbol_every)
        call qp_draw_symbol (curve%x_symb(i), curve%y_symb(i), clip = graph%clip, color = z_color())
    enddo    
    
  else
    call qp_draw_symbols (curve%x_symb, curve%y_symb, symbol_every = curve%symbol_every, clip = graph%clip)
  endif
endif

if (curve%draw_symbol_index .and. allocated(curve%ix_symb)) then
  if (size(curve%ix_symb) > 0) have_data = .true.
  do i = 1, size(curve%x_symb)
    if (mod(i-1, curve%symbol_every) /= 0) cycle
    if (graph%clip) then
      if (curve%x_symb(i) < graph%x%min .or. curve%x_symb(i) > graph%x%max)  cycle
      if (curve%y_symb(i) < graph%y%min .or. curve%y_symb(i) > graph%y%max) cycle
    endif
    write (num_str, '(i0)') curve%ix_symb(i)
    call qp_draw_text (num_str, curve%x_symb(i), curve%y_symb(i))
  enddo
endif

if (curve%draw_error_bars .and. allocated(curve%err_symb)) then
  if (size(curve%err_symb) > 0) have_data = .true.
  call qp_from_inch_rel(0.02_rp, 0.0_rp, dx, dummy)  ! Bar width at top & bottom is 0.02"
  do i = 1, size(curve%x_symb)
    if (mod(i-1, curve%symbol_every) /= 0) cycle
    if (graph%clip) then
      if (curve%x_symb(i) < graph%x%min .or. curve%x_symb(i) > graph%x%max)  cycle
      if (curve%y_symb(i) < graph%y%min .or. curve%y_symb(i) > graph%y%max) cycle
    endif
    call qp_draw_line (curve%x_symb(i), curve%x_symb(i), curve%y_symb(i) - curve%err_symb(i), &
                                                                   curve%y_symb(i) + curve%err_symb(i)) 
    call qp_draw_line (curve%x_symb(i)-dx, curve%x_symb(i)+dx, curve%y_symb(i) - curve%err_symb(i), &
                                                                   curve%y_symb(i) - curve%err_symb(i)) 
    call qp_draw_line (curve%x_symb(i)-dx, curve%x_symb(i)+dx, curve%y_symb(i) + curve%err_symb(i), &
                                                                   curve%y_symb(i) + curve%err_symb(i))
  enddo
endif

if (curve%draw_line .and. allocated(curve%x_line)) then
  if (size(curve%x_line) > 0) have_data = .true.
  call qp_set_line ('PLOT', curve%line)
  call qp_draw_polyline (curve%x_line, curve%y_line, clip = graph%clip, style = 'PLOT')
endif

call qp_use_axis (y = 'Y')  ! reset

!-----------------------------------------------
contains

function z_color() result (color)
character(16) color
real(rp) :: z
if (dz==0) then
  color = 'black'
else
  z = (curve%z_symb(i) - curve%z_color%min)/dz
  z = 1-z ! Make red -> purple with PLPlot
  color = qp_enum_to_string(qp_continuous_color(z), 'color')
endif
end function z_color

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
type (qp_line_struct) l

integer i
logical have_data
character(16) num_str

! Draw

if (.not. allocated(curve%x_line)) then
  have_data = .false.
  return
endif

l = curve%line
call qp_draw_histogram (curve%x_line, curve%y_line, line_color = l%color, fill_color = l%color, fill_pattern = l%pattern) 
have_data = .true.

end subroutine tao_draw_histogram_data

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_draw_graph_axes (plot, graph)
!
! Routine to draw a just the graph part of a data graph.
! The calling routine takes care of drawing any curves.
!
! Input:
!   plot  -- Tao_plot_struct: Plot containing the graph.
!   graph -- Tao_graph_struct: Graph to plot.
!-

subroutine tao_draw_graph_axes (plot, graph)

implicit none

type (tao_plot_struct) plot
type (tao_graph_struct), target :: graph
type (tao_curve_struct), pointer :: curve
type (qp_line_struct), allocatable :: line(:)
type (qp_symbol_struct), allocatable :: symbol(:)

integer i, j, k, nc, iu
real(rp) x, y, x1
character(100), allocatable :: text(:)

! Set scales, margens, etc

call qp_set_layout (box = graph%box, margin = graph%margin)
call qp_set_layout (x_axis = graph%x, x2_mirrors_x = .true.)
call qp_set_layout (y_axis = graph%y, y2_axis = graph%y2, y2_mirrors_y = graph%y2_mirrors_y)

if (graph%draw_title) then
  if (s%plot_page%draw_graph_title_suffix) then
    call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
  else
    call qp_set_graph (title = graph%title)
  endif
else
  call qp_set_graph (title = '')
endif

if (graph%draw_axes) call qp_draw_axes (draw_grid = graph%draw_grid)

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

if (graph%limited .and. graph%clip .and. s%global%draw_curve_off_scale_warn) &
  call qp_draw_text ('**Curve Off Scale**', -0.30_rp, -0.15_rp, '%/GRAPH/RT', color = 'red') 

! Draw the text legend if there is one

if (any(graph%text_legend_out /= ' ')) call qp_draw_text_legend (graph%text_legend_out, &
       graph%text_legend_origin%x, graph%text_legend_origin%y, graph%text_legend_origin%units)

! Draw the curve legend if needed

if (.not. allocated(graph%curve)) return
nc = size(graph%curve)
allocate (text(nc), symbol(nc), line(nc))

text = ''
symbol%type = ''
line%width = -1

do i = 1, nc
  ! Text
  curve => graph%curve(i)
  if (.not. curve%valid) then
    if (curve%why_invalid /= 'IGNORE') text(i) = curve%why_invalid
    cycle
  endif

  text(i) = curve%legend_text
  if (text(i) == '') text(i) = curve%data_type

  if (size(s%u) > 1 .and. .not. all(graph%curve%ix_universe == graph%curve(1)%ix_universe)) then
    iu = curve%ix_universe
    if (iu == -1) iu = graph%ix_universe
    if (iu == -1) iu = s%global%default_universe
    text(i) = int_str(iu) // '@' // text(i)
  endif

  if (.not. all(graph%curve%component == graph%curve(1)%component)) then
    text(i) = trim(text(i)) // ' ' // trim(curve%component)
  endif

  ! Symbol to display
  symbol(i) = curve%symbol
  line(i) = curve%line
  if (size(curve%x_line) == 0 .or. .not. curve%draw_line) line(i)%width = -1 ! Do not draw
  if (.not. curve%draw_symbols) symbol(i)%type = ''  ! Do not draw
enddo

if (graph%draw_curve_legend .and. nc > 1) then
  call qp_draw_curve_legend (graph%curve_legend_origin%x, graph%curve_legend_origin%y, &
            graph%curve_legend_origin%units, line, s%plot_page%curve_legend_line_len, &
            symbol, text, s%plot_page%curve_legend_text_offset)
endif


! Draw any curve info messages

if (graph%draw_curve_legend) then
  j = 0
  do i = 1, nc
    curve => graph%curve(i)
    if (curve%message_text == '') cycle
    j = j + 1
    text(j) = curve%message_text
  enddo

  if (j > 1) then
    call qp_draw_text_legend (text(1:j), 0.50_rp, 0.95_rp, '%GRAPH/LB')
  endif
endif

!

deallocate (text, symbol, line)

end subroutine tao_draw_graph_axes

end module
