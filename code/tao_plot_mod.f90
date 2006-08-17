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

real(rp) location(4), dx, dy

integer i, j, k

character(16) :: r_name = 'tao_plot_out'
character(3) view_str

logical found

! inits

if (.not. s%global%plot_on) return
call tao_create_plot_window () ! This routine knows not to create multiple windows.

call qp_clear_page

! print the title 
do i = 1, size(s%plot_page%title)
  if (s%plot_page%title(i)%draw_it)                                         &
    call qp_draw_text (s%plot_page%title(i)%string, s%plot_page%title(i)%x, &
                     s%plot_page%title(i)%y, s%plot_page%title(i)%units,    &
                     s%plot_page%title(i)%justify)
enddo

! Draw view universe

write (view_str, '(i3)') s%global%u_view
call qp_draw_text ('View Universe:' // view_str, &
              -2.0_rp, -2.0_rp, 'POINTS/PAGE/RT', 'RT')

! loop over all plots

do i = 1, size(s%plot_page%region)

  plot => s%plot_page%region(i)%plot
  if (.not. s%plot_page%region(i)%visible) cycle

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

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (tao_var_struct), pointer :: var
type (tao_keyboard_struct), pointer :: key

integer i, j, k, n, m, p, j_ele, j_att, ix_var
real(rp) :: y_here, norm, v, x1, x2, y1, y2
real(rp) :: dy_key, text_scale
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

call qp_get_qp_parameters (text_scale = text_scale)
dy_key = 12 * text_scale

do k = s%global%ix_key_bank+1, s%global%ix_key_bank+10
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
                             height = dy_key-1.0_rp, uniform_spacing = .true.)
  
write (fmt, '(a, i2.2, a, i2.2, a)') &
        '(i2, 2x, a', j_ele, ', 2x, a', j_att, ', 3a12, 2x, a, 3x, l)'

write (str, '(i2, a)') s%global%ix_key_bank/10, ':'
y_here = y_here - dy_key
call qp_draw_text (str, 5.0_rp, y_here, 'POINTS/GRAPH', &
                          height = dy_key-1.0_rp, uniform_spacing = .true.)

do i = 1, 10

  k = i + s%global%ix_key_bank
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
                          height = dy_key-1.0_rp, uniform_spacing = .true.)
  y_here = y_here - dy_key
enddo

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_plot_floor_plan (plot, graph)

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (ring_struct), pointer :: lat
type (ele_struct), pointer :: ele
type (floor_position_struct) screen1, screen2

integer i, j, icol, ix1, ix2, isu

real(rp) r, y

character(80) str
character(20) :: r_name = 'tao_plot_lat_layout'
character(16) shape

! Each graph is a separate lattice layout (presumably for different universes). 
! setup the placement of the graph on the plot page.

call qp_set_layout (x_axis = plot%x, y_axis = graph%y, box = graph%box, margin = graph%margin)
  
isu = graph%ix_universe
! if garph%ix_universe .eq. 0 then graph currently viewed universe
if (isu .eq. 0) then
  lat => s%u(s%global%u_view)%model%lat
else
  lat => s%u(isu)%model%lat
endif
  
! loop over all elements in the lattice. 

do i = 1, lat%n_ele_max

  ele => lat%ele_(i)
  if (ele%control_type == multipass_lord$) cycle
  if (ele%control_type == super_slave$) cycle

  call find_element_ends (lat, i, ix1, ix2)
  call floor_to_screen_coords (lat%ele_(ix1)%floor, screen1)
  call floor_to_screen_coords (lat%ele_(ix2)%floor, screen2)

  ! Only draw those element that are within bounds.

  if (min(screen1%x, screen2%x) > plot%x%max) cycle
  if (max(screen1%x, screen2%x) < plot%x%min) cycle

  if (min(screen1%y, screen2%y) > graph%y%max) cycle
  if (max(screen1%y, screen2%y) < graph%y%min) cycle

  ! Only those elements with ele%ix_pointer > 0 are to be drawn.
  ! All others are drawn with a line or are

  j = ele%ix_pointer
  if (j < 1) then
    cycle
  endif

  ! Here if element is to be drawn...

  call qp_translate_to_color_index (s%plot_page%ele_shape(j)%color, icol)

  y = s%plot_page%ele_shape(j)%dy_pix/2

  shape = s%plot_page%ele_shape(j)%shape

  if (shape == 'VAR_BOX' .or. shape == 'ASYM_VAR_BOX') then
    select case (ele%key)
    case (quadrupole$)
      r = ele%value(k1$)
    case (sextupole$)
      r = ele%value(k2$)
    case (octupole$)
      r = ele%value(k3$)
    case (solenoid$)
      r = ele%value(ks$)
    end select
  end if



  select case (s%plot_page%ele_shape(j)%shape)

  case ('BOX')

  case ('VAR_BOX')

  case ('ASYM_VAR_BOX')

  case ('XBOX')

  case default
    print *, 'ERROR: UNKNOWN SHAPE: ', s%plot_page%ele_shape(j)%shape
    call err_exit
  end select

enddo

!--------------------------------------------------------------------------
contains

subroutine floor_to_screen_coords (floor, screen)

type (floor_position_struct) floor, screen

!

! Mapping from floor coords to screen coords is:
!   Floor   Screen 
!    z   ->  -x
!    x   ->  -y

screen%x = -floor%z
screen%y = -floor%x
screen%theta = pi - floor%theta

end subroutine

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_plot_lat_layout (plot, graph)

type (tao_plot_struct) :: plot
type (tao_graph_struct) :: graph
type (ring_struct), pointer :: lat
type (ele_struct), pointer :: ele

real(rp) r, x1, x2, y1, y2, y, s_pos, y_off, y_bottom, y_top, height

integer i, j, k, kk, ix, ix1, ix2, isu
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
  
isu = graph%ix_universe
! if garph%ix_universe .eq. 0 then graph currently viewed universe
if (isu == 0) then
  lat => s%u(s%global%u_view)%model%lat
else
  lat => s%u(isu)%model%lat
endif
  
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

  ele => lat%ele_(i)
  if (ele%control_type == multipass_lord$) cycle
  if (ele%control_type == super_slave$) cycle

  if (plot%x_axis_type .eq. 's') then
    call find_element_ends (lat, i, ix1, ix2)
    x1 = lat%ele_(ix1)%s
    x2 = lat%ele_(ix2)%s
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

  j = ele%ix_pointer
  if (j < 1) then
    call qp_draw_line (x1, x2, 0.0_rp, 0.0_rp)
    cycle
  endif

  ! Here if element is to be drawn...

  call qp_translate_to_color_index (s%plot_page%ele_shape(j)%color, icol)

  y = s%plot_page%ele_shape(j)%dy_pix/2

  shape = s%plot_page%ele_shape(j)%shape

  if (shape == 'VAR_BOX' .or. shape == 'ASYM_VAR_BOX') then
    select case (ele%key)
    case (quadrupole$)
      r = ele%value(k1$)
    case (sextupole$)
      r = ele%value(k2$)
    case (octupole$)
      r = ele%value(k3$)
    case (solenoid$)
      r = ele%value(ks$)
    end select
  end if



  select case (s%plot_page%ele_shape(j)%shape)

  case ('BOX')
    call qp_draw_rectangle (x1, x2,  -y, y, color = icol)

  case ('VAR_BOX')
    call qp_draw_rectangle (x1, x2,  -y*r, y*r, color = icol)

  case ('ASYM_VAR_BOX')
    call qp_draw_rectangle (x1, x2,  0.0_rp, y*r, color = icol)

  case ('XBOX')
    call qp_draw_rectangle (x1, x2,  -y, y, color = icol)
    call qp_draw_line (x1, x2,  y, -y, color = icol)
    call qp_draw_line (x1, x2, -y,  y, color = icol)

  case default
    print *, 'ERROR: UNKNOWN SHAPE: ', s%plot_page%ele_shape(j)%shape
    call err_exit
  end select

  
  if (s%global%label_lattice_elements .and. s%plot_page%ele_shape(j)%plot_name) then
    j_label = j_label + 1
    if (j_label == s%global%n_lat_layout_label_rows) j_label = 0
    y_off = y_bottom   ! + 12.0_rp * j_label 
    height = 0.8 * s%plot_page%text_height 
    call qp_draw_text (ele%name, ele%s-ele%value(l$)/2, y_off, &
                                 height = height, justify = 'CB', ANGLE = 90.0_rp)
  endif

  call qp_draw_line (x1, x2, 0.0_rp, 0.0_rp)

enddo

! This is for drawing the key numbers under the appropriate elements

if (s%global%label_keys) then
  do k = s%global%ix_key_bank+1, s%global%ix_key_bank+10
    if (k > ubound(s%key, 1)) cycle
    ix_var = s%key(k)%ix_var
    if (ix_var < 1) cycle
    do ixv = 1, size(s%var(ix_var)%this)
      if (s%var(ix_var)%this(ixv)%ix_uni /= isu) cycle
      ix = s%var(ix_var)%this(ixv)%ix_ele
      kk = mod(k - s%global%ix_key_bank, 10)
      write (str, '(i1)') kk
      if (ix > lat%n_ele_use) then
        do j = lat%ele_(ix)%ix1_slave, lat%ele_(ix)%ix2_slave
          ix1 = lat%control_(j)%ix_slave
          s_pos = lat%ele_(ix1)%s - lat%ele_(ix1)%value(l$)/2
          call qp_draw_text (trim(str), s_pos, y_top, &
                              justify = 'CT', height = 10.0_rp)  
        enddo
      else
        s_pos = lat%ele_(ix)%s - lat%ele_(ix)%value(l$)/2
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
