!+
! Subroutine tao_plot_out ()
!
! Subroutine to draw the plots on the plot window.
!
! Input:
!-

subroutine tao_plot_out ()

use tao_mod
use quick_plot
use tao_single_mod
use tao_plot_window_mod

implicit none

type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (qp_rect_struct) border1, border2

real(rp) location(4), dx, dy
integer i, j, k
character(16) :: r_name = 'tao_plot_out'

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

    select case (graph%type)
    case ('data')
      call plot_data 
    case ('lat_layout')
      call plot_lat_layout  
    case ('key_table')
      call plot_key_table
    case default
      call out_io (s_fatal$, r_name, 'UNKNOWN PLOT TYPE: ' // graph%type)
    end select
  enddo

enddo

!--------------------------------------------------------------------------
contains

subroutine plot_data

call qp_set_layout (box = graph%box, margin = graph%margin)
call qp_set_layout (x_axis = plot%x, y_axis = graph%y, y2_axis = graph%y2)
call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
call qp_draw_axes

if (any(graph%curve%limited) .and. graph%clip) &
  call qp_draw_text ('**Limited**', -0.18_rp, -0.15_rp, '%/GRAPH/RT', color = red$) 

! loop over all the curves of the graph and draw them

do k = 1, size(graph%curve)
  curve => graph%curve(k)
  call qp_set_symbol (curve%symbol)
  call qp_set_line ('PLOT', curve%line) 
  if (curve%draw_symbols) call qp_draw_data (curve%x_symb, curve%y_symb, &
                                      .false., curve%symbol_every, graph%clip)
  if (curve%draw_line) call qp_draw_data (curve%x_line, curve%y_line, &
                                               curve%draw_line, 0, graph%clip)
enddo

! draw the legend if there is one

if (any(graph%legend /= ' ')) call qp_draw_legend (graph%legend, &
          graph%legend_origin%x, graph%legend_origin%y, graph%legend_origin%units)

end subroutine

!--------------------------------------------------------------------------
! contains

subroutine plot_key_table

integer i, j, k, j_ele, j_att, ix_var
real(rp) :: y, norm
real(rp) :: dy_key = 12
character(80) str, fmt, str2

!

j_ele = 4
j_att = 5

do k = s%global%ix_key_bank+1, s%global%ix_key_bank+10
  if (k > ubound(s%key, 1)) cycle
  ix_var = s%key(k)%ix_var
  if (ix_var == 0) cycle
  j_ele = max(j_ele, len_trim(s%var(ix_var)%ele_name))
  j_att = max(j_att, len_trim(s%var(ix_var)%attrib_name))
enddo

write (fmt, '(a, i5, a, i2, a)') &
                        '(3x, a, ', j_ele-2, 'x, a, ', j_att-2, 'x, a)'
write (str, fmt) 'ix  Name', 'Attrib', &
                              'Value   Value0    Delta  Uni useit_opt'
y = 10 * dy_key + 5.0
call qp_draw_text (str, 5.0_rp, y, 'POINTS/PAGE', &
                             height = dy_key-1.0_rp, uniform_spacing = .true.)
  
write (fmt, '(a, i2.2, a, i2.2, a)') &
        '(i2, 2x, a', j_ele, ', 2x, a', j_att, ', 3f9.4, 2x, a, 3x, l)'

write (str, '(i2, a)') s%global%ix_key_bank/10, ':'
y = 9 * dy_key + 5
call qp_draw_text (str, 5.0_rp, y, 'POINTS/PAGE', &
                          height = dy_key-1.0_rp, uniform_spacing = .true.)

do i = 1, 10

  k = i + s%global%ix_key_bank
  if (k > ubound(s%key, 1)) cycle
  ix_var = s%key(k)%ix_var
  j = mod(i, 10)

  if (s%key(k)%ix_var == 0) then
    write (str, '(i2)') j
  else
    norm = s%key(k)%normalizer
    str2 = tao_var_uni_string(s%var(ix_var))
    write (str, fmt) j, s%var(ix_var)%ele_name, &
      s%var(ix_var)%attrib_name, s%var(ix_var)%model_value, &
      s%key(k)%val0/norm, s%key(k)%delta/norm, &
      trim(str2), s%var(ix_var)%useit_opt
  endif

  y = (10-i) * dy_key + 5
  call qp_draw_text (str, 25.0_rp, y, 'POINTS/PAGE', &
                          height = dy_key-1.0_rp, uniform_spacing = .true.)
enddo

end subroutine

!--------------------------------------------------------------------------
! contains

subroutine plot_lat_layout

type (ring_struct), pointer :: lat
type (ele_struct), pointer :: ele

real(rp) x1, x2, y1, y2, y, s_pos, y_off, y_bottom, y_top

integer i, j, k, kk, ix, ix1, ix2, isu
integer icol, ix_var, ixv, j_label

character(80) str

! Each graph is a separate lattice layout (presumably for different universes). 
! setup the placement of the graph on the plot page.

call qp_set_layout (x_axis = plot%x, margin = graph%margin)
isu = graph%ix_universe
! if garph%ix_universe .eq. 0 then graph currently viewed universe
if (isu .eq. 0) then
  lat => s%u(s%global%u_view)%model
else
  lat => s%u(isu)%model
endif
  
call qp_set_layout (box = graph%box, margin = graph%margin)
call qp_get_layout_attrib ('GRAPH', x1, x2, y1, y2, 'POINTS/GRAPH')
y_bottom = -(y2 + 12 * (s%global%n_lat_layout_label_rows - 1)) / 2
y_top = y_bottom + y2
call qp_set_axis ('Y', y_bottom, y_top, 1, 0)
  
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

  y =  s%plot_page%ele_shape(j)%dy_pix/2

  select case (s%plot_page%ele_shape(j)%shape)

  case ('BOX')
    call qp_draw_rectangle (x1, x2,  -y, y, color = icol)

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
    y_off = y_bottom + 12.0_rp * j_label 
    call qp_draw_text (ele%name, ele%s-ele%value(l$)/2, y_off, justify = 'CB')
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

end subroutine
