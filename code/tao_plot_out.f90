!+
! Subroutine tao_plot_out (s)
!
! Subroutine to draw the plots on the plot window.
!
! Input:
!   s -- Tao_super_universe_struct:
!-

subroutine tao_plot_out (s)

use tao_mod
use quick_plot
use tao_single_mod

implicit none

type (tao_super_universe_struct), target :: s
type (tao_plot_struct), pointer :: plot
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve
type (qp_rect_struct) border1, border2

real(rp) region(4), dx, dy
integer i, j, k, graph_box(4)
character(16) :: r_name = 'tao_plot_out'

! inits

if (.not. s%global%plot_on) return
call qp_clear_page

! loop over all plots

do i = 1, size(s%plot_page%plot)

  plot => s%plot_page%plot(i)
  if (.not. plot%visible) cycle

! set the s%plot_page border for this particular region

  region = s%plot_page%plot(i)%region%location
  border1%units = '%PAGE'
  call qp_convert_rectangle_rel (s%plot_page%border, border1)
  dx = 1 - (border1%x2 - border1%x1)
  dy = 1 - (border1%y2 - border1%y1)
  border2%x1 = border1%x1 + dx * region(1)
  border2%x2 = border1%x2 + dx * (1 - region(2))
  border2%y1 = border1%y1 + dy * region(3)
  border2%y2 = border1%y2 + dy * (1 - region(4))
  border2%units = '%PAGE'
  call qp_set_layout (page_border = border2)

  select case (plot%type)
  case ('data')
    call plot_data 
  case ('lat_layout')
    call plot_lat_layout
  case ('key_table')
    call plot_key_table
  case default
    call out_io (s_fatal$, r_name, 'UNKNOWN PLOT TYPE: ' // plot%type)
  end select

enddo

!--------------------------------------------------------------------------
contains

subroutine plot_data

! loop over all the graphs of the plot and draw the axes.

do j = 1, size(plot%graph)
  graph => plot%graph(j)
  graph_box = (/ graph%this_box, plot%box_layout /)
  call qp_set_layout (box = graph_box, margin = graph%margin)
  call qp_set_layout (x_axis = plot%x, y_axis = graph%y, y2_axis = graph%y2)
  call qp_set_graph (title = trim(graph%title) // ' ' // graph%title_suffix)
  call qp_draw_axes

! loop over all the curves of the graph and draw them

  do k = 1, size(graph%curve)
    curve => graph%curve(k)
    call qp_set_symbol (curve%symbol)
    call qp_set_line ('PLOT', curve%line) 
    call qp_draw_data (curve%x_symb, curve%y_symb, curve%draw_line, &
                                               curve%symbol_every, graph%clip)
    call qp_draw_data (curve%x_line, curve%y_line, curve%draw_line, &
                                                           0, graph%clip)
  enddo
enddo

end subroutine

!--------------------------------------------------------------------------
! contains

subroutine plot_key_table

integer i, j, k, j_ele, j_att, ix_var
real(rp) :: y, norm
real(rp) :: dy_key = 12
character(3) str3
character(80) str, fmt

!

j_ele = 4
j_att = 5

do k = s%global%ix_key_bank+1, s%global%ix_key_bank+10
  ix_var = s%key(k)%ix_var
  if (ix_var == 0) cycle
  j_ele = max(j_ele, len_trim(s%var(ix_var)%ele_name))
  j_att = max(j_att, len_trim(s%var(ix_var)%attrib_name))
enddo

write (str, '(3x, a, <j_ele-2>x, a, <j_att-2>x, a)') 'ix  Name',  &
           'Attrib',  'Value   Value0    Delta  Ring#'
y = 10 * dy_key + 5.0
call qp_draw_text (str, 5.0_rp, y, 'POINTS/PAGE', &
                             height = dy_key-1.0_rp, uniform_spacing = .true.)
  
write (fmt, '(a, i2.2, a, i2.2, a)') &
                   '(a, i2, 2x, a', j_ele, ', 2x, a', j_att, ', 3f9.4, 2x, a)'

do i = 1, 10
  k = i + s%global%ix_key_bank
  ix_var = s%key(k)%ix_var
  j = mod(i, 10)
  str3 = '  '
  if (i == 1) write (str3, '(i2, a)') s%global%ix_key_bank/10, ':'

  if (s%key(k)%ix_var == 0) then
    write (str, '(a, i2)') str3, j
  else
    norm = s%key(k)%normalizer
    write (str, fmt) str3, j, s%var(ix_var)%ele_name, &
      s%var(ix_var)%attrib_name, s%var(ix_var)%model_value, &
      s%key(k)%val0/norm, s%key(k)%delta/norm, &
      trim(tao_var_uni_string(s%var(ix_var)))
  endif

  y = (10-i) * dy_key + 5
  call qp_draw_text (str, 5.0_rp, y, 'POINTS/PAGE', &
                          height = dy_key-1.0_rp, uniform_spacing = .true.)
enddo

end subroutine

!--------------------------------------------------------------------------
! contains

subroutine plot_lat_layout

type (ring_struct), pointer :: lat
type (ele_struct), pointer :: ele

real(rp) x1, x2, y, s_pos

integer i, j, k, kk, ix, ix1, ix2, isu
integer n_plot_tot, n_plot, icol, ix_var, ixv

character(80) str

!

n_plot_tot = count (s%u(:)%draw_lat_layout)
n_plot = n_plot_tot

do isu = 1, size(s%u)

  if (.not. s%u(isu)%draw_lat_layout) cycle
  lat => s%u(isu)%model

  call qp_set_box (1, n_plot, 1, n_plot_tot)
  n_plot = n_plot - 1
  call qp_set_axis ('Y', -70.0_rp, 30.0_rp, 1, 0)
  
  do i = 1, lat%n_ele_max

    ele => lat%ele_(i)
    call find_element_ends (lat, i, ix1, ix2)
    x1 = lat%ele_(ix1)%s
    x2 = lat%ele_(ix2)%s

    if (x1 > plot%x%max) cycle
    if (x2 < plot%x%min) cycle

    j = ele%ix_pointer
    if (j < 1) then
      call qp_draw_line (x1, x2, 0.0_rp, 0.0_rp)
      cycle
    endif

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

    if (s%global%label_lattice_elements .and. s%plot_page%ele_shape(j)%plot_name) &
                call qp_draw_text (ele%name, ele%s-ele%value(l$)/2, &
                -25.0_rp, justify = 'CT', height = 10.0_rp)

    call qp_draw_line (x1, x2, 0.0_rp, 0.0_rp)

  enddo

!

  if (s%global%label_keys) then
    do k = s%global%ix_key_bank+1, s%global%ix_key_bank+10
      ix_var = s%key(k)%ix_var
      do ixv = 1, size(s%var(ix_var)%this)
        if (s%var(ix_var)%this(ixv)%ix_uni /= isu) cycle
        ix = s%var(ix_var)%this(ixv)%ix_ele
        kk = mod(k - s%global%ix_key_bank, 10)
        write (str, '(i1)') kk
        if (ix > lat%n_ele_use) then
          do j = lat%ele_(ix)%ix1_slave, lat%ele_(ix)%ix2_slave
            ix1 = lat%control_(j)%ix_slave
            s_pos = lat%ele_(ix1)%s - lat%ele_(ix1)%value(l$)/2
            call qp_draw_text (trim(str), s_pos, -50.0_rp, &
                              justify = 'CT', height = 10.0_rp)  
          enddo
        else
          s_pos = lat%ele_(ix)%s - lat%ele_(ix)%value(l$)/2
          call qp_draw_text (trim(str), s_pos, -50.0_rp, &
                                justify = 'CT', height = 10.0_rp)  
        endif
      enddo
    enddo
  endif

enddo


end subroutine

end subroutine
