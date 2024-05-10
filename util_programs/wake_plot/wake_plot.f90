program wake_plot

use bmad
use wake_mod
use quick_plot

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: wake_ele
type (ele_pointer_struct), allocatable, target :: eles(:)
type (wake_sr_struct), pointer :: sr
type (wake_lr_struct), pointer :: lr
type (wake_lr_mode_struct), allocatable :: lr_mode(:)
type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p1, p2

real(rp) :: plot_min, plot_max, plot_limit, plot_size(2), text_scale, zero6(6) = 0, vec(6)
real(rp) xy_leading(2), xy_trailing(2), leading_charge, z
real(rp), allocatable :: time(:), wx(:), wy(:), wz(:), z_knot(:), wx_knot(:), wy_knot(:), wz_knot(:)

integer n, n_points, ix_wake, m_order, i, im, iz, ik, n_loc

logical make_plot, err, err_flag, ok, draw_knot_points

character(40) wake_ele_name, ans, x_axis_label, y_axis_label, who, plot_type
character(200) init_file, lat_file, postscript_file

namelist / params / plot_limit, make_plot, plot_size, who, plot_type, &
                    wake_ele_name, lat_file, n_points, postscript_file, &
                    draw_knot_points, text_scale, x_axis_label, y_axis_label, &
                    xy_leading, xy_trailing, ix_wake, m_order, leading_charge

! Read input

init_file = 'wake_plot.init'
if (command_argument_count() > 0) call get_command_argument(1, init_file)
print '(2a)', 'Init file: ', trim(init_file)

leading_charge = 1
n_points = 1001
make_plot = .true.
postscript_file = ''
draw_knot_points = .true.
text_scale = 1.2
ix_wake = 0
m_order = 0

open (1, file = init_file, status = 'old')
read (1, nml = params)
close (1)

allocate (time(n_points), wx(n_points), wy(n_points), wz(n_points))

! Init

call bmad_parser(lat_file, lat, err_flag = err)
if (err) stop

!

call lat_ele_locator (wake_ele_name, lat, eles, n_loc)

if (n_loc == 0) then
  print '(2a)', 'CANNOT FIND WAKE_ELE: ', trim(wake_ele_name)
  stop
endif

if (n_loc > 1) then
  print '(2a)', 'More than one wake_ele found with name: ', trim(wake_ele_name)
  print '(a)',  'Using first one...'
endif

wake_ele => eles(1)%ele
if (.not. associated(wake_ele%wake)) then
  print *, 'NO WAKES IN ELEMENT: ' // wake_ele%name
  stop
endif

!

call reallocate_bunch(bunch, 2)
p1 => bunch%particle(1)
p2 => bunch%particle(2)
call init_coord(p1, zero6, wake_ele, upstream_end$)
call init_coord(p2, zero6, wake_ele, upstream_end$)

p1%vec(1:3:2) = xy_leading
p2%vec(1:3:2) = xy_trailing
vec = p2%vec

if (plot_type == 'wake') then
  p1%charge = 1
else
  p1%charge = leading_charge
endif

! 

select case (who(1:2))

case ('lr')
  lr => wake_ele%wake%lr
  bmad_com%lr_wakes_on = .true.
  plot_max = plot_limit
  plot_min = 0

  if (size(lr%mode) > 1) then
    if (ix_wake == 0) then
      lr%mode = [lr%mode(ix_wake)]
    else
      n = count(lr%mode%m == m_order)
      if (n == 0) then
        print *, 'Number of LR wakes is zero with m_order = ' // int_str(m_order)
        stop
      endif
      allocate(lr_mode(n))
      n = 0
      do im = 1, size(lr%mode)
        if (lr%mode(im)%m /= m_order) cycle
        n = n + 1
        lr_mode(n) = lr%mode(im)
      enddo
      lr%mode = lr_mode
    endif
  endif

  p2%vec(5) = -1
  call order_particles_in_z(bunch)

  do i = 1, n_points
    z = plot_min + (plot_max - plot_min) * (i - 1.0_rp) / (n_points - 1.0_rp)
    p2%vec = vec
    p2%vec(5) = z
    time(i) = z
    call zero_lr_wakes_in_lat(lat)
    call track1_lr_wake(bunch, wake_ele)
    wx(i) = p2%vec(2)
    wy(i) = p2%vec(4)
    wz(i) = p2%vec(6)
  enddo

!

case ('sr')
  sr => wake_ele%wake%sr
  bmad_com%sr_wakes_on = .true.
  plot_max = 0
  plot_min = -abs(plot_limit)

  select case (who)
  case ('sr-mode')
    allocate(sr%z(0))
  case ('sr-long')
    allocate(sr%z(0))
    allocate(sr%trans(0))
    if (ix_wake /= 0) sr%long = [sr%long(ix_wake)]
  case ('sr-trans')
    allocate(sr%z(0))
    allocate(sr%long(0))
    if (ix_wake /= 0) sr%trans = [sr%trans(ix_wake)]
  case ('sr-z')
    allocate(sr%trans(0))
    allocate(sr%long(0))
    if (ix_wake /= 0) sr%z = [sr%z(ix_wake)]
  end select

  p2%vec(5) = -1
  call order_particles_in_z(bunch)

  do i = 1, n_points
    z = plot_min + (plot_max - plot_min) * (i - 1.0_rp) / (n_points - 1.0_rp)
    p2%vec = vec
    p2%vec(5) = z
    time(i) = z
    call track1_sr_wake(bunch, wake_ele)
    wx(i) = p2%vec(2)
    wy(i) = p2%vec(4)
    wz(i) = p2%vec(6)
  enddo

  if (draw_knot_points) then
    n = 0
    allocate(z_knot(0), wx_knot(0), wy_knot(0), wz_knot(0))
    do iz = 1, size(sr%z)
      do ik = 1, size(sr%z(iz)%w)
        z = sr%z(iz)%w(ik)%x0
        if (any(z_knot == z)) cycle
        call track1_sr_wake(bunch, wake_ele)
        n = n + 1
        call re_allocate(z_knot, n)
        call re_allocate(wx_knot, n)
        call re_allocate(wy_knot, n)
        call re_allocate(wz_knot, n)
        z_knot(n) = z
        wx_knot(n) = p2%vec(2)
        wy_knot(n) = p2%vec(4)
        wz_knot(n) = p2%vec(6)
      enddo
    enddo
  endif

end select

! Plotting

if (make_plot) call make_this_plot ('X')
if (postscript_file /= '') call make_this_plot ('PS')

!------------------------------------------------------
contains

subroutine make_this_plot(who)

type (qp_line_struct), allocatable :: line(:)
type (qp_symbol_struct), allocatable :: symbol(:)

real(rp) height
character(*) who
character(16) :: color(6) = [character(16):: 'blue', 'red', 'green', 'cyan', 'magenta', 'yellow']
integer j, k, n, id

!

if (who == 'X') then
  call qp_open_page ('X', id, plot_size(1), plot_size(2), 'POINTS')
elseif (plot_size(1) > plot_size(2)) then
  call qp_open_page ('PS-L', id, plot_size(1), plot_size(2), 'POINTS', plot_file = postscript_file)
else
  call qp_open_page ('PS', id, plot_size(1), plot_size(2), 'POINTS', plot_file = postscript_file)
endif

call qp_get_text_attrib ('AXIS_NUMBERS', height)
call qp_set_text_attrib ('AXIS_NUMBERS', text_scale*height)

call qp_get_text_attrib ('AXIS_LABEL', height)
call qp_set_text_attrib ('AXIS_LABEL', text_scale*height)

call qp_set_page_border (0.01_rp, 0.01_rp, 0.01_rp, 0.01_rp, '%PAGE')
call qp_set_margin (0.10_rp, 0.02_rp, 0.07_rp, 0.01_rp, '%PAGE')

call qp_calc_and_set_axis ('X', plot_min, plot_max, 4, 8, 'ZERO_AT_END')

!

call qp_set_box (1, 3, 1, 3)
call qp_calc_and_set_axis ('Y', minval(wx), maxval(wx), 4, 8, 'GENERAL')
call qp_draw_axes (x_axis_label, y_axis_label)

call qp_draw_data (time, wx, .true., 0)
if (draw_knot_points) then
  call qp_draw_symbols (z_knot, wx_knot)
endif

!

call qp_set_box (1, 2, 1, 3)
call qp_calc_and_set_axis ('Y', minval(wy), maxval(wy), 4, 8, 'GENERAL')
call qp_draw_axes (x_axis_label, y_axis_label)

call qp_draw_data (time, wy, .true., 0)
if (draw_knot_points) then
  call qp_draw_symbols (z_knot, wy_knot)
endif

!

call qp_set_box (1, 1, 1, 3)
call qp_calc_and_set_axis ('Y', minval(wz), maxval(wz), 4, 8, 'GENERAL')
call qp_draw_axes (x_axis_label, y_axis_label)

call qp_draw_data (time, wz, .true., 0)
if (draw_knot_points) then
  call qp_draw_symbols (z_knot, wz_knot)
endif

!

if (who == 'X') then
  write (*, '(a)', advance = 'NO') 'Hit any key to continue.'
  read (*, '(a)') ans
endif

call qp_close_page()

end subroutine

end program
