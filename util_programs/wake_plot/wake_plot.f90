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
type (wake_sr_z_long_struct), pointer :: srz
type (wake_lr_mode_struct), allocatable :: lr_mode(:)
type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p1, p2

real(rp) :: x_min, x_max, plot_size(2), text_scale, zero6(6) = 0
real(rp) xy_leading(2), xy_trailing(2), vec2(6), z
real(rp), allocatable :: x_axis(:), wx(:), wy(:), wz(:)

integer n, n_points, ix_wake, m_order, i, im, iz, ik, n_loc, ix, nw, n1, n2

logical make_plot, err, err_flag, ok

character(40) wake_ele_name, ans, who, plot_type
character(200) init_file, lat_file, postscript_file

namelist / params / x_min, x_max, make_plot, plot_size, who, plot_type, &
                    wake_ele_name, lat_file, n_points, postscript_file, &
                    text_scale, xy_leading, xy_trailing, ix_wake, m_order

! Read input

init_file = 'wake_plot.init'
if (command_argument_count() > 0) call get_command_argument(1, init_file)
print '(2a)', 'Init file: ', trim(init_file)

n_points = 1001
make_plot = .true.
postscript_file = ''
text_scale = 1.0
ix_wake = 0
m_order = 0
x_min = 0
x_max = 0

open (1, file = init_file, status = 'old')
read (1, nml = params)
close (1)

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
vec2 = p2%vec

p1%charge = 1
p2%charge = 1e-10
bunch%charge_live = 1

!

if (who == 'sr-z-long' .and. plot_type == 'wake') then
  srz => wake_ele%wake%sr%z_long
  if (size(srz%w) == 0) then
    print *, 'No SR Z_long wake.'
    stop
  endif

  if (x_min < -srz%z0) x_min = -srz%z0
  if (x_max > srz%z0) x_max = srz%z0

  if (x_min == 0 .and. x_max == 0) then
    x_min = -srz%z0
    x_max = srz%z0
  endif

  nw = (size(srz%w) - 1) / 2
  n1 = nint(x_min / srz%dz)
  n2 = nint(x_max / srz%dz)
  allocate(wz(n1:n2), x_axis(n1:n2))

  do im = n1, n2
    x_axis(im) = im * srz%dz
    wz(im) = srz%w(im + nw + 1)
  enddo

  if (make_plot) call make_this_plot ('X')
  if (postscript_file /= '') call make_this_plot ('PS')
  stop
endif

! 

allocate (x_axis(n_points), wx(n_points), wy(n_points), wz(n_points))

select case (who(1:2))
case ('lr')
  lr => wake_ele%wake%lr
  bmad_com%lr_wakes_on = .true.

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
    z = x_min + (x_max - x_min) * (i - 1.0_rp) / (n_points - 1.0_rp)
    p2%vec = vec2
    p2%vec(5) = z
    x_axis(i) = z
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

  if (x_min == 0 .and. x_max == 0) then
    if (sr%z_long%z0 /= 0) then
      x_min = -sr%z_long%z0
      x_max =  sr%z_long%z0
    else
      print *, 'X_min nor X_max set! I do not know what to do!'
      stop
    endif
  endif

  if (who /= 'sr-z-long') then
    srz => sr%z_long
    deallocate (srz%w, srz%fw, srz%fbunch, srz%w_out)
    allocate (srz%w(0), srz%fw(0), srz%fbunch(0), srz%w_out(0))
  endif

  if (who /= 'sr-mode' .and. who /= 'sr-trans') then
    deallocate (sr%trans)
    allocate (sr%trans(0))
  endif

  if (who /= 'sr_mode' .and. who /= 'sr_long') then
    deallocate (sr%long)
    allocate (sr%long(0))
  endif


  do i = 1, n_points
    z = x_min + (x_max - x_min) * (i - 1.0_rp) / (n_points - 1.0_rp)
    p2%vec = vec2
    p2%vec(5) = z
    x_axis(i) = z
    call order_particles_in_z(bunch)
    call track1_sr_wake(bunch, wake_ele)
    wx(i) = p2%vec(2)
    wy(i) = p2%vec(4)
    wz(i) = p2%vec(6)
  enddo

end select

! Plotting

if (make_plot) call make_this_plot ('X')
if (postscript_file /= '') call make_this_plot ('PS')

!------------------------------------------------------
contains

subroutine make_this_plot(plot_window)

type (qp_line_struct), allocatable :: line(:)
type (qp_symbol_struct), allocatable :: symbol(:)

real(rp) height
character(*) plot_window
character(16) :: y_label, x_axis_str
integer j, k, n, id

!

if (plot_window == 'X') then
  call qp_open_page ('X', id, plot_size(1), plot_size(2), 'POINTS')
elseif (plot_size(1) > plot_size(2)) then
  call qp_open_page ('PS-L', id, plot_size(1), plot_size(2), 'POINTS', plot_file = postscript_file)
else
  call qp_open_page ('PS', id, plot_size(1), plot_size(2), 'POINTS', plot_file = postscript_file)
endif

select case (plot_type)
case ('wake'); y_label = 'Wake'
case ('kick'); y_label = 'Kick'
case default; y_label = 'Wake'
end select

select case (plot_type)
case ('lr');  x_axis_str = 'Time'
case default; x_axis_str = 'Z'
end select

call qp_get_text_attrib ('AXIS_NUMBERS', height)
call qp_set_text_attrib ('AXIS_NUMBERS', text_scale*height)

call qp_get_text_attrib ('AXIS_LABEL', height)
call qp_set_text_attrib ('AXIS_LABEL', text_scale*height)

call qp_set_page_border (0.01_rp, 0.01_rp, 0.01_rp, 0.01_rp, '%PAGE')
call qp_set_margin (0.15_rp, 0.02_rp, 0.07_rp, 0.01_rp, '%PAGE')

! 

if (who == 'sr-z-long' .and. plot_type == 'wake') then
  call qp_calc_and_set_axis ('X', x_min, x_max, 4, 8, 'GENERAL')
  call qp_set_box (1, 1, 1, 1)
  call qp_calc_and_set_axis ('Y', minval(wz), maxval(wz), 4, 8, 'GENERAL')
  call qp_draw_axes (x_axis_str, trim(y_label) // '_z')

  call qp_draw_data (x_axis, wz, .true., 0)

!

else
  call qp_calc_and_set_axis ('X', x_min, x_max, 4, 8, 'GENERAL')

  call qp_set_box (1, 3, 1, 3)
  call qp_calc_and_set_axis ('Y', minval(wx), maxval(wx), 4, 8, 'GENERAL')
  call qp_draw_axes (x_axis_str, trim(y_label) // '_x')

  call qp_draw_data (x_axis, wx, .true., 0)

  !

  call qp_set_box (1, 2, 1, 3)
  call qp_calc_and_set_axis ('Y', minval(wy), maxval(wy), 4, 8, 'GENERAL')
  call qp_draw_axes (x_axis_str, trim(y_label) // '_y')

  call qp_draw_data (x_axis, wy, .true., 0)

  !

  call qp_set_box (1, 1, 1, 3)
  call qp_calc_and_set_axis ('Y', minval(wz), maxval(wz), 4, 8, 'GENERAL')
  call qp_draw_axes (x_axis_str, trim(y_label) // '_z')
  call qp_draw_data (x_axis, wz, .true., 0)
endif

!

if (plot_window == 'X') then
  write (*, '(a)', advance = 'NO') 'Hit any key to continue.'
  read (*, '(a)') ans
endif

call qp_close_page()

end subroutine

end program
