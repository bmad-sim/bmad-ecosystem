
program controller_function_plot

use bmad
use quick_plot

implicit none

type (lat_struct) lat
type (ele_pointer_struct), allocatable, target :: eles(:), rampers(:)
type (ele_struct), pointer :: controller, slave
type (ele_struct) ramp_slave
type (control_struct), pointer :: ctl

real(rp) var, var_min, var_max, plot_size(2), text_scale, height
real(rp), allocatable :: table(:,:), knots(:,:)

integer, allocatable :: ramper_list(:), slave_list2(:)
integer n_var_points, slave_list(100)
integer i, j, n, id, n_loc, ix_var, n_slave, n_slave_max

logical make_plot, err, err_flag, ok, draw_knot_points

character(200) init_file, lat_file, table_file_name, err_str, postscript_file
character(40) controller_name, control_var, ans
character(80), allocatable :: attrib_name(:)

namelist / params / var_min, var_max, slave_list, table_file_name, make_plot, plot_size, &
                    controller_name, control_var, lat_file, n_var_points, postscript_file, &
                    draw_knot_points, text_scale

! Read input

init_file = 'controller_function_plot.init'
if (cesr_iargc() > 0) call cesr_getarg(1, init_file)
print '(2a)', 'Init file: ', trim(init_file)

slave_list = -1                   ! Set some defaults
n_var_points = 200
table_file_name = 'table.dat'
make_plot = .true.
postscript_file = ''
draw_knot_points = .true.
text_scale = 1.2

open (1, file = init_file, status = 'old')
read (1, nml = params)
close (1)

! Init

call bmad_parser(lat_file, lat, err_flag = err)
if (err) stop

!

call lat_ele_locator (controller_name, lat, eles, n_loc)

if (n_loc == 0) then
  print '(2a)', 'CANNOT FIND CONTROLLER: ', trim(controller_name)
  stop
endif

if (n_loc > 1) then
  print '(2a)', 'More than one controller found with name: ', trim(controller_name)
  print '(a)',  'Using first one...'
endif

controller => eles(1)%ele
if (controller%key /= overlay$ .and. controller%key /= group$ .and. controller%key /= ramper$) then
  print '(2a)', 'ELEMENT IS NOT A CONTROLLER: ', trim(controller%name)
  stop
endif

if (control_var == '') then
  ix_var = 1
else
  ix_var = -1
  do i = 1, size(controller%control%var)
    if (control_var == controller%control%var(i)%name) ix_var = i
  enddo
  if (ix_var == -1) then
    print '(2a)', 'NO CONTROL VARIABLE FOUND WITH THIS NAME: ', trim(control_var)
    stop
  endif
endif

if (controller%key == ramper$) then
  n_slave_max = size(controller%control%ramp)
else
  n_slave_max = controller%n_slave
endif

do j = 1, size(slave_list)
  if (slave_list(j) < 0) exit
  if (slave_list(j) < 1 .or. slave_list(j) > n_slave_max) then
    print '(a, i0)', 'INDEX IN SLAVE_LIST IS OUT OF BOUNDS: ', slave_list(j)
    stop
  endif
enddo

n_slave = j - 1
if (n_slave == 0) then
  allocate (slave_list2(n_slave_max))
  n_slave = n_slave_max
  do j = 1, n_slave
    slave_list2(j) = j
  enddo
else
  allocate (slave_list2(n_slave))
  slave_list2 = slave_list(1:n_slave)
endif

allocate (table(0:n_var_points, 0:n_slave))

if (var_min == var_max) then
  if (controller%control%type == expression$) then
    print '(a)', 'CANNOT FIND CONTROL BOUNDS SINCE CONTROLLER USES EXPRESSIONS AND NOT A LIST OF KNOT POINTS.'
    print '(a)', 'IN THIS CASE PLEASE SET VAR_MIN AND VAR_MAX SO THAT THERE IS A FINITE RANGE THAT CAN BE USED.'
    stop
  endif

  n = size(controller%control%x_knot)
  var_min = controller%control%x_knot(1)
  var_max = controller%control%x_knot(n)
endif

! Start knot table

if (controller%control%type == expression$) then
  draw_knot_points = .false.
elseif (draw_knot_points) then
  allocate (knots(size(controller%control%x_knot), 0:n_slave))
  knots(:,0) = controller%control%x_knot
endif

! Make table and plot

ramp_slave%key = int_garbage$   ! Signals routines to skip some inconvienient bookkeeping.
call lat_ele_locator ('RAMPER::*', lat, rampers, n_loc)

do i = 0, n_var_points
  var = var_min + i * (var_max - var_min) / n_var_points
  controller%control%var(ix_var)%value = var
  table(i,0) = var

  select case (controller%key)
  case (ramper$)
    call apply_ramper (ramp_slave, rampers, err_flag)
  case (overlay$)
  do j = 1, n_slave
    slave => pointer_to_slave(controller, slave_list2(j))
    call makeup_control_slave(lat, slave, err_flag)
  enddo
  case (group$)
    call makeup_group_lord (lat, controller, err_flag)
  end select

  call re_allocate(attrib_name, n_slave)

  do j = 1, n_slave
    if (controller%key == ramper$) then
      ctl => controller%control%ramp(slave_list2(j))
    else
      slave => pointer_to_slave(controller, slave_list2(j), ctl)
    endif

    if (draw_knot_points .and. i == 0) knots(:,j) = ctl%y_knot
    table(i,j) = ctl%value
    attrib_name(j) = trim(ctl%slave_name) // '[' // trim(ctl%attribute) // ']'
  enddo
enddo

! Write table

if (table_file_name /= '') then
  open (1, file = table_file_name)
  write (1, '(a, 20a16)') '#  Ix   Var    ',  (attrib_name(j), j = 1, n_slave)
  do i = 0, n_var_points
    write (1, '(i5, 21es16.8)') i, (table(i,j), j = 0, n_slave)
  enddo
endif

! Plotting

if (make_plot) call make_this_plot ('X')
if (postscript_file /= '') call make_this_plot ('PS')

!------------------------------------------------------
contains

subroutine make_this_plot(who)

type (qp_line_struct), allocatable :: line(:)
type (qp_symbol_struct), allocatable :: symbol(:)

character(*) who
character(16) :: color(6) = [character(16):: 'blue', 'red', 'green', 'cyan', 'magenta', 'yellow']
integer j, k, n

!

allocate (line(n_slave), symbol(n_slave))

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
call qp_set_box (1, 1, 1, 1)
call qp_calc_and_set_axis ('X', var_min, var_max, 4, 8, 'GENERAL')
call qp_calc_and_set_axis ('Y', minval(table(:,1:)), maxval(table(:,1:)), 4, 8, 'GENERAL')
if (n_slave == 1) then
  call qp_draw_axes (controller%control%var(ix_var)%name, attrib_name(1))
else
  call qp_draw_axes (controller%control%var(ix_var)%name, 'Slave Value')
endif

do j = 1, n_slave
  k = modulo(j-1,6) + 1
  line(j)%color = color(k)
  symbol(j)%color = color(k)
  call qp_set_line_attrib ('PLOT', color = line(j)%color)
  call qp_draw_data (table(:,0), table(:,j), .true., 0)
  if (draw_knot_points) then
    call qp_set_symbol(symbol(j))
    call qp_draw_symbols (knots(:,0), knots(:,j))
  endif
enddo

n = min(n_slave, 20)
call qp_draw_curve_legend (0.1_rp, -0.1_rp, '%/GRAPH/LT', line(1:n), text = attrib_name(1:n))

if (who == 'X') then
  write (*, '(a)', advance = 'NO') 'Hit any key to continue.'
  read (*, '(a)') ans
endif

call qp_close_page()

end subroutine

end program
