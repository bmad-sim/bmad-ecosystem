
program controller_function_plot

use bmad
use quick_plot

implicit none

type (lat_struct) lat
type (ele_pointer_struct), allocatable, target :: eles(:)
type (ele_struct), pointer :: ele, slave
type (control_struct), pointer :: ctl
type (qp_line_struct) line(20)

real(rp) var, var_min, var_max, plot_size(2)
real(rp), allocatable :: table(:,:)

integer n_curve_points, slave_list(20)
integer i, j, n, id, n_loc, ix_var, n_slave, n_slave_max

logical make_plot, make_postscript, err, err_flag, ok

character(200) init_file, lat_file, table_file_name, err_str
character(40) controller_name, control_var
character(40) attrib_name(20)

namelist / params / var_min, var_max, slave_list, table_file_name, make_plot, plot_size, &
                    controller_name, control_var, lat_file, n_curve_points, make_postscript

! Read input

init_file = 'controller_function_plot'
if (cesr_iargc() > 0) call cesr_getarg(1, init_file)
print '(a)', 'Init file: ', trim(init_file)

slave_list = -1
n_curve_points = 200
table_file_name = 'table.dat'
make_plot = .true.
make_postscript = .false.

open (1, file = init_file)
read (1, nml = params)
close (1)

! Init

call bmad_parser(lat_file, lat, err_flag = err)
if (err) stop

call lat_ele_locator (controller_name, lat, eles, n_loc)

if (n_loc == 0) then
  print '(2a)', 'CANNOT FIND CONTROLLER: ', trim(controller_name)
  stop
endif

if (n_loc > 1) then
  print '(2a)', 'More than one controller found with name: ', trim(controller_name)
  print '(a)',  'Using first one...'
endif

ele => eles(1)%ele
if (ele%key /= overlay$ .and. ele%key /= group$ .and. ele%key /= ramper$) then
  print '(2a)', 'ELEMENT IS NOT A CONTROLLER: ', trim(ele%name)
  stop
endif

if (control_var == '') then
  ix_var = 1
else
  ix_var = -1
  do i = 1, size(ele%control%var)
    if (control_var == ele%control%var(i)%name) ix_var = i
  enddo
  if (ix_var == -1) then
    print '(2a)', 'NO CONTROL VARIABLE FOUND WITH THIS NAME: ', trim(control_var)
    stop
  endif
endif

if (ele%key == ramper$) then
  n_slave_max = size(ele%control%ramp)
else
  n_slave_max = ele%n_slave
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
  n_slave = n_slave_max
  do j = 1, n_slave
    slave_list(j) = j
  enddo
endif

allocate (table(0:n_slave, 0:n_curve_points))

if (var_min == var_max) then
  if (ele%control%type == expression$) then
    print '(a)', 'CANNOT FIND CONTROL BOUNDS SINCE CONTROLLER USES EXPRESSIONS AND NOT A LIST OF KNOT POINTS.'
    print '(a)', 'IN THIS CASE PLEASE SET VAR_MIN AND VAR_MAX SO THAT THERE IS A FINITE RANGE THAT CAN BE USED.'
    stop
  endif

  n = size(ele%control%x_knot)
  var_min = ele%control%x_knot(1)
  var_max = ele%control%x_knot(n)
endif

! Make table



do i = 0, n_curve_points
  var = var_min + i * (var_max - var_min) / n_curve_points
  ele%control%var(ix_var)%value = var
  table(i,0) = var

  do j = 1, n_slave
    if (ele%key == ramper$) then
      ctl => ele%control%ramp(slave_list(j))
    else
      slave => pointer_to_slave(ele, slave_list(j), ctl)
    endif

    if (ele%control%type == expression$) then
      call evaluate_expression_stack(ctl%stack, ctl%value, err_flag, err_str, ele%control%var, .false.)
    else
      call spline_akima_interpolate (ele%control%x_knot, ctl%y_knot, ele%control%var(1)%value, ok, ctl%value)
    endif

    table(i,j) = ctl%value
    attrib_name(j) = ctl%attribute
  enddo
enddo

! Write table

if (table_file_name /= '') then
  open (1, file = table_file_name, status = 'NEW')
  write (1, '(a, 20a16)') '#  Ix   Var    ',  (attrib_name(j), j = 1, n_slave)
  do i = 0, n_curve_points
    write (1, '(i5, 21es16.8)') i, (table(i,j), j = 0, n_slave)
  enddo
endif

! Plotting

if (make_plot) then
  line(1:6)%color = [character(16):: 'blue', 'red', 'green', 'cyan', 'magenta', 'yellow']
  call qp_open_page ('X', id, plot_size(1), plot_size(2), 'POINTS')
  call qp_set_page_border (0.02_rp, 0.02_rp, 0.02_rp, 0.02_rp, '%PAGE')
  call qp_set_box (1, 1, 1, 1)
  call qp_calc_and_set_axis ('X', var_min, var_max, 4, 8, 'GENERAL')
  call qp_calc_and_set_axis ('Y', minval(table(:,1:)), maxval(table(:,1:)), 4, 8, 'GENERAL')
  call qp_draw_axes (ele%control%var(ix_var)%name, 'Slave Value')
  do j = 1, n_slave
    call qp_set_line_attrib ('PLOT', color = line(j)%color)
    call qp_draw_data (table(:,0), table(:,j), .true., 0)
  enddo
  call qp_draw_curve_legend (0.1_rp, 0.1_rp, '%/GRAPH/LT', line(1:n_slave), text = attrib_name(1:n_slave))
endif

end program
