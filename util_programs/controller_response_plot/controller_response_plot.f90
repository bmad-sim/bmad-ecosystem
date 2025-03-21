
program controller_response_plot

use bmad
use quick_plot

implicit none

type (lat_struct) lat
type (ele_pointer_struct), allocatable, target :: eles(:), rf_eles(:)
type (ele_struct), pointer :: controller, slave
type (control_struct), pointer :: ctl
type (control_ramp1_struct), pointer :: rmp
type (coord_struct) orbit
type (qp_legend_struct) :: legend = qp_legend_struct()

real(rp) :: control_var, input_var_min, input_var_max, plot_size(2), text_scale, height, vzero(6) = 0
real(rp) curve_min, curve_max, z, dz, volt_old, pot, value
real(rp), allocatable :: table(:,:), knots(:,:)

integer, allocatable :: ramper_list(:), curve_list2(:)
integer curve_n_points, input_var_n_points
integer i, ii, j, jj, n, id, n_ramp_loc, n_control_loc, ix_var, n_curve, n_slave_max, n_rf

logical make_plot, err, err_flag, ok, draw_knot_points

character(200) init_file, lat_file, table_file_name, err_str, postscript_file, curve_list
character(40) controller_name, input_var, ans, x_axis_label, y_axis_label, input_var_name
character(60), allocatable :: attrib_name(:), time(:)

namelist / params / input_var_min, input_var_max, curve_list, table_file_name, make_plot, plot_size, &
                    controller_name, input_var, lat_file, input_var_n_points, postscript_file, &
                    draw_knot_points, text_scale, curve_min, curve_max, curve_n_points

! Read input

init_file = 'controller_response_plot.init'
if (command_argument_count() > 0) call get_command_argument(1, init_file)
print '(2a)', 'Init file: ', trim(init_file)

curve_list = '*'                   ! Set some defaults
input_var_n_points = 201
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

call lat_ele_locator (controller_name, lat, eles, n_control_loc)

if (n_control_loc == 0) then
  print '(2a)', 'CANNOT FIND CONTROLLER: ', trim(controller_name)
  stop
endif

! Vary RF

if (input_var == 'z') then
  input_var_name = 'z_pos'
  x_axis_label = 'z'
  y_axis_label = curve_list
  draw_knot_points = .false.
  n_curve = curve_n_points
  allocate (table(0:input_var_n_points-1, 0:n_curve))
  call lat_ele_locator ('rfcavity::*', lat, rf_eles, n_rf)
  call re_allocate(attrib_name, n_curve)

  do i = 1, n_curve
    control_var = curve_min + (i-1) * (curve_max - curve_min) / max(1, curve_n_points - 1)
    write (attrib_name(i), '(a, es10.3)') 'time:', control_var

    do ii = 1, n_control_loc
      controller => eles(ii)%ele
      controller%control%var(1)%value = control_var    ! 1 is hard coded. Need to change this.
    enddo

    call apply_all_rampers(lat, err)
    call lattice_bookkeeper(lat)

    do j = 0, input_var_n_points-1
      call init_coord (orbit, vzero, rf_eles(1)%ele, upstream_end$)
      dz = (input_var_max - input_var_min) / max(1, input_var_n_points-1)
      z = input_var_min + j * dz
      table(j,0) = z
      orbit%vec(5) = z

      do jj = 1, n_rf
        call track1(orbit, rf_eles(jj)%ele, lat%param, orbit)
      enddo
      table(j, i) = orbit%vec(6) * orbit%p0c
    enddo

    if (curve_list == 'rf::potential') then   ! Integrate voltage
      volt_old = table(0, i)
      table(0, i) = 0
      do j = 1, input_var_n_points-1
        pot = table(j-1, i) + (volt_old + table(j, i)) * dz
        volt_old = table(j, i)
        table(j, i) = pot
      enddo

      pot = minval(table(:, i))
      table(:, i) = table(:, i) - pot
    endif
  enddo

! Vary a controller var

else
  input_var_name = 'Cntl Var'

  if (n_control_loc > 1) then
    print '(2a)', 'More than one controller found with name: ', trim(controller_name)
    print '(a)',  'Using first one...'
  endif

  controller => eles(1)%ele
  if (controller%key /= overlay$ .and. controller%key /= group$ .and. controller%key /= ramper$) then
    print '(2a)', 'ELEMENT IS NOT A CONTROLLER: ', trim(controller%name)
    stop
  endif

  if (input_var == '') then
    ix_var = 1
  else
    ix_var = -1
    do i = 1, size(controller%control%var)
      if (input_var == controller%control%var(i)%name) ix_var = i
    enddo
    if (ix_var == -1) then
      print '(2a)', 'NO CONTROL VARIABLE FOUND WITH THIS NAME: ', trim(input_var)
      stop
    endif
  endif

  if (controller%key == ramper$) then
    n_slave_max = size(controller%control%ramp)
  else
    n_slave_max = controller%n_slave
  endif

  if (curve_list == 'slave::*') then
    n_curve = n_slave_max
    allocate (curve_list2(n_curve))
    do j = 1, n_curve
      curve_list2(j) = j
    enddo
  else
    n_curve = 1
    allocate (curve_list2(n_curve))
    read(curve_list, *) curve_list2(1)
  endif

  allocate (table(0:input_var_n_points-1, 0:n_curve))

  if (input_var_min == input_var_max) then
    if (.not. allocated(controller%control%x_knot)) then
      print '(a)', 'CANNOT FIND CONTROL BOUNDS SINCE CONTROLLER USES EXPRESSIONS AND NOT A LIST OF KNOT POINTS.'
      print '(a)', 'IN THIS CASE PLEASE SET INPUT_VAR_MIN AND INPUT_VAR_MAX SO THAT THERE IS A FINITE RANGE THAT CAN BE USED.'
      stop
    endif

    n = size(controller%control%x_knot)
    input_var_min = controller%control%x_knot(1)
    input_var_max = controller%control%x_knot(n)
  endif

  ! Start knot table

  if (allocated(controller%control%x_knot)) then
    allocate (knots(size(controller%control%x_knot), 0:n_curve))
    knots(:,0) = controller%control%x_knot
  else
    draw_knot_points = .false.
  endif

  ! Make table and plot

  call re_allocate(attrib_name, n_curve)

  do i = 0, input_var_n_points-1
    control_var = input_var_min + i * (input_var_max - input_var_min) / max(1, input_var_n_points-1)
    controller%control%var(ix_var)%value = control_var
    table(i,0) = control_var

    select case (controller%key)

    case (overlay$)
      do j = 1, n_curve
        slave => pointer_to_slave(controller, curve_list2(j))
        call makeup_control_slave(lat, slave, err_flag)
      enddo

    case (group$)
      call makeup_group_lord (lat, controller, err_flag)
    end select

    do j = 1, n_curve
      if (controller%key == ramper$) then
        rmp => controller%control%ramp(curve_list2(j))
        value = ramper_value(controller, rmp, err_flag)
        if (draw_knot_points .and. i == 0) knots(:,j) = rmp%y_knot
        table(i,j) = value
        attrib_name(j) = trim(rmp%slave_name) // '[' // trim(rmp%attribute) // ']'

      else
        slave => pointer_to_slave(controller, curve_list2(j), ctl)
        if (draw_knot_points .and. i == 0) knots(:,j) = ctl%y_knot
        table(i,j) = ctl%value
        attrib_name(j) = trim(ctl%slave_name) // '[' // trim(ctl%attribute) // ']'
      endif
    enddo
  enddo

  !

  x_axis_label = controller%control%var(ix_var)%name
  if (n_curve == 1) then
    y_axis_label = attrib_name(1)
  else
    y_axis_label = 'Slave Value'
  endif
endif

! Write table

if (table_file_name /= '') then
  open (1, file = table_file_name)
  write (1, '(2a, t24, 20a16)') '#  Ix   ', trim(input_var_name),  (attrib_name(j), j = 1, n_curve)
  do i = 0, input_var_n_points-1
    write (1, '(i5, 21es16.8)') i, (table(i,j), j = 0, n_curve)
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

allocate (line(n_curve), symbol(n_curve))

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
call qp_calc_and_set_axis ('X', input_var_min, input_var_max, 4, 8, 'GENERAL')
call qp_calc_and_set_axis ('Y', minval(table(:,1:)), maxval(table(:,1:)), 4, 8, 'GENERAL')
call qp_draw_axes (x_axis_label, y_axis_label)

do j = 1, n_curve
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

n = min(n_curve, 20)
call qp_draw_curve_legend (qp_point_struct(0.1_rp, -0.1_rp, '%/GRAPH/LT'), legend, line(1:n), text = attrib_name(1:n))

if (who == 'X') then
  write (*, '(a)', advance = 'NO') 'Hit any key to continue.'
  read (*, '(a)') ans
endif

call qp_close_page()

end subroutine

end program
