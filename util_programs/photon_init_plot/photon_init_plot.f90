program photon_init_plot

use bmad
use quick_plot

implicit none

type (lat_struct) lat
type (coord_struct) p_start, p_end

real(rp), allocatable :: x_dat(:), y_dat(:)
real(rp) x, x_min, x_max, plot_size(2), height, text_scale
integer i, n, ip, n_bin, n_photons, ix_var, iy_var, i_chan
logical ignore_out_of_bounds

character(200) init_file, lat_name
character(40) x_var, y_var, plot_type, ans

namelist / params / lat_name, n_bin, n_photons, x_min, x_max, x_var, y_var, plot_type, ignore_out_of_bounds, &
              plot_size, text_scale

!


text_scale = 1.2

call get_command_argument(1, init_file)
if (init_file == '') init_file = 'photon_init_plot.init'
open (1, file = init_file, status = 'old')
read (1, nml = params)
close (1)

!

select case (x_var)
case ('x');     ix_var = 1
case ('vx');    ix_var = 2
case ('y');     ix_var = 3
case ('vy');    ix_var = 4
case ('z');     ix_var = 5
case ('vz');    ix_var = 6
case ('E');     ix_var = 7
case default
  print *, 'Unknown x_var: ' // x_var
  stop
end select

if (plot_type == '2D') then
  select case (y_var)
  case ('x');     iy_var = 1
  case ('vx');    iy_var = 2
  case ('y');     iy_var = 3
  case ('vy');    iy_var = 4
  case ('z');     iy_var = 5
  case ('vz');    iy_var = 6
  case ('E');     iy_var = 7
  case default
    print *, 'Unknown y_var: ' // y_var
    stop
  end select
endif

!

call bmad_parser(lat_name, lat)
allocate(y_dat(0:n_bin-1), x_dat(0:n_bin-1))

do i = 0, n_bin-1
  x_dat(i) = x_min + (i + 0.5_rp) * (x_max - x_min) / n_bin 
enddo

y_dat = 0
do ip = 1, n_photons
  call init_photon_from_a_photon_init_ele (lat%ele(1), lat%param, p_end)
  if (ix_var == 7) then
    x = p_end%p0c
  else
    x = p_end%vec(ix_var)
  endif
  n = floor(n_bin * (x - x_min) / (x_max - x_min))
  if (ignore_out_of_bounds .and. (n < 0 .or. n >= n_bin)) cycle
  n = max(min(n, n_bin-1), 0)
  y_dat(n) = y_dat(n) + 1
enddo

y_dat = y_dat * n_bin / ((x_max - x_min) * n_photons)

! Plot

call qp_open_page ('X', i_chan, plot_size(1), plot_size(2), 'POINTS')
call qp_get_text_attrib ('AXIS_NUMBERS', height)
call qp_set_text_attrib ('AXIS_NUMBERS', text_scale*height)

call qp_get_text_attrib ('AXIS_LABEL', height)
call qp_set_text_attrib ('AXIS_LABEL', text_scale*height)

call qp_set_page_border (0.01_rp, 0.01_rp, 0.01_rp, 0.01_rp, '%PAGE')
call qp_set_margin (0.10_rp, 0.02_rp, 0.07_rp, 0.01_rp, '%PAGE')
call qp_set_box (1, 1, 1, 1)
call qp_calc_and_set_axis ('X', x_min, x_max, 4, 8, 'GENERAL')
call qp_calc_and_set_axis ('Y', 0.0_rp, maxval(y_dat), 4, 8, 'GENERAL')

call qp_draw_axes(x_var, 'Probability', draw_grid = .true.)
call qp_draw_histogram(x_dat, y_dat)

  write (*, '(a)', advance = 'NO') 'Hit any key to continue.'
  read (*, '(a)') ans


end program
