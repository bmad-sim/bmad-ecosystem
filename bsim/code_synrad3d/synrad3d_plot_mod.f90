module synrad3d_plot_mod

use synrad3d_track_mod
use quick_plot
use input_mod

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_plot_diffuse_probability (plot_param, lat)
!
! Routine to plot diffuse reflection probability curves
!
! Input;
!   plot_param -- sr3d_plot_param_struct: Plot parameters.
!-

subroutine sr3d_plot_diffuse_probability (plot_param, lat)

implicit none

type (lat_struct), target :: lat
type (sr3d_plot_param_struct) plot_param
type (photon_reflect_surface_struct), pointer :: surface

real(rp) energy, angle, theta_out, phi_out, x(100), y1(100), y2(100)
integer i, ix, i_chan

character(80) ans

!

energy = 34.0
angle = 1.45
surface => sr3d_com%surface(1)

do i = 1, 100
  x(i) = i / 100.0_rp
  y1(i) = 0
  y2(i) = 0
enddo

!

call qp_open_page ('X', i_chan, plot_param%window_width, plot_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
call qp_set_margin (0.07_rp, 0.05_rp, 0.10_rp, 0.05_rp, '%PAGE')

!

do

  call qp_calc_and_set_axis ('X', 0.0_rp, 1.0_rp, 10, 16, 'GENERAL')
  call qp_calc_and_set_axis ('Y', min(minval(y1), minval(y2)), max(maxval(y1), maxval(y2)), 6, 10, 'GENERAL')

  call qp_clear_page

  call qp_draw_graph (x, y1, 'x', 'cum_x', 'cheb', .true., 0)
  call qp_draw_polyline (x, y2, line_pattern = 'solid')

  print *, 'surface_roughness_rms:    ', surface%surface_roughness_rms
  print *, 'roughness_correlation_len:', surface%roughness_correlation_len

  call read_a_line ('Input: ', ans)
  call string_trim(ans, ans, ix)
  if (ix == 0) cycle

enddo

end subroutine sr3d_plot_diffuse_probability

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_plot_reflection_probability (plot_param, lat)
!
! Routine to plot reflection probability curves
!
! Input;
!   plot_param -- sr3d_plot_param_struct: Plot parameters.
!   lat        -- lat_struct: Lattice
!-

subroutine sr3d_plot_reflection_probability (plot_param, lat)

implicit none

type (lat_struct), target :: lat
type (sr3d_plot_param_struct) plot_param
type (photon_reflect_surface_struct), pointer :: surface
type (qp_legend_struct) legend

type yplot
  real(rp), allocatable :: y_reflect(:)
  real(rp), allocatable :: y_rel_specular(:)
  real(rp), allocatable :: y(:)
  character(40) label
  type (qp_line_struct) line
end type

type (yplot), allocatable :: ny(:)

real(rp), target :: angle_min, angle_max, ev_min, ev_max, angle, ev
real(rp) value, value1, value2, y_min, y_max
real(rp), allocatable :: x(:)

integer i, j, n, ix, ios, n_lines, i_chan, n_max, d_max, ix_surface

logical fixed_energy, logic

character(80) ans, head_lab, descrip, fmt, text
character(16) x_lab, y_lab, param, reflection_type

! init

angle_min = 0
angle_max = 40

ev_min = 50
ev_max = 150

fixed_energy = .true.
y_lab = 'Reflectivity'
reflection_type = 'total'
head_lab = 'Reflectivity (Total)'

ix_surface = 1
surface => sr3d_com%surface(ix_surface)

n_lines = 3
allocate (ny(n_lines))

do i = 1, n_lines
  allocate (ny(i)%y_reflect(plot_param%n_pt))
  allocate (ny(i)%y_rel_specular(plot_param%n_pt))
  allocate (ny(i)%y(plot_param%n_pt))
enddo

allocate (x(plot_param%n_pt))

call qp_open_page ('X', i_chan, plot_param%window_width, plot_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
call qp_set_margin (0.07_rp, 0.05_rp, 0.10_rp, 0.05_rp, '%PAGE')

! Endless plotting loop:

do 

  ! Get data

  do i = 1, n_lines
    do j = 1, plot_param%n_pt
      if (fixed_energy) then
        angle = angle_min + (j - 1) * (angle_max - angle_min) / (plot_param%n_pt - 1)
        ev = ev_min + (i - 1) * (ev_max - ev_min) / max(1, (n_lines - 1))
        write (ny(i)%label, '(a, f0.1)') 'Energy (eV) = ', ev
        x(j) = angle
        x_lab = 'Angle'
      else
        ev = ev_min + (j - 1) * (ev_max - ev_min) / (plot_param%n_pt - 1)
        angle = angle_min + (i - 1) * (angle_max - angle_min) / max(1, (n_lines - 1))
        write (ny(i)%label, '(a, f0.1)') 'Angle = ', angle
        x(j) = ev
        x_lab = 'Energy (eV)'
      endif
      call photon_reflectivity (angle*pi/180, ev, surface, ny(i)%y_reflect(j), ny(i)%y_rel_specular(j))
    enddo

    select case (reflection_type)
    case ('total')
      ny(i)%y = ny(i)%y_reflect
    case ('specular')
      ny(i)%y = ny(i)%y_rel_specular * ny(i)%y_reflect
    case ('diffuse')
      ny(i)%y = (1 - ny(i)%y_rel_specular) * ny(i)%y_reflect
    case ('%specular')
      ny(i)%y = ny(i)%y_rel_specular
    case default
      call err_exit
    end select

  enddo

  ! plot

  y_min = ny(1)%y(1)
  y_max = ny(1)%y(1)

  do i = 1, n_lines
    y_min = min(y_min, minval(ny(i)%y))
    y_max = max(y_max, maxval(ny(i)%y))
  enddo

  call qp_calc_and_set_axis ('X', x(1), x(plot_param%n_pt), 10, 16, 'GENERAL')
  call qp_calc_and_set_axis ('Y', y_min, y_max, 6, 10, 'GENERAL')

  call qp_clear_page

  call qp_draw_graph (x, ny(1)%y, x_lab, y_lab, head_lab, .true., 0)

  do i = 1, n_lines
    call qp_draw_polyline (x, ny(i)%y, line_pattern = qp_line_pattern_name(i))
    ny(i)%line%pattern = qp_line_pattern_name(i)
  enddo

  legend = qp_legend_struct()
  legend%line_length = 40.0_rp  
  call qp_draw_curve_legend (qp_point_struct(0.5_rp, 0.0_rp, '%GRAPH/LT'), legend, ny(:)%line, text = ny(:)%label)
  write (text, '(a, i0)') 'ix_surface = ', ix_surface
  call qp_draw_text (text, 0.98_rp, 0.95_rp, '%GRAPH', 'RT')

  ! Get input:

  print *

  print '(a)', 'Surfaces Defined:'
  n_max = maxval(len_trim(sr3d_com%surface%name))
  d_max = maxval(len_trim(sr3d_com%surface%description))
  write (fmt, '(a, i0, a, i0, a)') '(3x, i3, 2x, a', n_max+5, ', a', d_max+7, ', 2a)'
  do i = 1, size (sr3d_com%surface)
    descrip = sr3d_com%surface(i)%description
    if (descrip /= '') descrip = '"' // trim(descrip) // '"'
    print fmt, i, sr3d_com%surface(i)%name, descrip, 'File: ', trim(sr3d_com%surface(i)%reflectivity_file)
  enddo
  print *
  print '(a)', 'Commands:'
  print '(a)', '   energy   <ev_min> <ev_max>         ! energies to plot at'
  print '(a)', '   angle    <angle_min> <angle_max>   ! angles to plot at'
  print '(a)', '   n_lines  <num_lines_to_draw>'
  print '(a)', '   fixed_energy <T|F>                 ! Lines of constant energy or angle?'
  print '(a)', '   ix_surface <n>                     ! Surface index'
  print '(a)', '   write                              ! Write plot points to file'
  print '(a)', '   type <total|specular|%specular|diffuse>  ! %specular = specular / total'
  call read_a_line ('Input: ', ans)
  call string_trim(ans, ans, ix)
  if (ix == 0) cycle
  call match_word (ans(1:ix), ['ix_surface  ', 'energy      ', 'angle       ', 'type        ', &
                               'n_lines     ', 'fixed_energy', 'write       '], n, matched_name = param)
  if (n < 1) then
    print *, 'CANNOT PARSE THIS.'
    cycle
  endif

  call string_trim(ans(ix+1:), ans, ix)

  select case (param)

  case ('write')
    open (10, file = 'reflection_probability.dat')

    if (fixed_energy) then
      write (10, '(a)') 'X_axis: Angle'
      do i = 1, n_lines
        ev = ev_min + (i - 1) * (ev_max - ev_min) / max(1, (n_lines - 1))
        write (10, '(i3, 2x, a, f10.1)') i, 'Energy (eV) =', ev
      enddo

    else
      write (10, '(a)') 'X_axis: Energy'
      do i = 1, n_lines
        angle = angle_min + (i - 1) * (angle_max - angle_min) / max(1, (n_lines - 1))
        write (10, '(i3, 2x, a, f10.4)') i, 'Angle =', angle
      enddo
    endif

    write (10, '(a)') '              x        y1        y2        y3'

    do i = 1, size(x)
      write (10, '(i5, 100f10.5)') i, x(i), (ny(j)%y(i), j = 1, n_lines)
    enddo

    close (10)

    print *, 'Plot file: reflection_probability.dat'

  case ('ix_surface')

    read (ans, *, iostat = ios) ix
    if (ios /= 0) then
      print *, 'BAD INTEGER'
      cycle
    endif

    if (ix < 1 .or. ix > size(sr3d_com%surface)) then
      print *, 'SURFACE INDEX OUT OF RANGE.'
      cycle
    endif

    ix_surface = ix
    surface => sr3d_com%surface(ix_surface)

  case ('type')

    call match_word (ans(1:ix), ['total    ', 'specular ', 'diffuse  ', '%specular'], n, matched_name = ans)
    if (n < 1) then
      print *, 'CANNOT PARSE THIS.'
      cycle
    endif
    select case (ans)
    case ('total')
      head_lab = 'Total Reflectivity Probability'
    case ('specular')
      head_lab = 'Specular Reflectivity Probability'
    case ('%specular')
      head_lab = 'Specular/Total Relative Reflectivity Probability'
    case ('diffuse')
      head_lab = 'Diffuse Reflectivity Probability'
    end select
    reflection_type = ans

  case ('energy', 'angle')
    
    read (ans, *, iostat = ios) value1, value2
    
    if (ans(ix+1:) == '' .or. ios /= 0) then
      print *, 'CANNOT READ VALUES'
      cycle
    endif

    if (param == 'energy') then
      ev_min = value1
      ev_max = value2
    else
      angle_min = max(0.0_rp, value1)
      angle_max = min(90.0_rp, value2)
    endif

  case ('n_lines')
    read (ans, *, iostat = ios) n_lines
    if (ans == '' .or. ios /= 0) then
      print *, 'CANNOT READ VALUE'
      cycle
    endif
    n_lines = max(1, n_lines)
    deallocate (ny)
    allocate (ny(n_lines))

  case ('fixed_energy')
    read (ans, *, iostat = ios) logic
    if (ans == '' .or. ios /= 0) then
      print *, 'CANNOT READ VALUE'
      cycle
    endif
    fixed_energy = logic

  end select
enddo

end subroutine sr3d_plot_reflection_probability

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine sr3d_plot_wall_vs_s (plot_param, lat, plane, wall_hit_file)
!
! Routine to interactively plot (x, s) .or. (y, s) section of the wall.
! Note: This routine never returns to the main program.
!
! Input:
!   plot_param    -- sr3d_plot_param_struct: Plot parameters.
!   lat           -- lat_struct: Lattice with wall.
!   plane         -- character(*): section. 'xs' or. 'ys'
!   wall_hit_file -- character(*): Photon trajectory file for plotting the trajectory.
!-

subroutine sr3d_plot_wall_vs_s (plot_param, lat, plane)

implicit none

type (sr3d_plot_param_struct) plot_param
type (sr3d_photon_track_struct), target :: photon
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch, branch2

real(rp), target :: xy_min, xy_max, s_min, s_max, r_max, x_wall, y_wall
real(rp) dummy, r1(3), r2(3)
real(rp), allocatable :: s(:), xy_in(:), xy_out(:)
real(rp), pointer :: photon_xy, wall_xy

integer i, n, ix, iw, i_chan, ios, ix_branch
integer n_phot1, n_phot2, n_hit1, n_hit2

character(*) plane
character(16) plane_str, color, line_pattern
character(40) :: ans

logical xy_user_good, s_user_good, no_wall_here, found_wall, good_wall_hit, good_photon_track
logical, allocatable :: no_wall(:), phantom(:)

! Open plotting window

branch => lat%branch(0)
photon%now%ix_branch = 0

call qp_open_page ('X', i_chan, plot_param%window_width, plot_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

xy_user_good = .false.
s_user_good = .false.
r_max = 100
n = plot_param%n_pt
allocate(s(n), xy_in(n), xy_out(n), no_wall(n), phantom(n))

phantom = .false.

if (plane == 'xs') then
  plane_str = 'X (cm)'
  photon_xy => photon%now%orb%vec(1)
  wall_xy => x_wall
elseif (plane == 'ys') then
  plane_str = 'Y (cm)'
  photon_xy => photon%now%orb%vec(3)
  wall_xy => y_wall
else
  call err_exit
endif

! Print wall info

! Command Loop

do

  ! Determine s min/max

  if (.not. s_user_good) then
    s_min = 0
    s_max = branch%ele(branch%n_ele_track)%s
  endif

  call qp_calc_and_set_axis ('X', s_min, s_max, 10, 16, 'GENERAL')

  ! find transverse min/max

  if (.not. xy_user_good) then
    xy_min = 0
    xy_max = 0

    do iw = 1, size(branch%wall3d)
      do i = 1, size(s)
        s(i) = s_min + (i - 1) * (s_max - s_min) / (size(s) - 1)

        photon%now%orb%vec = 0
        photon%now%orb%s = s(i)
        photon%now%ix_wall3d = iw

        photon_xy = -r_max
        call sr3d_find_wall_point (photon, branch, x_wall, y_wall, no_wall_here)
        if (no_wall_here) cycle

        xy_min = min(xy_min, 101*wall_xy)
        xy_max = max(xy_max, 101*wall_xy)

        photon_xy = r_max
        call sr3d_find_wall_point (photon, branch, x_wall, y_wall, no_wall_here)

        xy_min = min(xy_min, 101*wall_xy)
        xy_max = max(xy_max, 101*wall_xy)
      enddo
    enddo
  endif

  ! Get xy data points

  call qp_clear_page
  call qp_calc_and_set_axis ('Y', xy_min, xy_max, 6, 10, 'GENERAL')
  call qp_set_margin (0.07_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
  call qp_draw_graph (s, xy_in, 'S (m)', plane_str, '', .false., 0)

  do iw = 1, size(branch%wall3d)
    do i = 1, size(s)

      s(i) = s_min + (i - 1) * (s_max - s_min) / (size(s) - 1)

      photon%now%orb%vec = 0
      photon%now%orb%s = s(i)
      photon%now%ix_wall3d = iw

      photon_xy = -r_max
      call sr3d_find_wall_point (photon, branch, x_wall, y_wall, no_wall_here)
      xy_in(i) = 100 * wall_xy

      photon_xy = r_max
      call sr3d_find_wall_point (photon, branch, x_wall, y_wall, no_wall_here)
      xy_out(i) = 100 * wall_xy

      no_wall(i) = no_wall_here

      if (.not. no_wall_here) then
        phantom(i) = (branch%wall3d(iw)%section(photon%now%ix_wall_section)%surface%name == 'PHANTOM')
      endif
    enddo

    ! Now plot

    found_wall = .false.
    do i = 1, size(s)

      if (no_wall(i)) then
        if (found_wall) then
          if (phantom(i-1)) then
            color = 'red'
            line_pattern = 'dashed'
          else
            color = 'black'
            line_pattern = 'solid'
          endif
          if (i > 2) then
            call qp_draw_polyline (s(i-2:i-1), xy_in(i-2:i-1), color = color, line_pattern = line_pattern)
            call qp_draw_polyline (s(i-2:i-1), xy_out(i-2:i-1), color = color, line_pattern = line_pattern)
          endif
          call qp_draw_line(s(i-1), s(i-1), xy_in(i-1), xy_out(i-1), color = color, line_pattern = line_pattern)
        endif
        found_wall = .false.

      else  ! have a wall
        if (phantom(i)) then
          color = 'red'
          line_pattern = 'dashed'
        else
          color = 'black'
          line_pattern = 'solid'
        endif

        if (.not. found_wall .or. i == size(s)) then
          call qp_draw_line(s(i), s(i), xy_in(i), xy_out(i), color = color, line_pattern = line_pattern)
        endif
        if (found_wall) then
          call qp_draw_polyline (s(i-1:i), xy_in(i-1:i), color = color, line_pattern = line_pattern)
          call qp_draw_polyline (s(i-1:i), xy_out(i-1:i), color = color, line_pattern = line_pattern)
        endif
        found_wall = .true.
      endif
    enddo
  enddo

  ! Plot photon trajectories if trajectory file exists

  good_photon_track = .false.
  if (sr3d_params%photon_track_file /= '') inquire (file = sr3d_params%photon_track_file, exist = good_photon_track)

  good_wall_hit = .false.
  if (sr3d_params%wall_hit_file /= '') inquire (file = sr3d_params%wall_hit_file, exist = good_wall_hit)

  if (good_photon_track) then
    open (10, file = sr3d_params%photon_track_file, iostat = ios)
    read (10, *, iostat = ios) n_phot2, dummy, r2
    do 
      n_phot1 = n_phot2; r1 = r2
      read (10, *, iostat = ios) n_phot2, dummy, r2, ix_branch
      if (ios /= 0) exit
      if (n_phot2 /= n_phot1) cycle
      ! If photon has traveled more than half the lattice length assume photon has "wrapped around" and do not plot
      if (2 * norm2(r2 - r1) > lat%branch(ix_branch)%param%total_length) cycle
      ! Plot
      if (plane == 'xs') then
        call qp_draw_line(r1(3), r2(3), 100*r1(1), 100*r2(1))
      else
        call qp_draw_line(r1(3), r2(3), 100*r1(2), 100*r2(2))
      endif
    enddo
    close (10)

  elseif (good_wall_hit) then
    open (10, file = sr3d_params%wall_hit_file, iostat = ios)
    read (10, *, iostat = ios) n_phot2, n_hit2, dummy, r2
    do 
      n_phot1 = n_phot2; n_hit1 = n_hit2; r1 = r2
      read (10, *, iostat = ios) n_phot2, n_hit2, dummy, r2, ix_branch
      if (ios /= 0) exit
      if (n_phot2 /= n_phot1) cycle
      ! If photon has traveled more than half the lattice length assume photon has "wrapped around" and do not plot
      if (2 * norm2(r2 - r1) > lat%branch(ix_branch)%param%total_length) cycle
      ! Plot
      if (plane == 'xs') then
        call qp_draw_line(r1(3), r2(3), 100*r1(1), 100*r2(1))
      else
        call qp_draw_line(r1(3), r2(3), 100*r1(2), 100*r2(2))
      endif
    enddo
    close (10)
  endif

  ! Query

  print '(a)', 'Syntax: "x", "y", or "s" followed by <min> <max> values.'
  print '(a)', '        or "b" followed by branch name or index. Branch index ranges from 0 to', ubound(lat%branch, 1)
  print '(a)', '[<min> = "auto" --> autoscale] Examples: "x auto", "s 10 60"'
  call read_a_line ('Input: ', ans)

  call string_trim (ans, ans, ix)
  if (ans(1:2) == 's ') then
    call string_trim(ans(2:), ans, ix)
    if (ans == 'auto') then
      s_user_good = .false.
    else
      read (ans, *, iostat = ios) s_min, s_max
      if (ios /= 0) then
        print *, 'CANNOT DECODE MIN/MAX VALUES'
      else
        s_user_good = .true.
      endif
    endif

  elseif (ans(1:2) == 'x ' .or. ans(1:2) == 'y ') then
    call string_trim(ans(2:), ans, ix)
    if (ans == 'auto') then
      xy_user_good = .false.
    else
      read (ans, *, iostat = ios) xy_min, xy_max
      if (ios /= 0) then
        print *, 'CANNOT DECODE MIN/MAX VALUES'
      else
        xy_user_good = .true.
      endif
    endif

  elseif (ans(1:2) == 'b ') then
    call string_trim(ans(2:), ans, ix)
    branch2 => pointer_to_branch(ans, lat)
    if (.not. associated(branch2)) then
      print *, 'BAD BRANCH NAME OR INDEX'
      cycle
    endif
    branch => branch2
    photon%now%ix_branch = branch%ix_branch
    s_min = 0
    s_max = branch%ele(branch%n_ele_track)%s

  else
    print *, 'I DO NOT UNDERSTAND THIS...'
  endif

enddo

end subroutine sr3d_plot_wall_vs_s 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_plot_wall_cross_sections (plot_param, lat)
!
! Routine to interactively plot wall (x,y) cross-sections at constant s.
! Note: This routine never returns to the main program.
!
! Input:
!   plot_param  -- sr3d_plot_param_struct: Plotting parameters.
!   lat         -- Lat_struct: lattice
!-

subroutine sr3d_plot_wall_cross_sections (plot_param, lat)

implicit none

type (sr3d_plot_param_struct) plot_param
type (wall3d_section_struct), pointer :: section, sec
type (sr3d_photon_track_struct) photon
type (lat_struct), target :: lat
type (wall3d_struct), pointer :: wall3d, wall3d_select
type (wall3d_vertex_struct), pointer :: v(:)
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch, branch2

real(rp), allocatable :: x(:), y(:)
real(rp), allocatable :: x_norm(:), y_norm(:)
real(rp) s_pos, x_min, x_max, y_min, y_max, x_min_user, x_max_user, y_min_user, y_max_user
real(rp) theta, r, r1, r2, r_max, s_pos_old, rr(2), minn, maxx

integer i, n, j, iw, ix, ix_section, i_in, ios, i_chan, iu, ie_max
integer ix0, iw_wall, n_sec_max

character(100) :: ans, label, label2, sub_label
character(8) v_str

logical at_section, draw_norm, reverse_x_axis, no_wall_here, write_cross_section, set_x, set_y

! Open plotting window

branch => lat%branch(0)
photon%now%ix_branch = 0

call qp_open_page ('X', i_chan, plot_param%window_width, plot_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

write_cross_section = .false.
draw_norm = .false.
reverse_x_axis = .false.
x_min_user = real_garbage$
x_max_user = real_garbage$
y_min_user = real_garbage$
y_max_user = real_garbage$
r_max = 100
n = plot_param%n_pt
ie_max = branch%n_ele_track

allocate (x(n), y(n))
allocate (x_norm(n), y_norm(n))

ix_section = 1
iw_wall = 1
wall3d_select => branch%wall3d(iw_wall)
n_sec_max = ubound(wall3d_select%section, 1)
s_pos = wall3d_select%section(ix_section)%s
s_pos_old = s_pos
at_section = .true.

! Loop

do

  ! Find the wall cross-section at the given s value.
  ! We characterize the cross-section by an array of points with straight lines drawn between the points.
  ! The i^th point can be characterized by (r_i, theta_i) with theta_i being linear in i.
  ! This is an approximation to the true shape but it is good enough for plotting and serves as
  ! an independent check on the routines used to detect intersections of the photon with the wall.

  photon%now%orb%s = s_pos

  ! Now plot

  if (at_section) then
    write (sub_label, '(2(a, i0), 2x, 2a)') 'Wall:Section #: ', wall3d_select%ix_wall3d, ':', ix_section, &
                                                                 'Name: ', wall3d_select%section(ix_section)%name
    if (s_pos_old == s_pos) then
      write (label, '(a, f0.3, 3x, a)') 'S: ', s_pos, trim(sub_label)
                
    else
      print '(2(2x, a, f0.3), 2x, a)', 'S: ', s_pos, 'S-S_prev: ', s_pos-s_pos_old, sub_label
      write (label, '(2(2x, a, f0.3), 2x, a)') 'S: ', s_pos, 'S-S_prev: ', s_pos-s_pos_old, trim(sub_label)
    endif
    label2 = 'Surface: ' // wall3d_select%section(ix_section)%surface%name
  else
    write (label, '(a, f0.3)') 'S: ', s_pos
    ! %species used for section index.
    n = modulo(photon%now%ix_wall_section, size(wall3d_select%section)) + 1
    label2 = 'Surface: ' // wall3d_select%section(n)%surface%name
  endif

  call qp_clear_page

  x_min = 0;   x_max = 0
  y_min = 0;   y_max = 0
  do iw = 1, size(branch%wall3d)
    wall3d => branch%wall3d(iw)
    photon%now%ix_wall3d = iw
    call calc_this_outline (wall3d, x, y, no_wall_here)
    if (no_wall_here) cycle
    x_min = min (x_min, 1.01*minval(x))
    y_min = min (y_min, 1.01*minval(y))
    x_max = max (x_max, 1.01*maxval(x))
    y_max = max (y_max, 1.01*maxval(y))
  enddo

  if (x_min_user /= real_garbage$) x_min = x_min_user
  if (x_max_user /= real_garbage$) x_max = x_max_user
  if (y_min_user /= real_garbage$) y_min = y_min_user
  if (y_max_user /= real_garbage$) y_max = y_max_user

  call qp_calc_and_set_axis ('X', x_min, x_max, 10, 16, 'GENERAL')
  call qp_calc_and_set_axis ('Y', y_min, y_max, 6, 10, 'GENERAL')

  call qp_set_margin (0.07_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

  set_x = (x_min_user /= real_garbage$ .or. x_max_user /= real_garbage$)
  set_y = (y_min_user /= real_garbage$ .or. y_max_user /= real_garbage$)

  if (set_x .and. .not. set_y) then
    call qp_eliminate_xy_distortion('Y')
  elseif (.not. set_x .and. set_y) then
    call qp_eliminate_xy_distortion('X')
  elseif (.not. set_x .and. .not. set_y) then
    call qp_eliminate_xy_distortion()
  endif

  if (reverse_x_axis) then
    call qp_get_axis_attrib('X', minn, maxx)
    call qp_set_axis ('X', maxx, minn) 
  endif

  call qp_draw_graph (x, y, 'X (cm)', 'Y (cm)', label, .false., 0)

  do iw = 1, size(branch%wall3d)
    wall3d => branch%wall3d(iw)
    photon%now%ix_wall3d = iw
    call calc_this_outline (wall3d, x, y, no_wall_here)
    if (no_wall_here) cycle
    call qp_draw_polyline (x, y)

    if (write_cross_section) then
      write (iu, *) '# Sub_chamber index =', iw
      write (iu, *) '#       x           y'
      do j = 1, size(x)
        write (iu, '(2f12.6)') x(j), y(j)
      enddo
    endif

    if (draw_norm) then
      do j = 1, size(x), 4
        call qp_draw_line(x(j), x(j) + x_norm(j), y(j), y(j) + y_norm(j))
      enddo
    endif

    if (at_section .and. iw == wall3d_select%ix_wall3d) then
      sec => wall3d%section(ix_section)
      v => sec%v
      do i = 1, size(v)
        call qp_draw_symbol (100 * (v(i)%x + sec%r0(1)), 100 * (v(i)%y + sec%r0(2)))
        write (v_str, '(a, i0, a)') 'v(', i, ')'
        call qp_draw_text (v_str, 100 * (v(i)%x + sec%r0(1)), 100 * (v(i)%y + sec%r0(2)))
      enddo
      rr = sec%r0
      call qp_draw_symbol (100 * rr(1), 100 * rr(2))
      call qp_draw_text ('r0', 100 * rr(1), 100 * rr(2))
    endif
  enddo

  call qp_draw_text (label2, 0.5_rp, 0.98_rp, '%/GRAPH/LB', 'CT')

  ! What element at this s-position?

  ele => pointer_to_element_at_s (branch, s_pos, .true.)
  call qp_draw_text ('Lattice Element: ' // trim(ele%name), 0.5_rp, 0.9_rp, '%/GRAPH/LB', 'CT')

  ! Query

  print *
  print '(a)',     'Commands:'
  print '(a)',     '   iw <sub-ch>        ! Select sub-chamber to use with other commands'
  print '(a)',     '   <CR>               ! Next section of selected sub-chamber (increment viewed section index by 1).'
  print '(a)',     '   b                  ! Back section of selected sub-chamber (decrement viewed section index by 1).'
  print '(a)',     '   <Section #>        ! Index of section to view'
  print '(a)',     '   s <s-value>        ! Plot section at <s-value>.'
  print '(a)',     '   x  <x-min> <x-max> ! Set horizontal plot scale. If only one number given, set min = -number, max = number.'
  print '(a)',     '   x auto             ! Auto scale plot.'
  print '(a)',     '   y  <y-min> <y-max> ! Set horizontal plot scale. If only one number given, set min = -number, max = number.'
  print '(a)',     '   y auto             ! Auto scale plot.'
  print '(a)',     '   write              ! Write (x,y) points to a file.'
  print '(a, i0)', '   branch <name>      ! Name or index of branch. Branch indexes for this lattice range from 0 to ', ubound(lat%branch, 1)
  print '(a)',     '   normal             ! Toggle drawing of a set of vectors normal to the wall'
  print '(a)',     '   reverse            ! Toggle reversing the x-axis to point left for +x'
  print '(a)',     '   list               ! List sections for current lattice branch.'

  if (write_cross_section) then
    close (iu)
    print *
    print *, 'Written: cross_section.dat'
  endif

  call read_a_line ('Input: ', ans)

  call string_trim (ans, ans, ix)

  if (ans(1:1) == 's') then
    read (ans(2:), *, iostat = ios) s_pos
    if (ios /= 0 .or. s_pos < branch%ele(0)%s .or. s_pos > branch%ele(ie_max)%s) then
      print *, 'Cannot read s-position or s-position out of range.'
      cycle
    endif
    at_section = .false.

  elseif (ans(1:2) == 'iw') then
    read (ans(3:), *, iostat = ios) ix
    if (ios /= 0 .or. ix < 1 .or. ix > size(branch%wall3d)) then
      print *, 'Cannot read sub-chamber index or sub-chamber index out of range.'
      cycle
    endif
    wall3d_select => branch%wall3d(ix)
    n_sec_max = ubound(wall3d_select%section, 1)

    ix_section = 1
    s_pos = wall3d_select%section(1)%s
    at_section = .true.

  elseif (ans(1:1) == 'x') then
    call string_trim(ans(2:), ans, ix)
    if (ans == 'auto') then
      x_min_user = real_garbage$
      x_max_user = real_garbage$
    else
      read (ans, *, iostat = ios) r1
      if (ios /= 0) then
        print *, 'Cannot read x-scale'
        cycle
      endif
      call string_trim(ans(ix+1:), ans, ix)
      if (ix == 0) then
        x_min_user = -r1
        x_max_user = r1
      else
        read (ans, *, iostat = ios) r2
        if (ios /= 0) then
          print *, 'Cannot read x-scale'
          cycle
        endif
        x_min_user = r1
        x_max_user = r2
      endif
    endif

  elseif (ans(1:1) == 'y') then
    call string_trim(ans(2:), ans, ix)
    if (ans == 'auto') then
      y_min_user = real_garbage$
      y_max_user = real_garbage$
    else
      read (ans, *, iostat = ios) r1
      if (ios /= 0) then
        print *, 'Cannot read y-scale'
        cycle
      endif
      call string_trim(ans(ix+1:), ans, ix)
      if (ix == 0) then
        y_min_user = -r1
        y_max_user = r1
      else
        read (ans, *, iostat = ios) r2
        if (ios /= 0) then
          print *, 'Cannot read x-scale'
          cycle
        endif
        y_min_user = r1
        y_max_user = r2
      endif
    endif

  elseif (ans == '') then
    ix_section = modulo(ix_section, n_sec_max) + 1
    s_pos_old = s_pos
    s_pos = wall3d_select%section(ix_section)%s
    at_section = .true.

  elseif (index('normal', ans(1:ix)) == 1) then
    draw_norm = .not. draw_norm

  elseif (index('reverse', ans(1:ix)) == 1) then
    reverse_x_axis = .not. reverse_x_axis

  elseif (index('write', ans(1:ix)) == 1) then
    iu = lunget()
    open (iu, file = 'cross_section.dat')
    write (iu, *) '# s = ', s_pos
    write_cross_section = .true.

  elseif (ans == 'b') then
    ix_section = ix_section - 1
    if (ix_section < 1) ix_section = ix_section + n_sec_max
    s_pos_old = s_pos
    s_pos = wall3d_select%section(ix_section)%s
    at_section = .true.

  elseif (index('branch', ans(1:ix)) == 1) then
    call string_trim (ans(ix+1:), ans, ix)
    branch2 => pointer_to_branch(ans, lat)
    if (.not. associated(branch2)) then
      print *, 'BAD BRANCH NAME OR INDEX'
      cycle
    endif
    branch => branch2
    photon%now%ix_branch = branch%ix_branch
    ix_section = 1
    iw_wall = 1
    wall3d_select => branch%wall3d(iw_wall)
    n_sec_max = ubound(wall3d_select%section, 1)
    s_pos = wall3d_select%section(ix_section)%s
    s_pos_old = s_pos
    at_section = .true.

  elseif (index('list', ans(1:ix)) == 1) then
    do iw = 1, size(branch%wall3d)
      wall3d => branch%wall3d(iw)
      print '(a, i6, 2x, a)', 'Sub-section:', iw, wall3d%name
      do i = 1, min(1000, ubound(wall3d%section, 1))
        section => wall3d%section(i)
        if (associated(section%surface)) then
          print '(2i8, f14.6, 2x, a20, 2x, a16, a)', iw, i, section%s, section%name, &
              wall3d_section_type_name(section%type), trim(section%surface%name)
        else
          print '(2i8, f14.6, 2x, a20, 2x, a16)', iw, i, section%s, section%name, wall3d_section_type_name(section%type)
        endif
      enddo
    enddo

  else
    read (ans, *, iostat = ios) i_in
    if (ios /= 0) then
      print *, 'Cannot read section index number'
      cycle
    endif
    if (i_in < 0 .or. i_in > n_sec_max) then
      print '(a, i0, a)', 'Number is out of range! (maximum = ', n_sec_max, ')'
      cycle
    endif
    ix_section = i_in
    s_pos_old = s_pos
    s_pos = wall3d_select%section(ix_section)%s
    at_section = .true.

  endif

enddo

!-------------------------------------------------------------------------
contains

subroutine calc_this_outline (wall3d, x, y, no_wall_here)

type (wall3d_struct) wall3d
real(rp) :: x(:), y(:)
real(rp) theta, dw_perp(3)
integer i, j
logical no_wall_here

! Idea is to see where photon path from photon%old (at the wall origin) to %now intersects the wall.
! photon%now is at 1 meter radius which is assumed to be outside the wall.

do i = 1, size(x)

  theta = (i-1) * twopi / (size(x) - 1)
  photon%now%orb%vec(1) = r_max * cos(theta)  ! Relative to wall origin.
  photon%now%orb%vec(3) = r_max * sin(theta)  ! Relative to wall origin.

  call sr3d_find_wall_point (photon, branch, x(i), y(i), no_wall_here, dw_perp)
  x_norm(i) = dw_perp(1)
  y_norm(i) = dw_perp(2)

  if (no_wall_here) return

enddo

x = x * 100; y = y * 100

end subroutine calc_this_outline

end subroutine sr3d_plot_wall_cross_sections

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_find_wall_point (photon, branch, x_wall, y_wall, no_wall_here, dw_perp)
!
! Routine to find the point on the specified sub-chamber wall where the half-line
! drawn from the sub-chamber origin point through the (x, y) photon point intersects
! the wall.
!
! Input:
!   photon    -- sr3d_photon_track_struct
!     %now%vec(1)     -- x-position relative to wall origin
!     %now%vec(3)     -- y-position relative to wall origin
!     %now%s          -- s-position
!     %now%ix_wall3d  -- sub-chamber wall index.
!     %now%ix_wall_section -- Section of wall index.
!
! Output:
!   x_wall        -- real(rp): x wall point. Not relative to wall origin. Return zero if no wall here..
!   y_wall        -- real(rp): y wall point. Not relative to wall origin. Return zero if no wall here.
!   no_wall_here  -- logical: Set True if the subchamber does not extend longitudinally to the given s-position.
!   dw_perp(3)    -- real(rp), optional: Vector which is outward perpendicular to the wall.
!-

subroutine sr3d_find_wall_point (photon, branch, x_wall, y_wall, no_wall_here, dw_perp)

implicit none

type (sr3d_photon_track_struct) photon
type (branch_struct) branch

real(rp) x_wall, y_wall
real(rp), optional :: dw_perp(3)
real(rp) r_wall, r_part, origin(3)

logical no_wall_here ! No wall at this s-position?

!

x_wall = 0
y_wall = 0

photon%now%orb%ix_ele = element_at_s (branch%lat, photon%now%orb%s, .true., branch%ix_branch)
call sr3d_photon_d_radius (photon%now, branch, no_wall_here, dw_perp, origin)
if (no_wall_here) return

if (photon%now%d_radius < 0) then
  print *, 'INTERNAL COMPUTATION ERROR!'
  call err_exit
endif

r_part = sqrt((photon%now%orb%vec(1) - origin(1))**2 + (photon%now%orb%vec(3) - origin(2))**2)
r_wall = r_part - photon%now%d_radius

x_wall = origin(1) + (photon%now%orb%vec(1) - origin(1)) * r_wall / r_part
y_wall = origin(2) + (photon%now%orb%vec(3) - origin(2)) * r_wall / r_part

end subroutine sr3d_find_wall_point

end module
