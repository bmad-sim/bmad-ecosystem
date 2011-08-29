module synrad3d_plot_mod

use synrad3d_track_mod
use quick_plot
use input_mod

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_plot_reflection_probability (plot_param)
!
! Routine to plot reflection probability curves
!
! Input;
!   plot_param -- sr3d_plot_param_struct: Plot parameters.
!-

subroutine sr3d_plot_reflection_probability (plot_param)

use photon_reflection_mod

implicit none

type yplot
  real(rp), allocatable :: y_rough(:)
  real(rp), allocatable :: y_smooth(:)
  real(rp), allocatable :: y(:)
  character(40) label
  type (qp_line_struct) line
end type

type (sr3d_plot_param_struct) plot_param
type (yplot), allocatable :: ny(:)

real(rp), target :: angle_min, angle_max, ev_min, ev_max, angle, ev
real(rp) value, value1, value2, y_min, y_max
real(rp), allocatable :: x(:)

integer i, j, n, ix, ios, n_lines, i_chan

logical fixed_energy, logic

character(80) ans, head_lab
character(16) x_lab, y_lab, param, reflection_type

! init

angle_min = 0
angle_max = 40

ev_min = 50
ev_max = 150

fixed_energy = .true.
y_lab = 'Reflectivity'
reflection_type = 'rough'
head_lab = 'Reflectivity (ROUGH)'

n_lines = 3
allocate (ny(n_lines))

do i = 1, n_lines
  allocate (ny(i)%y_rough(plot_param%n_pt))
  allocate (ny(i)%y_smooth(plot_param%n_pt))
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
      call photon_reflectivity (angle*pi/180, ev, ny(i)%y_rough(j), ny(i)%y_smooth(j))
    enddo

    select case (reflection_type)
    case ('rough')
      ny(i)%y = ny(i)%y_rough
    case ('smooth')
      ny(i)%y = ny(i)%y_smooth
    case ('diff')
      ny(i)%y = ny(i)%y_smooth - ny(i)%y_rough
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
    call qp_draw_polyline (x, ny(i)%y, line_pattern = i)
    ny(i)%line%pattern = i
  enddo

  call qp_draw_curve_legend (0.5_rp, 0.0_rp, '%GRAPH/LT', ny(:)%line, 40.0_rp, text = ny(:)%label, text_offset = 10.0_rp)

  ! Get input:

  print *
  print '(a)', 'Commands:'
  print '(a)', '   energy   <ev_min> <ev_max>         ! energies to plot at'
  print '(a)', '   angle    <angle_min> <angle_max>   ! angles to plot at'
  print '(a)', '   n_lines  <num_lines_to_draw>'
  print '(a)', '   fixed_energy <t|f>  '
  print '(a)', '   type <rough|smooth|diff>           ! diff = smooth - rough'
  call read_a_line ('Input: ', ans)
  call string_trim(ans, ans, ix)
  if (ix == 0) cycle
  call match_word (ans(1:ix), ['energy      ', 'angle       ', 'type     ', &
                               'n_lines     ', 'fixed_energy'], n, matched_name = param)
  if (n < 1) then
    print *, 'CANNOT PARSE THIS.'
    cycle
  endif

  call string_trim(ans(ix+1:), ans, ix)

  select case (param)

  case ('type')

    call match_word (ans(1:ix), ['rough ', 'smooth', 'diff  '], n, matched_name = ans)
    if (n < 1) then
      print *, 'CANNOT PARSE THIS.'
      cycle
    endif
    select case (ans)
    case ('rough')
      head_lab = 'Reflectivity (ROUGH)'
    case ('smooth')
      head_lab = 'Reflectivity (SMOOTH)'
    case ('diff')
      head_lab = 'Reflectivity (SMOOTH - ROUGH)'
    end select
    reflection_type = ans

  case ('energy', 'angle')
    
    read (ans, *, iostat = ios) value1, value2
    
    if (ans(ix+1:) == '' .or. ios /= 0) then
      print *, 'CANNOT READ VALUES'
      cycle
    endif

    if (param == 'energy') then
      ev_min = max(value1, ev_min)
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
! subroutine sr3d_plot_wall_vs_s (plot_param, wall, lat, plane)
!
! Routine to interactively plot (x, s) .or. (y, s) section of the wall.
! Note: This routine never returns to the main program.
!
! Input:
!   plot_param -- sr3d_plot_param_struct: Plot parameters.
!   wall       -- sr3d_wall_struct: Wall structure.
!   lat        -- Lat_struct: lattice
!   plane      -- Character(*): section. 'x' or. 'y'
!-

subroutine sr3d_plot_wall_vs_s (plot_param, wall, lat, plane)

implicit none

type (sr3d_plot_param_struct) plot_param
type (sr3d_wall_struct), target :: wall
type (sr3d_photon_track_struct), target :: photon
type (lat_struct) lat

real(rp), target :: xy_min, xy_max, s_min, s_max, r_max, x_wall, y_wall
real(rp), allocatable :: s(:), xy_in(:), xy_out(:)
real(rp), pointer :: photon_xy, wall_xy

integer i, ix, i_chan, ios

character(*) plane
character(16) plane_str
character(40) :: ans

logical xy_user_good, s_user_good

! Open plotting window

call qp_open_page ('X', i_chan, plot_param%window_width, plot_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

xy_user_good = .false.
s_user_good = .false.
r_max = 100
allocate(s(plot_param%n_pt), xy_in(plot_param%n_pt), xy_out(plot_param%n_pt))

if (plane == 'x') then
  plane_str = 'X (cm)'
  photon_xy => photon%now%vec(1)
  wall_xy => x_wall
elseif (plane == 'y') then
  plane_str = 'Y (cm)'
  photon_xy => photon%now%vec(3)
  wall_xy => y_wall
else
  call err_exit
endif

! Print wall info

do i = 0, wall%n_pt_max
  print '(i4, 2x, a, f12.2)', i, wall%pt(i)%basic_shape(1:12), wall%pt(i)%s
enddo

! Loop

do

  ! Determine s min/max

  if (.not. s_user_good) then
    s_min = wall%pt(0)%s
    s_max = wall%pt(wall%n_pt_max)%s
  endif

  call qp_calc_and_set_axis ('X', s_min, s_max, 10, 16, 'GENERAL')

  ! Get xy data points

  do i = 1, size(s)

    s(i) = s_min + (i - 1) * (s_max - s_min) / (size(s) - 1)

    photon%now%vec    = 0
    photon%now%vec(5) = s(i)

    photon_xy = -r_max
    call sr3d_find_wall_point (wall, lat, photon, x_wall, y_wall)
    xy_in(i) = wall_xy

    photon_xy = r_max
    call sr3d_find_wall_point (wall, lat, photon, x_wall, y_wall)
    xy_out(i) = wall_xy

  enddo

  xy_in = xy_in * 100; xy_out = xy_out * 100

  ! Now plot

  call qp_clear_page
  if (.not. xy_user_good) then
    xy_min = 1.01 * minval(xy_in)
    xy_max = 1.01 * maxval(xy_out)
  endif

  call qp_calc_and_set_axis ('Y', xy_min, xy_max, 6, 10, 'GENERAL')
  call qp_set_margin (0.07_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
  call qp_draw_graph (s, xy_in, 'S (m)', plane_str, '', .true., 0)
  call qp_draw_polyline (s, xy_out)

  ! Query

  print *, 'Syntax: "x", "y", or "s" followed by <min> <max> values.'
  print *, '[<min> = "auto" --> autoscale] Example: "x auto", "s 10 60"'
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

  else
    print *, 'I DO NOT UNDERSTAND THIS...'
  endif

enddo

end subroutine sr3d_plot_wall_vs_s 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_plot_wall_cross_sections (plot_param, wall, lat, extra)
!
! Routine to interactively plot wall (x,y) cross-sections at constant s.
! Note: This routine never returns to the main program.
!
! Input:
!   wall    -- sr3d_wall_struct: Wall structure.
!   lat     -- Lat_struct: lattice
!   extra   -- Character(*): 
!                   '-norm' -> Also plot a set of lines normal to the wall.
!                             This is used as a check to the wall normal calculation.
!                   '-rcross' -> flip x-axis
!-

subroutine sr3d_plot_wall_cross_sections (plot_param, wall, lat, extra)

implicit none

type (sr3d_plot_param_struct) plot_param
type (sr3d_wall_struct), target :: wall
type (sr3d_wall_pt_struct), pointer :: pt
type (sr3d_photon_track_struct) photon
type (lat_struct) lat

real(rp), allocatable :: x(:), y(:)
real(rp) s_pos, x_max, y_max, theta, r, x_max_user, r_max, s_pos_old
real(rp), allocatable :: x1_norm(:), y1_norm(:), x2_norm(:), y2_norm(:)
real(rp) minn, maxx

integer i, j, ix, ix_section, i_in, ios, i_chan, n, iu

character(100) :: ans, label
character(*) extra

logical at_section

! Open plotting window

call qp_open_page ('X', i_chan, plot_param%window_width, plot_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

x_max_user = -1
r_max = 100
n = plot_param%n_pt
allocate (x(n), y(n))
allocate (x1_norm(n), y1_norm(n), x2_norm(n), y2_norm(n))

! Print wall info

do i = 0, wall%n_pt_max
  pt => wall%pt(i)
  print '(i4, 2x, a, f12.2)', i, pt%basic_shape(1:12), pt%s
enddo

ix_section = 0
s_pos = wall%pt(ix_section)%s
s_pos_old = s_pos
at_section = .true.

! Loop

do

  ! Find the wall cross-section at the given s value.
  ! We characterize the cross-section by an array of points with straight lines drawn between the points.
  ! The i^th point can be characterized by (r_i, theta_i) with theta_i being linear in i.
  ! This is an approximation to the true shape but it is good enough for plotting and serves as
  ! an independent check on the routines used to detect intersections of the photon with the wall.

  photon%now%vec(5) = s_pos

  do i = 1, size(x)

    ! Idea is to see where photon path from photon%old (at the origin) to %now intersects the wall.
    ! photon%now is at 1 meter radius which is assumed to be outside the wall.

    theta = (i-1) * twopi / (size(x) - 1)
    photon%now%vec(1) = r_max * cos(theta)  
    photon%now%vec(3) = r_max * sin(theta)

    if (extra == '-norm' .and. modulo(i, 4) == 0) then
      j = (i / 4)
      call sr3d_find_wall_point (wall, lat, photon, x(i), y(i), x1_norm(j), x2_norm(j), y1_norm(j), y2_norm(j))
    else
      call sr3d_find_wall_point (wall, lat, photon, x(i), y(i))
    endif

  enddo

  x = x * 100; y = y * 100
  x1_norm = x1_norm * 100; x2_norm = x2_norm * 100
  y1_norm = y1_norm * 100; y2_norm = y2_norm * 100

  ! Now plot

  if (at_section) then
    if (s_pos_old == s_pos) then
      write (label, '(a, f0.3, a, i0, 2a)') 'S: ', s_pos, '   Section #: ', ix_section, '  Name: ', wall%pt(ix_section)%name
    else
      print '(2(a, f0.3), a, i0, 2a)', 'S: ', s_pos, '  dS: ', s_pos-s_pos_old, &
                                '   Section #: ', ix_section, '  Name: ', wall%pt(ix_section)%name
      write (label, '(2(a, f0.3), a, i0, 2a)') 'S: ', s_pos, '  dS: ', s_pos-s_pos_old, &
                                '   Section #: ', ix_section, '  Name: ', wall%pt(ix_section)%name
    endif
  else
    write (label, '(a, f0.3)') 'S: ', s_pos
  endif
  call qp_clear_page
  x_max = 1.01 * maxval(abs(x)); y_max = 1.01 * maxval(abs(y))
  if (x_max_user > 0) x_max = x_max_user
  call qp_calc_and_set_axis ('X', -x_max, x_max, 10, 16, 'ZERO_SYMMETRIC')
  call qp_calc_and_set_axis ('Y', -y_max, y_max, 6, 10, 'ZERO_SYMMETRIC')
  call qp_set_margin (0.07_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

  if (x_max_user > 0) then
    call qp_eliminate_xy_distortion('Y')
  else
    call qp_eliminate_xy_distortion()
  endif

  if (extra == '-rcross') then
    call qp_get_axis_attrib('X', minn, maxx)
    call qp_set_axis ('X', maxx, minn) 
  endif

  call qp_draw_graph (x, y, 'X', 'Y', label, .true., 0)

  if (extra == '-norm') then
    do j = 1, size(x1_norm)
      call qp_draw_line(x1_norm(j), x2_norm(j), y1_norm(j), y2_norm(j))
    enddo
  endif

  ! Query
  print *
  print '(a)', 'Commands:'
  print '(a)', '   <CR>             ! Next section (increment viewed section index by 1).'
  print '(a)', '   b                ! Back section (decrement viewed section index by 1).'
  print '(a)', '   <Section #>      ! Index of section to view'
  print '(a)', '   s <s-value>      ! Plot section at <s-value>.'
  print '(a)', '   x <x-max>        ! Set horizontal plot scale. Vertical will be scaled to match.'
  print '(a)', '   x auto           ! Auto scale plot.'
  print '(a)', '   write            ! Write (x,y) points to a file.'
  call read_a_line ('Input: ', ans)

  call string_trim (ans, ans, ix)

  if (ans(1:1) == 's') then
    read (ans(2:), *, iostat = ios) s_pos
    if (ios /= 0 .or. s_pos < wall%pt(0)%s .or. s_pos > wall%pt(wall%n_pt_max)%s) then
      print *, 'Cannot read s-position or s-position out of range.'
      cycle
    endif
    at_section = .false.

  elseif (ans(1:1) == 'x') then
    call string_trim(ans(2:), ans, ix)
    if (ans == 'auto') then
      x_max_user = -1
    else
      read (ans, *, iostat = ios) r
      if (ios /= 0) then
        print *, 'Cannot read x-scale'
        cycle
      endif
      x_max_user = r
    endif

  elseif (ans == '') then
    ix_section = modulo(ix_section + 1, wall%n_pt_max + 1)
    s_pos_old = s_pos
    s_pos = wall%pt(ix_section)%s
    at_section = .true.

  elseif (index('write', ans(1:ix)) == 1) then
    iu = lunget()
    open (iu, file = 'cross_section.dat')
    write (iu, *) '  s = ', s_pos
    write (iu, *) '        x           y'
    do j = 1, size(x)
      write (iu, '(2f12.6)') x(j), y(j)
    enddo
    close (iu)
    print *, 'Writen: cross_section.dat'

  elseif (ans == 'b') then
    ix_section = modulo(ix_section - 1, wall%n_pt_max + 1)
    s_pos_old = s_pos
    s_pos = wall%pt(ix_section)%s
    at_section = .true.

  else
    read (ans, *, iostat = ios) i_in
    if (ios /= 0) then
      print *, 'Cannot read section index number'
      cycle
    endif
    if (i_in < 0 .or. i_in > wall%n_pt_max) then
      print *, 'Number is out of range!'
      cycle
    endif
    ix_section = i_in
    s_pos_old = s_pos
    s_pos = wall%pt(ix_section)%s
    at_section = .true.

  endif

enddo

end subroutine sr3d_plot_wall_cross_sections

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine sr3d_find_wall_point (wall, lat, photon, x_wall, y_wall, x1_norm, x2_norm, y1_norm, y2_norm)

implicit none

type (sr3d_wall_struct), target :: wall
type (sr3d_photon_track_struct) photon
type (lat_struct) lat

real(rp) x_wall, y_wall
real(rp), optional :: x1_norm, x2_norm, y1_norm, y2_norm

real(rp) tri_vert0(3), tri_vert1(3), tri_vert2(3)
real(rp) dtrack, d_radius, r_old, dw_perp(3)

integer j, ixp

logical is_through

!

call sr3d_get_wall_index (photon%now, wall, ixp)

photon%old%vec = 0
photon%old%vec(5) = photon%now%vec(5)
r_old = sqrt(photon%now%vec(1)**2 + photon%now%vec(3)**2)

if (wall%pt(ixp+1)%basic_shape == 'gen_shape_mesh') then
  do j = 1, 2*size(wall%pt(ixp)%gen_shape%v)
    call sr3d_get_mesh_wall_triangle_pts (wall%pt(ixp), wall%pt(ixp+1), j, tri_vert0, tri_vert1, tri_vert2)
    call sr3d_mesh_triangle_intersect (photon, tri_vert0, tri_vert1, tri_vert2, is_through, dtrack)
    if (is_through) then
      x_wall = dtrack * photon%now%vec(1) / r_old
      y_wall = dtrack * photon%now%vec(3) / r_old
      exit
    endif
  enddo
  if (.not. is_through) then  ! did not find intersection with meshes
    print *, 'INTERNAL COMPUTATION ERROR!'
    call err_exit
  endif

else
  call sr3d_photon_d_radius (photon%now, wall, d_radius)
  if (d_radius < 0) then
    print *, 'INTERNAL COMPUTATION ERROR!'
    call err_exit
  endif
  x_wall = (r_old - d_radius) * photon%now%vec(1) / r_old
  y_wall = (r_old - d_radius) * photon%now%vec(3) / r_old

  ! The length of the normal vector is 1 cm.

  if (present(x1_norm)) then
    photon%now%vec(1) = x_wall; photon%now%vec(3) = y_wall
    call sr3d_photon_d_radius (photon%now, wall, d_radius, lat, dw_perp)
    x1_norm = x_wall;                  y1_norm = y_wall
    x2_norm = x_wall + dw_perp(1)/100; y2_norm = y_wall + dw_perp(2)/100
  endif

endif

end subroutine sr3d_find_wall_point

end module
