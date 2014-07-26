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

subroutine sr3d_plot_reflection_probability (plot_param, wall, sr3d_com)

implicit none

type yplot
  real(rp), allocatable :: y_reflect(:)
  real(rp), allocatable :: y_rel_specular(:)
  real(rp), allocatable :: y(:)
  character(40) label
  type (qp_line_struct) line
end type

type (sr3d_plot_param_struct) plot_param
type (sr3d_wall_struct), target :: wall
type (yplot), allocatable :: ny(:)
type (photon_reflect_surface_struct), pointer :: surface
type (sr3d_common_struct), target :: sr3d_com

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
reflection_type = 'total'
head_lab = 'Reflectivity (Total)'

surface => sr3d_com%surface(1)

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
    call qp_draw_polyline (x, ny(i)%y, line_pattern = i)
    ny(i)%line%pattern = i
  enddo

  call qp_draw_curve_legend (0.5_rp, 0.0_rp, '%GRAPH/LT', ny(:)%line, 40.0_rp, text = ny(:)%label, text_offset = 10.0_rp)

  ! Get input:

  print *
  print '(a)', 'Surfaces Defined:'
  do i = 1, size (sr3d_com%surface)
    print '(3x, i3, 2x, a)', i, trim(sr3d_com%surface(i)%descrip)
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

    surface => sr3d_com%surface(ix)

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
!   plane      -- Character(*): section. 'xs' or. 'ys'
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

if (plane == 'xs') then
  plane_str = 'X (cm)'
  photon_xy => photon%now%vec(1)
  wall_xy => x_wall
elseif (plane == 'ys') then
  plane_str = 'Y (cm)'
  photon_xy => photon%now%vec(3)
  wall_xy => y_wall
else
  call err_exit
endif

! Print wall info

do i = 0, wall%n_section_max
  print '(i4, 2x, a, f12.2)', i, wall%section(i)%basic_shape(1:16), wall%section(i)%s
enddo

! Loop

do

  ! Determine s min/max

  if (.not. s_user_good) then
    s_min = wall%section(0)%s
    s_max = wall%section(wall%n_section_max)%s
  endif

  call qp_calc_and_set_axis ('X', s_min, s_max, 10, 16, 'GENERAL')

  ! Get xy data points

  do i = 1, size(s)

    s(i) = s_min + (i - 1) * (s_max - s_min) / (size(s) - 1)

    photon%now%vec    = 0
    photon%now%s = s(i)

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
! Subroutine sr3d_plot_wall_cross_sections (plot_param, wall, lat)
!
! Routine to interactively plot wall (x,y) cross-sections at constant s.
! Note: This routine never returns to the main program.
!
! Input:
!   wall    -- sr3d_wall_struct: Wall structure.
!   lat     -- Lat_struct: lattice
!-

subroutine sr3d_plot_wall_cross_sections (plot_param, wall, lat)

implicit none

type (sr3d_plot_param_struct) plot_param
type (sr3d_wall_struct), target :: wall
type (sr3d_wall_section_struct), pointer :: section
type (sr3d_photon_track_struct) photon
type (lat_struct) lat
type (wall3d_vertex_struct), pointer :: v(:)

real(rp), allocatable :: x(:), y(:)
real(rp) s_pos, x_max, y_max, theta, r, x_max_user, r_max, s_pos_old
real(rp), allocatable :: x1_norm(:), y1_norm(:), x2_norm(:), y2_norm(:)
real(rp) minn, maxx

integer i, j, ix, ix_section, i_in, ios, i_chan, n, iu, n_norm_max, ix0

character(100) :: ans, label, label2
character(8) v_str

logical at_section, draw_norm, reverse_x_axis
logical, allocatable :: in_ante(:)

! Open plotting window

call qp_open_page ('X', i_chan, plot_param%window_width, plot_param%window_height, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

draw_norm = .false.
reverse_x_axis = .false.
x_max_user = -1
r_max = 100
n = plot_param%n_pt
allocate (x(n), y(n), in_ante(n))
allocate (x1_norm(n), y1_norm(n), x2_norm(n), y2_norm(n))

! Print wall info

do i = 0, min(1000, wall%n_section_max)
  section => wall%section(i)
  print '(i8, 2x, a, f14.6, 2x, a)', i, section%basic_shape(1:12), section%s, section%name
enddo

ix_section = 0
s_pos = wall%section(ix_section)%s
s_pos_old = s_pos
at_section = .true.

! Loop

do

  ! Find the wall cross-section at the given s value.
  ! We characterize the cross-section by an array of points with straight lines drawn between the points.
  ! The i^th point can be characterized by (r_i, theta_i) with theta_i being linear in i.
  ! This is an approximation to the true shape but it is good enough for plotting and serves as
  ! an independent check on the routines used to detect intersections of the photon with the wall.

  photon%now%s = s_pos

  do i = 1, size(x)

    ! Idea is to see where photon path from photon%old (at the origin) to %now intersects the wall.
    ! photon%now is at 1 meter radius which is assumed to be outside the wall.

    theta = (i-1) * twopi / (size(x) - 1)
    photon%now%vec(1) = r_max * cos(theta)  
    photon%now%vec(3) = r_max * sin(theta)

    if (draw_norm .and. modulo(i, 4) == 0) then
      j = (i / 4)
      call sr3d_find_wall_point (wall, lat, photon, x(i), y(i), x1_norm(j), x2_norm(j), y1_norm(j), y2_norm(j), in_ante = in_ante(i))
      n_norm_max = j
    else
      call sr3d_find_wall_point (wall, lat, photon, x(i), y(i), in_ante = in_ante(i))
    endif

  enddo

  x = x * 100; y = y * 100

  ! Now plot

  if (at_section) then
    if (s_pos_old == s_pos) then
      write (label, '(a, f0.3, a, i0, 2a)') 'S: ', s_pos, '   Section #: ', ix_section, '  Name: ', wall%section(ix_section)%name
    else
      print '(2(a, f0.3), a, i0, 2a)', 'S: ', s_pos, '  dS: ', s_pos-s_pos_old, &
                                '   Section #: ', ix_section, '  Name: ', wall%section(ix_section)%name
      write (label, '(2(a, f0.3), a, i0, 2a)') 'S: ', s_pos, '  dS: ', s_pos-s_pos_old, &
                                '   Section #: ', ix_section, '  Name: ', wall%section(ix_section)%name
    endif
    label2 = 'Surface: ' // wall%section(ix_section)%surface%descrip
  else
    write (label, '(a, f0.3)') 'S: ', s_pos
    label2 = 'Surface: ' // wall%section(photon%now%ix_wall+1)%surface%descrip
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

  if (reverse_x_axis) then
    call qp_get_axis_attrib('X', minn, maxx)
    call qp_set_axis ('X', maxx, minn) 
  endif

  call qp_draw_graph (x, y, 'X (cm)', 'Y (cm)', label, .false., 0)

  ix0 = 1
  do ix = 2, size(x)
    if (ix /= size(x)) then
      if (in_ante(ix+1) == in_ante(ix0)) cycle
    endif

    if (in_ante(ix0)) then
      call qp_draw_polyline (x(ix0:ix), y(ix0:ix), color = red$)
    else
      call qp_draw_polyline (x(ix0:ix), y(ix0:ix))
    endif
    ix0 = ix+1
  enddo

  call qp_draw_text (label2, 0.5_rp, 0.98_rp, '%/GRAPH/LB', 'CT')

  if (draw_norm) then
    do j = 1, n_norm_max
      call qp_draw_line(100*x1_norm(j), 100*x2_norm(j), 100*y1_norm(j), 100*y2_norm(j))
    enddo
  endif

  if (at_section .and. wall%section(ix_section)%basic_shape == "gen_shape") then
    v => wall%section(ix_section)%gen_shape%wall3d_section%v
    do i = 1, size(v)
      call qp_draw_symbol (100 * v(i)%x, 100 * v(i)%y)
      write (v_str, '(a, i0, a)') 'v(', i, ')'
      call qp_draw_text (v_str, 100 * v(i)%x, 100 * v(i)%y)
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
  print '(a)', '   normal           ! Toggle drawing of a set of vectors normal to the wall'
  print '(a)', '   reverse          ! Toggle reversing the x-axis to point left for +x'

  call read_a_line ('Input: ', ans)

  call string_trim (ans, ans, ix)

  if (ans(1:1) == 's') then
    read (ans(2:), *, iostat = ios) s_pos
    if (ios /= 0 .or. s_pos < wall%section(0)%s .or. s_pos > wall%section(wall%n_section_max)%s) then
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
    ix_section = modulo(ix_section + 1, wall%n_section_max + 1)
    s_pos_old = s_pos
    s_pos = wall%section(ix_section)%s
    at_section = .true.

  elseif (index('normal', ans(1:ix)) == 1) then
    draw_norm = .not. draw_norm

  elseif (index('reverse', ans(1:ix)) == 1) then
    reverse_x_axis = .not. reverse_x_axis

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
    ix_section = modulo(ix_section - 1, wall%n_section_max + 1)
    s_pos_old = s_pos
    s_pos = wall%section(ix_section)%s
    at_section = .true.

  else
    read (ans, *, iostat = ios) i_in
    if (ios /= 0) then
      print *, 'Cannot read section index number'
      cycle
    endif
    if (i_in < 0 .or. i_in > wall%n_section_max) then
      print *, 'Number is out of range!'
      cycle
    endif
    ix_section = i_in
    s_pos_old = s_pos
    s_pos = wall%section(ix_section)%s
    at_section = .true.

  endif

enddo

end subroutine sr3d_plot_wall_cross_sections

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine sr3d_find_wall_point (wall, lat, photon, x_wall, y_wall, x1_norm, x2_norm, y1_norm, y2_norm, in_ante)

implicit none

type (sr3d_wall_struct), target :: wall
type (sr3d_photon_track_struct) photon
type (lat_struct) lat

real(rp) x_wall, y_wall
real(rp), optional :: x1_norm, x2_norm, y1_norm, y2_norm

real(rp) tri_vert0(3), tri_vert1(3), tri_vert2(3)
real(rp) dtrack, d_radius, r_old, dw_perp(3)

integer j, ixp

logical, optional :: in_ante
logical is_through

!

call sr3d_get_wall_index (photon%now, wall, ixp)

photon%old%vec = 0
photon%old%s = photon%now%s
r_old = sqrt(photon%now%vec(1)**2 + photon%now%vec(3)**2)
photon%now%track_len = photon%old%track_len + r_old

if (wall%section(ixp+1)%basic_shape == 'triangular') then
  do j = 1, 2*size(wall%section(ixp)%gen_shape%wall3d_section%v)
    call sr3d_get_mesh_wall_triangle_pts (wall%section(ixp), wall%section(ixp+1), j, tri_vert0, tri_vert1, tri_vert2)
    call sr3d_mesh_triangle_intersect (photon, tri_vert0, tri_vert1, tri_vert2, is_through, dtrack)
    if (is_through) then
      x_wall = dtrack * photon%now%vec(1) / r_old
      y_wall = dtrack * photon%now%vec(3) / r_old
      photon%now%ix_triangle = j
      exit
    endif
  enddo
  if (.not. is_through) then  ! did not find intersection with meshes
    print *, 'INTERNAL COMPUTATION ERROR!'
    call err_exit
  endif

else
  call sr3d_photon_d_radius (photon%now, wall, d_radius, in_antechamber = in_ante)
  if (d_radius < 0) then
    print *, 'INTERNAL COMPUTATION ERROR!'
    call err_exit
  endif
  x_wall = (r_old - d_radius) * photon%now%vec(1) / r_old
  y_wall = (r_old - d_radius) * photon%now%vec(3) / r_old
endif

! The length of the normal vector is 1 cm.

if (present(x1_norm)) then
  photon%now%vec(1) = x_wall; photon%now%vec(3) = y_wall
  call sr3d_photon_d_radius (photon%now, wall, d_radius, lat, dw_perp, in_ante)
  x1_norm = x_wall;                  y1_norm = y_wall
  x2_norm = x_wall + dw_perp(1)/100; y2_norm = y_wall + dw_perp(2)/100
endif

end subroutine sr3d_find_wall_point

end module
