module synrad3d_plot_mod

use synrad3d_track_mod
use quick_plot
use input_mod

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_plot_wall_vs_s (wall, plane)
!
! Routine to interactively plot (x, s) .or. (y, s) section of the wall.
! Note: This routine never returns to the main program.
!
! Input:
!   wall    -- sr3d_wall_struct: Wall structure.
!   plane   -- Character(*): section. 'x' or. 'y'
!-

subroutine sr3d_plot_wall_vs_s (wall, plane)

implicit none

type (sr3d_wall_struct), target :: wall
type (sr3d_photon_track_struct), target :: photon

real(rp), target :: xy_min, xy_max, s_min, s_max, r_max, x_wall, y_wall
real(rp) s(200), xy_in(200), xy_out(200)
real(rp), pointer :: photon_xy, wall_xy

integer i, ix, i_chan, ios

character(*) plane
character(16) plane_str
character(40) :: ans = 'first'

logical xy_user_good, s_user_good

! Open plotting window

call qp_open_page ('X', i_chan, 800.0_rp, 400.0_rp, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

xy_user_good = .false.
s_user_good = .false.
r_max = 100

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

  ! Query

  if (ans /= 'first') then
    print *, 'Syntax: "x", "y", or "s" followed by <min> <max> values.'
    print *, '[<min> = "auto" --> autoscale] Example: "x auto", "s 10 60"'
    call read_a_line ('Input: ', ans)
  endif

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

  elseif (ans == 'first') then
    ans = ''

  else
    print *, 'I DO NOT UNDERSTAND THIS...'
  endif

  ! Determine s min/max

  if (.not. s_user_good) then
    s_min = wall%pt(0)%s
    s_max = wall%pt(wall%n_pt_max)%s
  endif

  call qp_calc_and_set_axis ('X', s_min, s_max, 6, 10, 'GENERAL')

  ! Get xy data points

  do i = 1, size(s)

    s(i) = s_min + (i - 1) * (s_max - s_min) / (size(s) - 1)

    photon%now%vec    = 0
    photon%now%vec(5) = s(i)

    photon_xy = -r_max
    call sr3d_find_wall_point (wall, photon, x_wall, y_wall)
    xy_in(i) = wall_xy

    photon_xy = r_max
    call sr3d_find_wall_point (wall, photon, x_wall, y_wall)
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

enddo

end subroutine sr3d_plot_wall_vs_s 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_plot_wall_cross_sections (wall, plot_norm)
!
! Routine to interactively plot wall cross-sections
! Note: This routine never returns to the main program.
!
! Input:
!   wall      -- sr3d_wall_struct: Wall structure.
!   plot_norm -- logical: If True then also plot a set of lines normal to the wall.
!                  This is used as a check to the wall normal calculation.
!-

subroutine sr3d_plot_wall_cross_sections (wall, plot_norm)

implicit none

type (sr3d_wall_struct), target :: wall
type (sr3d_wall_pt_struct), pointer :: pt
type (sr3d_photon_track_struct) photon

real(rp) s_pos, x(400), y(400), x_max, y_max, theta, r, x_max_user, r_max
real(rp) x1_norm(100), y1_norm(100), x2_norm(100), y2_norm(100)

integer i, j, ix, ix_section, i_in, ios, i_chan

character(80) :: ans = 'first', label

logical plot_norm, at_section

! Open plotting window

call qp_open_page ('X', i_chan, 800.0_rp, 400.0_rp, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

x_max_user = -1
r_max = 100

! Print wall info

do i = 0, wall%n_pt_max
  pt => wall%pt(i)
  print '(i4, 2x, a, f12.2)', i, pt%basic_shape(1:12), pt%s
enddo

ix_section = wall%n_pt_max

! Loop

do

  ! Query

  if (ans /= 'first') then
    call read_a_line ('Input: "<Section #>", "<CR>" (Next sec), "b" (Back sec), "s <s_value>", "x <x_max>" ("x auto" -> autoscale)', ans)
  endif

  call string_trim (ans, ans, ix)
  if (ans == 'first') ans = ''

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
      read (ans(2:), *, iostat = ios) r
      if (ios /= 0) then
        print *, 'Cannot read x-scale'
        cycle
      endif
      x_max_user = r
    endif

  elseif (ans == '') then
    ix_section = modulo(ix_section + 1, wall%n_pt_max + 1)
    s_pos = wall%pt(ix_section)%s
    at_section = .true.

  elseif (ans == 'b') then
    ix_section = modulo(ix_section - 1, wall%n_pt_max + 1)
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
    s_pos = wall%pt(ix_section)%s
    at_section = .true.

  endif

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

    if (plot_norm .and. modulo(i, 4) == 0) then
      j = (i / 4)
      call sr3d_find_wall_point (wall, photon, x(i), y(i), x1_norm(j), x2_norm(j), y1_norm(j), y2_norm(j))
    else
      call sr3d_find_wall_point (wall, photon, x(i), y(i))
    endif

  enddo

  x = x * 100; y = y * 100
  x1_norm = x1_norm * 100; x2_norm = x2_norm * 100
  y1_norm = y1_norm * 100; y2_norm = y2_norm * 100

  ! Now plot

  if (at_section) then
    write (label, '(a, f0.3, a, i0, 2a)') 'S: ', s_pos, '   Section #: ', ix_section, '  Name: ', wall%pt(ix_section)%name
  else
    write (label, '(a, f0.3)') 'S: ', s_pos
  endif
  call qp_clear_page
  x_max = 1.01 * maxval(abs(x)); y_max = 1.01 * maxval(abs(y))
  if (x_max_user > 0) x_max = x_max_user
  call qp_calc_and_set_axis ('X', -x_max, x_max, 6, 10, 'ZERO_SYMMETRIC')
  call qp_calc_and_set_axis ('Y', -y_max, y_max, 6, 10, 'ZERO_SYMMETRIC')
  call qp_set_margin (0.07_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
  call qp_eliminate_xy_distortion()
  call qp_draw_graph (x, y, 'X', 'Y', label, .true., 0)

  if (plot_norm) then
    do j = 1, size(x1_norm)
      call qp_draw_line(x1_norm(j), x2_norm(j), y1_norm(j), y2_norm(j))
    enddo
  endif

enddo

end subroutine sr3d_plot_wall_cross_sections

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine sr3d_find_wall_point (wall, photon, x_wall, y_wall, x1_norm, x2_norm, y1_norm, y2_norm)

implicit none

type (sr3d_wall_struct), target :: wall
type (sr3d_photon_track_struct) photon

real(rp) x_wall, y_wall
real(rp), optional :: x1_norm, x2_norm, y1_norm, y2_norm

real(rp) tri_vert0(3), tri_vert1(3), tri_vert2(3)
real(rp) dtrack, d_radius, r_old, dw_perp(3)

integer j, ixp

logical plot_norm, is_through

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
    call sr3d_photon_d_radius (photon%now, wall, d_radius, dw_perp)
    x1_norm = x_wall;                  y1_norm = y_wall
    x2_norm = x_wall + dw_perp(1)/100; y2_norm = y_wall + dw_perp(2)/100
  endif

endif

end subroutine sr3d_find_wall_point

end module
