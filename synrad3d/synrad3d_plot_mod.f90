module synrad3d_plot_mod

use synrad3d_track_mod
use quick_plot
use input_mod

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_plot_wall_cross_sections (wall)
!
! Routine to interactively plot wall cross-sections
! Note: This routine never returns to the main program.
!
! Input:
!   wall -- sr3d_wall_struct: Wall structure.
!-

subroutine sr3d_plot_wall_cross_sections (wall)

implicit none

type (sr3d_wall_struct), target :: wall
type (sr3d_wall_pt_struct), pointer :: pt
type (sr3d_photon_track_struct) photon

real(rp) s_pos, dtrack, x(400), y(400), x_max, y_max, theta, d_radius, r
real(rp) tri_vert0(3), tri_vert1(3), tri_vert2(3), x_max_user

integer i, j, ix, ix_section, i_in, ios, i_chan, ixp

character(80) ans, label

logical is_through, first, at_section

! Open plotting window

call qp_open_page ('X', i_chan, 800.0_rp, 400.0_rp, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

x_max_user = -1
first = .true.

! Print wall info

do i = 0, wall%n_pt_max
  pt => wall%pt(i)
  print '(i4, 2x, a, f12.2)', i, pt%basic_shape(1:12), pt%s
enddo

ix_section = wall%n_pt_max

! Loop

do

  ! Query

  if (first) then
    first = .false.
    ans = ''
  else
    call read_a_line ('Input: "<Section #>", "<CR>" (Next sec), "b" (Back sec), "s <s_value>", "x <x_max>" (neg -> autoscale)', ans)
  endif

  call string_trim (ans, ans, ix)
  if (ans(1:1) == 's') then
    read (ans(2:), *, iostat = ios) s_pos
    if (ios /= 0 .or. s_pos < wall%pt(0)%s .or. s_pos > wall%pt(wall%n_pt_max)%s) then
      print *, 'Cannot read s-position or s-position out of range.'
      cycle
    endif
    at_section = .false.

  elseif (ans(1:1) == 'x') then
    read (ans(2:), *, iostat = ios) r
    if (ios /= 0) then
      print *, 'Cannot read x-scale'
      cycle
    endif
    x_max_user = r

  elseif (ans == '') then
    ix_section = modulo(ix_section + 1, wall%n_pt_max + 1)
    s_pos = wall%pt(ix_section)%s
    at_section = .true.

  elseif (ans == '') then
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

  photon%old%vec = 0
  photon%old%vec(5) = s_pos
  photon%now%vec(5) = s_pos

  call sr3d_get_wall_index (photon%now, wall, ixp)

  do i = 1, size(x)

    ! Idea is to see where photon path from photon%old (at the origin) to %now intersects the wall.
    ! photon%now is at 1 meter radius which is assumed to be outside the wall.

    theta = (i-1) * twopi / (size(x) - 1)
    photon%now%vec(1) = cos(theta)  
    photon%now%vec(3) = sin(theta)

    if (wall%pt(ixp+1)%basic_shape == 'gen_shape_mesh') then
      do j = 1, 2*size(wall%pt(ixp)%gen_shape%v)
        call sr3d_get_mesh_wall_triangle_pts (wall%pt(ixp), wall%pt(ixp+1), j, tri_vert0, tri_vert1, tri_vert2)
        call sr3d_mesh_triangle_intersect (photon, tri_vert0, tri_vert1, tri_vert2, is_through, dtrack)
        if (is_through) then
          x(i) = dtrack * photon%now%vec(1)
          y(i) = dtrack * photon%now%vec(3)
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
      x(i) = (1 - d_radius) * photon%now%vec(1)
      y(i) = (1 - d_radius) * photon%now%vec(3)

    endif
  enddo

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


enddo

end subroutine sr3d_plot_wall_cross_sections

end module
