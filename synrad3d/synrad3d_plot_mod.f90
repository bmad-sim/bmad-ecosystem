module sr3d_plot_mod

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
!   wall -- wall3d_struct: Wall structure.
!-

subroutine sr3d_plot_wall_cross_sections (wall)

implicit none

type (wall3d_struct), target :: wall
type (wall3d_pt_struct), pointer :: pt
type (photon3d_track_struct) photon

real(rp) s_pos, dtrack, x(400), y(400), x_max, y_max, theta, d_radius
real(rp) tri_vert0(3), tri_vert1(3), tri_vert2(3)

integer i, j, ix, i_last, i_in, ios, i_chan, ixp

character(40) ans, label

logical is_through

! Open plotting window

call qp_open_page ('X', i_chan, 500.0_rp, 400.0_rp, 'POINTS')
call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')

! Print wall info

do i = 0, wall%n_pt_max
  pt => wall%pt(i)
  print '(i4, 2x, a, f12.2)', i, pt%basic_shape(1:12), pt%s
enddo

i_last = wall%n_pt_max

! Loop

do

  ! Get s-position

  call read_a_line ('Input: "<Number>" (# of section), "s <s_value>", <CR> = Next section', ans)
  call string_trim (ans, ans, ix)
  if (ans(1:1) == 's') then
    read (ans(2:), *, iostat = ios) s_pos
    if (ios /= 0 .or. s_pos < wall%pt(0)%s .or. s_pos > wall%pt(wall%n_pt_max)%s) then
      print *, 'Cannot read s-position or s-position out of range.'
      cycle
    endif

  else
    if (ans == '') then
      i_last = modulo(i_last + 1, wall%n_pt_max + 1)
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
      i_last = i_in
    endif
    s_pos = wall%pt(i_last)%s

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

  write (label, '(a, f0.2)') 'S: ', s_pos
  call qp_clear_page
  x_max = maxval(abs(x)); y_max = maxval(abs(y))
  call qp_calc_and_set_axis ('X', -x_max, x_max, 6, 10, 'ZERO_SYMMETRIC')
  call qp_calc_and_set_axis ('Y', -y_max, y_max, 6, 10, 'ZERO_SYMMETRIC')
  call qp_set_margin (0.07_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
  call qp_eliminate_xy_distortion()
  call qp_draw_graph (x, y, 'X', 'Y', label, .true., 0)


enddo

end subroutine sr3d_plot_wall_cross_sections

end module
