!+
! Function distance_to_aperture (orbit, particle_at, ele, no_aperture_here) result (dist)
!
! Routine to calculate the distance of the particle from the wall aperture normalized by the distance
! from the origin to the wall.
! Distances are negative if the particle is inside the wall and positive if outside.
!
! This routine is similar to wall3d_d_radius except that this routine works with apertures set by the
! element parameters x1_limit, y1_limit, x2_limit, y2_limit.
!
! If in a patch element this routine assues the the orbit is with respect to the downstream coordinates
! if the particle is moving forward and with respect to the upstream coords if moving backwards.
!
! Input:
!   orbit           -- coord_struct: Particle position.
!   particle_at     -- integer: first_track_edge$, second_track_edge$, or in_between$
!   ele             -- ele_struct: Element containing aperture.
!
! Output:
!   no_aperture_here  -- logical: True if aperture does not exist at the longitudinal location of the particle.
!   dist              -- real(rp): Normalized distance of the particle from the aperture.
!-

recursive function distance_to_aperture (orbit, particle_at, ele, no_aperture_here) result (dist)

use wall3d_mod, dummy => distance_to_aperture

implicit none

type (coord_struct) orbit
type (ele_struct) ele
type (ele_struct), pointer :: lord

real(rp) dist, x_particle, y_particle, r_wall, x_lim, y_lim, x0, y0, r, position(6), d_radius, lord_dist
integer particle_at, physical_end, i
logical no_aperture_here, no_ap

! Init

dist = real_garbage$
no_aperture_here = .true.
if (.not. bmad_com%aperture_limit_on) return

! super_slave elements have the aperture info stored in the lord(s).
! In this case return the most positive number

if (ele%slave_status == super_slave$) then
  physical_end = physical_ele_end (particle_at, orbit, ele%orientation)
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%lord_status /= super_lord$) cycle
    if (.not. lord_edge_aligned (ele, physical_end, lord) .and. lord%aperture_at /= continuous$) cycle
    lord_dist = distance_to_aperture (orbit, particle_at, lord, no_ap)
    if (dist == real_garbage$ .or. lord_dist > dist) then
      dist = lord_dist
      no_aperture_here = no_ap
    endif
  enddo
  return
endif

! Custom calc

if (ele%aperture_type == custom_aperture$) then
  dist = distance_to_aperture_custom (orbit, particle_at, ele, no_aperture_here)
  return
endif

! Aperture here?

physical_end = physical_ele_end (particle_at, orbit, ele%orientation)
if (.not. at_this_ele_end (physical_end, ele%aperture_at)) return

! Calc dist

no_aperture_here = .false.

select case (ele%aperture_type)

case (elliptical$, rectangular$, auto_aperture$)
  x_lim = (ele%value(x1_limit$) + ele%value(x2_limit$)) / 2
  y_lim = (ele%value(y1_limit$) + ele%value(y2_limit$)) / 2
  x0 = (ele%value(x2_limit$) - ele%value(x1_limit$)) / 2
  y0 = (ele%value(y2_limit$) - ele%value(y1_limit$)) / 2
  x_particle = orbit%vec(1) - x0
  y_particle = orbit%vec(3) - y0

  if (x_lim == 0) x_lim = bmad_com%max_aperture_limit
  if (y_lim == 0) y_lim = bmad_com%max_aperture_limit

  if (ele%aperture_type == elliptical$) then
    r = (x_particle / x_lim)**2 + (y_particle / y_lim)**2
    dist = sqrt(r) - 1
  else      ! rectangular$, auto_aperture$
    if (abs(x_particle/x_lim) > abs(y_particle/y_lim)) then
      dist = abs(x_particle/x_lim) - 1
    else
      dist = abs(y_particle/y_lim) - 1
    endif
  endif

case (wall3d$)
  position = wall3d_to_position(orbit, ele)
  d_radius = wall3d_d_radius (position, ele, no_wall_here = no_aperture_here, radius_wall = r_wall)
  dist = d_radius / r_wall
end select

end function distance_to_aperture
