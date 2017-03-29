!+
! Function distance_to_aperture (orbit, particle_at, ele, no_aperture_here) result (dist)
!
! Routine to calculate the percentage distance from the particle to the wall aperture.
! Distances are negative if the particle is inside the wall.
!
! This routine is similar to wall3d_d_radius except that this routine works with apertures set by the
! element parameters x1_limit, y1_limit, x2_limit, y2_limit.
!
! Input:
!   orbit           -- coord_struct: Particle position.
!   particle_at     -- integer: 
!   ele             -- ele_struct: Element containing aperture.
!
! Output:
!   no_aperture_here  -- logical, optional: True if aperture does not exist at the
!                          longitudinal location of the particle.
!   dist              -- real(rp): Percentage distance to the aperture.
!-

function distance_to_aperture (orbit, particle_at, ele, no_aperture_here) result (dist)

use wall3d_mod, dummy => distance_to_aperture

implicit none

type (coord_struct) orbit
type (ele_struct) ele

real(rp) dist, x_particle, y_particle, r_wall, x_lim, y_lim, x0, y0, r, position(6), d_radius
integer particle_at, physical_end
logical no_aperture_here

!

dist = real_garbage$
no_aperture_here = .true.

if (ele%aperture_type == custom_aperture$) then
  dist = distance_to_aperture_custom (orbit, particle_at, ele, no_aperture_here)
  return
endif

!

physical_end = physical_ele_end (particle_at, orbit%direction, ele%orientation)
if (.not. at_this_ele_end (physical_end, ele%aperture_at)) return

x_lim = (ele%value(x1_limit$) + ele%value(x2_limit$)) / 2
y_lim = (ele%value(y1_limit$) + ele%value(y2_limit$)) / 2
x0 = (ele%value(x2_limit$) - ele%value(x1_limit$)) / 2
y0 = (ele%value(y2_limit$) - ele%value(y1_limit$)) / 2
x_particle = orbit%vec(1) - x0
y_particle = orbit%vec(3) - y0

if (.not. bmad_com%aperture_limit_on .or. x_lim == 0) x_lim = bmad_com%max_aperture_limit
if (.not. bmad_com%aperture_limit_on .or. y_lim == 0) y_lim = bmad_com%max_aperture_limit

select case (ele%aperture_type)

case (elliptical$)
  r = (x_particle / x_lim)**2 + (y_particle / y_lim)**2
  dist = sqrt(r) - 1

case (rectangular$, auto_aperture$)
  if (abs(x_particle/x_lim) > abs(y_particle/y_lim)) then
    dist = abs(x_particle/x_lim) - 1
  else
    dist = abs(y_particle/y_lim) - 1
  endif

case (wall3d$)
  position = [orbit%vec(1:4), orbit%s-ele%s_start, 1.0_rp]
  d_radius = wall3d_d_radius (position, ele, no_wall_here = no_aperture_here, radius_wall = r_wall)

  dist = d_radius / r_wall

end select

end function distance_to_aperture
