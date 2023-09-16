!+
! Function ray_is_outside_wall (ray, walls) result (is_outside)
!
! Routine to determine if a ray position is outside the vacuum chamber wall.
!
! Input:
!   ray   -- ray_struct: photon position.
!   walls -- walls_struct: Vacuum chamber wall.
!
! Output:
!   is_outside  -- logical: True if outside wall
!-

function ray_is_outside_wall (ray, walls) result (is_outside)

use synrad_interface, except => ray_is_outside_wall

implicit none

type (ray_struct) ray
type (walls_struct), target :: walls

real(rp) x_ray, s_ray, wind_p, wind_m
logical is_outside

! Exceptional cases at ends of lattice

if (ray%now%s == 0) then
  is_outside = (outside_at_end (walls%positive_x_wall, 0) .or. outside_at_end (walls%negative_x_wall, 0))
  return
endif

if (ray%now%s == walls%s_max) then
  is_outside = (outside_at_end (walls%positive_x_wall, walls%positive_x_wall%n_pt_max) .or. &
                outside_at_end (walls%negative_x_wall, walls%negative_x_wall%n_pt_max))
  return
endif

! The particle is inside the wall if the winding angle with respect to the particle of a local 
! section of the wall is negative for the positive_x side and positive on the other side.

x_ray = ray%now%vec(1)
s_ray = ray%now%s

wind_p = calc_winding_angle (walls%positive_x_wall)
wind_m = calc_winding_angle (walls%negative_x_wall)

if (wind_p < 0 .and. wind_m > 0) then
  is_outside = .false.
else
  is_outside = .true.
endif

!----------------------------------------------------------------------------
contains

function calc_winding_angle (side) result (theta)

type (wall_struct) side
type (ray_struct) ray2
real(rp) theta, dtheta
integer i, j, ix_pt

!

ray2 = ray
ray2%direction = 1
call get_initial_wall_pt (ray2, side, ix_pt)

theta = 0

do 
  dtheta = atan2(side%pt(ix_pt)%x - x_ray, side%pt(ix_pt)%s - s_ray) - & 
           atan2(side%pt(ix_pt-1)%x - x_ray, side%pt(ix_pt-1)%s - s_ray)
  dtheta = modulo2(dtheta, pi)

  theta = theta + dtheta

  if (side%pt(ix_pt)%s >= ray%now%s .and. .not. side%pt(ix_pt)%next_to_alley) exit
  ix_pt = ix_pt + 1
enddo

end function

!----------------------------------------------------------------------------
! contains

function outside_at_end (side, ix_pt) result (is_outside)

type (wall_struct) side
integer ix_pt
logical is_outside

!

if ((ray%now%vec(1) - side%pt(ix_pt)%x) * side%side > 0) then
  is_outside = .true.
else
  is_outside = .false.
endif

end function

end function
