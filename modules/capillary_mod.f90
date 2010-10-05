module capillary_mod

use cross_section_mod

type photon_coord_struct
  type (coord_struct) orb       ! Phase space: orb%vec = (x, vx/c, y, vy/c, s, vs/c)
  real(rp) energy
  real(rp) track_len
  integer ix_cross              ! Cross section index
end type

type photon_track_struct
  type (photon_coord_struct) old, now
end type

! This is for passing info to the photon_hit_func routine used by zbrent

type (photon_track_struct), pointer, private :: photon_com
type (photon_track_struct), private :: photon1_com
type (ele_struct), private, pointer :: ele_com

private photon_hit_func

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_a_capillary (orb, ele, lost)
!
! Routine to track through a capillary.
!
! Input:
!   orb -- Coord_struct: Input coordinates.
!   ele -- ele_struct: Capillary element
!
! Output:
!   orb  -- Coord_struct: Output coordinates.
!   lost -- Logical: Set true if particle is lost. False otherwise.
!-

subroutine track_a_capillary (orb, ele, lost)

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (photon_track_struct) photon

real(rp) pz

integer at_end

logical lost

! Calculate maximum step size

! Setup photon coords

photon%now%orb = orb
photon%now%orb%vec(5) = 0
photon%now%orb%vec(6) = sqrt(1 - orb%vec(2)**2 - orb%vec(4)**2)
photon%now%energy = ele%value(e_tot$) * (1 + orb%vec(6))
photon%now%track_len = 0
photon%now%ix_cross = 1

! Loop over all bounces

lost = .true.
at_end = 0

do
  call capillary_track_photon_to_wall (photon, ele, at_end)
  if (at_end == l_start$) return  ! Reflected backwards
  if (at_end == l_end$) exit      ! And done
  call capillary_reflect_photon (photon, ele, lost)
  if (lost) return
enddo

! Correct z by distance traveled by reference particle.

pz = orb%vec(6)
orb = photon%now%orb
orb%vec(5) =  ele%value(l$) - photon%now%track_len
orb%vec(6) = pz

lost = .false.

end subroutine track_a_capillary

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine capillary_track_photon_to_wall (photon, ele, at_end)
!
! Routine to track through a capillary.
!
! Input:
!   photon -- photon_track_struct: Input coordinates.
!   ele    -- ele_struct: Capillary element
!
! Output:
!   photon -- Photon_struct: Output coordinates.
!   at_end -- Integer: Set to l_start$, or l_end$ if at an end.
!-

subroutine capillary_track_photon_to_wall (photon, ele, at_end)

implicit none

type (ele_struct), target :: ele
type (photon_track_struct), target :: photon
type (cross_section_struct), pointer :: cross

real(rp) dr_max, ds_max, p_max, dlen
real(rp), pointer :: vec(:)

integer at_end

character(40) :: r_name = 'capillary_track_photon_to_wall'

! 

do

  cross => ele%cross_section(photon%now%ix_cross)
  vec => photon%now%orb%vec

  ! Calculate a resonable step size

  if (size(cross%v) == 1) then
    dr_max = 2 * max(cross%v(1)%radius_x, cross%v(1)%radius_y)
  else
    dr_max = 2 * max(maxval(abs(cross%v%x)), maxval(abs(cross%v%y)))
  endif

  ds_max = 2 * ele%s
  p_max = max(abs(vec(2)), abs(vec(4)))
 
  if (dr_max * abs(vec(6)) > ds_max * p_max) then
    dlen = ds_max / abs(vec(6))
  else
    dlen = dr_max / p_max
  endif

  ! Propagate and see if there is a hit

  call capillary_propagate_photon_a_step (photon, ele, dlen, .true.)

  if (capillary_photon_radius (photon%now, ele) > 1) then
    call capillary_photon_hit_spot_calc (photon, ele)
    return
  endif

  if (vec(5) == 0) then
    at_end = l_start$
    return
  endif

  if (vec(5) == ele%s) then
    at_end = l_end$
    return
  endif

enddo

end subroutine capillary_track_photon_to_wall 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine capillary_propagate_photon_a_step (photon, ele, dlen, stop_at_boundary)
!
! Routine to track a photon a step of a given length
!
! Input:
!   photon -- photon_track_struct: Input coordinates.
!   ele    -- ele_struct: Capillary element
!   dlen   -- Real(rp): Length to propagate a photon. The actual propagation length
!               may be less if stop_at_boundary is True.
!
! Output:
!   photon -- Photon_struct: Output coordinates.
!   stop_at_boundary -- Logical: If True then stop at cross-section boundary.
!-

subroutine capillary_propagate_photon_a_step (photon, ele, dlen, stop_at_boundary)

implicit none

type (ele_struct), target :: ele
type (photon_track_struct), target :: photon

real(rp) dlen, dl, s_stop 
real(rp), pointer :: vec(:)

integer ixc

logical stop_at_boundary

!

photon%old = photon%now
vec => photon%now%orb%vec

! no stop is easy

if (.not. stop_at_boundary) then
  vec(1) = vec(1) + dlen * vec(2)
  vec(3) = vec(3) + dlen * vec(4)
  vec(5) = vec(5) + dlen * vec(6)
  photon%now%track_len = photon%now%track_len + dlen
  return
endif

! Stop at boundary

call bracket_index(ele%cross_section%s, 1, size(ele%cross_section), vec(5), ixc)

if (vec(6) > 0) then
  s_stop = ele%cross_section(ixc+1)%s
else
  s_stop = ele%cross_section(ixc)%s
endif

if (abs(vec(6)) * dlen > abs(s_stop - vec(5))) then
  dl = (s_stop - vec(5)) / vec(6)
else
  dl = dlen
endif

vec(1) = vec(1) + dl * vec(2)
vec(3) = vec(3) + dl * vec(4)
vec(5) = vec(5) + dl * vec(6)
photon%now%track_len = photon%now%track_len + dl

end subroutine capillary_propagate_photon_a_step

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine capillary_photon_hit_spot_calc (photon, ele) result (radius)
!
! Routine to interpolate to where the photon has hit the capillary.
!
! Input:
!   photon -- photon_track_struct: Input coordinates.
!   ele    -- ele_struct: Capillary element
!
! Output:
!   photon -- photon_track_struct: Photon at capillary wall
!-

subroutine capillary_photon_hit_spot_calc (photon, ele)

use nr, only: zbrent

implicit none

type (ele_struct), target :: ele
type (photon_track_struct), target :: photon

real(rp) d_rad, track_len0, track_len

integer i

character(40) :: r_name = 'capillary_photon_hit_spot_calc'

! Bracket the hit point. 
! Note: After the first reflection, the photon will start at the wall so
! if photon%old is at the wall we must avoid bracketing this point.

photon1_com = photon
photon_com => photon
ele_com => ele

track_len0 = (photon%now%track_len + photon%old%track_len) / 2
do i = 1, 30
  d_rad = photon_hit_func(track_len0)
  if (d_rad < 0) exit
  track_len0 = (track_len0 + photon%old%track_len) / 2
  if (i == 30) then
    call out_io (s_abort$, r_name, 'ERROR: CANNOT FIND HIT SPOT REGION LOWER BOUND!')
    call err_exit
  endif
enddo

! Find where the photon hits.

track_len = zbrent (photon_hit_func, track_len0, photon%now%track_len, 1d-10)

! Cleanup

photon%now = photon%old
call capillary_propagate_photon_a_step (photon, ele, track_len-photon%now%track_len, .false.)

end subroutine capillary_photon_hit_spot_calc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function photon_hit_func (track_len) result (d_radius)
! 
! Routine to be used as an argument in zbrent in the capillary_photon_hit_spot_calc.
!
! Input:
!   track_len -- Real(rp): Place to position the photon.
!
! Output:
!   d_radius -- Real(rp): 

function photon_hit_func (track_len) result (d_radius)

real(rp), intent(in) :: track_len
real(rp) :: d_radius
real(rp) radius, d_track

! Easy case

if (track_len == photon_com%now%track_len) then
  d_radius = capillary_photon_radius (photon_com%now, ele_com) - 1
  return
endif

! Track starting from the present position (photon1_com%now) if track_length > photon1_com%now%track_len.
! Otherwise, track starting from the beginning of the region (photon%old).

if (track_len < photon1_com%now%track_len) then
  photon1_com = photon_com
  photon1_com%now = photon_com%old
endif

! And track

d_track = track_len - photon1_com%now%track_len
call capillary_propagate_photon_a_step (photon1_com, ele_com, d_track, .false.)
d_radius = capillary_photon_radius (photon1_com%now, ele_com) - 1

end function photon_hit_func

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine capillary_reflect_photon (photon, ele, absorbed)
!
! Routine to reflect a photon from the capillary wall.
!
! Input:
!   photon -- photon_track_struct: Input coordinates.
!   ele    -- ele_struct: Capillary element
!
! Output:
!   photon   -- Photon_struct: Output coordinates.
!   absorbed -- Logical: If true the photon has been absorbed.
!-

subroutine capillary_reflect_photon (photon, ele, absorbed)

implicit none

type (ele_struct), target :: ele
type (photon_track_struct), target :: photon

real(rp) perp(3), r, graze_angle, cos_perp
real(rp), pointer :: vec(:)

logical absorbed

!

photon%old = photon%now
r = capillary_photon_radius (photon%now, ele, perp)

vec => photon%now%orb%vec
cos_perp = dot_product (vec(2:6:2), perp)
vec(2:6:2) = vec(2:6:2) - 2 * cos_perp * perp

absorbed = .false.
graze_angle = pi/2 - acos(cos_perp)
if (graze_angle > ele%value(critical_angle_factor$) / photon%now%energy) absorbed = .true.

end subroutine capillary_reflect_photon

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function capillary_photon_radius (p_orb, ele, perp) result (radius)
!
! Routine to calculate the normalized radius = photon_radius / wall_radius.
!
! Input:
!   p_orb   -- photon_coord_struct: Input coordinates.
!   ele     -- ele_struct: Capillary element
!
! Output:
!   radius  -- Real(rp), Normalized radius: r_photon / r_wall
!   perp(3) -- Real(rp), optional: Perpendicular normal to the wall.
!-

function capillary_photon_radius (p_orb, ele, perp) result (radius)

implicit none

type (ele_struct), target :: ele
type (photon_coord_struct), target :: p_orb
type (cross_section_struct), pointer :: cc0, cc1

real(rp) radius, s00, s11, f, spline, r00_wall, r11_wall
real(rp) r0_wall, r1_wall, gradient0(2), gradient1(2), f_eff, ds_spline, theta_photon
real(rp), optional :: perp(3)
real(rp), pointer :: vec(:)

integer ix, n_slice

! There is a sigularity in the calculation when the photon is at the origin.
! To avoid this, just return radius = 0 for small radii.

vec => p_orb%orb%vec

if (abs(vec(1)) < 1e-10 .and. abs(vec(3)) < 1e-10) then
  radius = 0
  if (present (perp)) perp = 0
  return
endif

!

call bracket_index (ele%cross_section%s, 1, size(ele%cross_section), vec(5), ix)
p_orb%ix_cross = ix

if (ix == size(ele%cross_section)) ix = size(ele%cross_section) - 1

! The outward normal vector is discontinuous at the wall points.
! If at a wall point, use the correct part of the wall.

if (vec(5) == ele%cross_section(ix)%s .and. vec(6) > 0) then
  if (ix /= 0) then
    ix = ix - 1
  endif
endif

!

cc0 => ele%cross_section(ix)
cc1 => ele%cross_section(ix+1)

theta_photon = atan2(vec(3), vec(1))
call calc_wall_radius (cc0%v, theta_photon, r0_wall, gradient0)
call calc_wall_radius (cc1%v, theta_photon, r1_wall, gradient1)

! Phanom slices

n_slice = nint(cc0%n_slice_spline)
if (n_slice > 1) then
  ix = min(int(f * n_slice), n_slice - 1)

  s00 = float(ix) / n_slice
  spline = cc0%s_spline(1) * s00 + cc0%s_spline(2) * s00**2 + cc0%s_spline(3) * s00**3
  r00_wall = (1 - spline) * r0_wall + spline * r1_wall

  s11 = float(ix+1) / n_slice
  spline = cc0%s_spline(1) * s11 + cc0%s_spline(2) * s11**2 + cc0%s_spline(3) * s11**3
  r11_wall = (1 - spline) * r0_wall + spline * r1_wall

else
  s00 = cc0%s
  s11 = cc1%s
  r00_wall = r0_wall
  r11_wall = r1_wall
endif

!

f = (vec(5) - s00) / (s11 - s00)

if (n_slice == 0) then
  f_eff = cc0%s_spline(1) * f + cc0%s_spline(2) * f**2 + cc0%s_spline(3) * f**3
else
  f_eff = f
  ds_spline = s11 - s00
endif

radius = sqrt(vec(1)**2 + vec(3)**2) / ((1 - f_eff) * r00_wall + f_eff * r11_wall)

if (present (perp)) then
  perp(1:2) = (1 - f_eff) * gradient0 + f_eff * gradient1
  perp(3)   = -(r11_wall - r00_wall) / ds_spline
  perp = perp / sqrt(sum(perp**2))  ! Normalize
endif

end function capillary_photon_radius

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!+
! Subroutine calc_wall_radius (v, theta, r_wall, gradient)
!
! Routine to calculate the wall radius at a given angle for a given cross-section
! Additionally, the transverse directional derivative is calculated.
!
! Module needed:
!   use cross_setction_mod
!
! Input:
!   v(:)        -- cross_section_vertex_struct: Array of vertices.
!   theta       -- Real(rp): Angle to evaluate at.
!
! Output:
!   r_wall      -- Real(rp): Wall radius at given angle
!   gradient(2) -- Real(rp): (dx, dy) directional derivative of f = r - r_w(theta)
!                              evaluated at r = r_wall.
!-

subroutine calc_wall_radius (v, theta, r_wall, gradient)

implicit none

type (cross_section_vertex_struct), target :: v(:)
type (cross_section_vertex_struct), pointer :: v1, v2


real(rp) theta, r_wall, gradient(2)
real(rp) angle, numer, denom, ct, st, x0, y0, gx, gy, a, b, c
real(rp) cos_ang, sin_ang, radx, f

integer ix

! Bracket index if there is more than one vertex
! If there is only one vertex then must be an ellipse or circle

angle = theta

if (size(v) == 1) then
  v2 => v(1)
else
  if (angle < v(1)%angle) angle = ceiling((v(1)%angle-angle)/twopi) * twopi + angle
  call bracket_index (v%angle, 1, size(v), angle, ix)

  v1 => v(ix)
  v2 => v(ix+1)
endif

! Straight line case

if (v2%radius_x == 0) then
  cos_ang = cos(theta)
  sin_ang = sin(theta)
  numer = (v1%x * v2%y - v1%y * v2%x)
  denom = (cos_ang * (v2%y - v1%y) - sin_ang * (v2%x - v1%x))
  r_wall = numer / denom
  gradient =  [cos_ang, sin_ang] + [-sin_ang, cos_ang] * &
                      (sin_ang * (v2%y - v1%y) + cos_ang * (v2%x - v1%x)) * numer / denom**2
  return
endif

! Convert into unrotated frame if tilted ellipse

x0 = v2%x0; y0 = v2%y0

if (v2%tilt /= 0) then
  angle = angle - v2%tilt
  ct = cos(v2%tilt); st = sin(v2%tilt)
  x0 =  ct * v2%x0 + st * v2%y0
  y0 = -st * v2%x0 + ct * v2%y0
endif

! If ellipse then shrink along y-axis

if (v2%radius_y /= 0) then
  y0 = y0 * v2%radius_x / v2%radius_y
  angle = atan2(sin(angle) * v2%radius_x / v2%radius_y, cos(angle))
endif

cos_ang = cos(angle)
sin_ang = sin(angle)

! Find wall point and derivative

a = 1
b = 2 * (cos_ang * x0 + sin_ang * y0)
c = x0**2 + y0**2 - v2%radius_x**2
radx = sqrt(b**2 - 4 * a * c)
if (v2%radius_x < 0) then
  r_wall = (b - radx) / (2 * a)
  f = (-1 - b/radx) / (2 * a)
else
  r_wall = (b + radx) / (2 * a)
  f = (-1 + b/radx) / (2 * a)
endif
gradient = [cos_ang, sin_ang] - [sin_ang, -cos_ang] * f * 2 * (sin_ang * x0 - cos_ang * y0)

! If an ellipse then expand along y-axis

if (v2%radius_y /= 0) then
  r_wall = r_wall * sqrt(cos_ang**2 + (sin_ang * v2%radius_y / v2%radius_x)**2)
  gradient(2) = gradient(2) * v2%radius_y / v2%radius_x
  gradient = gradient / sqrt(1 + (v2%radius_y / v2%radius_x)**2)
endif

! If a tilted ellipse then rotate

if (v2%tilt /= 0) then
  ct = cos(v2%tilt); st = sin(v2%tilt)
  gradient =  [ct * gradient(1) -  st * gradient(2), &
               st * gradient(1) +  ct * gradient(2)]
endif

end subroutine calc_wall_radius

end module
