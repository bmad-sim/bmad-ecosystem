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
    dr_max = 2 * max(cross%v(1)%radius_x, cross%v(2)%radius_y)
  else
    dr_max = 2 * max(maxval(abs(cross%v%x)), maxval(abs(cross%v%y)))
  endif

  ds_max = ele%s
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
  vec(4) = vec(3) + dlen * vec(4)
  vec(6) = vec(5) + dlen * vec(6)
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
vec(4) = vec(3) + dl * vec(4)
vec(6) = vec(5) + dl * vec(6)
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
type (photon_track_struct), target :: photon, photon1

real(rp) d_rad, track_len0, track_len

integer i

character(40) :: r_name = 'capillary_photon_hit_spot_calc'

! Bracket the hit point. 
! Note: After the first reflection, the photon will start at the wall so
! if photon%old is at the wall we must avoid bracketing this point.

photon1 = photon

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


!-----------------------------------------------------------------------
contains

function photon_hit_func (track_len) result (d_radius)

real(rp) track_len, d_radius, radius, d_track

! Easy case

if (track_len == photon%now%track_len) then
  d_radius = capillary_photon_radius (photon%now, ele) - 1
  return
endif

! Track starting from the present position (photon1%now) if track_length > photon1%now%track_len.
! Otherwise, track starting from the beginning of the region (photon%old).

if (track_len < photon1%now%track_len) then
  photon1 = photon
  photon1%now = photon%old
endif

! And track

d_track = track_len - photon1%now%track_len
call capillary_propagate_photon_a_step (photon1, ele, d_track, .false.)
d_radius = capillary_photon_radius (photon1%now, ele) - 1

end function

end subroutine capillary_photon_hit_spot_calc

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

real(rp) radius, g0, g1, g00, g11, s00, s11, dw_x0, dw_x1, dw_y0, dw_y1, f, spline
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

call capillary_cross_section_params (cc0, vec, g0, dw_x0, dw_y0)
call capillary_cross_section_params (cc1, vec, g1, dw_x1, dw_y1)

! Phanom slices

n_slice = nint(cc0%n_slice_spline)
if (n_slice > 1) then
  ix = min(int(f * n_slice), n_slice - 1)

  s00 = float(ix) / n_slice
  spline = cc0%s_spline(1) * s00 + cc0%s_spline(2) * s00**2 + cc0%s_spline(3) * s00**3
  g00 = (1 - spline) * g0 + spline * g1

  s11 = float(ix+1) / n_slice
  spline = cc0%s_spline(1) * s11 + cc0%s_spline(2) * s11**2 + cc0%s_spline(3) * s11**3
  g11 = (1 - spline) * g0 + spline * g1

else
  s00 = cc0%s
  s11 = cc1%s
  g00 = g0
  g11 = g1
endif

!

f = (vec(5) - s00) / (s11 - s00)

radius = 1 / ((1 - f) * g00 + f * g11)

if (present (perp)) then
  perp(1) = (1 - f) * dw_x0 + f * dw_x1
  perp(2) = (1 - f) * dw_y0 + f * dw_y1
  perp(3) = (g00 - g11) / (s11 - s00)
  perp = perp / sqrt(sum(perp**2))  ! Normalize
endif

end function capillary_photon_radius

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine cross_section_params (cross, vec, g, dw_x, dw_y)
!
! Routine to compute parameters needed by capillary_photon_radius routine.
!
! Input:
!   cross   -- cross_section_struct: cross-section at a particular longitudinal location.
!   vec(6)  -- Real(rp): Photon phase space coords. 
!
! Output:
!   g              -- Real(rp): Radius of the wall / radius of the photon.
!   [dw_x, dw_y]   -- Real(rp): Transverse directional derivatives of -g.
!-

subroutine capillary_cross_section_params (cross, vec, g, dw_x, dw_y)

implicit none

type (cross_section_struct) cross, pt
type (cross_section_vertex_struct), pointer :: v(:)

real(rp) g, dw_x, dw_y, vec(6)
real(rp) r_wall, r_photon, gradient(2)

integer ix

! 

call calc_wall_radius (cross%v, atan2(vec(3), vec(1)), r_wall, gradient)
r_photon = sqrt(vec(1)**2 + vec(3)**2)
g = r_wall / r_photon
dw_x = gradient(1) * g / r_photon
dw_y = gradient(2) * g / r_photon

end subroutine capillary_cross_section_params

end module
