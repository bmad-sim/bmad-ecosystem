module capillary_mod

use bmad_struct
use bmad_interface

type photon_coord_struct
  type (coord_struct) orb       ! Phase space: orb%vec = (x, vx/c, y, vy/c, s, vs/c)
  real(rp) energy               ! in eV.
  real(rp) track_len            ! Total track length from the start of the element.
  integer ix_section            ! Cross section index
end type

type photon_track_struct
  type (photon_coord_struct) old, now
end type

! This is for passing info to the photon_hit_func routine used by zbrent

type (photon_track_struct), pointer, private, save :: photon_com
type (photon_track_struct), private, save :: photon1_com
type (ele_struct), private, pointer, save :: ele_com

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

integer status

logical lost

! Calculate maximum step size

! Setup photon coords

photon%now%orb = orb
photon%now%orb%vec(5) = 0
photon%now%orb%vec(6) = sqrt(1 - orb%vec(2)**2 - orb%vec(4)**2)
photon%now%energy = ele%value(e_tot$) * (1 + orb%vec(6))
photon%now%track_len = 0
photon%now%ix_section = 1

! Loop over all bounces

lost = .true.
status = not_lost$

do
  call capillary_track_photon_to_wall (photon, ele, status)
  if (status == entrance_end$ .or. status == lost$) return  ! Reflected backwards or lost
  if (status == exit_end$) exit      ! And done
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
! Subroutine capillary_track_photon_to_wall (photon, ele, status)
!
! Routine to track through a capillary.
!
! Input:
!   photon -- photon_track_struct: Input coordinates.
!   ele    -- ele_struct: Capillary element
!
! Output:
!   photon -- Photon_track_struct: Output coordinates.
!   status -- Integer: Set to entrance_end$, exit_end$ if at an end, 
!               or lost$ if absorbed.
!-

subroutine capillary_track_photon_to_wall (photon, ele, status)

implicit none

type (ele_struct), target :: ele
type (photon_track_struct), target :: photon
type (wall3d_section_struct), pointer :: section

real(rp) dr_max, ds_max, p_max, dlen
real(rp), pointer :: vec(:)

integer status

character(40) :: r_name = 'capillary_track_photon_to_wall'

! Check if outside wall

if (capillary_photon_d_radius (photon%now, ele) >= 0) then
  status = lost$
  return
endif

! propagate

do

  section => ele%wall3d%section(photon%now%ix_section)
  vec => photon%now%orb%vec

  ! Calculate a resonable step size

  if (size(section%v) == 1) then
    dr_max = 2 * max(section%v(1)%radius_x, section%v(1)%radius_y)
  else
    dr_max = 2 * max(maxval(abs(section%v%x)), maxval(abs(section%v%y)))
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

  if (capillary_photon_d_radius (photon%now, ele) > 0) then
    call capillary_photon_hit_spot_calc (photon, ele)
    return
  endif

  if (vec(5) == 0) then
    status = entrance_end$
    return
  endif

  if (vec(5) == ele%s) then
    status = exit_end$
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
!   photon           -- Photon_track_struct: Output coordinates.
!   stop_at_boundary -- Logical: If True then stop at cross-section boundary.
!-

subroutine capillary_propagate_photon_a_step (photon, ele, dlen, stop_at_boundary)

implicit none

type (ele_struct), target :: ele
type (photon_track_struct), target :: photon
type (wall3d_section_struct), pointer :: section(:)

real(rp) dlen, dl, s_stop 
real(rp), pointer :: vec(:)

integer ix, ixc

logical stop_at_boundary

!

photon%old = photon%now
vec => photon%now%orb%vec

! If we do not have to stop at a boundary then propagation is easy.

if (.not. stop_at_boundary) then
  vec(1) = vec(1) + dlen * vec(2)
  vec(3) = vec(3) + dlen * vec(4)
  vec(5) = vec(5) + dlen * vec(6)
  photon%now%track_len = photon%now%track_len + dlen
  return
endif

! Here if stopping at a boundary plane is to be done...
! First see where we need to stop.

section => ele%wall3d%section
call bracket_index(section%s, 1, size(section), vec(5), ixc)

if (vec(6) > 0) then   ! Forward going photon
  s_stop = section(ixc+1)%s
else   ! Backward going photon
  s_stop = section(ixc)%s
endif

! Now calculate the distance to track

if (abs(vec(6)) * dlen > abs(s_stop - vec(5))) then
  dl = (s_stop - vec(5)) / vec(6)
else
  dl = dlen
endif

! And track to the stopping point.

vec(1) = vec(1) + dl * vec(2)
vec(3) = vec(3) + dl * vec(4)
vec(5) = vec(5) + dl * vec(6)
photon%now%track_len = photon%now%track_len + dl

end subroutine capillary_propagate_photon_a_step

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine capillary_photon_hit_spot_calc (photon, ele)
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
    if (bmad_status%exit_on_error) call err_exit
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
!   d_radius -- Real(rp): r_photon - r_wall

function photon_hit_func (track_len) result (d_radius)

real(rp), intent(in) :: track_len
real(rp) :: d_radius
real(rp) radius, d_track

! Easy case

if (track_len == photon_com%now%track_len) then
  d_radius = capillary_photon_d_radius (photon_com%now, ele_com)
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
d_radius = capillary_photon_d_radius (photon1_com%now, ele_com) 

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
!   photon   -- Photon_track_struct: Output coordinates.
!   absorbed -- Logical: If true the photon has been absorbed.
!-

subroutine capillary_reflect_photon (photon, ele, absorbed)

implicit none

type (ele_struct), target :: ele
type (photon_track_struct), target :: photon

real(rp) perp(3), r, graze_angle, cos_perp
real(rp), pointer :: vec(:)

logical absorbed

! perp is a vector perpendicular to the surface tangent plane

photon%old = photon%now
r = capillary_photon_d_radius (photon%now, ele, perp)

! The component of the photon velocity that is perpendicular to the surface 
! tangent gets reflected.

vec => photon%now%orb%vec
cos_perp = dot_product (vec(2:6:2), perp)
vec(2:6:2) = vec(2:6:2) - 2 * cos_perp * perp

! Check for absorbtion if the graze angle is too large.

absorbed = .false.
graze_angle = pi/2 - acos(cos_perp)
if (graze_angle > ele%value(critical_angle_factor$) / photon%now%energy) absorbed = .true.

end subroutine capillary_reflect_photon

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function capillary_photon_d_radius (p_orb, ele, perp, err_flag) result (d_radius)
!
! Routine to calculate the normalized radius = photon_radius - wall_radius.
!
! Input:
!   p_orb   -- photon_coord_struct: Input coordinates.
!   ele     -- ele_struct: Capillary element
!
! Output:
!   d_radius  -- Real(rp), Normalized radius: r_photon - r_wall
!   perp(3)   -- Real(rp), optional: Perpendicular normal to the wall.
!   err_flag  -- Logical, optional: Set True if error, false otherwise.
!-

function capillary_photon_d_radius (p_orb, ele, perp, err_flag) result (d_radius)

implicit none

type (ele_struct), target :: ele
type (photon_coord_struct), target :: p_orb
type (wall3d_section_struct), pointer :: sec1, sec2

real(rp) d_radius, r_photon, s_rel, spline, cos_theta, sin_theta
real(rp) r1_wall, r2_wall, dr1_dtheta, dr2_dtheta, f_eff, ds
real(rp) p1, p2, dp1, dp2
real(rp), optional :: perp(3)
real(rp), pointer :: vec(:)

integer ix, n_slice, n_sec
logical, optional :: err_flag

character(32), parameter :: r_name = 'capillary_photon_d_radius' 

! The outward normal vector is discontinuous at the wall points.
! If at a wall point, use the correct part of the wall.

vec => p_orb%orb%vec

if (present(err_flag)) err_flag = .true.

n_sec = size(ele%wall3d%section)
if (vec(5) < ele%wall3d%section(1)%s .or. vec(5) > ele%wall3d%section(n_sec)%s) then
  call out_io (s_error$, r_name, 'PHOTON S-POSITION BEYOUND WALL: \ES14.3\ ', vec(5))
  d_radius = -1
  return
endif

call bracket_index (ele%wall3d%section%s, 1, size(ele%wall3d%section), vec(5), ix)
if (ix == n_sec) ix = n_sec - 1
if (vec(5) == ele%wall3d%section(ix)%s .and. vec(6) > 0 .and. ix /= lbound(ele%wall3d%section, 1)) ix = ix - 1
p_orb%ix_section = ix

! sec1 and sec2 are the cross-sections to either side of the photon.

sec1 => ele%wall3d%section(ix)
sec2 => ele%wall3d%section(ix+1)

if (vec(1) == 0 .and. vec(3) == 0) then
  r_photon = 0
  cos_theta = 1
  sin_theta = 0
else
  r_photon = sqrt(vec(1)**2 + vec(3)**2)
  cos_theta = vec(1) / r_photon
  sin_theta = vec(3) / r_photon
endif

! Calculate the radius values at the cross-sections.

call calc_wall_radius (sec1%v, cos_theta, sin_theta, r1_wall, dr1_dtheta)
call calc_wall_radius (sec2%v, cos_theta, sin_theta, r2_wall, dr2_dtheta)

! Interpolate to get d_radius

ds = sec2%s - sec1%s
s_rel = (vec(5) - sec1%s) / ds
p1 = 1 - s_rel + sec1%p1_coef(1)*s_rel + sec1%p1_coef(2)*s_rel**2 + sec1%p1_coef(3)*s_rel**3
p2 =     s_rel + sec1%p2_coef(1)*s_rel + sec1%p2_coef(2)*s_rel**2 + sec1%p2_coef(3)*s_rel**3

d_radius = r_photon - (p1 * r1_wall + p2 * r2_wall)

! Calculate the surface normal vector

if (present (perp)) then
  perp(1:2) = [cos_theta, sin_theta] - [-sin_theta, cos_theta] * &
                        (p1 * dr1_dtheta + p2 * dr2_dtheta) / r_photon
  dp1 = -1 + sec1%p1_coef(1) + 2 * sec1%p1_coef(2)*s_rel + 3 * sec1%p1_coef(3)*s_rel**2
  dp2 =  1 + sec1%p2_coef(1) + 2 * sec1%p2_coef(2)*s_rel + 3 * sec1%p2_coef(3)*s_rel**2
  perp(3)   = -(dp1 * r1_wall + dp2 * r2_wall) / ds
  perp = perp / sqrt(sum(perp**2))  ! Normalize vector length to 1.
endif

if (present(err_flag)) err_flag = .false.

end function capillary_photon_d_radius

end module
