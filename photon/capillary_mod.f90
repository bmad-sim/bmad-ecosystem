module capillary_mod

use wall3d_mod

type photon_coord_struct
  type (coord_struct) orb       ! Phase space: orb%vec = (x, vx/c, y, vy/c, s, vs/c)
  real(rp) track_len            ! Total track length from the start of the element.
  integer ix_section            ! Cross section index
end type

type photon_track_struct
  type (photon_coord_struct) old, now
end type

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_a_capillary (orb, ele)
!
! Routine to track through a capillary.
!
! Input:
!   orb -- Coord_struct: Input photon coordinates.
!   ele -- ele_struct: Capillary element
!
! Output:
!   orb  -- Coord_struct: Output photon coordinates.
!-

subroutine track_a_capillary (orb, ele)

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (photon_coord_struct), pointer :: now
type (photon_track_struct), target :: photon

real(rp) pz

! Calculate maximum step size

! Setup photon coords

now => photon%now
now%orb = orb
now%orb%vec(5) = 0
now%orb%location = inside$
now%track_len = 0
now%ix_section = 1

! Loop over all bounces

do
  call capillary_track_photon_to_wall (photon, ele)
  if (now%orb%location /= inside$ .or. now%orb%state /= alive$) exit  ! At end or lost
  call capillary_reflect_photon (photon, ele)
  if (now%orb%state /= alive$) exit
enddo

! Cleanup

orb = now%orb

end subroutine track_a_capillary

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine capillary_track_photon_to_wall (photon, ele)
!
! Routine to track through a capillary.
!
! Input:
!   photon -- photon_track_struct: Input coordinates.
!   ele    -- ele_struct: Capillary element
!
! Output:
!   photon -- Photon_track_struct: Output coordinates.
!-

subroutine capillary_track_photon_to_wall (photon, ele)

implicit none

type (ele_struct), target :: ele
type (photon_track_struct), target :: photon
type (wall3d_section_struct), pointer :: section

real(rp) dr_max, ds_max, p_max, dlen
real(rp), pointer :: vec(:)

character(40) :: r_name = 'capillary_track_photon_to_wall'

! Check if outside wall

if (wall3d_d_radius (photon%now%orb%vec, ele, 1, ix_section = photon%now%ix_section) >= 0) then
  photon%now%orb%state = lost$
  return
endif

! propagate

do

  section => ele%wall3d(1)%section(photon%now%ix_section)
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

  if (wall3d_d_radius (photon%now%orb%vec, ele, 1, ix_section = photon%now%ix_section) > 0) then
    call capillary_photon_hit_spot_calc (photon, ele)
    return
  endif

  if (vec(5) == 0) then
    photon%now%orb%location = upstream_end$
    return
  endif

  if (vec(5) == ele%value(l$)) then
    photon%now%orb%location = downstream_end$
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

section => ele%wall3d(1)%section
ixc = bracket_index(vec(5), section%s, 1)

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

use super_recipes_mod, only: super_zbrent

implicit none

type (ele_struct), target :: ele
type (photon_track_struct), target :: photon
type (photon_track_struct) :: photon1

real(rp) d_rad, track_len0, track_len

integer status, i

character(40) :: r_name = 'capillary_photon_hit_spot_calc'

! Bracket the hit point. 
! Note: After the first reflection, the photon will start at the wall so
! if photon%old is at the wall we must avoid bracketing this point.

photon1 = photon

track_len0 = (photon%now%track_len + photon%old%track_len) / 2
do i = 1, 30
  d_rad = photon_hit_func(track_len0, status)
  if (d_rad < 0) exit
  track_len0 = (track_len0 + photon%old%track_len) / 2
  if (i == 30) then
    call out_io (s_abort$, r_name, 'ERROR: CANNOT FIND HIT SPOT REGION LOWER BOUND IN ELE: '//trim(ele%name) )
    if (global_com%exit_on_error) call err_exit
  endif
enddo

! Find where the photon hits.

track_len = super_zbrent (photon_hit_func, track_len0, photon%now%track_len, 1e-15_rp, 1d-10, status)

! Cleanup

photon%now = photon%old
call capillary_propagate_photon_a_step (photon, ele, track_len-photon%now%track_len, .false.)



contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function photon_hit_func (track_len, status) result (d_radius)
! 
! Routine to be used as an argument in zbrent in the capillary_photon_hit_spot_calc.
!
! Input:
!   track_len -- Real(rp): Place to position the photon.
!
! Output:
!   d_radius  -- Real(rp): r_photon - r_wall
!   status    -- integer: Not set.
!-
function photon_hit_func (track_len, status) result (d_radius)

real(rp), intent(in) :: track_len
real(rp) :: d_radius
real(rp) radius, d_track
integer status

! Easy case

if (track_len == photon%now%track_len) then
  d_radius = wall3d_d_radius (photon%now%orb%vec, ele)
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
d_radius = wall3d_d_radius (photon1%now%orb%vec, ele) 

end function photon_hit_func

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
!   photon   -- Photon_track_struct: Output coordinates.
!   photon%now%orb%state -- Set to lost$ if absorbed.
!-

subroutine capillary_reflect_photon (photon, ele)

implicit none

type (ele_struct), target :: ele
type (photon_track_struct), target :: photon

real(rp) perp(3), r, graze_angle, cos_perp, energy
real(rp), pointer :: vec(:)

! perp is a vector perpendicular to the surface tangent plane

photon%old = photon%now
r = wall3d_d_radius (photon%now%orb%vec, ele, 1, perp)

! The component of the photon velocity that is perpendicular to the surface 
! tangent gets reflected.

vec => photon%now%orb%vec
cos_perp = dot_product (vec(2:6:2), perp)
vec(2:6:2) = vec(2:6:2) - 2 * cos_perp * perp

! Check for absorption if the graze angle is too large.

graze_angle = pi/2 - acos(cos_perp)
energy = photon%now%orb%p0c 
if (graze_angle > ele%value(critical_angle_factor$) / energy) photon%now%orb%state = lost$

end subroutine capillary_reflect_photon

end module
