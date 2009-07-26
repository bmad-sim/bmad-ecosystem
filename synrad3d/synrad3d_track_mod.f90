module synrad3d_track_mod

use photon_utils

type stop_pt_struct
  real(rp) s      ! Longitudinal position.
  integer ix_wall ! wall point index who's s position is at or just less than s_pos.
  integer ix_ele  ! element index in lat%ele(:) array at s_pos.
                  ! If at boundary, ix_ele points to the downstream element.
end type

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine photon_wall_reflect (p_orb, wall_pt)
!
! Routine to reflect a photon off of the wall.
!-

subroutine photon_wall_reflect (p_orb)

implicit none

type (photon_coord_struct) p_orb
type (wall_3d_pt_struct) wall_pt

real(rp) dx_parallel, dy_parallel, dx_perp0, dy_perp0
real(rp) dx_perp1, dy_perp1, dx_perp, dy_perp, denom

! Get the wall index for this section of the lattice

vec => p_orb%vec
call bracket_index (wall%pt%s, 0, wall%n_pt_max, vec(5), ix)
if (ix == wall%n_pt_max) ix = wall%n_pt_max - 1

wall0 => wall%pt(ix)
wall1 => wall%pt(ix+1)

! (dx_perp, dy_perp) is the normalized vector perpendicular to the wall
! at the photon hit point.

if (wall0%type == rectangular$) then
  if (abs(vec(1)/wall0%width2) > abs(vec(3)/wall0%height2)) then
    dx_perp0 = p_orb%now%vec(1)
    dy_perp0 = 0
  else
    dx_perp0 = 0
    dy_perp0 = p_orb%now%vec(3)
  endif
else
  dx_perp0 = wall1%height2**2 * p_orb%now%vec(1)
  dy_perp0 = wall1%width2**2 * p_orb%now%vec(3)
endif

denom = sqrt(dx_perp0**2 + dy_perp0**2)
dx_perp0 = dx_perp0 / denom
dy_perp0 = dy_perp0 / denom

if (wall1%type == rectangular$) then
  if (abs(vec(1)/wall1%width2) > abs(vec(3)/wall1%height2)) then
    dx_perp1 = p_orb%now%vec(1)
    dy_perp1 = 0
  else
    dx_perp1 = 0
    dy_perp1 = p_orb%now%vec(3)
  endif
else
  dx_perp1 = wall1%height2**2 * p_orb%now%vec(1)
  dy_perp1 = wall1%width2**2 * p_orb%now%vec(3)
endif

denom = sqrt(dx_perp1**2 + dy_perp1**2)
dx_perp1 = dx_perp1 / denom
dy_perp1 = dy_perp1 / denom

! Average the vectors to get the vector at the photon point.

f = (vec(5) - wall0%s) / (wall1%s - wall0%s)
dx_perp = (1 - f) * dx_perp0 + f * dx_perp1
dy_perp = (1 - f) * dy_perp0 + f * dy_perp1

denom = sqrt(dx_perp**2 + dy_perp**2)
dx_perp = dx_perp / denom
dy_perp = dy_perp / denom


! (dx_parallel, dy_parallel) is the normalized vector parallel to the wall
! at the photon hit point.

dx_parallel =  dy_perp 
dy_parallel = -dx_perp

! The perpendicular component gets reflected and the parallel component is invarient.

p_orb%now%vec(2) = (dx_parallel - dx_perp) * p_orb%now%vec(2)
p_orb%now%vec(4) = (dy_parallel - dy_perp) * p_orb%now%vec(4)

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine track_photon_to_wall (photon, lat, wall)
!
! Routine to propagate a synch radiation photon until it hits a wall.
!
! Modules needed:
!   use photon_mod
!
! Input:
!   photon    -- photon_coord_struct: photon with starting parameters set
!   lat       -- lat_struct: with twiss propagated and mat6s made
!   wall      -- wall_3d_struct: Beam chamber walls
!
! Output:
!   photon    -- photon_coord_struct: synch radiation photon propagated to wall
!-

subroutine track_photon_to_wall (photon, lat, wall)

implicit none

type (lat_struct), target :: lat
type (photon_coord_struct), target :: photon
type (wall_3d_struct), target :: wall
type (stop_pt_struct) stop_pt

real(rp) s_next

! The photon is tracked in a series of steps.

vec => photon%now%vec

do

  vr = max(abs(photon%now%vec(2)), abs(photon%now%vec(4))
  if (synrad3d_params%dr_track_step_max * vec(6) > &
      synrad3d_params%ds_track_step_max * vr)) then
    dlen = synrad3d_params%ds_track_step_max / vec(6)
  else
    dlen = synrad3d_params%ds_track_step_max / vr
  endif

  call propagate_photon (photon, dlen, lat, .true.)

  ! See if the photon has hit the wall.
  ! If so we calculate the exact hit spot where the photon crossed the
  ! wall boundry and return

  call photon_radius (photon%now, wall)
  if (photon%now%radius > 1) then
    call photon_hit_spot_calc (photon, wall, lat)
    return
  endif

enddo

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine photon_hit_spot_calc (photon, wall, lat)
!
! Routine to calculate where the photon has hit the wall.
!
! Modules needed:
!   use photon_mod
!
! Input:
!   photon  -- photon_coord_struct:
!   wall    -- wall_3d_struct: 
!   lat     -- lat_struct: Lattice
!   stop_pt -- stop_pt_struct: Present point.
!
! Output:
!   photon  -- photon_coord_struct: 
!			%now       -- If the photon has hit, the photon position is adjusted accordingly.
!-

subroutine photon_hit_spot_calc (photon, wall, lat, has_hit)

implicit none

type (lat_struct) lat
type (photon_coord_struct) :: photon, photon0, photon1, photon2
type (wall_3d_struct), target :: wall
type (wall_3d_pt_struct), pointer :: wall_pt

integer ix_wall, ix0, ix1, ix2, i

real(rp) del_s, s1, s_now, s_old
real(rp) del0, del1, del2

! Find where the photon hits.
! we need to iterate in a bend since the wall is actually curved.

photon0 = photon
photon1 = photon
photon2 = photon

if (photon%now%vec(5) < photon%old%vec(5)) then
  photon2%now = photon%old
elseif (photon%now%vec(5) > photon%old%vec(5)) then
  photon0%now = photon%old
endif

do i = 1, 20

  if (abs(photon0%now%radius - 1) < 1.0e-4) then
    photon1 = photon0
    exit
  elseif (abs(photon2%now%radius-1) < 1.0e-4) then
    photon1 = photon2
    exit
  endif

  if (i == 20) then
    print *, 'ERROR IN PHOTON_HIT_SPOT_CALC: CALCULATION IS NOT CONVERGING'
    call err_exit
  endif

	del0 = sqrt(photon0%now%radius) - 1
	del1 = sqrt(photon1%now%radius) - 1
	del2 = sqrt(photon2%now%radius) - 1

  s1 = (del2 * photon0%now%vec(5) - del0 * photon2%now%vec(5)) / (del2 - del0)

  if (s1 < photon1%now%vec(5)) then
    photon1%direction = -1
  else
    photon1%direction = +1
  endif
  call propagate_photon (photon1, s1, lat, .false.)

  call photon_radius (photon%now, wall)
  del1 = sqrt(photon1%now%radius) - 1

  if (s1 < photon0%now%vec(5)) then
    photon2 = photon0; del2 = del0
    photon0 = photon1; del0 = del1
  elseif (s1 > photon2%now%vec(5)) then
    photon0 = photon2; del0 = del2
    photon2 = photon1; del2 = del1
  elseif (sign(1.0_rp, del0) == sign(1.0_rp, del1)) then
    photon0 = photon1; del0 = del1
  else
    photon2 = photon1; del2 = del1
  endif

enddo

! cleanup...
! We assume that the travel length cannot be greater 
! then half the circumference.

photon%now = photon1%now

del_s = photon%now%vec(5) - photon%start%vec(5)
if (del_s*photon%direction < 0) then
  photon%track_len = lat%param%total_length - abs(del_s)
  photon%crossed_end = .true.
else
  photon%track_len = abs(del_s)
  photon%crossed_end = .false.
endif

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine propagate_photon (photon, dl_step, lat, stop_at_extremum)
!
! Routine to propagate a photon to a given spot
!
! Modules needed:
!   use photon_mod
!
! Input:
!   photon  -- photon_coord_struct: Photon to track.
!   dl_step -- Real(rp): Distance to track. Note: the propagation distance may not be exact
!               when going long distances.
!   lat     -- lat_struct: Lattice to track through.
!   stop_at_extremum 
!           -- Logical: If True then stop at transverse extremum in a bend
!                or at the ends of the lattice
!
! Output:
!   photon  -- photon_coord_struct: 
!			%now       -- If the photon has hit, the photon position is adjusted accordingly.

subroutine propagate_photon (photon, dl_step, lat, stop_at_extremum)

implicit none

type (lat_struct), target :: lat
type (photon_coord_struct), target :: photon

real(rp) dl_step, dl_left
real(rp) g, new_x, theta0, theta1, c_t0, c_t1
real(rp), pointer :: vec(6)

logical stop_at_extremum


! update old 

photon%old = photon%now
dl_left = dl_step
vec => photon%now%vec

! propagate the photon until we have gone a distance dl_step

propagation_loop: do

  ! If we are crossing over to a new element then update photon%ix_ele.

  if (vec(6) > 0) then
    do
      if (vec(5) >= lat%ele(photon%ix_ele)%s) then
        photon%ix_ele = photon%ix_ele + 1
        if (photon%ix_ele > lat%n_ele_track) then
          photon%ix_ele = 1
          s_target = s_target - lat%param%total_length
          photon%crossed_end = .not. photon%crossed_end
        endif
      elseif (vec(5) < lat%ele(photon%ix_ele-1)%s) then
        photon%ix_ele = photon%ix_ele - 1
        if (photon%ix_ele == 0) then
          print *, 'ERROR IN PROPAGATE_PHOTON: INTERNAL +ERROR'
          call err_exit
        endif
      else
        exit
      endif
    enddo

  else   ! direction = -1
    do
      if (vec(5) <= lat%ele(photon%ix_ele-1)%s) then
        photon%ix_ele = photon%ix_ele - 1
        if (photon%ix_ele .le. 0) then
          photon%ix_ele = lat%n_ele_track
          s_target = s_target + lat%param%total_length
          photon%crossed_end = .not. photon%crossed_end
        endif
      elseif (vec(5) > lat%ele(photon%ix_ele)%s) then
        photon%ix_ele = photon%ix_ele + 1
        if (photon%ix_ele == lat%n_ele_track+1) then
          print *, 'ERROR IN PROPAGATE_PHOTON: INTERNAL -ERROR'
          call err_exit
        endif
      else
        exit
      endif
    enddo
  endif

  ! In a bend

  if (lat%ele(photon%ix_ele)%key == sbend$ .and. lat%ele(photon%ix_ele)%value(g$) /= 0) then

    g = lat%ele(photon%ix_ele)%value(g$)
    radius = 1 / g

    ! Next position is determined by whether the distance to the element edge is 
    ! shorder than the distance left to travel.

    if (vec(6) > 0) then
      theta = (lat%ele(photon%ix_ele)%s - vec(5)) * g
    else
      theta = (lat%ele(photon%ix_ele-1)%s - vec(5)) * g
    endif

    tan_t = tan(theta)
    dl = tan_t * (radius + vec(2)) / (vec(6) - tan_t * vec(2))
    if (abs(tan_t * (radius + vec(2))) > dl_left * abs(vec(6) - tan_t * vec(2))) then
      dl = dl_left
      tan_t = (dl * vec(6)) / (radius + vec(1) + dl * vec(2))
      theta = atan(tan_t)
    else
      dl = tan_t * (radius + vec(2)) / (vec(6) - tan_t * vec(2))
    endif

    ! Check if we should actually be stopping at the extremum (minimal x)

    if (stop_at_extremum .and. vec(2) < 0) then 
      dl2 = -vec(2) * (radius + vec(1)) / (vec(2)**2 + vec(6)**2)
      if (dl2 < dl) then
        dl = dl2
        tan_t = (dl * vec(6)) / (radius + vec(1) + dl * vec(2))
        theta = atan(tan_t)
      endif
    endif

    ! Move to the stop point

    denom = sqrt((radius + vec(1) + dl * vec(2))**2 + (dl * vec(6))**2) 
    sin_t = (dl * vec(6)) / denom
    cos_t = (radius + vec(1) + dl * vec(2)) / denom
    v_x = vec(2); v_s = vec(6)
    vec(1) = denom - radius
    vec(2) = v_s * sin_t + v_x * cos_t
    vec(3) = vec(3) + dl * vec(4)
    vec(5) = radius * theta
    vec(6) = v_s * cos_t - v_x * sin_t

  ! Else we are not in a bend

  else

    ! Next position

    if (vec(6) > 0) then
      if (vec(6) * dl_left > lat%ele(photon%ix_ele)%s - vec(5)) then
        dl = (lat%ele(photon%ix_ele)%s - vec(5)) / vec(6)
      else
        dl = dl_left
      endif
    else
      if (-vec(6) * dl_left > vec(5) - lat%ele(photon%ix_ele-1)%s) then
        dl = -(lat%ele(photon%ix_ele)%s - vec(5)) / vec(6)
      else
        dl = dl_left
      endif
    endif
  
    ! And move to the next position

    vec(1) = vec(1) + dl * vec(2)
    vec(3) = vec(3) + dl * vec(4)
    vec(5) = vec(5) + dl * vec(6)
  endif

  !

  photon%track_len = photon%track_len + abs(dl)
  dl_left = dl_left - dl

  if (photon%crossed_end .and. lat%param%lattice_type == linear_lattice$) return

  if (s_next == s_target) return

enddo propagation_loop

end subroutine


end module
