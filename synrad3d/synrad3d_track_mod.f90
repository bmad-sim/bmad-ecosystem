module synrad3d_track_mod

use synrad3d_utils
use photon_reflection_mod

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine track_photon (photon, lat, wall)
!
! Routine to propagate a synch radiation photon until it gets absorbed by a wall.
!
! Modules needed:
!   use synrad3d_track_mod
!
! Input:
!   photon    -- photon3d_coord_struct: photon with starting parameters set.
!     %start    -- Starting coords.
!   lat       -- lat_struct: with twiss propagated and mat6s made
!   wall      -- wall3d_struct: Beam chamber walls
!
! Output:
!   photon    -- photon3d_coord_struct: synch radiation photon propagated until absorbtion.
!-

subroutine track_photon (photon, lat, wall)

implicit none

type (lat_struct), target :: lat
type (photon3d_track_struct), target :: photon
type (wall3d_struct), target :: wall

logical absorbed

!

photon%start%track_len = 0
photon%now = photon%start

do
  call track_photon_to_wall (photon, lat, wall)
  call reflect_photon (photon, wall, absorbed)
  if (absorbed) return
  photon%n_reflect = photon%n_reflect + 1
enddo

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
!   use synrad3d_track_mod
!
! Input:
!   photon    -- photon3d_coord_struct: photon with starting parameters set
!   lat       -- lat_struct: with twiss propagated and mat6s made
!   wall      -- wall3d_struct: Beam chamber walls
!
! Output:
!   photon    -- photon3d_coord_struct: synch radiation photon propagated to wall
!-

subroutine track_photon_to_wall (photon, lat, wall)

implicit none

type (lat_struct), target :: lat
type (photon3d_track_struct), target :: photon
type (wall3d_struct), target :: wall

real(rp) v_rad_max, dlen, radius
real(rp), pointer :: vec(:)

! The photon is tracked in a series of steps.

vec => photon%now%vec

do

  v_rad_max = max(abs(vec(2)), abs(vec(4)))
  if (synrad3d_params%dr_track_step_max * abs(vec(6)) > &
      synrad3d_params%ds_track_step_max * v_rad_max) then
    dlen = synrad3d_params%ds_track_step_max / abs(vec(6))
  else
    dlen = synrad3d_params%dr_track_step_max / v_rad_max
  endif

  call propagate_photon_a_step (photon, dlen, lat, wall, .true.)

  ! See if the photon has hit the wall.
  ! If so we calculate the exact hit spot where the photon crossed the
  ! wall boundry and return

  call photon_radius (photon%now, wall, radius)
  if (radius > 1) then
    call photon_hit_spot_calc (photon, wall, lat)
    return
  endif

enddo

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine propagate_photon_a_step (photon, dl_step, lat, wall, stop_at_check_pt)
!
! Routine to propagate a photon to a given spot
!
! Modules needed:
!   use synrad3d_track_mod
!
! Input:
!   photon  -- photon3d_coord_struct: Photon to track.
!   dl_step -- Real(rp): Distance to track. Note: the propagation distance may not be exact
!               when going long distances.
!   lat     -- lat_struct: Lattice to track through.
!   stop_at_check_pt 
!           -- Logical: If True, stop at a check point which is defined to be:
!                a) minimum x extremum in a bend, or
!                b) At wall point boundries.
!              Note: (b) guarantees that there will be check points at the ends of the lattice.
!
! Output:
!   photon  -- photon3d_coord_struct: 
!			%now       -- If the photon has hit, the photon position is adjusted accordingly.
!-

subroutine propagate_photon_a_step (photon, dl_step, lat, wall, stop_at_check_pt)

implicit none

type (lat_struct), target :: lat
type (photon3d_track_struct), target :: photon
type (wall3d_struct) wall
type (photon3d_coord_struct), pointer :: now

real(rp) dl_step, dl_left, s_stop, denom, v_x, v_s, sin_t, cos_t
real(rp) g, new_x, radius, theta, tan_t, dl, dl2
real(rp), pointer :: vec(:)

integer ixw

logical stop_at_check_pt, s_stop_is_check_pt, will_stop_at_s_stop


! update old 

photon%old = photon%now  ! Save for hit spot calc
now => photon%now
dl_left = dl_step

! propagate the photon until we have gone a distance dl_step

propagation_loop: do

  s_stop_is_check_pt = .false.
  if (stop_at_check_pt) call bracket_index (wall%pt%s, 0, wall%n_pt_max, now%vec(5), ixw)

  ! If we are crossing over to a new element then update now%ix_ele.

  if (now%vec(6) > 0) then
    do
      if (now%vec(5) >= lat%ele(now%ix_ele)%s) then
        if (now%ix_ele == lat%n_ele_track) then
          now%vec(5) = now%vec(5) - lat%param%total_length
          now%ix_ele = 1
          photon%crossed_end = .not. photon%crossed_end
          exit
        endif
        now%ix_ele = now%ix_ele + 1
      elseif (now%vec(5) < lat%ele(now%ix_ele-1)%s) then
        now%ix_ele = now%ix_ele - 1
        if (now%ix_ele == 0) then
          print *, 'ERROR IN PROPAGATE_PHOTON: INTERNAL +ERROR'
          call err_exit
        endif
      else
        exit
      endif
    enddo

    s_stop = lat%ele(now%ix_ele)%s

    if (stop_at_check_pt .and. ixw < wall%n_pt_max) then
      if (wall%pt(ixw+1)%s < s_stop) then
        s_stop = wall%pt(ixw+1)%s
        s_stop_is_check_pt = .true.
      endif
    endif

  else   ! direction = -1
    do
      if (now%vec(5) <= lat%ele(now%ix_ele-1)%s) then
        if (now%ix_ele <= 0) then
          now%vec(5) = now%vec(5) + lat%param%total_length
          now%ix_ele = lat%n_ele_track
          photon%crossed_end = .not. photon%crossed_end
          exit
        endif
        now%ix_ele = now%ix_ele - 1
      elseif (now%vec(5) > lat%ele(now%ix_ele)%s) then
        now%ix_ele = now%ix_ele + 1
        if (now%ix_ele == lat%n_ele_track+1) then
          print *, 'ERROR IN PROPAGATE_PHOTON: INTERNAL -ERROR'
          call err_exit
        endif
      else
        exit
      endif
    enddo

    s_stop = lat%ele(now%ix_ele-1)%s

    if (stop_at_check_pt .and. ixw > 0) then
      if (wall%pt(ixw)%s == now%vec(5)) ixw = ixw - 1
      if (wall%pt(ixw)%s > s_stop) then
        s_stop = wall%pt(ixw)%s
        s_stop_is_check_pt = .true.
      endif
    endif

  endif

  ! Propagate the photon a step.

  ! In a bend...

  if (lat%ele(now%ix_ele)%key == sbend$ .and. lat%ele(now%ix_ele)%value(g$) /= 0) then

    ! Next position is determined by whether the distance to the element edge is 
    ! shorder than the distance left to travel.

    g = lat%ele(now%ix_ele)%value(g$)
    radius = 1 / g
    theta = (s_stop - now%vec(5)) * g
    tan_t = tan(theta)
    dl = tan_t * (radius + now%vec(2)) / (now%vec(6) - tan_t * now%vec(2))

    if (abs(tan_t * (radius + now%vec(2))) > dl_left * abs(now%vec(6) - tan_t * now%vec(2))) then
      dl = dl_left
      tan_t = (dl * now%vec(6)) / (radius + now%vec(1) + dl * now%vec(2))
      theta = atan(tan_t)
      will_stop_at_s_stop = .false.
    else
      dl = tan_t * (radius + now%vec(2)) / (now%vec(6) - tan_t * now%vec(2))
      will_stop_at_s_stop = .true.
    endif

    ! Check if we should actually be stopping at the extremum (minimal x)

    if (stop_at_check_pt .and. now%vec(2) < 0) then 
      dl2 = -now%vec(2) * (radius + now%vec(1)) / (now%vec(2)**2 + now%vec(6)**2)
      if (dl2 < dl) then
        dl = 1.000001 * dl2 ! Add extra to make sure we are not short due to roundoff.
        tan_t = (dl * now%vec(6)) / (radius + now%vec(1) + dl * now%vec(2))
        theta = atan(tan_t)
        s_stop = now%vec(5) + radius * theta
        will_stop_at_s_stop = .true.
        s_stop_is_check_pt = .true.
      endif
    endif

    ! Move to the stop point

    denom = sqrt((radius + now%vec(1) + dl * now%vec(2))**2 + (dl * now%vec(6))**2) 
    sin_t = (dl * now%vec(6)) / denom
    cos_t = (radius + now%vec(1) + dl * now%vec(2)) / denom
    v_x = now%vec(2); v_s = now%vec(6)
    now%vec(1) = denom - radius
    now%vec(2) = v_s * sin_t + v_x * cos_t
    now%vec(3) = now%vec(3) + dl * now%vec(4)
    now%vec(5) = now%vec(5) + radius * theta
    now%vec(6) = v_s * cos_t - v_x * sin_t

  ! Else we are not in a bend

  else

    ! Next position

    if (abs(now%vec(6)) * dl_left > abs(s_stop - now%vec(5))) then
      dl = (s_stop - now%vec(5)) / now%vec(6)
      will_stop_at_s_stop = .true.
    else
      dl = dl_left
      will_stop_at_s_stop = .false.
    endif

    ! And move to the next position

    now%vec(1) = now%vec(1) + dl * now%vec(2)
    now%vec(3) = now%vec(3) + dl * now%vec(4)
    now%vec(5) = now%vec(5) + dl * now%vec(6)
  endif

  if (will_stop_at_s_stop) now%vec(5) = s_stop ! Needed to avoid roundoff errors.

  !

  photon%now%track_len = photon%now%track_len + dl
  dl_left = dl_left - dl

  if (dl_left == 0) return
  if (stop_at_check_pt .and. will_stop_at_s_stop .and. s_stop_is_check_pt) return

enddo propagation_loop

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
!   use synrad3d_track_mod
!
! Input:
!   photon  -- photon3d_coord_struct:
!   wall    -- wall3d_struct: 
!   lat     -- lat_struct: Lattice
!
! Output:
!   photon  -- photon3d_coord_struct: 
!			%now       -- If the photon has hit, the photon position is adjusted accordingly.
!-

subroutine photon_hit_spot_calc (photon, wall, lat)

implicit none

type (lat_struct) lat
type (photon3d_track_struct) :: photon, photon0, photon1, photon2
type (wall3d_struct), target :: wall
type (wall3d_pt_struct), pointer :: wall_pt

integer ix_wall, ix0, ix1, ix2, i

real(rp) del0, del1, del2, dl, radius

! Find where the photon hits.
! we need to iterate in a bend since the wall is actually curved.

photon0 = photon
photon0%now = photon%old
photon1 = photon0
photon2 = photon

call photon_radius (photon0%now, wall, radius)
del0 = radius - 1 ! Must be negative

call photon_radius (photon2%now, wall, radius)
del2 = radius - 1 ! Must be positive

do i = 1, 20

  if (abs(del0) < 1.0e-4) then
    photon1 = photon0
    exit
  elseif (abs(del2) < 1.0e-4) then
    photon1 = photon2
    exit
  endif

  if (i == 20) then
    print *, 'ERROR IN PHOTON_HIT_SPOT_CALC: CALCULATION IS NOT CONVERGING'
    call err_exit
  endif

  ! Try linear interpolation

  dl = -del0 * (photon2%now%track_len - photon0%now%track_len) / (del2 - del0)
  call propagate_photon_a_step (photon1, dl, lat, wall, .false.)

  call photon_radius (photon1%now, wall, radius)
  del1 = radius - 1

  if (del1 < 0) then
    photon0 = photon1; del0 = del1
  elseif (del1 > 0) then
    photon2 = photon1; del2 = del1
    photon1 = photon0
  endif

  ! Linear interpolation will fail badly if the photon is doing from one side of the 
  ! chamber to another. So also do a divide in half.

  dl = (photon2%now%track_len - photon0%now%track_len) / 2
  call propagate_photon_a_step (photon1, dl, lat, wall, .false.)

  call photon_radius (photon1%now, wall, radius)
  del1 = radius - 1

  if (del1 < 0) then
    photon0 = photon1; del0 = del1
  elseif (del1 > 0) then
    photon2 = photon1; del2 = del1
    photon1 = photon0
  endif

enddo

! cleanup...

photon = photon1

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine reflect_photon (photon, wall, absorbed)
!
! Routine to reflect a photon off of the wall.
!-

subroutine reflect_photon (photon, wall, absorbed)

implicit none

type (photon3d_track_struct), target :: photon
type (wall3d_struct), target :: wall
type (wall3d_pt_struct), pointer :: wall0, wall1

real(rp) dx_parallel, dy_parallel, dx_perp0, dy_perp0
real(rp) dx_perp1, dy_perp1, dx_perp, dy_perp, denom, f
real(rp) dot_parallel, dot_perp, r, graze_angle, reflectivity
real(rp), pointer :: vec(:)

integer ix

logical absorbed

! Check if reflections allowed

if (.not. synrad3d_params%allow_reflections) then
  absorbed = .true.
  return
endif

! Get the wall index for this section of the lattice

vec => photon%now%vec
call bracket_index (wall%pt%s, 0, wall%n_pt_max, vec(5), ix)
if (ix == wall%n_pt_max) ix = wall%n_pt_max - 1

wall0 => wall%pt(ix)
wall1 => wall%pt(ix+1)

! (dx_perp, dy_perp) is the normalized vector perpendicular to the wall
! at the photon hit point.

if (wall0%type == 'rectangular') then
  if (abs(vec(1)/wall0%width2) > abs(vec(3)/wall0%height2)) then
    dx_perp0 = vec(1)
    dy_perp0 = 0
  else
    dx_perp0 = 0
    dy_perp0 = vec(3)
  endif
elseif (wall0%type == 'elliptical') then
  dx_perp0 = wall1%height2**2 * vec(1)
  dy_perp0 = wall1%width2**2 * vec(3)
else
  print *, 'BAD WALL%TYPE: ' // wall0%type, ix
  call err_exit
endif

denom = sqrt(dx_perp0**2 + dy_perp0**2)
dx_perp0 = dx_perp0 / denom
dy_perp0 = dy_perp0 / denom

if (wall1%type == 'rectangular') then
  if (abs(vec(1)/wall1%width2) > abs(vec(3)/wall1%height2)) then
    dx_perp1 = vec(1)
    dy_perp1 = 0
  else
    dx_perp1 = 0
    dy_perp1 = vec(3)
  endif
elseif (wall1%type == 'elliptical') then
  dx_perp1 = wall1%height2**2 * vec(1)
  dy_perp1 = wall1%width2**2  * vec(3)
else
  print *, 'BAD WALL%TYPE: ' // wall1%type, ix+1
  call err_exit
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

! absorbtion

graze_angle = pi/2 - acos(abs(vec(2) * dx_perp + vec(4) * dy_perp))
call photon_reflectivity (graze_angle, photon%now%energy, reflectivity)
call ran_uniform(r)
absorbed = (r > reflectivity)
if (absorbed) return  ! Do not reflect if absorbed

! Reflect the ray.
! The perpendicular component gets reflected and the parallel component is invarient.

dot_parallel = dx_parallel * vec(2) + dy_parallel * vec(4)
dot_perp     = dx_perp     * vec(2) + dy_perp     * vec(4)

vec(2) = dot_parallel * dx_parallel - dot_perp * dx_perp
vec(4) = dot_parallel * dy_parallel - dot_perp * dy_perp

end subroutine

end module
