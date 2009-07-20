module photon_track_mod

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
! Subroutine photon_specular_reflect (photon, wall_pt)
!
! Routine to specularly reflect a photon off of the wall.
!-

subroutine photon_specular_reflect (photon, wall_pt)

implicit none

type (photon_coord_struct) photon
type (wall_2d_pt_struct) wall_pt

real(rp) dx_parallel, dy_parallel, dx_perp, dy_perp, denom

! The wall is assumed elliptical and that the photon is at the wall
! (dx_parallel, dy_parallel) is the normalized vector parallel to the wall
! at the photon hit point.

dx_parallel =  wall_pt%width2**2 * photon%now%vec(3)
dy_parallel = -wall_pt%height2**2 * photon%now%vec(1)

denom = sqrt(dx_parallel**2 + dy_parallel**2)
dx_parallel = dx_parallel / denom
dy_parallel = dy_parallel / denom

! (dx_perp, dy_perp) is the normalized vector perpendicular to the wall
! at the photon hit point.

dx_perp = -dy_parallel
dy_perp =  dx_parallel

! The perpendicular component gets reflected and the parallel component is invarient.

photon%now%vec(1) = (dx_parallel - dx_perp) * photon%now%vec(1)
photon%now%vec(3) = (dy_parallel - dy_perp) * photon%now%vec(3)

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
!   wall      -- wall_2d_struct: Beam chamber walls
!
! Output:
!   photon    -- photon_coord_struct: synch radiation photon propagated to wall
!-

subroutine track_photon_to_wall (photon, lat, wall)

implicit none

type (lat_struct), target :: lat
type (photon_coord_struct), target :: photon
type (wall_2d_struct), target :: wall
type (stop_pt_struct) stop_pt

real(rp) s_next

logical is_hit

! The photon is tracked in a series of steps.
! Get the first stopping point.

call get_initial_stop_pt (photon, wall, lat, stop_pt)

! propagation loop:

do

  call propagate_photon (photon, stop_pt%s, lat, .true.)

  ! See if the photon has hit the wall.
  ! If so we calculate the exact hit spot where the photon crossed the
  ! wall boundry and return

  call photon_hit_spot_calc (photon, wall, lat, is_hit)
  if (is_hit) return

  call get_next_stop_pt (photon, wall, lat, stop_pt)

enddo

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! subroutine get_initial_stop_pt (photon, wall, lat, stop_pt)
!
! subroutine to get the first point to stop at when tracking a photon.
!
! Modules needed:
!   use photon_mod
!
! Input:
!   photon  -- photon_coord_struct:
!   wall    -- wall_2d_struct: Beam chamber wall.
!   lat     -- lat_struct: Lattice.
!
! Output:
!   stop_pt -- stop_pt_struct: Point to stop at
!-

subroutine get_initial_wall_pt (photon, wall, lat, stop_pt)

implicit none

type (photon_coord_struct) photon
type (wall_2d_struct) wall
type (lat_struct) lat
type (stop_pt_struct) stop_pt

integer ix0, ix1, ix2

! Error check

if (wall%n_pt_max == 0) then
  print *, 'THE WALL HAS NOT BEEN DEFINED!'
  call err_exit
endif

! The stop point will be at the next downstream wall point or the next element boundary
! whichever is closest.

! edge cases

if (photon%now%vec(5) == lat%param%total_length) then
  stop_pt%s = lat%param%total_length
  stop_pt%ix_wall = wall%n_pt_max
  if (photon%direction == 1) then
    stop_pt%ix_ele = 1
  else
    stop_pt%ix_ele = lat%n_ele_track 
  endif
  return
endif

if (photon%now%vec(5) == 0) then
  stop_pt%s = 0
  stop_pt%ix_wall = 0
  if (photon%direction == 1) then
    stop_pt%ix_ele = 1
  else
    stop_pt%ix_ele = lat%n_ele_track 
  endif
  return
endif

! Calc next element s.

call ele_at_s (lat, photon%now%vec(5), stop_pt%ix_ele)

if (lat%ele(stop_pt%ix_ele)%s == photon%now%vec(5)) then
  if (photon%direction == -1) then
    stop_pt%ix_ele = stop_pt%ix_ele - 1 
  endif
endif

! Calc next wall s 

ix0 = 0
ix2 = wall%n_pt_max

do
  ix1 = (ix2 + ix0) / 2
  if (wall%pt(ix1)%s < photon%now%vec(5)) then
    ix0 = ix1
  elseif (wall%pt(ix1)%s > photon%now%vec(5)) then
    ix2 = ix1
  elseif (photon%direction == 1) then   ! here wall%pt(ix)%s == photon%now%vec(5)
    ix2 = ix1
  else
    ix0 = ix1
  endif
  if (ix2 - ix0 == 1) then
    if (photon%direction == 1) then
      stop_pt%ix_wall = ix2
    else
      stop_pt%ix_wall = ix0
    endif
    exit
  endif
enddo

! and compare the two points.

if (photon%direction == 1) then
  stop_pt%s = min(lat%ele(stop_pt%ix_ele)%s, wall%pt(stop_pt%ix_wall)%s)
else
  stop_pt%s = max(lat%ele(stop_pt%ix_ele-1)%s, wall%pt(stop_pt%ix_wall)%s)
endif

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine get_next_stop_pt (photon, wall, lat, stop_pt)
!
! Routine to get the next point in increasing or decreasing s depending 
! upon photon%direction.
!
! Modules needed:
!   use photon_mod
!
! Input:
!   photon  -- photon_coord_struct:
!   wall    -- wall_2d_struct: 
!   lat     -- lat_struct: Lattice
!   stop_pt -- stop_pt_struct: Present point
!
! Output:
!   stop_pt -- stop_pt_struct: Point to stop at
!-

subroutine get_next_wall_pt (photon, wall, lat, stop_pt)

implicit none

type (photon_coord_struct) photon
type (wall_2d_struct) wall
type (stop_pt_struct) stop_pt
type (lat_struct) lat

integer direct
real(rp) s_now

! wrap around cases

direct = photon%direction
s_now = photon%now%vec(5)

if (stop_pt%s == 0 .and. direct == -1) then
  stop_pt%s = lat%param%total_length
  stop_pt%ix_wall = wall%n_pt_max
  stop_pt%ix_ele = lat%n_ele_track 
  return
endif

if (stop_pt%s == lat%param%total_length .and. direct == 1) then
  stop_pt%s = 0
  stop_pt%ix_wall = 0
  stop_pt%ix_ele = 1
  return
endif

! Normal case
! Next stop point is nearest element boundary or wall point.

do
  if (direct == -1) then
    if (lat%ele(stop_pt%ix_ele-1)%s < s_now) exit
    stop_pt%ix_ele = stop_pt%ix_ele - 1
  else
    if (lat%ele(stop_pt%ix_ele)%s > s_now) exit
    stop_pt%ix_ele = stop_pt%ix_ele + 1
  endif
enddo

stop_pt%ix_wall = stop_pt%ix_wall + direct

if (photon%direction == 1) then
  stop_pt%s = min(lat%ele(stop_pt%ix_ele)%s, wall%pt(stop_pt%ix_wall)%s)
else
  stop_pt%s = max(lat%ele(stop_pt%ix_ele-1)%s, wall%pt(stop_pt%ix_wall)%s)
endif

end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine photon_hit_spot_calc (photon, wall, lat, has_hit)
!
! Routine to calculate if the photon has hit the wall and, if so, where.
!
! Modules needed:
!   use photon_mod
!
! Input:
!   photon  -- photon_coord_struct:
!   wall    -- wall_2d_struct: 
!   lat     -- lat_struct: Lattice
!   stop_pt -- stop_pt_struct: Present point.
!
! Output:
!		has_hit -- Logical: Set True if photon has hit the wall.
!   photon  -- photon_coord_struct: 
!			%now       -- If the photon has hit, the photon position is adjusted accordingly.
!-

subroutine photon_hit_spot_calc (photon, wall, lat, has_hit)

implicit none

type (lat_struct) lat
type (photon_coord_struct) :: photon, photon0, photon1, photon2
type (wall_2d_struct), target :: wall
type (wall_2d_pt_struct), pointer :: wall_pt

integer ix_wall, ix0, ix1, ix2, i

real(rp) r_now, r_old, r0, r1, r2
real(rp) del_s, s1, s_now, s_old
real(rp) del0, del1, del2

logical has_hit, circular

! Find chamber size at present point and see if there is a hit.

has_hit = .false.  ! assume no hit

s_now = photon%now%vec(5)
call wall_at_s (wall, s_now, wall_pt)
r_now = (photon%now%vec(1)/wall_pt%width2)**2 + (photon%now%vec(3)/wall_pt%height2)**2
if (r_now < 1) return

! Find where the photon hits.
! we need to iterate in a bend since the wall is actually curved.

s_old = photon%old%vec(5)
call wall_at_s (wall, s_old, wall_pt)
r_old = (photon%old%vec(1)/wall_pt%width2)**2 + (photon%old%vec(3)/wall_pt%height2)**2

photon0 = photon; r0 = r_now
photon1 = photon; r1 = r_now
photon2 = photon; r2 = r_now

if (photon%now%vec(5) < photon%old%vec(5)) then
  photon2%now = photon%old
	r2 = r_old
elseif (photon%now%vec(5) > photon%old%vec(5)) then
  photon0%now = photon%old
	r0 = r_old
endif

do i = 1, 20

  if (abs(r0 - 1) < 1.0e-4) then
    photon1 = photon0
    exit
  elseif (abs(r2-1) < 1.0e-4) then
    photon1 = photon2
    exit
  endif

  if (i == 20) then
    print *, 'ERROR IN PHOTON_HIT_SPOT_CALC: CALCULATION IS NOT CONVERGING'
    call err_exit
  endif

	del0 = sqrt(r0) - 1
	del1 = sqrt(r1) - 1
	del2 = sqrt(r2) - 1

  s1 = (del2 * photon0%now%vec(5) - del0 * photon2%now%vec(5)) / (del2 - del0)

  if (s1 < photon1%now%vec(5)) then
    photon1%direction = -1
  else
    photon1%direction = +1
  endif
  call propagate_photon (photon1, s1, lat, .false.)

  call wall_at_s (wall, s1, wall_pt)
	r1 = (photon1%now%vec(1)/wall_pt%width2)**2 + (photon1%now%vec(3)/wall_pt%height2)**2
  del1 = sqrt(r1) - 1

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
! Subroutine propagate_photon (photon, s_end, lat, stop_at_extremum)
!
! Routine to propagate a photon to a given spot
!
! Modules needed:
!   use photon_mod
!
! Input:
!   photon  -- photon_coord_struct: Photon to track
!   s_end   -- Real(rp): Longitu
!   lat     -- lat_struct: Lattice to track through

!
! Output:
!		has_hit -- Logical: Set True if photon has hit the wall.
!   photon  -- photon_coord_struct: 
!			%now       -- If the photon has hit, the photon position is adjusted accordingly.

subroutine propagate_photon (photon, s_end, lat, stop_at_extremum)

implicit none

type (lat_struct), target :: lat
type (photon_coord_struct), target :: photon

real(rp) s_end, s_next, del_s, s_target
real(rp) g, new_x, theta0, theta1, c_t0, c_t1

logical stop_at_extremum

! find the target

s_target = s_end

if ((s_target - photon%now%vec(5)) * photon%direction < 0) s_target = &
        s_target + photon%direction * lat%param%total_length

! update old (but only if we have moved or gone through the IP)

if (s_target == photon%now%vec(5) .and. &
          abs(s_end - photon%now%vec(5)) /= lat%param%total_length) return

photon%old = photon%now

! propagate the photon until we get to s_end

propagation_loop: do

  ! If we are crossing over to a new element then update photon%ix_ele.
  ! Additionally, if we cross the lattice end we need to 
  ! reset photon%now%vec(5) and photon%ix_ele.

  if (photon%direction == 1) then
    do
      if (photon%now%vec(5) .ge. lat%ele(photon%ix_ele)%s) then
        photon%ix_ele = photon%ix_ele + 1
        if (photon%ix_ele > lat%n_ele_track) then
          photon%ix_ele = 1
          photon%now%vec(5) = photon%now%vec(5) - lat%param%total_length
          s_target = s_target - lat%param%total_length
          photon%crossed_end = .not. photon%crossed_end
        endif
      elseif (photon%now%vec(5) .lt. lat%ele(photon%ix_ele-1)%s) then
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
      if (photon%now%vec(5) .le. lat%ele(photon%ix_ele-1)%s) then
        photon%ix_ele = photon%ix_ele - 1
        if (photon%ix_ele .le. 0) then
          photon%ix_ele = lat%n_ele_track
          photon%now%vec(5) = photon%now%vec(5) + lat%param%total_length
          s_target = s_target + lat%param%total_length
          photon%crossed_end = .not. photon%crossed_end
        endif
      elseif (photon%now%vec(5) .gt. lat%ele(photon%ix_ele)%s) then
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

  ! Next position

  if (photon%direction == 1) then
    s_next = min (s_target, lat%ele(photon%ix_ele)%s)
  else
    s_next = max (s_target, lat%ele(photon%ix_ele-1)%s)
  endif

  del_s = s_next - photon%now%vec(5)

  ! In a bend: Exact formula is:
  ! new_x = (rho * (cos(theta0) - cos(theta1)) + photon%now%vec(1) * cos(theta0)) cos(theta1)
  ! We stop then theta = 0 since that is an extremum.

  if (lat%ele(photon%ix_ele)%key == sbend$ .and. lat%ele(photon%ix_ele)%value(g$) /= 0) then
    g = lat%ele(photon%ix_ele)%value(g$)
    theta0 = photon%now%vec(2)
    theta1 = photon%now%vec(2) + del_s * g
    ! if theta has changed sign then there is an extremum 
    if (stop_at_extremum .and. theta0 * theta1 < 0) then 
      del_s = -theta0 / g            ! step to extremum
      theta1 = 0
      s_next = photon%now%vec(5) + del_s
      s_target = s_next
    endif
    c_t0 = -(theta0**2)/2 + theta0**4/24
    c_t1 = -(theta1**2)/2 + theta1**4/24
    new_x = ((c_t0 - c_t1) / g + photon%now%vec(1) * cos(theta0)) / cos(theta1)
    photon%now%vec(1) = new_x
    photon%now%vec(2) = theta1
  else
    photon%now%vec(1) = photon%now%vec(1) + del_s * tan(photon%now%vec(2))
  endif

  photon%now%vec(3) = photon%now%vec(3) + del_s * tan(photon%now%vec(4))

  photon%track_len = photon%track_len + abs(del_s)
  photon%now%vec(5) = s_next

  if (photon%crossed_end .and. lat%param%lattice_type == linear_lattice$) return

  if (s_next == s_target) return

enddo propagation_loop

end subroutine


end module
