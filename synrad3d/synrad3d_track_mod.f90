module synrad3d_track_mod

use synrad3d_utils
use photon_reflection_mod

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_track_photon (photon, lat, wall)
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

subroutine sr3d_track_photon (photon, lat, wall)

implicit none

type (lat_struct), target :: lat
type (photon3d_track_struct), target :: photon
type (photon3d_coord_struct), allocatable :: p_temp(:)
type (wall3d_struct), target :: wall

integer n
logical absorbed

!

photon%start%track_len = 0
photon%now = photon%start

! if (.not. sr3d_params%debug_on) then
!   if (allocated(photon%reflect)) deallocate(photon%reflect)
!   allocate (photon%reflect(0:0))
!   photon%reflect(0) = photon%start
! endif

do

  call sr3d_track_photon_to_wall (photon, lat, wall)

  ! if (.not. sr3d_params%debug_on) then
  !   n = size(photon%reflect)
  !   allocate (p_temp(n))
  !   p_temp = photon%reflect
  !   deallocate (photon%reflect)
  !   allocate (photon%reflect(0:n))
  !   photon%reflect(0:n-1) = p_temp
  !   photon%reflect(n) = photon%now
  !   deallocate(p_temp)
  ! endif

  if (sr3d_params%stop_if_hit_antechamber .and. photon%hit_antechamber) return

  call sr3d_reflect_photon (photon, wall, absorbed)
  if (absorbed) return

  photon%n_reflect = photon%n_reflect + 1

enddo

end subroutine sr3d_track_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_track_photon_to_wall (photon, lat, wall)
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

subroutine sr3d_track_photon_to_wall (photon, lat, wall)

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
  if (sr3d_params%dr_track_step_max * abs(vec(6)) > &
      sr3d_params%ds_track_step_max * v_rad_max) then
    dlen = sr3d_params%ds_track_step_max / abs(vec(6))
  else
    dlen = sr3d_params%dr_track_step_max / v_rad_max
  endif

  call sr3d_propagate_photon_a_step (photon, dlen, lat, wall, .true.)

  ! See if the photon has hit the wall.
  ! If so we calculate the exact hit spot where the photon crossed the
  ! wall boundry and return

  call sr3d_photon_radius (photon%now, wall, radius)
  if (radius > 1) then
    call sr3d_photon_hit_spot_calc (photon, wall, lat)
    return
  endif

enddo

end subroutine sr3d_track_photon_to_wall

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_propagate_photon_a_step (photon, dl_step, lat, wall, stop_at_check_pt)
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

subroutine sr3d_propagate_photon_a_step (photon, dl_step, lat, wall, stop_at_check_pt)

implicit none

type (lat_struct), target :: lat
type (photon3d_track_struct), target :: photon
type (wall3d_struct) wall
type (photon3d_coord_struct), pointer :: now

real(rp) dl_step, dl_left, s_stop, denom, v_x, v_s, sin_t, cos_t
real(rp) g, new_x, radius, theta, tan_t, dl, dl2, ct, st
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
          photon%crossed_lat_end = .not. photon%crossed_lat_end
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
          photon%crossed_lat_end = .not. photon%crossed_lat_end
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
    ! shorter than the distance left to travel.

    g = lat%ele(now%ix_ele)%value(g$)
    radius = 1 / g
    theta = (s_stop - now%vec(5)) * g
    tan_t = tan(theta)
    dl = tan_t * (radius + now%vec(1)) / (now%vec(6) - tan_t * now%vec(2))

    if (abs(tan_t * (radius + now%vec(1))) > dl_left * abs(now%vec(6) - tan_t * now%vec(2))) then
      dl = dl_left
      tan_t = (dl * now%vec(6)) / (radius + now%vec(1) + dl * now%vec(2))
      theta = atan(tan_t)
      will_stop_at_s_stop = .false.
    else
      dl = tan_t * (radius + now%vec(1)) / (now%vec(6) - tan_t * now%vec(2))
      will_stop_at_s_stop = .true.
    endif

    ! Check if we should actually be stopping at the extremum (minimal x)

    if (stop_at_check_pt .and. now%vec(2) * g < 0) then 
      dl2 = -now%vec(2) * (radius + now%vec(1)) / (now%vec(2)**2 + now%vec(6)**2)
      if (dl2 < dl) then
        dl = 1.00000001 * dl2 ! Add extra to make sure we are not short due to roundoff.
        tan_t = (dl * now%vec(6)) / (radius + now%vec(1) + dl * now%vec(2))
        theta = atan(tan_t)
        s_stop = now%vec(5) + radius * theta
        will_stop_at_s_stop = .true.
        s_stop_is_check_pt = .true.
      endif
    endif

    ! Move to the stop point. 
    ! Need to remember that radius can be negative.
    st = dl * now%vec(6)
    ct = radius + now%vec(1) + dl * now%vec(2)
    if (abs(st) < 1e-3 * ct) then
      denom = sign (ct * (1 + (st/ct)**2/2 + (st/ct)**4/8), radius)
    else
      denom = sign (sqrt((radius + now%vec(1) + dl * now%vec(2))**2 + (dl * now%vec(6))**2), radius)
    endif
    sin_t = st / denom
    cos_t = ct / denom
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

end subroutine sr3d_propagate_photon_a_step

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_photon_hit_spot_calc (photon, wall, lat)
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

subroutine sr3d_photon_hit_spot_calc (photon, wall, lat)

use nr, only: zbrent

implicit none

type (lat_struct) lat
type (photon3d_track_struct) :: photon, photon1
type (wall3d_struct), target :: wall

real(rp) track_length, radius, tol, d_rad

integer i

! Find where the photon hits.
! Note: After the first reflection, the photon will start at the wall so
! we must avoid selecting this point as the next hit spot!

photon1 = photon

track_length = (photon%now%track_len + photon%old%track_len) / 2
do i = 1, 30
  d_rad = photon_hit_func(track_length)
  if (d_rad < 0) exit
  track_length = (track_length + photon%old%track_len) / 2
  if (i == 30) then
    print *, 'ERROR: CANNOT FIND HIT SPOT REGION LOWER BOUND!'
    call err_exit
  endif
enddo

tol = 1e-10 * (photon%now%track_len + photon%old%track_len)
track_length = zbrent (photon_hit_func, track_length, photon%now%track_len, tol)

! Cleanup

photon%now = photon%old
call sr3d_propagate_photon_a_step (photon, track_length-photon%now%track_len, lat, wall, .false.)
call sr3d_photon_radius (photon%now, wall, radius, hit_antechamber = photon%hit_antechamber)

!-----------------------------------------------------------------------
contains

function photon_hit_func (track_len) result (d_radius)

real(rp) track_len, d_radius, radius, d_track

! Easy case

if (track_len == photon%now%track_len) then
  call sr3d_photon_radius (photon%now, wall, radius)
  d_radius = radius - 1
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
call sr3d_propagate_photon_a_step (photon1, d_track, lat, wall, .false.)
call sr3d_photon_radius (photon1%now, wall, radius)
d_radius = radius - 1

end function

end subroutine sr3d_photon_hit_spot_calc 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_reflect_photon (photon, wall, absorbed)
!
! Routine to reflect a photon off of the wall.
!-

subroutine sr3d_reflect_photon (photon, wall, absorbed)

implicit none

type (photon3d_track_struct), target :: photon
type (wall3d_struct), target :: wall
type (wall3d_pt_struct), pointer :: wall0, wall1

real(rp) cos_perp, dw_perp(3), denom, f
real(rp) r, graze_angle, reflectivity
real(rp), pointer :: vec(:)

integer ix, iu

logical absorbed

! Check if reflections allowed

if (.not. sr3d_params%allow_reflections) then
  absorbed = .true.
  return
endif

! get the perpendicular outward normal to the wall

vec => photon%now%vec
photon%old = photon%now

call sr3d_photon_radius (photon%now, wall, r, dw_perp)

! cos_perp is the component of the photon velocity perpendicular to the wall.
! since the photon is striking the wall from the inside this must be positive.

cos_perp = dot_product (vec(2:6:2), dw_perp)
graze_angle = pi/2 - acos(cos_perp)
call photon_reflectivity (graze_angle, photon%now%energy, reflectivity)

iu = sr3d_params%iu_reflect_file 
if (iu /= 0) then
  write (iu, *) '*********************************************'
  write (iu, '(2i8, 3f10.4, 10x, 2f12.6)') photon%ix_photon, photon%n_reflect, &
                                 dw_perp, cos_perp, reflectivity
  write (iu, '(6f12.6)') photon%old%vec
endif

if (cos_perp < 0) then
  print *, 'ERROR: PHOTON AT WALL HAS VELOCITY DIRECTED INWARD!', cos_perp
  print *, '       START COORDS: ', photon%start%vec
  print *, '       ENERGY: ', photon%start%energy, photon%ix_photon
  print *, '       WILL EXIT HERE...'
  call err_exit
endif

! absorbtion

call ran_uniform(r)
absorbed = (r > reflectivity)
if (absorbed) return  ! Do not reflect if absorbed

! Reflect the ray.
! The perpendicular component gets reflected and the parallel component is invarient.

vec(2:6:2) = vec(2:6:2) - 2 * cos_perp * dw_perp

if (iu /= 0) then
  write (iu, '(3(12x, f12.6))') photon%now%vec(2:6:2)
endif

end subroutine sr3d_reflect_photon

end module
