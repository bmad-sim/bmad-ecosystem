module synrad3d_track_mod

use synrad3d_utils
use photon_reflection_mod

! These common variables are needed for the Num Rec routine zbrent.

private sr3d_photon_hit_func
type (photon3d_track_struct), pointer, private :: photon_com
type (photon3d_track_struct), private :: com, photon1_com
type (wall3d_struct), private, pointer :: wall_com
type (lat_struct), private, pointer :: lat_com

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_track_photon (photon, lat, wall, wall_hit)
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
!   photon      -- photon3d_coord_struct: synch radiation photon propagated until absorbtion.
!   wall_hit(:) -- photon3d_wall_hit_struct: Array of wall hit data.
!-

subroutine sr3d_track_photon (photon, lat, wall, wall_hit)

implicit none

type (lat_struct), target :: lat
type (photon3d_track_struct), target :: photon
type (wall3d_struct), target :: wall
type (photon3d_wall_hit_struct), allocatable :: wall_hit(:)

logical absorbed

!

photon%start%track_len = 0
photon%now = photon%start
wall_hit(0)%after_reflect = photon%start

!

do
  call sr3d_track_photon_to_wall (photon, lat, wall, wall_hit)
  call sr3d_reflect_photon (photon, wall, wall_hit, absorbed)
  if (absorbed .and. sr3d_params%allow_absorbtion) return
enddo

end subroutine sr3d_track_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_track_photon_to_wall (photon, lat, wall, wall_hit)
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

subroutine sr3d_track_photon_to_wall (photon, lat, wall, wall_hit)

implicit none

type (lat_struct), target :: lat
type (photon3d_track_struct), target :: photon
type (wall3d_struct), target :: wall
type (photon3d_wall_hit_struct), allocatable :: wall_hit(:)

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

  if (sr3d_photon_through_wall(photon, wall)) then
    call sr3d_photon_hit_spot_calc (photon, wall, lat, wall_hit)
    return
  endif

enddo

end subroutine sr3d_track_photon_to_wall

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function sr3d_photon_through_wall (photon, wall) result (is_through)
!
! Routine to determine if a photon has crossed through the wall
!
! Input:
!   photon  -- photon3d_track_struct
!   wall    -- wall3d_struct: Wall
!
! Output:
!   is_through -- Logical: True if through.
!-

function sr3d_photon_through_wall (photon, wall) result (is_through)

implicit none

type (photon3d_track_struct) photon
type (wall3d_struct) wall

real(rp) d_radius
real(rp) tri_vert0(3), tri_vert1(3), tri_vert2(3)

integer i, ix

logical is_through

! The present position is between wall%pt(ix)%s and wall%pt(ix+1)%s

call sr3d_get_wall_index (photon%now, wall, ix)

if (wall%pt(ix+1)%basic_shape == 'gen_shape_mesh') then
  do i = 1, 2*size(wall%pt(ix)%gen_shape%v)
    call sr3d_get_mesh_wall_triangle_pts (wall%pt(ix), wall%pt(ix+1), i, tri_vert0, tri_vert1, tri_vert2)
    call sr3d_mesh_triangle_intersect (photon, tri_vert0, tri_vert1, tri_vert2, is_through)
    if (is_through) return 
  enddo
else
  call sr3d_photon_d_radius (photon%now, wall, d_radius)
  is_through = (d_radius > 0)
endif

end function sr3d_photon_through_wall

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_mesh_triangle_intersect (photon, pt0, pt1, pt2, intersect, dtrack_len)
!
! Routine to find the intersection point between the photon trajectory
! and a triangular mesh surface between the points photon%old and photon%now. 
!
!
! Note: This calculation assumes a straight reference coordinates.
!
! Input:
!   photon    -- photon3d_track_struct:
!     %old        -- Original point
!     %now        -- current point                
!   pt0(3), pt1(3), pt2(3) 
!             -- Real(rp): Triangle vertex points.
!
! Output:
!   intersect  -- Logical: True if there is an intersection. between photon%old and photon%now.
!   dlen       -- Real(rp), optional: track length from photon%old to intersection point.
!                  Negative if no intersection.
!-

subroutine sr3d_mesh_triangle_intersect (photon, pt0, pt1, pt2, intersect, dtrack_len)

implicit none

type (photon3d_track_struct) photon

real(rp) pt0(3), pt1(3), pt2(3)
real(rp), optional :: dtrack_len
real(rp) ratio
real(rp) mat3(3,3), abc_vec(3)

logical intersect, ok, line_intersect

!

if (present(dtrack_len)) dtrack_len = -1

! Solve the matrix equation:
!   photon_old + a * (photon_now - photon_old) = pt0 + b * (pt1 - pt0) + c * (pt2 = pt0)

mat3(1:3,1) = photon%now%vec(1:5:2) - photon%old%vec(1:5:2)
mat3(1:3,2) = pt0 - pt1
mat3(1:3,3) = pt0 - pt2

call mat_inverse (mat3, mat3, ok)
if (.not. ok) then   ! EG: photon trajectory is parallel to triangle surface.
  intersect = .false.
  return
endif

abc_vec = matmul(mat3, pt0  - photon%old%vec(1:5:2))

! Intersection of photon trajectory between "old" and "now" positions with the triangle if:
!   0 < a < 1
!   0 < b
!   0 < c
!   b + c < 1

line_intersect = (abc_vec(1) >= 0 .and. abc_vec(2) >= 0 .and. &
                  abc_vec(3) >= 0 .and. abc_vec(2) + abc_vec(3) <= 1)

intersect = (line_intersect .and. abc_vec(1) <= 1)

if (present(dtrack_len) .and. line_intersect) then
  dtrack_len = abc_vec(1) * (photon%now%track_len - photon%old%track_len)
endif

end subroutine sr3d_mesh_triangle_intersect

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

! propagate the photon a number of sub-steps until we have gone a distance dl_step
! A sub-step might be to the next element boundary since tracking in a bend is
!  different from all other elements.

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
! Subroutine sr3d_photon_hit_spot_calc (photon, wall, lat, wall_hit)
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

subroutine sr3d_photon_hit_spot_calc (photon, wall, lat, wall_hit)

use nr, only: zbrent

implicit none

type (lat_struct), target :: lat
type (photon3d_track_struct), target :: photon
type (wall3d_struct), target :: wall
type (photon3d_wall_hit_struct), allocatable :: wall_hit(:)

real(rp) track_len0, radius, d_rad, r0, r1, track_len

integer i

! For debugging

if (photon%ix_photon_generated == sr3d_params%ix_generated_warn) then
  print *
  print *, '*************************************************************'
  print *, 'Hit:', photon%n_wall_hit
  call sr3d_photon_d_radius (photon%old, wall, r0)
  call sr3d_photon_d_radius (photon%now, wall, r1)
  print *, 'photon%old:', photon%old%vec, photon%old%track_len, r0
  print *, 'photon%now:', photon%now%vec, photon%now%track_len, r1
endif

! Bracket the hit point. 
! Note: After the first reflection, the photon will start at the wall so
! if photon%old is at the wall we must avoid bracketing this point.

photon1_com = photon
photon_com => photon
wall_com => wall
lat_com => lat
photon_com%now%ix_triangle = -1

if (wall_hit(photon%n_wall_hit)%after_reflect%track_len == photon%old%track_len) then

  track_len0 = (photon%now%track_len + photon%old%track_len) / 2
  do i = 1, 30
    d_rad = sr3d_photon_hit_func(track_len0)
    if (photon%ix_photon_generated == sr3d_params%ix_generated_warn) then
      print *
      print *, 'track_len, d_rad:', track_len0, d_rad
      print *, 'photon1_com%now:', i, photon1_com%now%vec, photon1_com%now%track_len
    endif
    if (d_rad < 0) exit
    track_len0 = (track_len0 + photon%old%track_len) / 2
    if (i == 30) then
      print *, 'ERROR: CANNOT FIND HIT SPOT REGION LOWER BOUND!'
      print *, '       Photon:', photon%ix_photon, photon%ix_photon_generated, photon%n_wall_hit, photon%start%energy
      print *, '       Start: ', photon%start%vec
      call print_hit_points (10, photon, wall_hit, '(6es25.15)')
      call err_exit
    endif
  enddo

else
  track_len0 = photon%old%track_len
endif

! Find where the photon hits.

track_len = zbrent (sr3d_photon_hit_func, track_len0, photon%now%track_len, 1d-10)

! Cleanup

photon%now = photon%old
call sr3d_propagate_photon_a_step (photon, track_len-photon%now%track_len, lat, wall, .false.)
call sr3d_photon_d_radius (photon%now, wall, in_antechamber = photon%hit_antechamber)
photon_com%now%ix_triangle = photon1_com%now%ix_triangle

end subroutine sr3d_photon_hit_spot_calc 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function sr3d_photon_hit_func (track_len) result (d_radius)
! 
! Routine to be used as an argument in zbrent in the sr3d_photon_hit_spot_calc.
!
! Input:
!   track_len -- Real(rp): Place to position the photon.
!
! Output:
!   d_radius -- Real(rp): 
!-

function sr3d_photon_hit_func (track_len) result (d_radius)

implicit none

real(rp), intent(in) :: track_len
real(rp) d_radius, radius, d_track

! Easy case at the end of the track

if (track_len == photon_com%now%track_len) then
  if (wall_com%pt(photon_com%now%ix_wall+1)%basic_shape == 'gen_shape_mesh') then
    call sr3d_mesh_d_radius (photon_com, wall_com, d_radius)
  else
    call sr3d_photon_d_radius (photon_com%now, wall_com, d_radius)
  endif
  return
endif

! At the beginning of the track the mesh calc has porblems with zero length steps.
! So just interpolate from the end

if (track_len == photon_com%old%track_len .and. &
      wall_com%pt(photon_com%now%ix_wall+1)%basic_shape == 'gen_shape_mesh') then
  call sr3d_mesh_d_radius (photon_com, wall_com, d_radius)
  d_radius = d_radius - (photon_com%now%track_len - photon_com%old%track_len)  
  return
endif

! Determine start of tracking.
! If track_length > photon1_com%now%track_len: 
!   Track starting from the present position (photon1_com%now).
! Otherwise:
!   Track starting from the beginning of the region (photon_com%old).

if (track_len < photon1_com%now%track_len) then
  photon1_com = photon_com
  photon1_com%now = photon_com%old
endif

! And track to track_len position.

d_track = track_len - photon1_com%now%track_len
call sr3d_propagate_photon_a_step (photon1_com, d_track, lat_com, wall_com, .false.)

if (wall_com%pt(photon_com%now%ix_wall+1)%basic_shape == 'gen_shape_mesh') then
  call sr3d_mesh_d_radius (photon1_com, wall_com, d_radius)
else
  call sr3d_photon_d_radius (photon1_com%now, wall_com, d_radius)
endif

end function sr3d_photon_hit_func

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_mesh_d_radius (photon, wall, d_radius)
!
! Routine to calculate an effective d_radius for a photon
!-

subroutine sr3d_mesh_d_radius (photon, wall, d_radius)

implicit none

type (photon3d_track_struct), target :: photon
type (wall3d_struct), target :: wall

real(rp) d_radius, dlen, dlen_min
real(rp) tv0(3), tv1(3), tv2(3)

integer i, ix

logical is_through

character(20) :: r_name = 'sr3d_mesh_d_radius'

! 

ix = photon%now%ix_wall

! To save time try old triangle

if (photon%now%ix_triangle > 0) then
  call sr3d_get_mesh_wall_triangle_pts (wall%pt(ix), wall%pt(ix+1), photon%now%ix_triangle, tv0, tv1, tv2)
  call sr3d_mesh_triangle_intersect (photon, tv0, tv1, tv2, is_through, dlen)
  if (dlen > 0) then
    d_radius = (photon%now%track_len - photon%old%track_len) - dlen 
    return
  endif
endif

! Now must look at all triangles and pick first tirangle through

photon%now%ix_triangle = -1
dlen_min = 1d30   ! Something large

do i = 1, 2*size(wall%pt(ix)%gen_shape%v)
  call sr3d_get_mesh_wall_triangle_pts (wall%pt(ix), wall%pt(ix+1), i, tv0, tv1, tv2)
  call sr3d_mesh_triangle_intersect (photon, tv0, tv1, tv2, is_through, dlen)
  if (dlen <= 0) cycle
  if (dlen < dlen_min) then
    photon%now%ix_triangle = i
    dlen_min = dlen
  endif
enddo

if (photon%now%ix_triangle == -1) then
  call out_io (s_fatal$, r_name, 'NO INTERSECTION WITH WALL FOUND!')
  call err_exit
endif

d_radius = (photon%now%track_len - photon%old%track_len) - dlen_min

end subroutine sr3d_mesh_d_radius

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_reflect_photon (photon, wall, wall_hit, absorbed)
!
! Routine to reflect a photon off of the wall.
!-

subroutine sr3d_reflect_photon (photon, wall, wall_hit, absorbed)

implicit none

type (photon3d_track_struct), target :: photon
type (wall3d_struct), target :: wall
type (wall3d_pt_struct), pointer :: wall0, wall1
type (photon3d_wall_hit_struct), allocatable :: wall_hit(:)
type (photon3d_wall_hit_struct), allocatable :: hit_temp(:)

real(rp) cos_perp, dw_perp(3), denom, f, r
real(rp) graze_angle, reflectivity, dvec(3)

integer ix, iu
integer n_old, n_wall_hit

logical absorbed

!

n_old = ubound(wall_hit, 1)
n_wall_hit = photon%n_wall_hit + 1
if (n_old < n_wall_hit) then
  allocate (hit_temp(0:n_old))
  hit_temp = wall_hit
  deallocate (wall_hit)
  allocate (wall_hit(0:2*n_wall_hit))
  wall_hit(0:n_old) = hit_temp
  deallocate(hit_temp)
endif

photon%n_wall_hit = n_wall_hit
wall_hit(n_wall_hit)%before_reflect = photon%now
wall_hit(n_wall_hit)%dw_perp = 0
wall_hit(n_wall_hit)%cos_perp = 0
wall_hit(n_wall_hit)%reflectivity = -1

! Check if reflections allowed or hit antechamber

if (.not. sr3d_params%allow_reflections .or. &
    (sr3d_params%stop_if_hit_antechamber .and. photon%hit_antechamber)) then
  wall_hit(n_wall_hit)%dw_perp = 0
  wall_hit(n_wall_hit)%cos_perp = 0
  wall_hit(n_wall_hit)%reflectivity = 0
  wall_hit(n_wall_hit)%after_reflect%vec = 0
  absorbed = .true.
  return
endif

! get the perpendicular outward normal to the wall

photon%old = photon%now

call sr3d_photon_d_radius (photon%now, wall, dw_perp = dw_perp)

! cos_perp is the component of the photon velocity perpendicular to the wall.
! since the photon is striking the wall from the inside this must be positive.

cos_perp = dot_product (photon%now%vec(2:6:2), dw_perp)
graze_angle = pi/2 - acos(cos_perp)
dvec = -2 * cos_perp * dw_perp

call photon_reflectivity (graze_angle, photon%now%energy, reflectivity)

if (cos_perp < 0) then
  print *, 'ERROR: PHOTON AT WALL HAS VELOCITY DIRECTED INWARD!', cos_perp
  print *, '       START COORDS: ', photon%start%vec
  print *, '       ENERGY: ', photon%start%energy, photon%ix_photon
  print *, '       WILL EXIT HERE...'
  call err_exit
endif

! Record

wall_hit(n_wall_hit)%dw_perp = dw_perp
wall_hit(n_wall_hit)%cos_perp = cos_perp
wall_hit(n_wall_hit)%reflectivity = reflectivity
wall_hit(n_wall_hit)%after_reflect = photon%now
wall_hit(n_wall_hit)%after_reflect%vec(2:6:2) = photon%now%vec(2:6:2) + dvec

! absorbtion

call ran_uniform(r)
absorbed = (r > reflectivity)
if (absorbed .and. sr3d_params%allow_absorbtion) return  ! Do not reflect if absorbed

! Reflect the ray.
! The perpendicular component gets reflected and the parallel component is invarient.

photon%now%vec(2:6:2) = photon%now%vec(2:6:2) + dvec

end subroutine sr3d_reflect_photon

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------

subroutine print_hit_points (iu_hit_file, photon, wall_hit, fmt)

implicit none

type (photon3d_track_struct), target :: photon
type (photon3d_wall_hit_struct), pointer :: hit
type (photon3d_wall_hit_struct), target :: wall_hit(:)

integer iu, n, iu_hit_file

character(20) fm
character(*), optional :: fmt
!


fm = '(6f12.6)'
if (present(fmt)) fm = fmt

iu = iu_hit_file 
if (iu == 0) return

write (iu, *) '*********************************************'
write (iu, '(2i8, f10.1)') photon%ix_photon, 0, photon%start%energy
write (iu, fm) photon%start%vec

do n = 1, photon%n_wall_hit
  hit => wall_hit(n)
  write (iu, *) '*********************************************'
  write (iu, '(2i8, f10.1)') photon%ix_photon, n, hit%before_reflect%energy
  write (iu, fm) hit%before_reflect%vec
  write (iu, '(3(12x, f12.6))') hit%after_reflect%vec(2:6:2)
  write (iu, '(3f10.4, 10x, 2f12.6)') hit%dw_perp, hit%cos_perp, hit%reflectivity
enddo

end subroutine

end module
