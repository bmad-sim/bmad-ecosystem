module synrad3d_track_mod

use synrad3d_utils
use synrad3d_output_mod

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_track_photon (photon, lat, wall, wall_hit, err)
!
! Routine to propagate a synch radiation photon until it gets absorbed by a wall.
!
! Modules needed:
!   use synrad3d_track_mod
!
! Input:
!   photon    -- sr3d_photon_coord_struct: photon with starting parameters set.
!     %start    -- Starting coords.
!   lat       -- lat_struct: with twiss propagated and mat6s made
!   wall      -- sr3d_wall_struct: Beam chamber walls
!   one_reflection_only
!               -- Logical, optional: If present and True then only one reflection is allowed
!
! Output:
!   photon      -- sr3d_photon_coord_struct: synch radiation photon propagated until absorption.
!   wall_hit(:) -- sr3d_photon_wall_hit_struct: Array of wall hit data.
!   err         -- Tracking calculation failed.
!-

subroutine sr3d_track_photon (photon, lat, wall, wall_hit, err, one_reflection_only)

implicit none

type (lat_struct), target :: lat
type (sr3d_photon_track_struct), target :: photon
type (sr3d_wall_struct), target :: wall
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)

logical absorbed, err
logical, optional :: one_reflection_only

!

photon%start%track_len = 0
photon%now = photon%start
wall_hit(0)%after_reflect = photon%start

call ran_default_state (get_state = sr3d_params%ran_state)  ! Save 

!

err = .false.

do
  call sr3d_track_photon_to_wall (photon, lat, wall, wall_hit, err)
  if (err) return
  call sr3d_reflect_photon (photon, wall, lat, wall_hit, absorbed, err)
  if (absorbed .or. err .or. logic_option(.false., one_reflection_only)) return
enddo

end subroutine sr3d_track_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_track_photon_to_wall (photon, lat, wall, wall_hit, err)
!
! Routine to propagate a synch radiation photon until it hits a wall.
!
! Modules needed:
!   use synrad3d_track_mod
!
! Input:
!   photon    -- sr3d_photon_coord_struct: photon with starting parameters set
!   lat       -- lat_struct: with twiss propagated and mat6s made
!   wall      -- sr3d_wall_struct: Beam chamber walls.
!
! Output:
!   photon    -- sr3d_photon_coord_struct: synch radiation photon propagated to wall
!   err       -- Tracking calculation failed.
!-

subroutine sr3d_track_photon_to_wall (photon, lat, wall, wall_hit, err)

implicit none

type (lat_struct), target :: lat
type (sr3d_photon_track_struct), target :: photon
type (sr3d_wall_struct), target :: wall
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)

real(rp) v_rad_max, dlen, radius
real(rp), pointer :: vec(:)

logical err

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

  call sr3d_photon_status_calc (photon, wall)
  if (photon%status == at_lat_end$) return
  if (photon%status == is_through_wall$) then
    call sr3d_photon_hit_spot_calc (photon, wall, lat, wall_hit, err)
    return
  endif

enddo

end subroutine sr3d_track_photon_to_wall

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
!-
subroutine sr3d_check_if_photon_init_coords_outside_wall (photon_start, wall, is_inside, num_ignore_generated_outside_wall)

type (sr3d_photon_coord_struct) photon_start
type (sr3d_photon_track_struct) photon
type (sr3d_wall_struct), target :: wall

real(rp) d_radius
integer num_ignore_generated_outside_wall
logical is_inside

! For triangular "outside wall" is defined here to be true if a ray from the 
! center to photon_start crosses a wall boundary.

photon%now = photon_start
photon%old = photon_start
photon%old%vec(1:3:2) = 0   ! 
call sr3d_photon_status_calc (photon, wall)

is_inside = .true.

if (photon%status /= inside_the_wall$ .and. photon%status /= at_lat_end$) then
  is_inside = .false.
  print *,              'ERROR: INITIALIZED PHOTON IS OUTSIDE THE WALL!', photon%ix_photon_generated
  print '(a, 6f10.4)', '        INITIALIZATION PT: ', photon_start%vec      

  num_ignore_generated_outside_wall = num_ignore_generated_outside_wall - 1
  if (num_ignore_generated_outside_wall < 0) then
    print '(a)', '       STOPPING SYNRAD3D DUE TO NUMBER OF PHOTONS GENERATED OUTSIDE'
    print '(a)', '       THE WALL EXCEEDING NUM_IGNORE_GENERATED_OUTSIDE_WALL VALUE!'
    stop
  endif

endif

end subroutine sr3d_check_if_photon_init_coords_outside_wall 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_photon_status_calc (photon, wall) 
!
! Routine to determine if a photon has crossed through the wall or
! is at the end of a linear lattice
!
! Input:
!   photon  -- sr3d_photon_track_struct
!   wall    -- sr3d_wall_struct: Wall
!
! Output:
!   photon%status -- Integer: is_through_wall$, at_lat_end$, or inside_the_wall$
!-

subroutine sr3d_photon_status_calc (photon, wall) 

implicit none

type (sr3d_photon_track_struct) photon
type (sr3d_wall_struct) wall

real(rp) d_radius
real(rp) tri_vert0(3), tri_vert1(3), tri_vert2(3)

integer i, ix

logical is_through, checked

! check for particle outside wall

photon%status = inside_the_wall$
checked = .false.

if (wall%has_triangular_sections) then
  call sr3d_get_wall_index (photon%now, wall, ix)
  if (wall%section(ix+1)%basic_shape == 'triangular') then
    do i = 1, 2*size(wall%section(ix)%gen_shape%wall3d_section%v)
      call sr3d_get_mesh_wall_triangle_pts (wall%section(ix), wall%section(ix+1), i, tri_vert0, tri_vert1, tri_vert2)
      call sr3d_mesh_triangle_intersect (photon, tri_vert0, tri_vert1, tri_vert2, is_through)
      if (is_through) then
        photon%status = is_through_wall$
        return 
      endif
    enddo
    checked = .true.
  endif
endif

if (.not. checked) then
  call sr3d_photon_d_radius (photon%now, wall, d_radius, check_safe = .true.)
  if (d_radius > 0) then
    photon%status = is_through_wall$
    return 
  endif    
endif

! Is through if at ends of a linear lattice

if (wall%geometry == open$) then
  if (photon%now%s == 0 .and. photon%now%vec(6) < 0) photon%status = at_lat_end$
  if (photon%now%s == wall%section(wall%n_section_max)%s .and. photon%now%vec(6) > 0) photon%status = at_lat_end$
endif

end subroutine sr3d_photon_status_calc

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
!   photon    -- sr3d_photon_track_struct:
!     %old        -- Original point
!     %now        -- current point                
!   pt0(3), pt1(3), pt2(3) 
!             -- Real(rp): Triangle vertex points.
!
! Output:
!   intersect  -- Logical: True if there is an intersection. between photon%old and photon%now.
!   dtrack_len -- Real(rp), optional: track length from photon%old to intersection point.
!                  Negative if no intersection.
!-

subroutine sr3d_mesh_triangle_intersect (photon, pt0, pt1, pt2, intersect, dtrack_len)

implicit none

type (sr3d_photon_track_struct) photon

real(rp) pt0(3), pt1(3), pt2(3), eps
real(rp), optional :: dtrack_len
real(rp) ratio
real(rp) mat3(3,3), abc_vec(3)

integer ip
logical intersect, ok, line_intersect

!

if (present(dtrack_len)) dtrack_len = -1

! Solve the matrix equation:
!   photon_old + a * (photon_now - photon_old) = pt0 + b * (pt1 - pt0) + c * (pt2 - pt0)

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
!   0 <= a <= 1
!   0 <= b
!   0 <= c
!   b + c <= 1
! Roundoff errors can be a problem here. Better to count something as intersecting
! that is not then miss something that is. So use eps.

ip = precision(1.0_rp)  ! precision
eps = 10.0_rp**(1-ip)   ! Something small used for preventing round off errors

line_intersect = (abc_vec(1) >= 0 .and. abc_vec(2) >= 0 .and. &
                  abc_vec(3) >= 0 .and. abc_vec(2) + abc_vec(3) <= 1+eps)

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
!   photon  -- sr3d_photon_coord_struct: Photon to track.
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
!   photon  -- sr3d_photon_coord_struct: 
!			%now       -- If the photon has hit, the photon position is adjusted accordingly.
!-

subroutine sr3d_propagate_photon_a_step (photon, dl_step, lat, wall, stop_at_check_pt)

implicit none

type (lat_struct), target :: lat
type (sr3d_photon_track_struct), target :: photon
type (sr3d_wall_struct) wall
type (sr3d_photon_coord_struct), pointer :: now
type (ele_struct), pointer :: ele

real(rp) dl_step, dl_left, s_stop, denom, v_x, v_s, sin_t, cos_t
real(rp) g, new_x, radius, theta, tan_t, dl, dl2, ct, st
real(rp), pointer :: vec(:)

integer ixw

logical stop_at_check_pt, check_section_here

! update old 

photon%old = photon%now  ! Save for hit spot calc
now => photon%now
dl_left = dl_step

! propagate the photon a number of sub-steps until we have gone a distance dl_step
! A sub-step might be to the next element boundary since tracking in a bend is
!  different from all other elements.

propagation_loop: do

  check_section_here = .false.
  if (stop_at_check_pt) then
    call bracket_index2 (wall%section%s, 0, wall%n_section_max, now%s, now%ix_wall, ixw)
    now%ix_wall = ixw
  endif

  ! If we are crossing over to a new element then update now%ix_ele.

  if (now%vec(6) > 0) then
    do
      if (now%s >= lat%ele(now%ix_ele)%s) then
        if (now%ix_ele == lat%n_ele_track) then
          if (lat%param%geometry == open$) return
          now%s = now%s - lat%param%total_length
          now%ix_ele = 1
          photon%crossed_lat_end = .not. photon%crossed_lat_end
          exit
        endif
        now%ix_ele = now%ix_ele + 1
      elseif (now%s < lat%ele(now%ix_ele-1)%s) then
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

    if (stop_at_check_pt .and. ixw < wall%n_section_max) then
      if (wall%section(ixw+1)%s < s_stop) then
        s_stop = wall%section(ixw+1)%s
        check_section_here = .true.
      endif
    endif

  else   ! direction = -1
    do
      if (now%s <= lat%ele(now%ix_ele-1)%s) then
        if (now%ix_ele <= 1) then
          if (lat%param%geometry == open$) return
          now%s = now%s + lat%param%total_length
          now%ix_ele = lat%n_ele_track
          photon%crossed_lat_end = .not. photon%crossed_lat_end
          exit
        endif
        now%ix_ele = now%ix_ele - 1
      elseif (now%s > lat%ele(now%ix_ele)%s) then
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
      if (wall%section(ixw)%s == now%s) ixw = ixw - 1
      if (wall%section(ixw)%s > s_stop) then
        s_stop = wall%section(ixw)%s
        check_section_here = .true.
      endif
    endif

  endif

  ! Propagate the photon a step.

  ! In a bend...

  ele => lat%ele(now%ix_ele)
  if (ele%key == sbend$ .and. ele%value(g$) /= 0) then

    ! Rotate to element reference frame (bend in x-plane) if bend is tilted.

    if (ele%value(ref_tilt_tot$) /= 0) call tilt_coords(ele%value(ref_tilt_tot$), now%vec)

    ! Next position is determined by whether the distance to the element edge is 
    ! shorter than the distance left to travel.

    g = ele%value(g$)
    radius = 1 / g
    theta = (s_stop - now%s) * g
    tan_t = tan(theta)

    if (abs(tan_t * (radius + now%vec(1))) > dl_left * abs(now%vec(6) - tan_t * now%vec(2))) then
      dl = dl_left
      tan_t = (dl * now%vec(6)) / (radius + now%vec(1) + dl * now%vec(2))
      theta = atan(tan_t)
      check_section_here = .false.
      s_stop = now%s + radius * theta
    else
      dl = tan_t * (radius + now%vec(1)) / (now%vec(6) - tan_t * now%vec(2))
    endif

    ! Check if we should actually be stopping at the extremum (minimal x)

    if (stop_at_check_pt .and. now%vec(2) * g < 0) then 
      dl2 = -now%vec(2) * (radius + now%vec(1)) / (now%vec(2)**2 + now%vec(6)**2)
      if (dl2 < dl) then
        dl = dl2 * (1 + sr3d_params%significant_length) ! Add extra to make sure we are not short due to roundoff.
        tan_t = (dl * now%vec(6)) / (radius + now%vec(1) + dl * now%vec(2))
        theta = atan(tan_t)
        s_stop = now%s + radius * theta
        check_section_here = .true.
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
    now%s = s_stop
    now%vec(6) = v_s * cos_t - v_x * sin_t

    if (ele%value(ref_tilt_tot$) /= 0) call tilt_coords(-ele%value(ref_tilt_tot$), now%vec)

  ! Else we are not in a bend

  else

    ! Next position

    if (abs(now%vec(6)) * dl_left > abs(s_stop - now%s)) then
      dl = (s_stop - now%s) / now%vec(6)
    else
      dl = dl_left
      check_section_here = .false.
      s_stop = now%s + dl * now%vec(6)
    endif

    ! And move to the next position

    now%vec(1) = now%vec(1) + dl * now%vec(2)
    now%vec(3) = now%vec(3) + dl * now%vec(4)
    now%s = s_stop
  endif

  !

  photon%now%track_len = photon%now%track_len + dl
  dl_left = dl_left - dl

  if (dl_left == 0) exit
  if (stop_at_check_pt .and. check_section_here) exit

enddo propagation_loop

end subroutine sr3d_propagate_photon_a_step

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_photon_hit_spot_calc (photon, wall, lat, wall_hit, err)
!
! Routine to calculate where the photon has hit the wall.
!
! Modules needed:
!   use synrad3d_track_mod
!
! Input:
!   photon    -- sr3d_photon_coord_struct:
!   wall      -- sr3d_wall_struct: 
!   lat       -- lat_struct: Lattice
!
! Output:
!   photon    -- sr3d_photon_coord_struct: 
!			%now       -- If the photon has hit, the photon position is adjusted accordingly.
!   err       -- Tracking calculation failed.
!-

subroutine sr3d_photon_hit_spot_calc (photon, wall, lat, wall_hit, err)

use super_recipes_mod

implicit none

type (lat_struct), target :: lat
type (sr3d_photon_track_struct), target :: photon
type (sr3d_wall_struct), target :: wall
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)
type (sr3d_photon_track_struct) :: photon1

real(rp) r0, r1, track_len
real(rp) track_len0, track_len1, d_rad0, d_rad1

integer i

logical err
logical :: in_zbrent

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

photon1 = photon
photon%now%ix_triangle = -1
track_len1 = photon%now%track_len
d_rad0 = real_garbage$
d_rad1 = real_garbage$
in_zbrent = .false.

if (wall_hit(photon%n_wall_hit)%after_reflect%track_len == photon%old%track_len) then

  track_len0 = (photon%now%track_len + 3*photon%old%track_len) / 4
  do i = 1, 30
    d_rad0 = sr3d_photon_hit_func(track_len0)
    if (photon%ix_photon_generated == sr3d_params%ix_generated_warn) then
      print *
      print *, 'track_len, d_rad0:', track_len0, d_rad0
      print *, 'photon1%now:', i, photon1%now%vec, photon1%now%track_len
    endif
    if (d_rad0 < 0) exit
    track_len1 = track_len0; d_rad1 = d_rad0
    track_len0 = (track_len0 + 3*photon%old%track_len) / 4
    if (i == 30) then
      print *, 'ERROR: CANNOT FIND HIT SPOT REGION LOWER BOUND!'
      print *, '       Photon:', photon%ix_photon, photon%ix_photon_generated, photon%n_wall_hit, photon%start%energy
      print *, '       Start: ', photon%start%vec
      print *, '       Now:   ', photon%now%vec
      print *, '       WILL IGNORE THIS PHOTON.'
      call print_hit_points (-1, photon, wall_hit, .true.)
      err = .true.
      return
    endif
  enddo

else
  track_len0 = photon%old%track_len
endif

! Find where the photon hits.

in_zbrent = .true.
track_len = super_zbrent (sr3d_photon_hit_func, track_len0, track_len1, sr3d_params%significant_length, err)
if (err) then
  call print_hit_points (-1, photon, wall_hit, .true.)
  print *, 'WILL IGNORE THIS PHOTON.'
  print *, '       Photon:', photon%ix_photon, photon%ix_photon_generated, photon%n_wall_hit, photon%start%energy
  print *, '       Start: ', photon%start%vec
  print *, '       Now:   ', photon%now%vec
  return
endif

! Cleanup

photon%now = photon%old
call sr3d_propagate_photon_a_step (photon, track_len-photon%now%track_len, lat, wall, .false.)
call sr3d_photon_d_radius (photon%now, wall, d_rad0, in_antechamber = photon%hit_antechamber)
photon%now%ix_triangle = photon1%now%ix_triangle

!---------------------------------------------------------------------------
contains

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

! Easy case at the ends of the track.
! The reason why we are carful about reusing d_rad0 and d_rad1 is that 
! roundoff can cause calculated radius at the end points to shift from positive 
! to negative which will case zbrent to crash.

if (in_zbrent) then
  if (track_len == track_len0 .and. d_rad0 /= real_garbage$) then
    d_radius = d_rad0
    return
  elseif (track_len == track_len1 .and. d_rad1 /= real_garbage$) then
    d_radius = d_rad1
    return
  endif
endif

! At the beginning of the track the mesh calc has problems with zero length steps.
! So just interpolate from the end

if (track_len == photon%old%track_len .and. &
      wall%section(photon%now%ix_wall+1)%basic_shape == 'triangular') then
  call sr3d_mesh_d_radius (photon, wall, d_radius)
  d_radius = d_radius - (photon%now%track_len - photon%old%track_len)  
  return
endif

! Determine start of tracking.
! If track_length > photon1%now%track_len: 
!   Track starting from the present position (photon1%now).
! Otherwise:
!   Track starting from the beginning of the region (photon%old).

if (track_len < photon1%now%track_len) then
  photon1 = photon
  photon1%now = photon%old
endif

! And track to track_len position.

d_track = track_len - photon1%now%track_len
call sr3d_propagate_photon_a_step (photon1, d_track, lat, wall, .false.)

if (wall%section(photon%now%ix_wall+1)%basic_shape == 'triangular') then
  call sr3d_mesh_d_radius (photon1, wall, d_radius)
else
  call sr3d_photon_d_radius (photon1%now, wall, d_radius)
endif

end function sr3d_photon_hit_func

end subroutine sr3d_photon_hit_spot_calc 

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

type (sr3d_photon_track_struct), target :: photon
type (sr3d_wall_struct), target :: wall

real(rp) d_radius, dlen, dlen_min
real(rp) tv0(3), tv1(3), tv2(3)

integer i, ix

logical is_through

character(20) :: r_name = 'sr3d_mesh_d_radius'

! 

ix = photon%now%ix_wall

! To save time try old triangle

if (photon%now%ix_triangle > 0) then
  call sr3d_get_mesh_wall_triangle_pts (wall%section(ix), wall%section(ix+1), photon%now%ix_triangle, tv0, tv1, tv2)
  call sr3d_mesh_triangle_intersect (photon, tv0, tv1, tv2, is_through, dlen)
  if (dlen > 0) then
    d_radius = (photon%now%track_len - photon%old%track_len) - dlen 
    return
  endif
endif

! Now must look at all triangles and pick first tirangle through

photon%now%ix_triangle = -1
dlen_min = 1d30   ! Something large

do i = 1, 2*size(wall%section(ix)%gen_shape%wall3d_section%v)
  call sr3d_get_mesh_wall_triangle_pts (wall%section(ix), wall%section(ix+1), i, tv0, tv1, tv2)
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
! Subroutine sr3d_reflect_photon (photon, wall, lat, wall_hit, absorbed, err_flag)
!
! Routine to reflect a photon off of the chamber wall.
!
! Additionally: this routine will calculate if the photon is to be absorbed or reflected.
! The absorption calculation involves calculating the reflection probability and then,
! using a random number generator, deciding if the photon is indeed absorbed.
!
! Input:
!   photon   -- sr3d_photon_track_struct: Photon position.
!   wall     -- sr3d_wall_struct: Chamber wall.
!   lat      -- lat_struct: Lattice
!
! Output:
!   wall_hit(:) -- sr3d_photon_wall_hit_struct: Array recording where the photon has hit the wall.
!   absorbed    -- Logical: Set True if photon is absorbed.
!   err_flag    -- Logical: Set True if an error found. Not touched otherwise.
!-

subroutine sr3d_reflect_photon (photon, wall, lat, wall_hit, absorbed, err_flag)

implicit none

type (sr3d_photon_track_struct), target :: photon
type (sr3d_wall_struct), target :: wall
type (lat_struct) lat
type (sr3d_wall_section_struct), pointer :: wall0, wall1
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)
type (sr3d_photon_wall_hit_struct), allocatable :: hit_temp(:)
type (photon_reflect_surface_struct), pointer :: surface

real(rp) cos_perp, dw_perp(3), denom, f, r, d_rad, theta_diffuse, phi_diffuse
real(rp) graze_angle, reflectivity, rel_reflect_specular, dvec(3)
real(rp) vec_in_plane(3), vec_out_plane(3)

integer ix, iu
integer n_old, n_wall_hit

logical absorbed, err_flag

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
wall_hit(n_wall_hit)%cos_perp_in = 0
wall_hit(n_wall_hit)%cos_perp_out = 0
wall_hit(n_wall_hit)%reflectivity = 0
wall_hit(n_wall_hit)%after_reflect%vec = 0

absorbed = .true.

! Check if reflections allowed or hit antechamber

if (.not. sr3d_params%allow_reflections .or. photon%status == at_lat_end$ .or. &
    (sr3d_params%stop_if_hit_antechamber .and. photon%hit_antechamber)) return

! get the perpendicular outward normal to the wall

photon%old = photon%now

call sr3d_photon_d_radius (photon%now, wall, d_rad, lat, dw_perp)

! cos_perp is the component of the photon velocity perpendicular to the wall.
! since the photon is striking the wall from the inside this must be positive.

cos_perp = dot_product (photon%now%vec(2:6:2), dw_perp)
graze_angle = pi/2 - acos(cos_perp)
dvec = -2 * cos_perp * dw_perp

surface => wall%section(photon%now%ix_wall+1)%gen_shape%wall3d_section%surface

call photon_reflectivity (graze_angle, photon%now%energy, surface, reflectivity, rel_reflect_specular)
wall_hit(n_wall_hit)%reflectivity = reflectivity

if (cos_perp < 0) then
  print *, 'ERROR: PHOTON AT WALL HAS VELOCITY DIRECTED INWARD!', cos_perp
  print *, '       dw_perp:', dw_perp
  print *, '       Photon: ', photon%ix_photon, photon%ix_photon_generated, photon%n_wall_hit, photon%start%energy
  print *, '       Start:  ', photon%start%vec
  print *, '       Now:    ', photon%now%vec
  print *, '       WILL IGNORE THIS PHOTON...'
  call print_hit_points (-1, photon, wall_hit, .true.)
  err_flag = .true.
  return
endif

! absorption or reflection...
! For specular reflection the perpendicular component gets reflected and the parallel component is invarient.

call ran_uniform(r)
if (.not. sr3d_params%allow_absorption) reflectivity = 1

if (r <= reflectivity) then
  absorbed = .false.

  if (sr3d_params%specular_reflection_only .or. r < reflectivity * rel_reflect_specular) then
    photon%now%vec(2:6:2) = photon%now%vec(2:6:2) + dvec

  else
    call photon_diffuse_scattering (graze_angle, photon%now%energy, surface, theta_diffuse, phi_diffuse)
    ! vec_in_plane is normalized vector perpendicular to dw_perp and in plane of photon & dw_perp.
    vec_in_plane = photon%now%vec(2:6:2) - dw_perp * cos_perp  
    vec_in_plane = vec_in_plane / sqrt(dot_product(vec_in_plane, vec_in_plane))  ! Normalize to 1.
    vec_out_plane = cross_product(dw_perp, vec_in_plane)
    photon%now%vec(2:6:2) = -cos(theta_diffuse) * dw_perp + sin(theta_diffuse) * &
                            (vec_in_plane * cos(phi_diffuse) + vec_out_plane * sin(phi_diffuse))
  endif
endif

! Record

wall_hit(n_wall_hit)%dw_perp = dw_perp
wall_hit(n_wall_hit)%cos_perp_in = cos_perp
wall_hit(n_wall_hit)%after_reflect = photon%now
wall_hit(n_wall_hit)%cos_perp_out = dot_product (photon%now%vec(2:6:2), dw_perp)

end subroutine sr3d_reflect_photon

end module
