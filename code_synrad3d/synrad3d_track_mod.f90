module synrad3d_track_mod

use synrad3d_utils
use synrad3d_output_mod

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_track_photon (photon, lat, wall_hit, err)
!
! Routine to propagate a synch radiation photon until it gets absorbed by a wall.
!
! Modules needed:
!   use synrad3d_track_mod
!
! Input:
!   photon      -- sr3d_photon_track_struct: photon with starting parameters set.
!     %start       -- Starting coords.
!   branch      -- branch_struct: Lattice branch with twiss propagated and mat6s made
!   one_reflection_only
!               -- Logical, optional: If present and True then only one reflection is allowed
!
! Output:
!   photon      -- sr3d_photon_track_struct: synch radiation photon propagated until absorption.
!   wall_hit(:) -- sr3d_photon_wall_hit_struct: Array of wall hit data.
!   err         -- Tracking calculation failed.
!-

subroutine sr3d_track_photon (photon, lat, wall_hit, err, one_reflection_only)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch, branch2
type (ele_struct), pointer :: ele0, ele1
type (sr3d_photon_track_struct), target :: photon, photon2
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)
type (coord_struct), pointer :: orb
type (sr3d_branch_overlap_struct), pointer :: bo
type (floor_position_struct) floor

real(rp) dw_perp(3), d_radius
integer iw, status, ij, ie_start, ie_end

logical absorbed, err, no_wall_here
logical, optional :: one_reflection_only

character(*), parameter :: r_name = 's43d_track_photon'

!

if (sign_of(photon%start%orb%vec(6)) /= photon%start%orb%direction) then
  call out_io (s_fatal$, r_name, '%vec(6) does not agree with %direction')
  call err_exit
endif

photon%start%orb%path_len = 0
photon%now = photon%start
photon%crossed_lat_end = .false.
wall_hit(0)%after_reflect = photon%start%orb

call ran_default_state (get_state = sr3d_params%ran_state)  ! Save 

!

branch => lat%branch(photon%now%ix_branch)

err = .false.

main_loop: do
  call sr3d_track_photon_to_wall (photon, lat, wall_hit, err)
  if (err) return

  ! Switch to another sub-chamber within same branch?

  do iw = 1, size(branch%wall3d)
    if (iw == photon%now%ix_wall3d) cycle
    status = sr3d_photon_status_calc (photon, branch, iw)

    if (status == inside_the_wall$) then
      photon%now%ix_wall3d = iw
      photon%status = status
      cycle main_loop
    endif

    if (photon%status == at_transverse_wall$) then
      call sr3d_photon_d_radius (photon%now, branch, no_wall_here, d_radius, dw_perp, ix_wall3d = iw)
      if (dot_product (photon%now%orb%vec(2:6:2), dw_perp) < 0) then   ! Traveling into sub-chamber
        photon%now%ix_wall3d = iw
        photon%status = status
        cycle main_loop
      endif
    endif  
  enddo

  ! Switch to a subchamber in a different branch?

  if (allocated(sr3d_com%branch_overlap)) then
    do ij = 1, size(sr3d_com%branch_overlap)
      bo => sr3d_com%branch_overlap(ij)
      if (branch%ix_branch == bo%ix_branch1) then
        branch2 => lat%branch(bo%ix_branch2)
        ie_start = bo%ix_ele2_start
        ie_end = bo%ix_ele2_end
      elseif (branch%ix_branch == bo%ix_branch2) then
        branch2 => lat%branch(bo%ix_branch1)
        ie_start = bo%ix_ele1_start
        ie_end = bo%ix_ele1_end
      else
        cycle
      endif

      photon2%now = photon%now
      photon2%now%ix_branch = branch2%ix_branch
      orb => photon2%now%orb
      ele0 => branch%ele(orb%ix_ele)
      floor%r = [orb%vec(1), orb%vec(3), orb%s - (ele0%s - ele0%value(l$))]
      call axis_angle_to_w_mat ([-orb%vec(4), orb%vec(2), 0.0_rp], asin(sqrt(orb%vec(2)**2 + orb%vec(4)**2)), floor%w)
      floor = coords_local_curvilinear_to_floor(floor, ele0, calculate_angles = .false.)
      floor = coords_floor_to_curvilinear (floor, branch2%ele((ie_start+ie_end)/2), ele1, status)
      if (status == outside$ .or. status == patch_problem$) cycle

      orb%vec(1:5:2) = floor%r
      orb%vec(2:6:2) = floor%w(1:3,3) * sign_of(photon%now%orb%vec(6))
      orb%ix_ele = ele1%ix_ele

      do iw = 1, size(branch2%wall3d)
        status = sr3d_photon_status_calc (photon2, branch2, iw)
        if (status /= inside_the_wall$) cycle
        photon%now = photon2%now
        photon%now%ix_wall3d = iw
        photon%status = status
        branch => branch2
        cycle main_loop
      enddo
    enddo
  endif

  ! Reflect photon

  if (sr3d_params%iu_photon_track > 0) call sr3d_record_photon_position('RECORD_TRACK_POINT', photon)
  call sr3d_reflect_photon (photon, branch, wall_hit, absorbed, err)
  if (absorbed .or. err .or. logic_option(.false., one_reflection_only)) return
enddo main_loop

end subroutine sr3d_track_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_track_photon_to_wall (photon, lat, wall_hit, err)
!
! Routine to propagate a synch radiation photon until it hits a wall.
!
! Modules needed:
!   use synrad3d_track_mod
!
! Input:
!   photon    -- sr3d_photon_track_struct: photon with starting parameters set
!   lat       -- lat_struct: Lattice.
!
! Output:
!   photon    -- sr3d_photon_track_struct: synch radiation photon propagated to wall
!   err       -- Tracking calculation failed.
!-

subroutine sr3d_track_photon_to_wall (photon, lat, wall_hit, err)

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (sr3d_photon_track_struct), target :: photon
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)

real(rp) v_rad_max, dlen, radius
real(rp), pointer :: vec(:)

integer i, status
integer, pointer :: ix_w(:)
logical err

! The photon is tracked in a series of steps.

vec => photon%now%orb%vec
branch => lat%branch(photon%now%ix_branch)

main_loop: do

  if (sr3d_params%iu_photon_track > 0) call sr3d_record_photon_position('RECORD_TRACK_POINT', photon)

  v_rad_max = max(abs(vec(2)), abs(vec(4)))
  if (sr3d_params%dr_track_step_max * abs(vec(6)) > &
      sr3d_params%ds_track_step_max * v_rad_max) then
    dlen = sr3d_params%ds_track_step_max / abs(vec(6))
  else
    dlen = sr3d_params%dr_track_step_max / v_rad_max
  endif

  call sr3d_propagate_photon_a_step (photon, branch, dlen, .true.)

  ! See if there is a fast sub-chamber to switch to

  ix_w => sr3d_com%branch(branch%ix_branch)%fast(photon%now%ix_wall3d)%ix_wall3d
  do i = 1, size(ix_w)
    status = sr3d_photon_status_calc (photon, branch, ix_w(i))
    if (status == inside_the_wall$) then
      photon%now%ix_wall3d = ix_w(i)
      photon%status = status
      cycle main_loop
    endif
  enddo

  ! See if the photon has hit the wall.
  ! If so we calculate the exact hit spot where the photon crossed the
  ! wall boundry and return

  photon%status = sr3d_photon_status_calc (photon, branch)
  if (photon%status == at_wall_end$) return
  if (photon%status == at_transverse_wall$) return
  if (photon%status == is_through_wall$) then
    call sr3d_photon_hit_spot_calc (photon, branch, wall_hit, err)
    return
  endif

enddo main_loop

end subroutine sr3d_track_photon_to_wall

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine sr3d_check_if_photon_init_coords_outside_wall (photon_start, lat, is_inside, bad_photon_counter)
!
! Routine to determine if a photon is initially created outside the vacuum chamber wall.
! Also this routine will stop the program if too many photons have been created outside.
! 
!
! Input:
!   photon_start        -- sr3d_coord_struct: Photon coords.
!   lat                 -- lat_struct: Lattice 
!   bad_photon_counter  -- integer: Counter 
!
! Output:
!   photon_start%ix_wall3d -- Set to index of sub-chamber photon is inside.
!   is_inside           -- logical: True if inside vacuum chamber wall. False otherwise.
!   bad_photon_counter  -- integer: Counter decremented by one. When counter is zero the program will be stopped.
!-

subroutine sr3d_check_if_photon_init_coords_outside_wall (photon_start, lat, is_inside, bad_photon_counter)

type (sr3d_coord_struct) photon_start
type (sr3d_photon_track_struct) photon
type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

real(rp) d_radius
integer iw, bad_photon_counter, ix, status
logical is_inside

! 

photon%now = photon_start
photon%old = photon_start
photon%old%orb%vec(1:3:2) = 0   ! 
branch => lat%branch(photon_start%ix_branch)

do iw = 1, size(branch%wall3d)
  status = sr3d_photon_status_calc (photon, branch, iw)
  if (status == inside_the_wall$ .or. status == at_wall_end$) exit
enddo

photon_start%ix_wall3d = iw
is_inside = .true.

if (status /= inside_the_wall$ .and. status /= at_wall_end$) then
  is_inside = .false.
  ix = photon_start%orb%ix_ele
  print *,                       'ERROR: INITIALIZED PHOTON IS OUTSIDE THE WALL!'
  print '(a, 4f10.4, a, f12.4)', '       PHOTON COORDS:', photon_start%orb%vec(1:4), ',  S =', photon_start%orb%s
  print '(a, i0, 2x, a)',        '       AT ELEMENT: ', ix, trim(branch%ele(ix)%name)

  ! bad_photon_counter is decremented by one. When counter is zero the program will be stopped.

  bad_photon_counter = bad_photon_counter - 1
  if (bad_photon_counter < 0) then
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
! Function sr3d_photon_status_calc (photon, branch, ix_wall3d) result (status)
!
! Routine to determine if a photon has crossed through the wall or
! is at the end of a linear lattice
!
! Input:
!   photon    -- sr3d_photon_track_struct
!   branch    -- branch_struct: Lattice branch with associated wall.
!   ix_wall3d -- integer, optional: If present then override phton%now%ix_wall3d.
!
! Output:
!   status -- Integer: is_through_wall$, at_wall_end$, at_transverse_wall$, or inside_the_wall$
!-

function sr3d_photon_status_calc (photon, branch, ix_wall3d) result (status)

implicit none

type (sr3d_photon_track_struct) photon
type (sr3d_coord_struct) now
type (branch_struct), target :: branch
type (wall3d_struct), pointer :: wall3d

real(rp) d_radius, s

integer status
integer, optional :: ix_wall3d
integer i, ix, ixs, end_geometry

logical is_through, no_wall_here

! check for particle outside wall

status = inside_the_wall$

now = photon%now
call sr3d_get_section_index(now, branch, ix_wall3d)
call sr3d_photon_d_radius (now, branch, no_wall_here, d_radius, ix_wall3d = ix_wall3d)

if (abs(d_radius) < sr3d_params%significant_length .and. .not. no_wall_here) then
  status = at_transverse_wall$
  return
endif

if (d_radius > 0 .or. no_wall_here) then
  status = is_through_wall$
  return
endif    

! Is at the end of a linear lattice or at the end of the current sub-section?

end_geometry = sr3d_params%chamber_end_geometry
if (branch%ix_branch /= 0) end_geometry = branch%param%geometry

s = now%orb%s
ixs = now%ix_wall_section
wall3d => branch%wall3d(integer_option(now%ix_wall3d, ix_wall3d))

if (now%orb%vec(6) < 0) then
  if (end_geometry == open$ .and. s == 0) then
    status = at_wall_end$
    return
  endif
  if (s == wall3d%section(ixs)%s .and. wall3d%section(ixs)%type == wall_start$) status = at_wall_end$
endif

if (now%orb%vec(6) > 0) then
  ix = branch%n_ele_track
  if (end_geometry == open$ .and. s == branch%ele(ix)%s) status = at_wall_end$
  if (s == wall3d%section(ixs+1)%s .and. wall3d%section(ixs+1)%type == wall_end$) status = at_wall_end$
endif

end function sr3d_photon_status_calc

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_propagate_photon_a_step (photon, branch, dl_step, stop_at_check_pt)
!
! Routine to propagate a photon to a given spot
!
! Modules needed:
!   use synrad3d_track_mod
!
! Input:
!   photon  -- sr3d_photon_track_struct: Photon to track.
!   branch  -- branch_struct: Lattice branch with associated wall.
!   dl_step -- Real(rp): Distance to track. Note: the propagation distance may not be exact
!               when going long distances.
!   stop_at_check_pt 
!           -- Logical: If True, stop at a check point which is defined to be:
!                a) minimum x extremum in a bend, or
!                b) At wall point boundries or
!                c) At patch boundries.
!              Note: (b) guarantees that there will be check points at the ends of the lattice.
!
! Output:
!   photon  -- sr3d_photon_track_struct: 
!			%now       -- If the photon has hit, the photon position is adjusted accordingly.
!-

subroutine sr3d_propagate_photon_a_step (photon, branch, dl_step, stop_at_check_pt)

use track1_photon_mod

implicit none

type (branch_struct), target :: branch
type (sr3d_photon_track_struct), target :: photon
type (wall3d_struct), pointer :: wall3d
type (coord_struct), pointer :: now
type (coord_struct) orb
type (ele_struct), pointer :: ele

real(rp) dl_step, dl_left, s_stop, denom, v_x, v_s, sin_t, cos_t, sl
real(rp) g, new_x, radius, theta, tan_t, dl, dl2, ct, st, s_ent, s0, ds
real(rp), pointer :: vec(:)

integer ixw, stop_location, end_geometry

logical stop_at_check_pt, check_section_here

! update old 

photon%old = photon%now  ! Save for hit spot calc
now => photon%now%orb
dl_left = dl_step

wall3d => branch%wall3d(photon%now%ix_wall3d)

ele => branch%ele(now%ix_ele)
s0 = ele%s - ele%value(l$)
sl = sr3d_params%significant_length 

end_geometry = sr3d_params%chamber_end_geometry
if (branch%ix_branch /= 0) end_geometry = branch%param%geometry

call sr3d_get_section_index(photon%now, branch, photon%now%ix_wall3d)

! Since patch elements have a complicated geometry, the photon s-position may not be
! within the interval from the s-position at the ends of the patch

if (ele%key /= patch$ .and. (now%s < min(s0, ele%s) - sl .or. now%s > max(s0, ele%s) + sl)) then
  print *, '%IX_ELE BOOKKEEPING ERROR IN SR3D_PROPAGATING_PHOTON_A_STEP!'
  call sr3d_print_photon_info (photon)
  call err_exit
endif

if (dl_step < 0) then
  print *, 'DL_STEP BOOKKEEPING ERROR IN SR3D_PROPAGATING_PHOTON_A_STEP!'
  call sr3d_print_photon_info (photon)
  call err_exit
endif

if (ele%key /= patch$ .and. now%location /= inside$) then
  if (abs(now%s - s0) > sl .and. abs(now%s - ele%s) > sl) then
    print *, 'PHOTON LOCATION BOOKKEEPING ERROR IN SR3D_PROPAGATING_PHOTON_A_STEP!'
    call err_exit
  endif
endif

! propagate the photon a number of sub-steps until we have gone a distance dl_step
! A sub-step might be to the next element boundary since tracking in a bend is
!  different from all other elements.

propagation_loop: do

  check_section_here = .false.

  if (now%direction == 1) then

    ! Move to next element if at downstream end.

    if (now%location == downstream_end$) then
      if (now%ix_ele == branch%n_ele_track) then
        if (end_geometry == open$) return
        now%ix_ele = 1
        now%s = now%s - branch%param%total_length
        photon%crossed_lat_end = .not. photon%crossed_lat_end
      else
        now%ix_ele = now%ix_ele + 1
      endif

      now%location = upstream_end$
      ele => branch%ele(now%ix_ele)
      now%vec(5) = 0

      ! Entering a patch: Transform coordinates to be with respect to the downstream end.

      if (ele%key == patch$) then
        if (ele%orientation == -1) then
          print *, 'CODE FOR PATCH WITH ORIENTATION = -1 NOT YET IMPLEMENTED!'
          call err_exit
        endif
        call track_a_patch_photon (ele, now, .false.)
        ! Due to finite pitches, it is possible that now the photon is in some element after the patch.
        if (now%s > ele%s) then 
          now%ix_ele = element_at_s(branch, now%s, .false.)
          ele => branch%ele(now%ix_ele)
          now%location = inside$
        endif
      endif
    endif

    ! Set stop position for this step

    s_stop = branch%ele(now%ix_ele)%s
    stop_location = downstream_end$

    if (stop_at_check_pt) then
      call bracket_index2 (wall3d%section%s, 1, ubound(wall3d%section, 1), now%s, photon%now%ix_wall_section, ixw)
      photon%now%ix_wall_section = ixw
    endif

    if (stop_at_check_pt .and. ixw < ubound(wall3d%section, 1)) then
      if (wall3d%section(ixw+1)%s < s_stop) then
        s_stop = wall3d%section(ixw+1)%s
        stop_location = inside$
        check_section_here = .true.
      endif
    endif

  !----------------------

  else   ! direction = -1

    ! Move to next element if at upstream end.

    if (now%location == upstream_end$) then

      if (now%ix_ele == 1) then
        if (end_geometry == open$) return
        now%s = now%s + branch%param%total_length
        now%ix_ele = branch%n_ele_track 
        photon%crossed_lat_end = .not. photon%crossed_lat_end
      else
        now%ix_ele = now%ix_ele - 1
      endif

      ele => branch%ele(now%ix_ele)
      now%vec(5) = ele%value(l$)

    endif

    ! Set stop position for this step

    stop_location = upstream_end$
    s_stop = branch%ele(now%ix_ele-1)%s

    if (stop_at_check_pt) then
      call bracket_index2 (wall3d%section%s, 1, ubound(wall3d%section, 1), now%s, photon%now%ix_wall_section, ixw)
      photon%now%ix_wall_section = ixw
    endif

    if (stop_at_check_pt) then
      do
        if (ixw == 1) exit
        if (wall3d%section(ixw)%s /= now%s) exit
        ixw = ixw - 1
      enddo

      if (wall3d%section(ixw)%s > s_stop) then
        s_stop = wall3d%section(ixw)%s
        stop_location = inside$
        check_section_here = .true.
      endif
    endif

  endif

  !------------------------------------------
  ! Propagate the photon a step.

  ! In a bend...

  ele => branch%ele(now%ix_ele)
  if (ele%key == sbend$ .and. ele%value(g$) /= 0) then

    if (ele%orientation == -1) then
      print *, 'BENDS WITH ORIENTATION == -1 NOT YET IMPLEMENTED!'
      call err_exit
    endif

    ! Rotate to element reference frame (bend in x-plane) if bend is tilted.

    if (ele%value(ref_tilt_tot$) /= 0) call tilt_coords(ele%value(ref_tilt_tot$), now%vec)

    ! Next position is determined by whether the distance to the element edge is 
    ! shorter than the distance left to travel.

    g = ele%value(g$)
    radius = 1 / g
    theta = (s_stop - now%s) * g * ele%orientation
    tan_t = tan(theta)

    if (abs(tan_t * (radius + now%vec(1))) > dl_left * abs(now%vec(6) - tan_t * now%vec(2))) then
      dl = dl_left
      tan_t = (dl * now%vec(6)) / (radius + now%vec(1) + dl * now%vec(2))
      theta = atan(tan_t)
      s_stop = now%s + radius * theta * ele%orientation
      stop_location = inside$
      check_section_here = .false.
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
        s_stop = now%s + radius * theta * ele%orientation
        stop_location = inside$
        check_section_here = .true.
      endif
    endif

    ! Move to the stop point. 
    ! Need to remember that radius can be negative.

    st = dl * now%vec(6)
    ct = radius + now%vec(1) + dl * now%vec(2)
    denom = sign (sqrt(st**2 + ct**2), radius)
    sin_t = st / denom
    cos_t = ct / denom

    v_x = now%vec(2); v_s = now%vec(6)

    now%vec(1) = denom - radius
    now%vec(2) = v_s * sin_t + v_x * cos_t
    now%vec(3) = now%vec(3) + dl * now%vec(4)
    now%s = s_stop
    now%vec(5) = now%s - (ele%s - ele%value(l$))
    now%vec(6) = v_s * cos_t - v_x * sin_t

    if (ele%value(ref_tilt_tot$) /= 0) call tilt_coords(-ele%value(ref_tilt_tot$), now%vec)

  !----
  ! Else we are not in a bend

  else

    ! Next position...
    ! In a patch going backwards: Compute distance to travel by rotating to the upstream coords.

    ds = s_stop - now%s
    if (ele%key == patch$ .and. now%direction == -1) then
      orb = now
      call track_a_patch_photon (ele, orb, .false., .true.)
      ds = s_stop - orb%s
    endif

    if (abs(now%vec(6)) * dl_left > abs(ds)) then
      dl = ds / now%vec(6)
    else
      dl = dl_left * sign_of(ele%value(l$))
      check_section_here = .false.
      s_stop = now%s + dl * now%vec(6)
      stop_location = inside$
    endif

    ! And move to the next position

    now%vec(1) = now%vec(1) + dl * now%vec(2)
    now%vec(3) = now%vec(3) + dl * now%vec(4)
    now%vec(5) = now%vec(5) + dl * now%vec(6)
    now%s = s_stop

  endif

  ! Adjust photon coords when exiting a patch with direction == -1

  if (ele%key == patch$ .and. now%direction == -1 .and. stop_location == upstream_end$) then
    if (ele%orientation == -1) then
      print *, 'CODE FOR PATCH WITH ORIENTATION = -1 NOT YET IMPLEMENTED!'
      call err_exit
    endif
    call track_a_patch_photon (ele, now, .false., .true.)
  endif

  ! Path length increases even when moving though an element with negative length.
  ! This is important when zbrent is called to find the hit spot.

  now%path_len = now%path_len + abs(dl)
  dl_left = dl_left - abs(dl)
  now%location = stop_location

  if (dl_left == 0) exit
  if (stop_at_check_pt .and. check_section_here) exit

enddo propagation_loop

end subroutine sr3d_propagate_photon_a_step

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_photon_hit_spot_calc (photon, branch, wall_hit, err)
!
! Routine to calculate where the photon has hit the wall.
!
! Modules needed:
!   use synrad3d_track_mod
!
! Input:
!   photon    -- sr3d_photon_track_struct:
!   branch  -- branch_struct: Lattice branch with associated wall.
!
! Output:
!   photon    -- sr3d_photon_track_struct: 
!			%now       -- If the photon has hit, the photon position is adjusted accordingly.
!   err       -- Tracking calculation failed.
!-

subroutine sr3d_photon_hit_spot_calc (photon, branch, wall_hit, err)

use super_recipes_mod

implicit none

type (branch_struct), target :: branch
type (sr3d_photon_track_struct), target :: photon
type (wall3d_struct), pointer :: wall3d
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)
type (sr3d_photon_track_struct) :: photon1

real(rp) r0, r1, path_len
real(rp) path_len0, path_len1, d_rad0, d_rad1

integer i

logical err, no_wall_here
logical :: in_zbrent

! For debugging

if (photon%ix_photon_generated == sr3d_params%ix_generated_warn) then
  print *
  print *, '*************************************************************'
  print *, 'Hit:', photon%n_wall_hit
  call sr3d_photon_d_radius (photon%old, branch, no_wall_here, r0)
  call sr3d_photon_d_radius (photon%now, branch, no_wall_here, r1)
  print *, 'photon%old:', photon%old%orb%vec, photon%old%orb%path_len, r0
  print *, 'photon%now:', photon%now%orb%vec, photon%now%orb%path_len, r1
endif

! Bracket the hit point. 
! Note: After the first reflection, the photon will start at the wall so
! if photon%old is at the wall we must avoid bracketing this point.

wall3d => branch%wall3d(photon%now%ix_wall3d)
photon1 = photon
path_len1 = photon%now%orb%path_len
d_rad0 = real_garbage$
d_rad1 = real_garbage$
in_zbrent = .false.

if (wall_hit(photon%n_wall_hit)%after_reflect%path_len == photon%old%orb%path_len) then

  path_len0 = (photon%now%orb%path_len + 3*photon%old%orb%path_len) / 4
  do i = 1, 30
    d_rad0 = sr3d_photon_hit_func(path_len0)
    if (photon%ix_photon_generated == sr3d_params%ix_generated_warn) then
      print *
      print *, 'path_len, d_rad0:', path_len0, d_rad0
      print *, 'photon1%now:', i, photon1%now%orb%vec, photon1%now%orb%path_len
    endif
    if (d_rad0 < 0) exit
    path_len1 = path_len0; d_rad1 = d_rad0
    path_len0 = (path_len0 + 3*photon%old%orb%path_len) / 4
    if (i == 30) then
      print *, 'ERROR: CANNOT FIND HIT SPOT REGION LOWER BOUND!'
      print '(8x, a)', 'WILL IGNORE THIS PHOTON.'
      call sr3d_print_photon_info (photon)
      call sr3d_print_hit_points (-1, photon, wall_hit, branch, .true.)
      err = .true.
      return
    endif
  enddo

else
  path_len0 = photon%old%orb%path_len
endif

! Find where the photon hits.

in_zbrent = .true.
path_len = super_zbrent (sr3d_photon_hit_func, path_len0, path_len1, sr3d_params%significant_length, err)
if (err) then
  print *, 'WILL IGNORE THIS PHOTON.'
  call sr3d_print_photon_info (photon)
  call sr3d_print_hit_points (-1, photon, wall_hit, branch, .true.)
  return
endif

! Cleanup

photon%now = photon%old
call sr3d_propagate_photon_a_step (photon, branch, path_len-photon%now%orb%path_len, .false.)
call sr3d_photon_d_radius (photon%now, branch, no_wall_here, d_rad0)

!---------------------------------------------------------------------------
contains

!+
! Function sr3d_photon_hit_func (path_len) result (d_radius)
! 
! Routine to be used as an argument in zbrent in the sr3d_photon_hit_spot_calc.
!
! Input:
!   path_len -- Real(rp): Place to position the photon.
!
! Output:
!   d_radius -- Real(rp): 
!-

function sr3d_photon_hit_func (path_len) result (d_radius)

implicit none

real(rp), intent(in) :: path_len
real(rp) d_radius, d_track

! Easy case at the ends of the track.
! The reason why we are carful about reusing d_rad0 and d_rad1 is that 
! roundoff can cause calculated radius at the end points to shift from positive 
! to negative which will case zbrent to crash.

if (in_zbrent) then
  if (path_len == path_len0 .and. d_rad0 /= real_garbage$) then
    d_radius = d_rad0
    return
  elseif (path_len == path_len1 .and. d_rad1 /= real_garbage$) then
    d_radius = d_rad1
    return
  endif
endif

! Determine start of tracking.
! If path_length > photon1%now%orb%path_len: 
!   Track starting from the present position (photon1%now).
! Otherwise:
!   Track starting from the beginning of the region (photon%old).

if (path_len < photon1%now%orb%path_len) then
  photon1 = photon
  photon1%now = photon%old
endif

! And track to path_len position.

d_track = path_len - photon1%now%orb%path_len
call sr3d_propagate_photon_a_step (photon1, branch, d_track, .false.)

call sr3d_photon_d_radius (photon1%now, branch, no_wall_here, d_radius)

end function sr3d_photon_hit_func

end subroutine sr3d_photon_hit_spot_calc 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_reflect_photon (photon, branch, wall_hit, absorbed, err_flag)
!
! Routine to reflect a photon off of the chamber wall.
!
! Additionally: this routine will calculate if the photon is to be absorbed or reflected.
! The absorption calculation involves calculating the reflection probability and then,
! using a random number generator, deciding if the photon is indeed absorbed.
!
! Input:
!   photon   -- sr3d_photon_track_struct: Photon position.
!   branch   -- branch_struct: Lattice branch with associated wall.
!
! Output:
!   wall_hit(:) -- sr3d_photon_wall_hit_struct: Array recording where the photon has hit the wall.
!   absorbed    -- Logical: Set True if photon is absorbed.
!   err_flag    -- Logical: Set True if an error found. Not touched otherwise.
!-

subroutine sr3d_reflect_photon (photon, branch, wall_hit, absorbed, err_flag)

implicit none

type (sr3d_photon_track_struct), target :: photon
type (wall3d_struct), pointer :: wall3d
type (branch_struct), target :: branch
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)
type (sr3d_photon_wall_hit_struct), allocatable :: hit_temp(:)
type (photon_reflect_surface_struct), pointer :: surface

real(rp) cos_perp, dw_perp(3), denom, f, r, d_rad, theta_diffuse, phi_diffuse
real(rp) graze_angle, reflectivity, rel_reflect_specular, dvec(3), a, b
real(rp) vec_in_plane(3), vec_out_plane(3)

integer ix, iu
integer n_old, n_wall_hit

logical absorbed, err_flag, no_wall_here

character(*), parameter :: r_name = 'sr3d_reflect_photon'

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
wall_hit(n_wall_hit)%before_reflect = photon%now%orb
wall_hit(n_wall_hit)%dw_perp = 0
wall_hit(n_wall_hit)%cos_perp_in = 0
wall_hit(n_wall_hit)%cos_perp_out = 0
wall_hit(n_wall_hit)%reflectivity = 0
wall_hit(n_wall_hit)%after_reflect%vec = 0
wall_hit(n_wall_hit)%ix_wall3d = photon%now%ix_wall3d
wall_hit(n_wall_hit)%ix_branch = photon%now%ix_branch

absorbed = .true.

! Check if reflections allowed

if (.not. sr3d_params%allow_reflections) return

photon%old = photon%now

! get the perpendicular outward normal to the wall

if (photon%status == at_wall_end$) then

  if (photon%now%orb%vec(6) > 0) then
    dw_perp = [0, 0, 1]
  else
    dw_perp = [0, 0, -1]
  endif

else
  call sr3d_photon_d_radius (photon%now, branch, no_wall_here, d_rad, dw_perp)
endif

! cos_perp is the component of the photon velocity perpendicular to the wall.
! since the photon is striking the wall from the inside this must be positive.

cos_perp = dot_product (photon%now%orb%vec(2:6:2), dw_perp)
if (cos_perp > 1) then
  if (cos_perp < 1.0d0 + 1d-10) then ! Allow for roundoff
    cos_perp = 1
  else
    call out_io (s_error$, r_name, 'ERROR IN REFLECTION CALCULATATION!')
    err_flag = .true.
    return
  endif
endif

graze_angle = pi/2 - acos(cos_perp)
dvec = -2 * cos_perp * dw_perp

call sr3d_get_section_index(photon%now, branch, photon%now%ix_wall3d)
if (photon%now%ix_wall_section == not_set$) call sr3d_get_section_index (photon%now, branch)
surface => branch%wall3d(photon%now%ix_wall3d)%section(photon%now%ix_wall_section+1)%surface

! Get reflectivity coef.

if (surface%name == 'ABSORBER') then
  reflectivity = 0
else
  call photon_reflectivity (graze_angle, photon%now%orb%p0c, surface, reflectivity, rel_reflect_specular)
endif
wall_hit(n_wall_hit)%reflectivity = reflectivity

if (cos_perp < 0) then
  call out_io (s_error$, r_name, &
  'ERROR: PHOTON AT WALL HAS VELOCITY DIRECTED INWARD! \f10.5\ ', & 
  'WILL IGNORE THIS PHOTON...', &
  'dw_perp: \3f10.5\ ', & 
  r_array = [cos_perp, dw_perp])
  call sr3d_print_photon_info (photon)
  call sr3d_print_hit_points (-1, photon, wall_hit, branch, .true.)
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
    photon%now%orb%vec(2:6:2) = photon%now%orb%vec(2:6:2) + dvec

  else
    call photon_diffuse_scattering (graze_angle, photon%now%orb%p0c, surface, theta_diffuse, phi_diffuse)
    ! vec_in_plane is normalized vector perpendicular to dw_perp and in plane of photon & dw_perp.
    vec_in_plane = photon%now%orb%vec(2:6:2) - dw_perp * cos_perp  
    vec_in_plane = vec_in_plane / sqrt(dot_product(vec_in_plane, vec_in_plane))  ! Normalize to 1.
    vec_out_plane = cross_product(dw_perp, vec_in_plane)
    photon%now%orb%vec(2:6:2) = -cos(theta_diffuse) * dw_perp + sin(theta_diffuse) * &
                            (vec_in_plane * cos(phi_diffuse) + vec_out_plane * sin(phi_diffuse))
  endif
endif

if (photon%now%orb%vec(6) < 0) then
  photon%now%orb%direction = -1
else
  photon%now%orb%direction = 1
endif

! Record

wall_hit(n_wall_hit)%dw_perp = dw_perp
wall_hit(n_wall_hit)%cos_perp_in = cos_perp
wall_hit(n_wall_hit)%after_reflect = photon%now%orb
wall_hit(n_wall_hit)%cos_perp_out = dot_product (photon%now%orb%vec(2:6:2), dw_perp)

end subroutine sr3d_reflect_photon

end module
