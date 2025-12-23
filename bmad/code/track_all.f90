!+                       
! Subroutine track_all (lat, orbit, ix_branch, track_state, err_flag, orbit0, init_lost)
!
! Routine to track through a lattice branch. Will start at the beginning unless there
! is an active fixer in which case tracking will be forwards and backwards from the fixer.
!
! Input:
!   lat         -- lat_struct: Lat to track through.
!   orbit(0:)   -- coord_struct, allocatable: orbit(0) is the starting coordinates for tracking.
!                   If not allocated, the zero orbit will be used.
!   ix_branch   -- integer, optional: Index of branch to track. Default is 0 (main branch).
!   init_lost   -- logical, option: Default if False. If True, initialize orbit(N) terms that
!                    are not tracked through due to particle loss.
!
! Output:
!   orbit(0:)    -- coord_struct, allocatable: Orbit array.
!   track_state  -- integer, optional: Set to moving_forward$ if everything is OK.
!                     Otherwise: set to index of element where particle was lost.
!   err_flag     -- logical, optional: Set true if particle lost or error. False otherwise
!   orbit0(0:)   -- coord_struct, allocatable, optional: Orbit array for branch 0. Used to fill 
!                     in the orbit at lord elemenets. 
!                     Only needed when orbit(:) is not the orbit for branch 0.
!-

subroutine track_all (lat, orbit, ix_branch, track_state, err_flag, orbit0, init_lost)

use bmad_interface, except_dummy => track_all

implicit none

type (lat_struct), target :: lat
type (coord_struct), allocatable, target :: orbit(:)
type (coord_struct), optional, allocatable, target :: orbit0(:)
type (coord_struct), pointer :: orb(:)
type (ele_struct), pointer :: lord, slave
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

integer n, i, nn, ix_br, ix_fix
integer, optional :: ix_branch, track_state

logical, optional :: err_flag, init_lost
logical err
logical :: debug = .false.

character(12), parameter :: r_name = 'track_all'

! init

if (present(err_flag)) err_flag = .true.
ix_br = integer_option (0, ix_branch)
branch => lat%branch(ix_br)
ix_fix = branch%ix_fixer
if (branch%param%geometry == closed$) ix_fix = 0
if (present(track_state)) track_state = moving_forward$

if (.not. allocated(orbit)) call reallocate_coord (orbit, branch%n_ele_max)
if (size(orbit) < branch%n_ele_max) call reallocate_coord (orbit, branch%n_ele_max)

if (orbit(ix_fix)%state == not_set$) call init_coord(orbit(ix_fix), orbit(ix_fix)%vec, branch%ele(ix_fix), downstream_end$)

orbit(ix_fix)%ix_ele   = ix_fix
orbit(ix_fix)%ix_branch = ix_br
orbit(ix_fix)%location = downstream_end$

if (orbit(ix_fix)%species /= photon$) then
  call convert_pc_to (branch%ele(ix_fix)%value(p0c$) * (1 + orbit(ix_fix)%vec(6)), orbit(ix_fix)%species, beta = orbit(ix_fix)%beta)
  orbit(ix_fix)%p0c = branch%ele(ix_fix)%value(p0c$)
endif

! Track forward through elements.

if (present(err_flag)) err_flag = .false.

do n = ix_fix+1, branch%n_ele_track
  ele => branch%ele(n)
  call track1 (orbit(n-1), ele, branch%param, orbit(n), err_flag = err)

  ! check for lost particles

  if (err .or. .not. particle_is_moving_forward(orbit(n))) then
    if (present(track_state)) track_state = n
    call set_orbit_to_zero (orbit, n+1, branch%n_ele_track)
    if (orbit(n)%location == upstream_end$) orbit(n)%vec = 0 ! But do not reset orbit(n)%state
    if (present(err_flag)) err_flag = .true.
    if (logic_option(.false., init_lost)) orbit(n+1:branch%n_ele_track) = coord_struct()
    exit
  endif

  if (debug) then
    call out_io (s_blank$, r_name, ele%name, '\6es16.6\ ', r_array = orbit(n)%vec)
  endif
enddo

! Track reverse through elements

if (ix_fix > 0) then
  orbit(ix_fix-1) = orbit(ix_fix)
  orbit(ix_fix-1)%time_dir = -1
  orbit(ix_fix-1)%ix_ele = ix_fix-1

  do n = ix_fix-1, 1, -1
    ele => branch%ele(n)
    call track1 (orbit(n), ele, branch%param, orbit(n-1), err_flag = err)

    ! check for lost particles

    if (err .or. .not. particle_is_moving_forward(orbit(n-1), -1)) then
      if (present(track_state)) track_state = n-1
      call set_orbit_to_zero (orbit, 0, n-1)
      if (orbit(n-1)%location == downstream_end$) orbit(n-1)%vec = 0 ! But do not reset orbit(n-1)%state
      if (present(err_flag)) err_flag = .true.
      if (logic_option(.false., init_lost)) orbit(0:n-2) = coord_struct()
      exit
    endif

    if (debug) then
      call out_io (s_blank$, r_name, ele%name, '\6es16.6\ ', r_array = orbit(n-1)%vec)
    endif
  enddo
endif

! Fill in orbits for lord elements.

if (ix_br == 0) then
  orb => orbit
elseif (present(orbit0)) then
  orb => orbit0
else
  return   ! Nothing can be done
endif

do n = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(n) 
  if (lord%lord_status /= super_lord$) cycle
  slave => pointer_to_slave(lord, lord%n_slave)
  if (slave%ix_branch /= ix_br) cycle
  orb(lord%ix_ele) = orbit(slave%ix_ele)
enddo

end subroutine
