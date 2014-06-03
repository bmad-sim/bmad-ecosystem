!+                       
! Subroutine track_all (lat, orbit, ix_branch, track_state, err_flag)
!
! Subroutine to track through the lat.
!
! Note: If x_limit (or y_limit) for an element is zero then track_all will take
!       x_limit (or y_limit) as infinite.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat         -- lat_struct: Lat to track through.
!   orbit(0)    -- Coord_struct: Coordinates at beginning of lat.
!   ix_branch   -- Integer, optional: Index of branch to track. Default is 0 (main lattice).
!
! Output:
!   orbit(0:*)   -- Coord_struct: Orbit array.
!   track_state  -- Integer, optional: Set to moving_forward$ if everything is OK.
!                     Otherwise: set to index of element where particle was lost.
!   err_flag     -- Logical, optional: Set true if particle lost or error. False otherwise
!-

subroutine track_all (lat, orbit, ix_branch, track_state, err_flag)

use bmad_interface, except_dummy => track_all
use bookkeeper_mod, only: control_bookkeeper

implicit none

type (lat_struct), target :: lat
type (coord_struct), allocatable :: orbit(:)
type (ele_struct), pointer :: lord, slave
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

integer n, i, nn, ix_br
integer, optional :: ix_branch, track_state

logical, optional :: err_flag
logical err
logical :: debug = .false.

character(12), parameter :: r_name = 'track_all'

! init

if (present(err_flag)) err_flag = .true.
ix_br = integer_option (0, ix_branch)
branch => lat%branch(ix_br)
if (present(track_state)) track_state = moving_forward$

if (.not. allocated(orbit)) call reallocate_coord (orbit, branch%n_ele_max)
if (size(orbit) < branch%n_ele_max) call reallocate_coord (orbit, branch%n_ele_max)

if (orbit(0)%state == not_set$) call init_coord(orbit(0), orbit(0)%vec, branch%ele(0), .true., branch%param%particle) 
if (bmad_com%auto_bookkeeper) call control_bookkeeper (lat)

orbit(0)%ix_ele   = 0
orbit(0)%location = downstream_end$

if (orbit(0)%species /= photon$) then
  call convert_pc_to (branch%ele(0)%value(p0c$) * (1 + orbit(0)%vec(6)), orbit(0)%species, beta = orbit(0)%beta)
  orbit(0)%p0c = branch%ele(0)%value(p0c$)
endif

! track through elements.

if (present(err_flag)) err_flag = .false.

do n = 1, branch%n_ele_track

  ele => branch%ele(n)
  call track1 (orbit(n-1), ele, branch%param, orbit(n), err_flag = err)

  ! check for lost particles

  if (err .or. .not. particle_is_moving_forward(orbit(n))) then
    if (present(track_state)) track_state = n
    do nn = n+1, branch%n_ele_track
      orbit(nn)%vec = 0
      orbit(nn)%state = not_set$
    enddo
    if (present(err_flag)) err_flag = .true.
    exit
  endif

  if (debug) then
    call out_io (s_blank$, r_name, ele%name, '\6es16.6\ ', r_array = orbit(n)%vec)
  endif

enddo

! Fill in super_lord info.

do n = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(n) 
  if (lord%lord_status /= super_lord$) cycle
  slave => pointer_to_slave(lord, lord%n_slave)
  if (slave%ix_branch /= ix_br) cycle
  orbit(lord%ix_ele) = orbit(slave%ix_ele)
enddo

end subroutine
