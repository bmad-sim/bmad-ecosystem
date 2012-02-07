!+                       
! Subroutine track_all (lat, orbit, ix_branch)
!
! Subroutine to track through the lat.
!
! Note: If x_limit (or y_limit) for an element is zero then track_all will take
!       x_limit (or y_limit) as infinite.
!
! Note: If a particle does not make it through an lcavity because of lack of
!       sufficient energy, then orbit(ix_lost)%vec(6) will be < -1. 
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat       -- lat_struct: Lat to track through.
!     %param%aperture_limit_on -- Logical: Sets whether track_all looks to
!                                 see whether a particle hits an aperture or not.
!   orbit(0)  -- Coord_struct: Coordinates at beginning of lat.
!   ix_branch -- Integer, optional: Branch to track. Default is 0 (main lattice).
!
! Output:
!   lat%branch(ix_branch)%param -- Structure holding the info if the particle is lost.
!       %lost          -- Logical: Set True when a particle cannot make it 
!                           through an element.
!       %ix_lost       -- Integer: Set to index of element where particle is lost.
!       %end_lost_at   -- entrance_end$ or exit_end$.
!       %plane_lost_at -- x_plane$, y_plane$ (for apertures), or 
!                           z_plane$ (turned around in an lcavity).
!   orbit(0:*)  -- Coord_struct: Orbit array.
!-

subroutine track_all (lat, orbit, ix_branch)

use bmad_struct
use bmad_interface, except_dummy => track_all
use bookkeeper_mod, only: control_bookkeeper

implicit none

type (lat_struct), target :: lat
type (coord_struct), allocatable :: orbit(:)
type (ele_struct), pointer :: lord, slave
type (branch_struct), pointer :: branch

integer n, i, nn, ix_br
integer, optional :: ix_branch

logical :: debug = .false.

character(12), parameter :: r_name = 'track_all'

! init

ix_br = integer_option (0, ix_branch)
branch => lat%branch(ix_br)
branch%param%ix_lost = not_lost$

if (size(orbit) < branch%n_ele_max) call reallocate_coord (orbit, branch%n_ele_max)

if (bmad_com%auto_bookkeeper) call control_bookkeeper (lat)

if (branch%ele(0)%tracking_method /= time_runge_kutta$ .or. orbit(0)%status == not_set$) then
  orbit(0)%status = outside$
  call convert_pc_to (branch%ele(0)%value(p0c$) * (1 + orbit(0)%vec(6)), branch%param%particle, beta = orbit(0)%beta)
  orbit(0)%p0c = branch%ele(0)%value(p0c$)
endif

! track through elements.

do n = 1, branch%n_ele_track

  call track1 (orbit(n-1), branch%ele(n), branch%param, orbit(n))

  ! check for lost particles

  if (branch%param%lost) then
    branch%param%ix_lost = n
    do nn = n+1, branch%n_ele_track
      orbit(nn)%vec = 0
    enddo
    exit
  endif

  if (debug) then
    call out_io (s_blank$, r_name, branch%ele(n)%name, '\6es16.6\ ', r_array = orbit(n)%vec)
  endif

enddo

! Fill in super_lord info.
! This only works on the main branch at present

if (ix_br /= 0) return

do n = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(n) 
  if (lord%lord_status /= super_lord$) cycle
  slave => pointer_to_slave(lord, lord%n_slave)
  if (slave%ix_branch /= ix_br) cycle
  orbit(lord%ix_ele) = orbit(slave%ix_ele)
enddo

end subroutine
