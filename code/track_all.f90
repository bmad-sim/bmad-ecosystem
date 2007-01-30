!+                       
! Subroutine track_all (lat, orbit)
!
! Subroutine to track through the lat.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat      -- lat_struct: Lat to track through.
!     %param%aperture_limit_on -- Logical: Sets whether track_all looks to
!                                 see whether a particle hits an aperture or not.
!   orbit(0)  -- Coord_struct: Coordinates at beginning of lat.
!
! Output:
!   lat
!     %param%lost    -- Logical: Set True when a particle cannot make it 
!                         through an element.
!     %param%ix_lost -- Integer: Set to index of element where particle is lost.
!   orbit(0:*)  -- Coord_struct: Orbit array
!
! Note: If x_limit (or y_limit) for an element is zero then TRACK_ALL will take
!       x_limit (or y_limit) as infinite.
!-

#include "CESR_platform.inc"

subroutine track_all (lat, orbit)

  use bmad_struct
  use bmad_interface, except => track_all
  use bookkeeper_mod, only: control_bookkeeper

  implicit none

  type (lat_struct)  lat
  type (coord_struct), allocatable :: orbit(:)

  integer n, i, nn

  logical :: debug = .false.

! init

  if (size(orbit) < lat%n_ele_max+1) &
                  call reallocate_coord (orbit, lat%n_ele_max)

  lat%param%ix_lost = -1

  if (bmad_com%auto_bookkeeper) call control_bookkeeper (lat)

! track through elements.

  do n = 1, lat%n_ele_track

    call track1 (orbit(n-1), lat%ele(n), lat%param, orbit(n))

! check for lost particles

    if (lat%param%lost) then
      lat%param%ix_lost = n
      do nn = n+1, lat%n_ele_track
        orbit(nn)%vec = 0
      enddo
      return
    endif

    if (debug) then
      print *, lat%ele(n)%name
      print *, (orbit(n)%vec(i), i = 1, 6)
    endif

  enddo

end subroutine
