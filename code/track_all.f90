!+                       
! Subroutine track_all (ring, orbit_)
!
! Subroutine to track through the ring.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ring      -- Ring_struct: Ring to track through.
!     %param%aperture_limit_on -- Logical: Sets whether track_all looks to
!                                 see whether a particle hits an aperture or not.
!   orbit_(0) -- Coord_struct: Coordinates at beginning of ring.
!
! Output:
!   ring
!     %param%lost    -- Logical: Set True when a particle cannot make it 
!                         through an element.
!     %param%ix_lost -- Integer: Set to index of element where particle is lost.
!   orbit_(0:*) -- Coord_struct: Orbit array
!
! Note: If x_limit (or y_limit) for an element is zero then TRACK_ALL will take
!       x_limit (or y_limit) as infinite.
!-

#include "CESR_platform.inc"

subroutine track_all (ring, orbit_)

  use bmad_struct
  use bmad_interface
  use bookkeeper_mod, only: control_lord_bookkeeper

  implicit none

  type (ring_struct)  ring
  type (coord_struct), allocatable :: orbit_(:)

  integer n, i, nn

  logical debug / .false. /

! init

  if (size(orbit_) < ring%n_ele_max+1) &
                  call reallocate_coord (orbit_, ring%n_ele_max)

  ring%param%ix_lost = -1

  call control_lord_bookkeeper (ring)

! track through elements.

  do n = 1, ring%n_ele_use

    call track1 (orbit_(n-1), ring%ele_(n), ring%param, orbit_(n))

! check for lost particles

    if (ring%param%lost) then
      ring%param%ix_lost = n
      do nn = n+1, ring%n_ele_use
        orbit_(nn)%vec = 0
      enddo
      return
    endif

    if (debug) then
      print *, ring%ele_(n)%name
      print *, (orbit_(n)%vec(i), i = 1, 6)
    endif

  enddo

end subroutine
