!+                       
! Subroutine track_all (ring, orbit_)
!
! Subroutine to track through the ring.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring      -- Ring_struct: Ring to track through.
!     %param%aperture_limit_on -- Logical: Sets whether TRACK_ALL looks to
!                                   see whether a particle is lost or not
!   orbit_(0) -- Coord_struct: Coordinates at beginning of ring.
!
! Output:
!   ring
!     %param%lost    -- Logical: Set when a particle is lost with the 
!                         aperture limit on.
!     %param%ix_lost -- Integer: set to index of element where particle is lost
!   orbit_(0:*) -- Coord_struct: Orbit array
!
! Note: If x_limit (or y_limit) for an element is zero then TRACK_ALL will take
!       x_limit (or y_limit) as infinite.
!-


subroutine track_all (ring, orbit_)

  use bmad_struct
  implicit none

  type (ring_struct)  ring
  type (coord_struct)  orbit_(0:*)

  integer n, i

  logical debug / .false. /

! track through elements.

  ring%param%lost = .false.

  do n = 1, ring%n_ele_use

    if (ring%ele_(n)%key == coil$) then
      call coil_track (orbit_(n-1), n, ring, orbit_(n))
    else
      call track1 (orbit_(n-1), ring%ele_(n), ring%param, orbit_(n))
    endif

! check for lost particles

    if (ring%param%lost) then
      ring%param%ix_lost = n
      return
    endif

    if (debug) then
      type *, ring%ele_(n)%name
      type *, (orbit_(n)%vec(i), i = 1, 6)
    endif

  enddo

end subroutine
