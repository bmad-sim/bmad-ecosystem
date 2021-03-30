!+
! Subroutine angle_to_canonical_coords (orbit)
!
! Routine to convert from angle (x, x' = dx/ds, y, y' = dy/ds, z, pz) coordinates to
! canonical (x, px, y, py, z, pz) coordinates.
!
! Input:
!   orbit -- coord_struct: Orbit in angular coordinates.
!
! Output:
!   orbit -- coord_struct: Orbit in canonical coordinates.
!-

subroutine angle_to_canonical_coords (orbit)

use bmad_struct
implicit none

type (coord_struct) orbit
real(rp) f

!

f = (1 + orbit%vec(6)) / sqrt(1 + orbit%vec(2)**2 + orbit%vec(4)**2)
orbit%vec(2) = orbit%vec(2) * f
orbit%vec(4) = orbit%vec(4) * f

end subroutine angle_to_canonical_coords
