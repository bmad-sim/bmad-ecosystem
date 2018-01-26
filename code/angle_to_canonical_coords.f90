!+
! Subroutine angle_to_canonical_coords (orbit)
!
! Routine to convert from angle (x, x', y, y', z, z') coordinates to
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

!

orbit%vec(2) = orbit%vec(2) * (1 + orbit%vec(6))
orbit%vec(4) = orbit%vec(4) * (1 + orbit%vec(6))

end subroutine angle_to_canonical_coords
