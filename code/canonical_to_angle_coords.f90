!+
! Subroutine canonical_to_angle_coords (orbit)
!
! Routine to convert from canonical (x, px, y, py, z, pz) coordinates to
! angle (x, x', y, y', z, z') coordinates.
!
! Input:
!   orbit -- coord_struct: Orbit in canonical coordinates.
!
! Output:
!   orbit -- coord_struct: Orbit in angular coordinates.
!-

subroutine canonical_to_angle_coords (orbit)

use bmad_struct
implicit none

type (coord_struct) orbit

!

orbit%vec(2) = orbit%vec(2) / (1 + orbit%vec(6))
orbit%vec(4) = orbit%vec(4) / (1 + orbit%vec(6))

end subroutine canonical_to_angle_coords

