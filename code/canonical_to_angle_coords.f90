!+
! Subroutine canonical_to_angle_coords (orbit)
!
! Routine to convert from canonical (x, px, y, py, z, pz) coordinates to
! angle (x, x' = dx/ds, y, y' = dy/ds, z, pz) coordinates.
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
real(rp) f

!

f = 1 / sqrt((1 + orbit%vec(6))**2 - orbit%vec(2)**2 - orbit%vec(4)**2)
orbit%vec(2) = orbit%vec(2) * f
orbit%vec(4) = orbit%vec(4) * f

end subroutine canonical_to_angle_coords

