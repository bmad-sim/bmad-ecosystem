!+
! Subroutine type_coord (coord)
!
! Subroutine to type out a coordinate.
!
! Modules Needed:
!   use bmad
!
! Input:
!   coord -- Coord_struct: Coordinate
!-

#include "CESR_platform.inc"

subroutine type_coord (coord)

  use bmad_struct
  use bmad_interface, except => type_coord

  implicit none
  type (coord_struct)  coord

!

  print *, '(X, X''):', coord%vec(1), coord%vec(2)
  print *, '(Y, Y''):', coord%vec(3), coord%vec(4)
  print *, '(Z, Z''):', coord%vec(5), coord%vec(6)

end subroutine

