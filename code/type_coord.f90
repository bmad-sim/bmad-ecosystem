!+
! Subroutine TYPE_COORD (COORD)
!
! Subroutine to type out a coordinate.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     COORD -- Coord_struct: Coordinate
!-


subroutine type_coord (coord)

  use bmad_struct
  implicit none
  type (coord_struct)  coord

!

  type *, '(X, X''):', coord%x%pos, coord%x%vel
  type *, '(Y, Y''):', coord%y%pos, coord%y%vel
  type *, '(Z, Z''):', coord%z%pos, coord%z%vel

  return
  end

