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

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:32:01  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



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

