!+
! Subroutine TYPE_COORD (COORD)
!
! Subroutine to type out a coordinate.
!
! Modules Needed:
!   use bmad
!
! Input:
!     COORD -- Coord_struct: Coordinate
!-

!$Id$
!$Log$
!Revision 1.5  2003/05/02 15:44:04  dcs
!F90 standard conforming changes.
!
!Revision 1.4  2003/01/27 14:40:47  dcs
!bmad_version = 56
!
!Revision 1.3  2002/02/23 20:32:30  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:32:01  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine type_coord (coord)

  use bmad_struct
  use bmad_interface

  implicit none
  type (coord_struct)  coord

!

  print *, '(X, X''):', coord%x%pos, coord%x%vel
  print *, '(Y, Y''):', coord%y%pos, coord%y%vel
  print *, '(Z, Z''):', coord%z%pos, coord%z%vel

  return
  end

