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
!Revision 1.6  2003/06/04 17:55:57  dcs
!Eliminated x%pos, x%vel, etc. from coord_struct.
!
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

  print *, '(X, X''):', coord%vec(1), coord%vec(2)
  print *, '(Y, Y''):', coord%vec(3), coord%vec(4)
  print *, '(Z, Z''):', coord%vec(5), coord%vec(6)

  return
  end

