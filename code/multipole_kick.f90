!+
! Subroutine multipole_kick (knl, tilt, n, coord)
!
! Subroutine to put in the kick due to a multipole.
!
! Modules Needed:
!   use bmad
!                          
! Input:
!   knl   -- Real(rdef): Multipole strength (mad units).
!   tilt  -- Real(rdef): Multipole tilt.
!   n     -- Real(rdef): Multipole order.
!   coord -- Coord_struct:
!     %x%pos -- X position.
!     %y%pos -- Y position.
!
! Output:
!   coord -- Coord_struct: 
!     %x%vel -- X kick.
!     %y%vel -- Y kick.
!-

!$Id$
!$Log$
!Revision 1.4  2002/06/13 14:54:27  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:21  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:54  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine multipole_kick (knl, tilt, n, coord)

  use bmad
  implicit none

  type (coord_struct)  coord

  real(rdef) knl, tilt, x, y, sin_ang, cos_ang
  real(rdef) x_vel, y_vel

  integer n, m

! simple case

  if (knl == 0) return

! normal case

  if (tilt == 0) then
    x = coord%x%pos
    y = coord%y%pos
  else
    sin_ang = sin(tilt)
    cos_ang = cos(tilt)
    x =  coord%x%pos * cos_ang + coord%y%pos * sin_ang
    y = -coord%x%pos * sin_ang + coord%y%pos * cos_ang
  endif

  x_vel = 0
  y_vel = 0

  do m = 0, n, 2
    x_vel = x_vel + knl * c_multi(n, m) * mexp(x, n-m) * mexp(y, m)
  enddo

  do m = 1, n, 2
    y_vel = y_vel + knl * c_multi(n, m) * mexp(x, n-m) * mexp(y, m)
  enddo

  if (tilt == 0) then
    coord%x%vel = coord%x%vel + x_vel
    coord%y%vel = coord%y%vel + y_vel
  else
    coord%x%vel = coord%x%vel + x_vel * cos_ang - y_vel * sin_ang
    coord%y%vel = coord%y%vel + x_vel * sin_ang + y_vel * cos_ang
  endif

end subroutine
