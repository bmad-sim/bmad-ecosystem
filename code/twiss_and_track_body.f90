!+
! Subroutine twiss_and_track_body (ele1, ele2, param, del_s, ele3, start, end)
!
! Subroutine to propagate partially through ELE2 the twiss parameters and the
! orbit. The only difference between this and TWISS_AND_TRACK_PARTIAL is that
! in this routine, for a bend, there is no entrence face focusing. That is,
! the tracking starts just inside the bend.
!
! Modules needed:
!   use bmad
!
! Input:
!     ele1    -- Ele_struct: Structure holding the starting Twiss parameters.
!     ele2    -- Ele_struct: Element to partially track through.
!     del_s   -- Real(rdef): length to propagate.
!     param   -- Param_struct:
!     start   -- Coord_struct (OPTIONAL): Starting position
!                  If not present then START is taken to be 0.
!     end     -- Coord_struct (OPTIONAL): End position at DEL_S.
!
! Output:
!     ele3 -- Ele_struct (OPTIONAL): Structure for the Twiss results at DEL_S.
!       %x             -- X Twiss parameters
!       %y             -- Y Twiss parametser
!       %value(l$)     -- Set to DEL_S
!-

!$Id$
!$Log$
!Revision 1.5  2003/01/27 14:40:46  dcs
!bmad_version = 56
!
!Revision 1.4  2002/06/13 14:54:30  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:27  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:32:00  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine twiss_and_track_body (ele1, ele2, param, del_s, ele3, &
                                                               start, end)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct), optional :: ele3
  type (ele_struct) ele1, ele2, ele
  type (coord_struct), optional :: start, end
  type (coord_struct) c0, c1
  type (param_struct) param

  real(rdef) del_s, del
                    
! The only real(rdef) complication comes with a dipole where we have to negate
! the focusing of the faces. 

  ele = ele2

  if (present(start)) then
    c0 = start
  else
    c0%vec = 0
  endif

  if (ele2%key == sbend$) then
    ele%value(e1$) = -atan(c0%x%vel) * ele%value(g$) / (1 + c0%z%vel)
    ele%value(e2$) = 0
  endif

  ele%value(l$) = del_s

  call track1 (c0, ele, param, c1)

  if (present(end)) then
    end = c1
    if (ele2%key == sbend$) then
      del = -end%x%vel * ele%value(g$) / (1 + end%z%vel)
      end%x%vel = end%x%vel - del * end%x%pos
      end%y%vel = end%y%vel + del * end%y%pos
    endif
  endif

  if (present(ele3)) then
    if (ele2%key == sbend$) c1%vec = 0
    call make_mat6 (ele, param, c0, c1)   
    call twiss_propagate1(ele1, ele)      
    ele3 = ele
  endif

end subroutine
