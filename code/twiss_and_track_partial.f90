!+
! Subroutine twiss_and_track_partial (ele1, ele2, param, del_s, ele3, &
!                                                                  start, end)
!
! Subroutine to propagate partially through ELE2 the twiss parameters and the
! orbit. See also twiss_and_track_body and twiss_and_track_at_s.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1    -- Ele_struct: Structure holding the starting Twiss parameters.
!   ele2    -- Ele_struct: Element to partially track through.
!   del_s   -- Real(rp): length to propagate.
!   param   -- Param_struct:
!   start   -- Coord_struct, optional: Starting position
!                If not present then START is taken to be 0.
!   end     -- Coord_struct, optional: End position at DEL_S.
!
! Output:
!   ele3 -- Ele_struct, optional: Structure for the Twiss results at DEL_S.
!     %x             -- X Twiss parameters
!     %y             -- Y Twiss parametser
!     %value(l$)     -- Set to DEL_S
!-

#include "CESR_platform.inc"

subroutine twiss_and_track_partial (ele1, ele2, param, del_s, ele3, &
                                                               start, end)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct), optional :: ele3
  type (ele_struct) ele1, ele2, ele
  type (coord_struct), optional :: start, end
  type (coord_struct) c0, c1
  type (param_struct) param

  real(rp) del_s, del
                    
! The only real(rp) complication comes with a dipole where we have to negate
! the focusing of the exit face (we never get to the exit face since we are
! only partially tracking through).

  ele = ele2

  ele%value(l$) = del_s
  if (ele%key == wiggler$) ele%value(n_pole$) = &
                              ele2%value(n_pole$) * del_s / ele2%value(l$)

  if (present(start)) then
    c0 = start
  else
    c0%vec = 0
  endif

  call track1 (c0, ele, param, c1)

  if (present(end)) then
    end = c1
    if (ele2%key == sbend$) then
      del = -end%vec(2) * ele%value(g$) / (1 + end%vec(6))
      end%vec(2) = end%vec(2) - del * end%vec(1)
      end%vec(4) = end%vec(4) + del * end%vec(3)
    endif
  endif

  if (present(ele3)) then
    if (ele2%key == sbend$) c1%vec = 0
    call make_mat6 (ele, param, c0, c1)   
    call twiss_propagate1(ele1, ele)      
    ele3 = ele
  endif

end subroutine
