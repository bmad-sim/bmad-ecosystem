!+
! Subroutine TWISS_AND_TRACK_PARTIAL (ELE1, ELE2, PARAM, DEL_S, ELE3, &
!                                                                  START, END)
!
! Subroutine to propagate partially through ELE2 the twiss parameters and the
! orbit. See also TWISS_AND_TRACK_BODY.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!     ELE1    -- Ele_struct: Structure holding the starting Twiss parameters.
!     ELE2    -- Ele_struct: Element to partially track through.
!     DEL_S   -- Real: length to propagate.
!     PARAM   -- Param_struct:
!     START   -- Coord_struct (OPTIONAL): Starting position
!                  If not present then START is taken to be 0.
!     END     -- Coord_struct (OPTIONAL): End position at DEL_S.
!
! Output:
!     ELE3 -- Ele_struct (OPTIONAL): Structure for the Twiss results at DEL_S.
!       %X             -- X Twiss parameters
!       %Y             -- Y Twiss parametser
!       %VALUE(L$)     -- Set to DEL_S
!       %VALUE(ANGLE$) -- Set appropriately for an sbend$
!-

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

  real del_s, del
                    
! The only real complication comes with a dipole where we have to negate
! the focusing of the exit face (we never get to the exit face since we are
! only partially tracking through).

  ele = ele2

  if (ele2%key == sbend$) then
    ele%value(angle$) =  ele2%value(angle$) * del_s / ele2%value(l$)
    ele%value(e2$) = 0
  endif

  ele%value(l$) = del_s

  if (present(start)) then
    c0 = start
  else
    c0%vec = 0
  endif

  call track1 (c0, ele, param, c1)

  if (present(end)) then
    end = c1
    if (ele2%key == sbend$) then
      del = -end%x%vel / (ele%value(rho$) * (1 + end%z%vel))
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
