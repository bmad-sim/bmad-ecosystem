!+
! Subroutine twiss_and_track_partial (ele1, ele2, param, del_s, ele3, &
!                                                     start, end, body_only)
!
! Subroutine to propagate partially through ELE2 the twiss parameters and the
! orbit. See also twiss_and_track_at_s.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele1      -- Ele_struct: Structure holding the starting Twiss parameters.
!   ele2      -- Ele_struct: Element to partially track through.
!   del_s     -- Real(rp): length to propagate.
!   param     -- Param_struct:
!   start     -- Coord_struct, optional: Starting position
!                If not present then START is taken to be 0.
!   body_only -- Logical, optional: If present and True then for a bend
!                   no entrence face focusing is applied. That is, the 
!                   tracking starts just inside the element.
!
! Output:
!   ele3 -- Ele_struct, optional: Structure for the Twiss results at DEL_S.
!     %x             -- X Twiss parameters
!     %y             -- Y Twiss parametser
!     %value(l$)     -- Set to DEL_S
!   end  -- Coord_struct, optional: End position at DEL_S.
!-

#include "CESR_platform.inc"

subroutine twiss_and_track_partial (ele1, ele2, param, del_s, ele3, &
                                                     start, end, body_only)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct), optional :: ele3
  type (ele_struct) ele1, ele2
  type (ele_struct), save :: ele
  type (coord_struct), optional :: start, end
  type (coord_struct) c0, c1
  type (param_struct) param

  real(rp) del_s, del, l_orig
  logical, optional :: body_only
  character(20) :: r_name = 'twiss_and_track_partial'
  character(80) line

! Error check

  l_orig = ele2%value(l$)

  if (del_s < 0 .or. del_s > l_orig+1e-6) then
    write (line, '(a, f10.4)') 'DEL_S NEGATIVE OR LARGER THAN ELEMENT LENGTH: ', del_s
    call out_io (s_abort$, r_name, line)
    call err_exit
  endif

! Easy case when ele2 has zero length

  if (l_orig == 0) then
    if (present(ele3)) ele3 = ele2
    if (present(end)) then
      if (present(start)) then
        end = start 
      else
        end%vec = 0
      endif
    endif
    return
  endif

! The only real(rp) complication comes with a dipole where we have to negate
! the focusing of the exit face (we never get to the exit face since we are
! only partially tracking through).

  ele = ele2
  ele%value(l$) = del_s

  if (present(start)) then
    c0 = start
  else
    c0%vec = 0
  endif

  select case (ele%key)
  case (wiggler$) 
    ele%value(n_pole$) = ele2%value(n_pole$) * del_s / l_orig
    if (ele%tracking_method == taylor$) ele%tracking_method = symp_lie_bmad$
    if (ele%mat6_calc_method == taylor$) ele%mat6_calc_method = symp_lie_bmad$
  case (sbend$)
    if (logic_option(.false., body_only)) ele%value(e1$) = 0
    ele%value(e2$) = 0  
  end select

  ele%num_steps = max(nint(ele%num_steps * del_s / l_orig), 1)

  call attribute_bookkeeper (ele, param)

  call track1 (c0, ele, param, c1)

  if (present(end)) end = c1

  if (present(ele3)) then
    if (ele2%key == sbend$) c1%vec = 0
    call make_mat6 (ele, param, c0, c1)   
    call twiss_propagate1(ele1, ele)      
    ele3 = ele
  endif

end subroutine
