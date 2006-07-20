!+
! Subroutine twiss_and_track_partial (ele1, ele2, param, del_s, ele3, &
!                                                     start, end, body_only)
!
! Subroutine to propagate partially through ELE2 the twiss parameters and the
! orbit. 
!
! This routine assumes that the end point is inside the element so
! for bends a finite e2 exit pole face angle is always ignored.
!
! See also twiss_and_track_at_s.
!
! Note: This routine works by passing to the tracking and Twiss calculation
! routines an element with %value(l$) = del_s. For custom elements, 
! the true length of ele2 is passed in %value(l_original$).
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
!     %mat6(6,6)     -- Transfer matrix (1st order part of xfer map).
!     %vec0(6)       -- 0th order part of the transfer map.
!     %value(l$)     -- Set to DEL_S
!   end  -- Coord_struct, optional: End position at DEL_S.
!-

#include "CESR_platform.inc"

subroutine twiss_and_track_partial (ele1, ele2, param, del_s, ele3, &
                                                     start, end, body_only)

  use bmad_struct
  use bmad_interface, except => twiss_and_track_partial
  use bookkeeper_mod, only: attribute_bookkeeper

  implicit none

  type (ele_struct), optional :: ele3
  type (ele_struct) ele1, ele2
  type (ele_struct), save :: ele
  type (coord_struct), optional :: start, end
  type (coord_struct) c0, c1
  type (param_struct) param

  real(rp) del_s, l_orig, ratio
  integer track, mat6
  logical, optional :: body_only
  character(24) :: r_name = 'twiss_and_track_partial'
  character(80) line

! Error check

  l_orig = ele2%value(l$)

  if (del_s*l_orig < 0 .or. abs(del_s) > abs(l_orig)+1e-6) then
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

! The only real complication comes with a dipole where we have to negate
! the focusing of the exit face (we never get to the exit face since we are
! only partially tracking through).

  ele = ele2
  ele%value(l$) = del_s
  ratio = del_s / l_orig
  ele%value(hkick$) = ele2%value(hkick$) * ratio
  ele%value(vkick$) = ele2%value(vkick$) * ratio
  if (ele%key == hkicker$ .or. ele%key == vkicker$) &
                        ele%value(kick$) = ele2%value(kick$) * ratio

  if (present(start)) then
    c0 = start
  else
    c0%vec = 0
  endif

  track = ele%tracking_method
  mat6 = ele%mat6_calc_method

  select case (ele%key)
  case (wiggler$) 
    if (track == taylor$ .or. track == symp_map$ .or. track == symp_lie_ptc$) &
                                                  ele%tracking_method = symp_lie_bmad$
    if (mat6 == taylor$ .or. track == symp_map$ .or. track == symp_lie_ptc$) &
                                                  ele%mat6_calc_method = symp_lie_bmad$
  case (sbend$)
    if (track == taylor$ .or. track == symp_map$) ele%tracking_method = bmad_standard$
    if (mat6 == taylor$ .or. track == symp_map$)  ele%mat6_calc_method = bmad_standard$
    if (logic_option(.false., body_only)) ele%value(e1$) = 0
    ele%value(e2$) = 0  
  case default
    if (track == taylor$ .or. track == symp_map$) ele%tracking_method = bmad_standard$
    if (mat6 == taylor$ .or. track == symp_map$)  ele%mat6_calc_method = bmad_standard$
  end select

  call track1 (c0, ele, param, c1)

  if (present(end)) end = c1

  if (present(ele3)) then
    if (ele2%key == sbend$) c1%vec = 0
    call make_mat6 (ele, param, c0, c1)   
    call twiss_propagate1(ele1, ele)      
    ele3 = ele
  endif

end subroutine
