!+
! Subroutine twiss_and_track_partial (ele0, ele_track, param, del_s, ele_end, 
!                                                              orb_start, orb_end, body_only, err)
!
! Subroutine to propagate partially through ele_track the twiss parameters and the
! orbit. 
!
! This routine assumes that the end point is inside the element so
! for bends a finite e2 exit pole face angle is always ignored.
!
! For a bend, for del_s = 0, and with body_only = F, the entrance pole face is tracked.
!
! See also: 
!   twiss_and_track_at_s
!   twiss_and_track_intra_ele
!
! Note: This routine works by passing to the tracking and Twiss calculation
! routines an element with %value(l$) = del_s. For custom elements, 
! the true length of ele_track is passed in %value(l_original$).
!
! Modules needed:
!   use bmad
!
! Input:
!   ele0      -- Ele_struct: Structure holding the starting Twiss parameters.
!   ele_track -- Ele_struct: Element to partially track through.
!   del_s     -- Real(rp): length to propagate.
!   param     -- lat_param_struct:
!   orb_start -- Coord_struct, optional: Starting position
!                  If not present then ORB_START is taken to be 0.
!   body_only -- Logical, optional: If present and True then for a bend
!                   no entrence face focusing is applied. That is, the 
!                   tracking starts just inside the element.
!
! Output:
!   ele_end -- Ele_struct, optional: Structure for the Twiss results at DEL_S.
!     %a             -- X Twiss parameters
!     %b             -- Y Twiss parametser
!     %mat6(6,6)     -- Transfer matrix (1st order part of xfer map).
!     %vec0(6)       -- 0th order part of the transfer map.
!     %value(l$)     -- Set to DEL_S
!   orb_end -- Coord_struct, optional: End position at DEL_S.
!   err     -- Logical, optional: Set True if there is a problem like 
!                the particle gets lost in tracking
!-

subroutine twiss_and_track_partial (ele0, ele_track, param, del_s, ele_end, &
                                                       orb_start, orb_end, body_only, err)

use bmad_struct
use bmad_interface, except_dummy => twiss_and_track_partial
use bookkeeper_mod, only: attribute_bookkeeper, create_element_slice

implicit none

type (ele_struct), optional :: ele_end
type (ele_struct) ele0, ele_track
type (ele_struct), save :: runt
type (coord_struct), optional :: orb_start, orb_end
type (coord_struct) c0, c1
type (lat_param_struct) param

real(rp) del_s, l_orig, ratio
integer track, mat6
character(24) :: r_name = 'twiss_and_track_partial'
character(80) line
logical, optional :: body_only, err
logical error, track_entrance

! Error check

if (present(err)) err = .true.
l_orig = ele_track%value(l$)
track_entrance = .not. logic_option(.false., body_only)

if (del_s*l_orig < 0 .or. abs(del_s) > abs(l_orig)+1e-6) then
  write (line, '(a, f10.4)') 'DEL_S NEGATIVE OR LARGER THAN ELEMENT LENGTH: ', del_s
  call out_io (s_abort$, r_name, line)
  call err_exit
endif

! Easy case when del_s is zero.

if (del_s == 0 .and. .not. track_entrance) then

  if (present(ele_end)) then
    ele_end = ele0
    call mat_make_unit(ele_end%mat6)
    ele_end%vec0 = 0
  endif

  if (present(orb_end)) then
    if (present(orb_start)) then
      orb_end = orb_start 
    else
      orb_end%vec = 0
    endif
  endif

  if (present(err)) err = .false.
  return

endif

! Create the runt element to track through.

runt = ele_track
call create_element_slice (runt, ele_track, del_s, 0.0_rp, param, track_entrance, .false.)

! If not easy case then do tracking...

if (present(orb_start)) then
  c0 = orb_start
else
  c0%vec = 0
endif

track = runt%tracking_method
mat6 = runt%mat6_calc_method

select case (runt%key)
case (wiggler$) 
  if (track == taylor$ .or. track == symp_map$ .or. track == symp_lie_ptc$) &
                                                runt%tracking_method = symp_lie_bmad$
  if (mat6 == taylor$ .or. track == symp_map$ .or. track == symp_lie_ptc$) &
                                                runt%mat6_calc_method = symp_lie_bmad$
case default
  if (track == taylor$ .or. track == symp_map$) runt%tracking_method = bmad_standard$
  if (mat6 == taylor$ .or. track == symp_map$)  runt%mat6_calc_method = bmad_standard$
end select

call track1 (c0, runt, param, c1)
if (param%lost) return

if (present(orb_end)) orb_end = c1

if (present(ele_end)) then
  if (ele_track%key == sbend$) c1%vec = 0
  call make_mat6 (runt, param, c0, c1, .true.)   
  call twiss_propagate1(ele0, runt, error)      
  if (error) return
  ele_end = runt
endif

if (present(err)) err = .false.

end subroutine
