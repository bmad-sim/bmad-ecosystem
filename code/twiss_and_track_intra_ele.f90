!+
! Subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_entrance, 
!                   track_exit, orbit_start, orbit_end, ele_start, ele_end, err)
!
! Routine to track a particle within an element.
!
! See also: 
!   twiss_and_track_from_s_to_s
!   twiss_and_track_at_s
!
! The track_entrance and track_exit arguments determine whether entrance and exit effects are 
! included. These effects are bend edge focusing and aperture checks. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele            -- Ele_struct: Element to track through.
!   param          -- lat_param_struct:
!     %aperture_limit_on -- If True check if particle is lost by going outside
!                of the element aperture. 
!   l_start        -- Real(rp): Start position measured from the beginning of the element.
!   l_end          -- Real(rp): Stop position measured from the beginning of the element.
!   track_entrance -- Logical: If True then entrance effects are included in the tracking.
!                       But only if l_start = 0.
!   track_exit     -- Logical: If True then exit effects are included in the tracking but 
!                       only if l_end = ele%value(l$) (within bmad_com%significant_length tol).
!   orbit_start    -- Coord_struct, optional: Starting phase space coordinates at l_start.
!   ele_start      -- Ele_struct, optional: Holds the starting Twiss parameters at l_start.
!
! Output:
!   orbit_end  -- Coord_struct, optional: End phase space coordinates. 
!             If present then the orbit_start argument must also be present.
!   ele_end   -- Ele_struct, optional: Holds the ending Twiss parameters at l_end.
!                  The map (ele_end%mat6, ele_end%vec0) is map from l_start to l_end.
!             If present then the ele_start argument must also be present.
!   err       -- Logical, optional: Set True if there is a problem like 
!                  the particle gets lost in tracking
!-   

subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_entrance, &
                          track_exit, orbit_start, orbit_end, ele_start, ele_end, err)

use bookkeeper_mod

implicit none

type (coord_struct), optional :: orbit_start, orbit_end
type (coord_struct) orb_at_end
type (ele_struct), optional :: ele_start, ele_end
type (ele_struct) ele
type (lat_param_struct) param
type (ele_struct) :: runt

real(rp) l_start, l_end, mat6(6,6), vec0(6)

logical track_entrance, track_exit, do_entrance, do_exit, err_flag
logical, optional :: err

! Easy case when l_end = l_start

if (l_end == l_start .and. .not. (l_start == 0 .and. track_entrance)) then

  if (present(ele_end)) then
    call transfer_ele(ele_start, ele_end, .true.)
    call mat_make_unit(ele_end%mat6)
    ele_end%vec0 = 0
  endif

  if (present(orbit_end)) then
    orbit_end = orbit_start
  endif

  if (present(err)) err = .false.
  return

endif

! Construct a "runt" element to track through.

if (present(err)) err = .true.
call transfer_ele(ele, runt, .true.)
do_entrance = (track_entrance .and. l_start == 0)
do_exit = (track_exit .and. abs(l_end - ele%value(l$)) < bmad_com%significant_length)
call create_element_slice (runt, ele, l_end - l_start, l_start, param, do_entrance, do_exit, err_flag, ele_start)
if (err_flag) return

! Now track. 
! Must take care if orbit_start and orbit_end are the same actual argument so use temporary orb_at_end.

if (present(orbit_start)) then
  call track1 (orbit_start, runt, param, orb_at_end)
  if (.not. particle_is_moving_forward(orb_at_end)) return
endif

if (present(ele_end)) then
  if (present(orbit_start)) then
    call make_mat6 (runt, param, orbit_start, orb_at_end, .true.)
  else
    call make_mat6 (runt, param)
  endif
  call twiss_propagate1 (ele_start, runt, err_flag)
  call transfer_ele(runt, ele_end, .true.)
  if (err_flag) return
endif

if (present(orbit_end)) then
  orbit_end = orb_at_end
  orbit_end%ix_ele = ele%ix_ele  ! Since runt%ix_ele gets set to -2 to indicate it is a slice.
  if (.not. do_exit) orbit_end%location = inside$
endif

if (present(err)) err = .false.

end subroutine
