!+
! Subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_upstream_end, track_downstream_end, 
!                                          orbit_start, orbit_end, ele_start, ele_end, err, compute_floor_coords)
!
! Routine to track a particle within an element.
!
! See also: 
!   twiss_and_track_from_s_to_s
!   twiss_and_track_at_s
!
! The track_upstream_end and track_downstream_end arguments determine whether edge effects are 
! included. These effects include bend edge focusing and aperture checks. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele                  -- Ele_struct: Element to track through.
!   param                -- lat_param_struct:
!     %aperture_limit_on     -- If True check if particle is lost by going outside
!                                 of the element aperture. 
!   l_start              -- Real(rp): Start position measured from the beginning of the element.
!   l_end                -- Real(rp): Stop position measured from the beginning of the element.
!   track_upstream_end   -- Logical: If True then entrance effects are included in the tracking.
!                            But only if l_start = 0.
!   track_downstream_end -- Logical: If True then exit effects are included in the tracking but 
!                            only if l_end = ele%value(l$) (within bmad_com%significant_length tol).
!   orbit_start          -- Coord_struct, optional: Starting phase space coordinates at l_start.
!   ele_start            -- Ele_struct, optional: Holds the starting Twiss parameters at l_start.
!   compute_floor_coords -- logical, optional: If present and True then the global "floor" coordinates 
!                             will be calculated and put in ele_end%floor.
!
! Output:
!   orbit_end  -- Coord_struct, optional: End phase space coordinates. 
!                   If present then the orbit_start argument must also be present.
!   ele_end   -- Ele_struct, optional: Holds the ending Twiss parameters at l_end.
!                  The map (ele_end%mat6, ele_end%vec0) is map from l_start to l_end.
!                  If present then the ele_start argument must also be present.
!   err       -- Logical, optional: Set True if there is a problem like 
!                  the particle gets lost in tracking
!-   

subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_upstream_end, track_downstream_end, &
                                          orbit_start, orbit_end, ele_start, ele_end, err, compute_floor_coords)

use bookkeeper_mod, dummy => twiss_and_track_intra_ele

implicit none

type (coord_struct), optional :: orbit_start, orbit_end
type (coord_struct) orb_at_end
type (ele_struct), optional :: ele_start, ele_end
type (ele_struct) ele
type (lat_param_struct) param
type (ele_struct) :: runt

real(rp) l_start, l_end, mat6(6,6), vec0(6)

logical track_upstream_end, track_downstream_end, do_upstream, do_downstream, err_flag
logical, optional :: err, compute_floor_coords

character(*), parameter :: r_name = 'twiss_and_track_intra_ele'

!

if (present(err)) err = .true.

! zero length element

if (ele%value(l$) == 0) then
  if (.not. track_upstream_end .or. .not. track_downstream_end) then
    call out_io (s_fatal$, r_name, 'Partial tracking through a zero length element does not make sense: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
  ! Track and return
  call track1 (orbit_start, ele, param, orbit_end)
  if (present(err)) err = .false.
  return
endif

! Construct a "runt" element to track through.

call transfer_ele(ele, runt, .true.)

do_upstream = (track_upstream_end .and. l_start == 0)
if (present(orbit_start)) do_upstream = (do_upstream .and. orbit_start%location /= inside$)
do_downstream = (track_downstream_end .and. abs(l_end - ele%value(l$)) < bmad_com%significant_length)

call create_element_slice (runt, ele, l_end - l_start, l_start, param, do_upstream, do_downstream, err_flag, ele_start)
if (err_flag) return

! Now track. 
! Must take care if orbit_start and orbit_end are the same actual argument so use temporary orb_at_end.

if (present(ele_end)) then
  if (present(orbit_start)) then
    call make_mat6 (runt, param, orbit_start, orb_at_end, err_flag = err_flag)
    if (present(orbit_end)) then
      orbit_end = orb_at_end
      orbit_end%ix_ele = ele%ix_ele  ! Since runt%ix_ele gets set to -2 to indicate it is a slice.
      if (.not. do_downstream) orbit_end%location = inside$
    endif
  else
    call make_mat6 (runt, param, err_flag = err_flag)
  endif
  if (err_flag) return
  call twiss_propagate1 (ele_start, runt, err_flag)
  if (logic_option(.false., compute_floor_coords)) call ele_geometry (ele_start%floor, runt, runt%floor)
  call transfer_ele(runt, ele_end, .true.)
  if (err_flag) return

elseif (present(orbit_end)) then  ! and not present(ele_start)
  call track1 (orbit_start, runt, param, orbit_end)
  orbit_end%ix_ele = ele%ix_ele  ! Since runt%ix_ele gets set to -2 to indicate it is a slice.
  if (.not. do_downstream) orbit_end%location = inside$
endif

if (present(err)) err = .false.

end subroutine
