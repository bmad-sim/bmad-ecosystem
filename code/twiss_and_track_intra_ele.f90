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
!   l_start              -- Real(rp): Start position measured from the beginning of the element.
!   l_end                -- Real(rp): Stop position measured from the beginning of the element.
!   track_upstream_end   -- Logical: If True then entrance effects are included in the tracking.
!                            But only if l_start = 0 and orbit_start%location /= inside$.
!   track_downstream_end -- Logical: If True then exit effects are included in the tracking but 
!                            only if l_end = ele%value(l$) (within bmad_com%significant_length tol)
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

recursive subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_upstream_end, &
                       track_downstream_end, orbit_start, orbit_end, ele_start, ele_end, err, compute_floor_coords)

use bookkeeper_mod, dummy => twiss_and_track_intra_ele

implicit none

type (coord_struct), optional :: orbit_start, orbit_end
type (coord_struct) orb_at_end
type (ele_struct), optional :: ele_start, ele_end
type (ele_struct), target :: ele
type (lat_param_struct) param
type (ele_struct), target :: runt
type (ele_struct), pointer :: ele_p, slave

real(rp) l_start, l_end, mat6(6,6), vec0(6), l0, l1, s_start, s_end
integer ie

logical track_upstream_end, track_downstream_end, do_upstream, do_downstream, err_flag
logical track_up, track_down
logical, optional :: err, compute_floor_coords

character(*), parameter :: r_name = 'twiss_and_track_intra_ele'

!

if (present(err)) err = .true.

! If a super_lord then must track through the slaves.

if (ele%lord_status == super_lord$) then
  if (present(orbit_end)) orbit_end = orbit_start
  if (present(ele_end)) ele_end = ele_start
  s_start = ele%s_start + l_start
  s_end   = ele%s_start + l_end
  do ie = 1, ele%n_slave
    slave => pointer_to_slave (ele, ie)
    if (slave%s < s_start) cycle

    l0 = max(0.0_rp, s_start - slave%s_start)
    l1 = min(slave%value(l$), s_end - slave%s_start)
    track_up = (track_upstream_end .or. s_start < slave%s_start)
    track_down = (track_downstream_end .or. s_end > slave%s)
    call twiss_and_track_intra_ele (slave, param, l0, l1, track_up, &
                       track_down, orbit_end, orbit_end, ele_end, ele_end, err, compute_floor_coords)
    if (s_end <= slave%s + bmad_com%significant_length) exit
  enddo
  return
endif

! zero length element:
! Must ignore track_upstream_end and track_downstream_end since they do not make sense in this case.

call transfer_ele(ele, runt, .true.)

if (ele%value(l$) == 0) then
  do_upstream = .true.
  do_downstream = .true.
else
  do_upstream = (track_upstream_end .and. abs(l_start) < bmad_com%significant_length)
  if (present(orbit_start)) do_upstream = (do_upstream .and. orbit_start%location /= inside$)
  do_downstream = (track_downstream_end .and. abs(l_end - ele%value(l$)) < bmad_com%significant_length)
endif

! Special case: tracking though the whole element.
! Otherwise: Construct a "runt" element to track through.

if (do_upstream .and. do_downstream) then
  ele_p => ele
else
  call create_element_slice (runt, ele, l_end - l_start, l_start, param, do_upstream, do_downstream, err_flag, ele_start)
  if (err_flag) return
  ele_p => runt
endif

! Now track. 
! Must take care if orbit_start and orbit_end are the same actual argument so use temporary orb_at_end.

if (present(ele_end)) then
  if (present(orbit_start)) then
    call make_mat6 (ele_p, param, orbit_start, orb_at_end, err_flag = err_flag)
    if (present(orbit_end)) then
      orbit_end = orb_at_end
      orbit_end%ix_ele = ele%ix_ele  ! Since ele_p%ix_ele gets set to -2 to indicate it is a slice.
      if (.not. do_downstream) orbit_end%location = inside$
    endif
  else
    call make_mat6 (ele_p, param, err_flag = err_flag)
  endif
  if (err_flag) return
  call twiss_propagate1 (ele_start, ele_p, err_flag)
  if (logic_option(.false., compute_floor_coords)) call ele_geometry (ele_start%floor, ele_p, ele_p%floor)
  call transfer_ele(ele_p, ele_end, .true.)
  if (err_flag) return

elseif (present(orbit_end)) then  ! and not present(ele_start)
  call track1 (orbit_start, ele_p, param, orbit_end)
  orbit_end%ix_ele = ele%ix_ele  ! Since ele_p%ix_ele gets set to -2 to indicate it is a slice.
  if (.not. do_downstream) orbit_end%location = inside$
endif

if (present(err)) err = .false.

end subroutine
