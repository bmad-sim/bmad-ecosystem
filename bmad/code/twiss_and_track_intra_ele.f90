!+
! Subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_upstream_end, track_downstream_end, 
!                            orbit_start, orbit_end, ele_start, ele_end, err, compute_floor_coords, reuse_ele_end)
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
! Input:
!   ele                  -- ele_struct: Element to track through.
!   param                -- lat_param_struct:
!   l_start              -- real(rp): Start position measured from the beginning of the element.
!   l_end                -- real(rp): Stop position measured from the beginning of the element.
!   track_upstream_end   -- logical: If True then entrance effects are included in the tracking.
!                            But only if l_start = 0 and orbit_start%location /= inside$.
!   track_downstream_end -- Logical: If True then exit effects are included in the tracking but 
!                            only if l_end = ele%value(l$) (within bmad_com%significant_length tol)
!   orbit_start          -- coord_struct, optional: Starting phase space coordinates at l_start.
!   ele_start            -- ele_struct, optional: Holds the starting Twiss parameters at l_start.
!   ele_end              -- ele_struct, optional: If reuse_ele_end is set True then reuse ele_end from trancking
!                             instead of recomputing ele_end from scratch. This can save time.
!   compute_floor_coords -- logical, optional: If present and True then the global "floor" coordinates 
!                             (without misalignments) will be calculated and put in ele_end%floor.
!   reuse_ele_end        -- logical, optional: If present and True, and if ele_end has the correct 
!                             lonigitudianal length and key type, reuse ele_end from trancking instead of 
!                             recomputing ele_end from scratch. This can save time.
!
! Output:
!   orbit_end  -- Coord_struct, optional: End phase space coordinates. 
!                   If present then the orbit_start argument must also be present.
!   ele_end    -- ele_struct, optional: Holds the ending Twiss parameters at l_end (except for photons).
!                  The map (ele_end%mat6, ele_end%vec0) is map from l_start to l_end.
!                  If present, the ele_start argument must also be present.
!                  NOTE! In the calling routine, deallocate_ele_pointers must not be called on ele_end since
!                   the pointers of ele_end point to the same locations as ele. 
!   err        -- logical, optional: Set True if there is a problem like 
!                  the particle gets lost in tracking
!-   

recursive subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_upstream_end, &
           track_downstream_end, orbit_start, orbit_end, ele_start, ele_end, err, compute_floor_coords, reuse_ele_end)

use bmad_interface, dummy => twiss_and_track_intra_ele

implicit none

type (coord_struct), optional :: orbit_start, orbit_end
type (coord_struct) orb_at_start, orb_at_end
type (ele_struct), optional, target :: ele_start, ele_end
type (ele_struct), target :: ele
type (lat_param_struct) param
type (ele_struct), target :: runt
type (ele_struct), pointer :: ele_p, slave

real(rp) l_start, l_end, mat6(6,6), vec0(6), l0, l1, s_start, s_end, dlength
integer ie, species, dir

logical track_upstream_end, track_downstream_end, do_upstream, do_downstream, err_flag
logical track_up, track_down, length_ok
logical, optional :: err, compute_floor_coords, reuse_ele_end

character(*), parameter :: r_name = 'twiss_and_track_intra_ele'

!

if (present(err)) err = .true.
dir = 1

if (present(orbit_start)) then
  if (orbit_start%state /= alive$) return
  dir = orbit_start%time_dir * orbit_start%direction  ! Direction particle is being propagated.
endif

! If a super_lord then must track through the slaves.

if (ele%lord_status == super_lord$) then
  if (present(orbit_end)) orbit_end = orbit_start
  if (present(ele_end)) call transfer_ele (ele_start, ele_end, .true.)
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

! mirror, multilayer_mirror, and crystal are exceptional if floor coords needed

if (ele%key == mirror$ .or. ele%key == multilayer_mirror$ .or. ele%key == crystal$) then
  do_downstream = (track_downstream_end .or. (orbit_start%location == downstream_end$ .and. orbit_start%ix_ele == ele%ix_ele))
  if (.not. do_downstream) then  ! Then must be center
    if (present(orbit_end)) orbit_end = orbit_start   ! Note: Orbit does not make sense here.
    if (present(ele_end)) then
      ele_end = ele
      if (logic_option(.false., compute_floor_coords)) call ele_geometry (ele_start%floor, ele, ele_end%floor, 0.5_rp)
    endif
    return
  endif
endif



! zero length element:
! Must ignore track_upstream_end and track_downstream_end since they do not make sense in this case.

call transfer_ele(ele, runt, .true.)

if (ele%value(l$) == 0) then
  do_upstream = .true.
  do_downstream = .true.
elseif (dir == 1) then
  do_upstream = (track_upstream_end .and. abs(l_start) < bmad_com%significant_length)
  if (present(orbit_start)) do_upstream = (do_upstream .and. orbit_start%location /= inside$)
  do_downstream = (track_downstream_end .and. abs(l_end - ele%value(l$)) < bmad_com%significant_length)
else
  do_downstream = (track_downstream_end .and. abs(l_start - ele%value(l$)) < bmad_com%significant_length)
  if (present(orbit_start)) do_downstream = (do_downstream .and. orbit_start%location /= inside$)
  do_upstream = (track_upstream_end .and. abs(l_end) < bmad_com%significant_length)
endif

! Special case: tracking though the whole element.
! Otherwise: Construct a "runt" element to track through.

dlength = dir * (l_end - l_start)

ele_p => null()
if (present(ele_end)) then
  if (logic_option(.false., reuse_ele_end) .and. abs(ele_end%value(l$)- dlength) < bmad_com%significant_length .and. ele_end%key == ele%key) then
    ele_p => ele_end
  endif
endif

if (.not. associated(ele_p)) then
  if (do_upstream .and. do_downstream) then
    ele_p => ele
  else
    call create_element_slice (runt, ele, dlength, min(l_start, l_end), param, &
                              do_upstream, do_downstream, err_flag, ele_start, orb_in = orbit_start)
    if (err_flag) then
      if (present(orbit_end)) orbit_end%state = lost$
      return
    endif
    ele_p => runt
  endif
endif

! Now track. 
! Must take care if orbit_start and orbit_end are the same actual argument so use temporary orb_at_end.
! Also: In an lcavity, and depending upon how orbit_start was created (EG: with RK tracking), orbit_start%p0c may correspond
! to ele%value(p0c$) and not runt%value(p0c_start$). If so, shift orbit_start%p0c. Only shift if orbit_start%p0c is
! equal to ele%value(p0c$) to minimize the possibility of papering over a bug in the calling code.

if (present(orbit_start)) then
  species = orbit_start%species
  orb_at_start = orbit_start
  if (ele_p%value(p0c_start$) /= ele%value(p0c$) .and. &
                      .not. significant_difference(orbit_start%p0c, ele%value(p0c$), rel_tol = small_rel_change$)) then
    select case (dir)
    case (1);    call reference_energy_correction(ele_p, orb_at_start, upstream_end$)
    case (-1);   call reference_energy_correction(ele_p, orb_at_start, downstream_end$)
    end select
  endif
else
  species = default_tracking_species(param)
endif

if (present(ele_end) .and. species /= photon$) then
  if (present(orbit_start)) then
    call make_mat6 (ele_p, param, orb_at_start, orb_at_end, err_flag = err_flag)
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
  if (.not. associated(ele_p, ele_end)) call transfer_ele(ele_p, ele_end, .true.)
  if (err_flag) return

elseif (present(orbit_end)) then  ! and not present(ele_start)
  orbit_end = orb_at_start
  select case (dir)
  case (1);   orbit_end%s = ele_p%s_start
  case (-1);  orbit_end%s = ele_p%s
  end select

  call track1 (orbit_end, ele_p, param, orbit_end)

  select case (dir)
  case (1)
    orbit_end%ix_ele = ele%ix_ele;   ! Since ele_p%ix_ele gets set to -2 to indicate it is a slice.
    if (.not. do_downstream) orbit_end%location = inside$
  case (-1)  
    orbit_end%ix_ele = ele%ix_ele;   ! Since ele_p%ix_ele gets set to -2 to indicate it is a slice.
    if (.not. do_upstream) orbit_end%location = inside$
  end select

  if (present(ele_end)) then
    if (logic_option(.false., compute_floor_coords)) call ele_geometry (ele_start%floor, ele_p, ele_p%floor)
    if (.not. associated(ele_p, ele_end)) call transfer_ele(ele_p, ele_end, .true.)
  endif
endif

if (present(err)) err = .false.

end subroutine
