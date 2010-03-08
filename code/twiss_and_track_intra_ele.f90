!+
! Subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_entrance, 
!                               track_exit, orbit_start, orbit_end, ele_start, ele_end, err)
!
! Routine to track a particle within an element.
!
! See also: 
!   twiss_and_track_partial
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
!                       only if l_end = ele%value(l$) (within bmad_com%significant_longitudinal_length tol).
!   orbit_start    -- Coord_struct, optional: Starting phase space coordinates at l_start.
!   ele_start      -- Ele_struct, optional: Holds the starting Twiss parameters at l_start.
!
! Output:
!   param      -- lat_param_struct: 
!     %lost          -- Set True If the particle cannot make it through an element.
!                         Set False otherwise.
!     %plane_lost_at -- x_plane$, y_plane$ (for apertures), or 
!                         z_plane$ (turned around in an lcavity).
!     %end_lost_at   -- entrance_end$ or exit_end$.
!   orbit_end  -- Coord_struct, optional: End phase space coordinates. 
!             If present then the orbit_start argument must also be present.
!   ele_end   -- Ele_struct, optional: Holds the ending Twiss parameters and the transfer matrix.
!             If present then the orbit_start argument must also be present.
!   err       -- Logical, optional: Set True if there is a problem like 
!                  the particle gets lost in tracking
!-   

subroutine twiss_and_track_intra_ele (ele, param, l_start, l_end, track_entrance, &
                          track_exit, orbit_start, orbit_end, ele_start, ele_end, err)

use bookkeeper_mod

implicit none

type (coord_struct), optional :: orbit_start, orbit_end
type (ele_struct), optional ::ele_start, ele_end
type (ele_struct) ele
type (lat_param_struct) param
type (ele_struct), save :: runt

real(rp) l_start, l_end

integer track, mat6

logical track_entrance, track_exit, do_entrance, do_exit
logical, optional :: err

! Easy case when l_end = l_start

if (l_end == l_start .and. .not. (l_start == 0 .and. track_entrance)) then

  if (present(ele_end)) then
    ele_end = ele_start
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
! A superposition slave is a section of the lord element so makeup_super_slave1 
! does the needed modification to the runt element.

runt = ele
runt%value(l$) = l_end - l_start
do_entrance = (track_entrance .and. l_start == 0)
do_exit = (track_exit .and. abs(l_end - ele%value(l$)) < bmad_com%significant_longitudinal_length)
call makeup_super_slave1 (runt, ele, l_start, param, do_entrance, do_exit)
call attribute_bookkeeper (runt, param)

track = runt%tracking_method
mat6  = runt%mat6_calc_method

select case (runt%key)
case (wiggler$) 
  if (track == taylor$ .or. track == symp_map$ .or. track == symp_lie_ptc$) &
                                                runt%tracking_method = symp_lie_bmad$
  if (mat6 == taylor$ .or. track == symp_map$ .or. track == symp_lie_ptc$) &
                                                runt%mat6_calc_method = symp_lie_bmad$
case default
  if (track == taylor$ .or. track == symp_map$) runt%tracking_method  = bmad_standard$
  if (mat6 == taylor$ .or. track == symp_map$)  runt%mat6_calc_method = bmad_standard$
end select

! Now track

if (present(orbit_end)) then
  call track1 (orbit_start, runt, param, orbit_end)
  if (param%lost) return
endif

if (present(ele_end)) then
  call make_mat6 (runt, param)
  call twiss_propagate1 (ele_start, runt, err)
  ele_end = runt
  return
endif

if (present(err)) err = .false.

end subroutine
