!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+ 
! Subroutine track1_time_runge_kutta(start_orb, ele, param, end_orb, err_flag, track)
!
! Routine to track a particle through an element using 
! Runge-Kutta time-based tracking. Converts to and from element
! coordinates before and after tracking.
!
! Note that the argument tref in is arbitrary in the function
! convert_particle_coordinates_t_to_s, so the value ele%ref_time is
! used. Consistency is what matters.
!
! Modules Needed:
!   use time_tracker_mod
!
! Input:
!   start_orb                  -- coord_struct: starting position, t-based global
!   ele                    -- ele_struct: element
!    %value                -- real(rp): attribute values
!    %ref_time             -- real(rp): time ref particle passes exit end
!    %s                    -- real(rp): longitudinal ref position at exit end
!   param                  -- lat_param_struct: lattice parameters
!
! Output:
!   end_orb     -- coord_struct: end position, t-based global
!   err_flag    -- Logical: Set True if there is an error. False otherwise
!   track       -- track_struct (optional): Contains array of the step-by-step particle
!                    trajectory along with the field at these positions.
!                    When tracking through multiple elements, the trajectory in an element
!                    is appended to the existing trajectory. To reset: Set track%n_pt = -1.
!-


subroutine track1_time_runge_kutta (start_orb, ele, param, end_orb, err_flag, track)

use time_tracker_mod, dummy => track1_time_runge_kutta
use track1_mod, dummy2 => track1_time_runge_kutta

implicit none

type (coord_struct) :: start_orb, end_orb
type (coord_struct) :: ele_origin
type (lat_param_struct), target :: param
type (ele_struct), target :: ele
type (ele_struct), pointer :: hard_ele
type (track_struct), optional :: track

real(rp) dt_step, ref_time, vec6
real(rp) s_rel, time, s1, s2, del_s, p0c_save, s_save
real(rp) s_edge_track, s_edge_hard

integer :: i, hard_end

character(30), parameter :: r_name = 'track1_time_runge_kutta'

logical :: local_ref_frame = .true.
logical :: abs_time, err_flag, err

!---------------------------------

err_flag = .true.

end_orb = start_orb
time = particle_time(start_orb, ele)

! If element has zero length, skip tracking

if (ele%value(l$) .eq. 0) then
  
  !If saving tracks, allocate track array and save one point
  if ( present(track) ) then
    !Convert to global-s to local-t coordinates
    call convert_particle_coordinates_s_to_t(end_orb)
    !Tracks use vec(5) = s_rel
    end_orb%vec(5) = 0.0_rp
    call save_a_step (track, ele, param, .false., end_orb%vec(5), end_orb, s_save, time)
  endif

  ! Reset particle to s-coordinates
  end_orb = start_orb
  if (end_orb%direction == +1) then
    end_orb%location = downstream_end$
  else
    end_orb%location = upstream_end$
  endif
  err_flag = .false.
  return
end if

!------
! Specify initial time step.

dt_step = bmad_com%init_ds_adaptive_tracking / c_light

! Convert particle to element coordinates
! Kicks and multipoles should be turned off in offset_particle

! Particle is moving forward towards the entrance
if (end_orb%direction == +1 .and. end_orb%location /= inside$) then
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., set_multipoles = .false., ds_pos = 0.0_rp ) 

! Interior start, reference momentum is at the end. No edge kicks are given
elseif (end_orb%location == inside$) then
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., set_multipoles = .false., &
                        ds_pos = end_orb%s - (ele%s - ele%value(l$)) )

! Particle is at the exit surface, should be moving backwards
elseif (end_orb%direction == -1 .and. end_orb%location /= inside$) then
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., set_multipoles = .false., &
                        ds_pos = end_orb%s - (ele%s - ele%value(l$)) )

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE ENTERING ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! ele(s-based) -> ele(t-based)

call convert_particle_coordinates_s_to_t(end_orb)
s_rel =  end_orb%s - (ele%s - ele%value(l$) )
if ( abs(s_rel)  < bmad_com%significant_length ) then
  s_rel = 0
else if ( abs(s_rel - ele%value(l$))  < bmad_com%significant_length ) then
  s_rel = ele%value(l$)
endif

end_orb%vec(5) = s_rel

!------
!Check wall

if ( wall3d_d_radius(end_orb%vec, ele) > 0 ) then

  call out_io (s_info$, r_name, "PARTICLE STARTED IN REGION OUTSIDE OF WALL: "//trim(ele%name)//", SETTING TO LOST")
  !Particle won't be tracked, so set end = start with the saved state
  end_orb = start_orb
  end_orb%state = lost$
  return

endif

if ( present(track) ) then
  call save_a_step (track, ele, param, .false., end_orb%vec(5), end_orb, s_save, time)
endif

! Track through element

call odeint_bmad_time(end_orb, ele, param, 0.0_rp, ele%value(l$), time, dt_step, local_ref_frame, err, track)
if (err) return

!------
!Convert back to s-based coordinates

if (end_orb%location /= inside$ .and. end_orb%vec(6) < 0) then
  !Particle left entrance end going backwards
  !set reference time and momentum
  end_orb%p0c = ele%value(p0c_start$)
  end_orb%direction = -1
  ref_time = ele%ref_time - ele%value(delta_ref_time$)

  !ele(t-based) -> ele(s-based)
  call convert_particle_coordinates_t_to_s(end_orb, ref_time)
  ! call apply_element_edge_kick (end_orb, ele, param, upstream_end$)
  !unset
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., set_multipoles = .false.)

elseif (end_orb%state /= alive$) then
  !Particle is lost in the interior of the element.
  !  The reference is a the end of the element
  if (end_orb%vec(6) < 0) then
    end_orb%p0c = ele%value(p0c$)
    end_orb%direction = -1
  else 
    end_orb%p0c = ele%value(p0c$)
    end_orb%direction = 1
  end if

  ref_time = ele%ref_time
  !ele(t-based) -> ele(s-based)
  call convert_particle_coordinates_t_to_s(end_orb, ref_time)
  !unset
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., set_multipoles = .false., &
                          ds_pos = end_orb%s - (ele%s - ele%value(l$)) )

elseif (end_orb%location /= inside$ .and. end_orb%vec(6) >= 0) then
  !Particle left exit end going forward
  end_orb%p0c = ele%value(p0c$)
  end_orb%direction = 1
  ref_time = ele%ref_time
  !ele(t-based) -> ele(s-based)
  call convert_particle_coordinates_t_to_s(end_orb, ref_time)
  !call apply_element_edge_kick (end_orb, ele, param, downstream_end$)
  !unset
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., set_multipoles = .false.)

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE LEAVING ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

!Set relativistic beta
call convert_pc_to (end_orb%p0c * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)
err_flag = .false.

end subroutine

