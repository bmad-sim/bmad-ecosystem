!+ 
! Subroutine track1_time_runge_kutta(start_orb, ele, param, end_orb, err_flag, track)
!
! Routine to track a particle through an element using 
! Runge-Kutta time-based tracking. Converts to and from element
! coordinates before and after tracking.
!
! Modules Needed:
!   use time_tracker_mod
!
! Input:
!   start_orb   -- coord_struct: starting position, t-based global
!   ele         -- ele_struct: element
!   param       -- lat_param_struct: lattice parameters
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
type (em_field_struct) :: saved_field

real(rp) vec(6), d_radius
real(rp) s_rel, s0, s1, s2, ds_ref, del_s, p0c_save
real(rp) s_edge_track, s_edge_hard, rf_time, beta_ref

integer :: i, hard_end, t_dir

logical :: local_ref_frame = .true.
logical :: abs_time, err_flag, err, set_spin

character(*), parameter :: r_name = 'track1_time_runge_kutta'

!---------------------------------

if (ele%key /= patch$ .and. ele%value(l$) == 0) then
  call track_a_zero_length_element (start_orb, ele, param, end_orb, err_flag, track)
  return
endif

!

err_flag = .true.

t_dir = 1
end_orb = start_orb
rf_time = particle_rf_time (end_orb, ele, .true., end_orb%s - ele%s_start)
set_spin = (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$ .and. &
            (ele%field_calc == bmad_standard$ .or. ele%field_calc == fieldmap$) .and. &
            is_true(ele%value(spin_fringe_on$)))

! Relative s. Adjust to match element edges if close enough
s_rel =  end_orb%s - ele%s_start
if (abs(s_rel) < bmad_com%significant_length) then
  s_rel = 0
else if (abs(s_rel - ele%value(l$))  < bmad_com%significant_length) then
  s_rel = ele%value(l$)
endif

!------
! Check wall

call check_aperture_limit (end_orb, ele, in_between$, param)
if (end_orb%state /= alive$) then
  call out_io (s_info$, r_name, "PARTICLE STARTED IN REGION OUTSIDE OF WALL: "//trim(ele%name))
  !Particle won't be tracked, so set end = start with the saved state
  end_orb = start_orb
  end_orb%state = lost$
  return
endif

!------
! Convert particle to element coordinates
! Kicks and multipoles should be turned off in offset_particle

! Interior start, reference momentum is at the end.
if (end_orb%location == inside$) then
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., ds_pos =s_rel, set_spin = set_spin)
  if (ele%value(l$) < 0) t_dir = -1

elseif (ele%key == patch$) then
  if (start_orb%location == inside$) then
    call out_io (s_error$, r_name, 'TIME-RUNGE-KUTTA TRACKING STARTING WITH A INSIDE A PATCH IS NOT PERMITED: ' // ele%name)
    return
  endif
  call track_a_patch (ele, end_orb, .false., s0, ds_ref)
  beta_ref = ele%value(p0c$) / ele%value(E_tot$)
  end_orb%vec(5) = end_orb%vec(5) + (ds_ref + s0 * end_orb%direction * ele%orientation) * end_orb%beta / beta_ref 
  if (s0*ele%orientation > 0) t_dir = -1
  s_rel = s0

! Particle is at an end.
else
  call offset_particle (ele, param, set$, end_orb, set_hvkicks = .false., set_spin = set_spin)
  if (ele%value(l$) < 0) t_dir = -1
endif

! ele(s-based) -> ele(t-based)

call convert_particle_coordinates_s_to_t(end_orb, ele, s_rel)

if ( present(track) ) then
  ! here local_ref_frame is false to avoid calling offset_particle, because we are in time coordinates
  ! This should be the same as done inside odeint_bmad_time 
  call save_a_step (track, ele, param, .false., end_orb, end_orb%vec(5), .true., rf_time = rf_time)
  call em_field_calc (ele, param, end_orb%vec(5), end_orb, .true., saved_field, .false., err_flag, rf_time = rf_time)
  if (err_flag) return
  track%field(track%n_pt) = saved_field
endif

! Track through element

if ((ele%key == lcavity$ .or. ele%key == rfcavity$) .and. ele%field_calc == bmad_standard$ .and. &
                                                         ele%value(l$) < ele%value(l_hard_edge$)) then
  call out_io (s_error$, r_name, 'TIME-RUNGE-KUTTA TRACKING THROUGH RF CAVITY: ' // ele%name, &
                          'WILL NOT BE ACCURATE SINCE THE LENGTH IS LESS THAN THE HARD EDGE MODEL LENGTH.')
endif

call odeint_bmad_time(end_orb, ele, param, t_dir, rf_time, local_ref_frame, err, track)

if (err) return

!------
! Convert back to s-based coordinates.
! The particle is either dead or is alive and at an element end.

if (end_orb%location == upstream_end$) then
  end_orb%p0c = ele%value(p0c_start$)
  call convert_particle_coordinates_t_to_s(end_orb, ele, ele%value(ref_time_start$))
  end_orb%direction = -1  ! In case t_to_s conversion confused by roundoff error
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., set_spin = set_spin)

elseif (end_orb%location == downstream_end$) then
  end_orb%p0c = ele%value(p0c$)
  call convert_particle_coordinates_t_to_s(end_orb, ele, ele%ref_time)
  end_orb%direction = 1  ! In case t_to_s conversion confused by roundoff error
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., set_spin = set_spin)

elseif (end_orb%state /= alive$) then
  ! Particle is lost in the interior of the element.
  ! The reference energy in the interior is equal to the ref energy at the exit end of the element.
  if (end_orb%vec(6) < 0) then
    end_orb%p0c = ele%value(p0c$)
    end_orb%direction = -1
  else 
    end_orb%p0c = ele%value(p0c$)
    end_orb%direction = 1
  end if

  call convert_particle_coordinates_t_to_s(end_orb, ele, ele%ref_time)
  call offset_particle (ele, param, unset$, end_orb, set_hvkicks = .false., ds_pos = end_orb%s - ele%s_start, set_spin = set_spin)

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE LEAVING ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! Set relativistic beta

call convert_pc_to (end_orb%p0c * (1 + end_orb%vec(6)), end_orb%species, beta = end_orb%beta)
err_flag = .false.

end subroutine

