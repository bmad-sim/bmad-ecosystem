!+ 
! Subroutine track1_time_runge_kutta(orbit, ele, param, err_flag, track, t_end, dt_step)
!
! Routine to track a particle through an element using 
! Runge-Kutta time-based tracking. Converts to and from element
! coordinates before and after tracking.
!
! Input:
!   orbit       -- coord_struct: starting position, z-based coords
!   ele         -- ele_struct: element
!   param       -- lat_param_struct: lattice parameters
!   t_end       -- real(rp), optional: If present, maximum time to which the particle will be tracked.
!                   Used for tracking with given time steps. The time orb%t at which tracking stops 
!                   may be less than this if the particle gets to the end of the element
!   dt_step     -- real(rp), optional: If positive, next RK time step to take. 
!                   This overrides bmad_com%init_ds_adaptive_tracking. Used by track_bunch_time.
!
! Output:
!   orbit       -- coord_struct: end position, z-based coords
!   err_flag    -- Logical: Set True if there is an error. False otherwise
!   track       -- track_struct (optional): Contains array of the step-by-step particle
!                    trajectory along with the field at these positions.
!                    When tracking through multiple elements, the trajectory in an element
!                    is appended to the existing trajectory. To reset: Set track%n_pt = -1.
!   dt_step     -- real(rp), optional: Next RK time step that this tracker would take based on the error tolerance.
!                   Used by track_bunch_time.
!-

subroutine track1_time_runge_kutta (orbit, ele, param, err_flag, track, t_end, dt_step)

use time_tracker_mod, dummy => track1_time_runge_kutta
use bmad_interface, dummy2 => track1_time_runge_kutta

implicit none

type (coord_struct) :: orbit
type (coord_struct) :: ele_origin
type (lat_param_struct), target :: param
type (ele_struct), target :: ele
type (ele_struct), pointer :: hard_ele
type (track_struct), optional :: track
type (em_field_struct) :: saved_field

real(rp), optional :: t_end, dt_step
real(rp) vec(6), d_radius
real(rp) s_lab, s0, s1, s2, ds_ref, del_s, p0c_save, s_body
real(rp) s_edge_track, s_edge_hard, rf_time, beta_ref, r, dref_time

integer :: i, hard_end, t_dir

logical :: abs_time, err_flag, err, set_spin

character(*), parameter :: r_name = 'track1_time_runge_kutta'

!---------------------------------

if (ele%key /= patch$ .and. ele%value(l$) == 0) then
  call track_a_zero_length_element (orbit, ele, param, err_flag, track)
  return
endif

if (present(track)) then
  call save_a_step (track, ele, param, .false., orbit, s_lab, .true., rf_time = rf_time)
endif

!

err_flag = .true.

rf_time = particle_rf_time (orbit, ele, .true., orbit%s - ele%s_start)
beta_ref = ele%value(p0c$) / ele%value(E_tot$)
set_spin = (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$)

! s_lab is longitudinal position relative to entrance end (element body coords).
! Remember that if an element has reversed orientation, +s direction in body coords is opposite to
! the reference trajectory +s direction.
! Adjust to match element edges if close enough.

s_lab =  orbit%s - ele%s_start

if (abs(s_lab) < bmad_com%significant_length) then
  s_lab = 0
else if (abs(s_lab - ele%value(l$))  < bmad_com%significant_length) then
  s_lab = ele%value(l$)
endif

!------
! Check wall

call check_aperture_limit (orbit, ele, in_between$, param)
if (orbit%state /= alive$) then
  call out_io (s_info$, r_name, "PARTICLE STARTED IN REGION OUTSIDE OF WALL: "//trim(ele%name))
  orbit%state = lost$
  return
endif

!------
! Convert particle to element coordinates
! Kicks and multipoles should be turned off in offset_particle

! Interior start, reference momentum is at the end.
if (orbit%location == inside$) then
  call offset_particle (ele, set$, orbit, set_hvkicks = .false., s_pos = s_lab, s_out = s_body, set_spin = set_spin)
  t_dir = sign_of(ele%value(l$)) * orbit%time_dir

elseif (ele%key == patch$) then
  if (orbit%location == inside$) then
    call out_io (s_error$, r_name, 'TIME-RUNGE-KUTTA TRACKING STARTING WITH A INSIDE A PATCH IS NOT PERMITED: ' // ele%name)
    return
  endif
  call track_a_patch (ele, orbit, .false., s0, ds_ref)
  orbit%vec(5) = orbit%vec(5) + orbit%time_dir * (ds_ref + s0 * orbit%direction*orbit%time_dir * ele%orientation) * orbit%beta / beta_ref 
  t_dir = -sign_of(s0*ele%orientation)
  s_body = s0

! Particle is at an end.
else
  call offset_particle (ele, set$, orbit, set_hvkicks = .false., s_out = s_body, set_spin = set_spin)
  t_dir = sign_of(ele%value(l$)) * orbit%time_dir
endif

! Convert orbit coords to time based.

call convert_particle_coordinates_s_to_t (orbit, s_body, ele%orientation)

! Track through element

if ((ele%key == lcavity$ .or. ele%key == rfcavity$) .and. &
                  ele%field_calc == bmad_standard$ .and. ele%value(l$) < ele%value(l_active$)) then
  call out_io (s_error$, r_name, 'TIME-RUNGE-KUTTA TRACKING THROUGH RF CAVITY: ' // ele%name, &
                          'WILL NOT BE ACCURATE SINCE THE LENGTH IS LESS THAN THE HARD EDGE MODEL LENGTH.')
endif

call odeint_bmad_time(orbit, ele, param, t_dir, rf_time, err, track, t_end, dt_step)

if (err) return

!------
! Convert back to s-based coordinates.
! The particle is either dead or is alive and at an element end.

! For patch elements keeping track of the s-position is problematical. So set %s if particle is at the element end.
! Setting %s is not done for other types of elements since this could potentially cover up code problems.

if (orbit%location == upstream_end$) then
  orbit%p0c = ele%value(p0c_start$)
  call convert_particle_coordinates_t_to_s(orbit, ele, s_body)
  orbit%direction = -orbit%time_dir  ! In case t_to_s conversion confused by roundoff error.
  call offset_particle (ele, unset$, orbit, set_hvkicks = .false., set_spin = set_spin)
  if (ele%key == patch$) orbit%s = ele%s_start

elseif (orbit%location == downstream_end$) then
  orbit%p0c = ele%value(p0c$)
  call convert_particle_coordinates_t_to_s(orbit, ele, s_body)
  orbit%direction = orbit%time_dir  ! In case t_to_s conversion confused by roundoff error
  call offset_particle (ele, unset$, orbit, set_hvkicks = .false., set_spin = set_spin)
  if (ele%key == patch$) orbit%s = ele%s

elseif (orbit%state /= alive$) then
  ! Particle is lost in the interior of the element.
  ! The reference energy in the interior is equal to the ref energy at the exit end of the element.
  if (orbit%vec(6) < 0) then
    orbit%p0c = ele%value(p0c$)
    orbit%direction = -orbit%time_dir
  else 
    orbit%p0c = ele%value(p0c$)
    orbit%direction = orbit%time_dir
  end if

  call convert_particle_coordinates_t_to_s(orbit, ele, s_body)
  call offset_particle (ele, unset$, orbit, set_hvkicks = .false., s_pos = s_body, set_spin = set_spin)

elseif (present(t_end)) then
  orbit%p0c = ele%value(p0c$)
  call convert_particle_coordinates_t_to_s(orbit, ele, s_body)
  call offset_particle (ele, unset$, orbit, set_hvkicks = .false., s_pos = s_body, set_spin = set_spin)

else
  call out_io (s_fatal$, r_name, 'CONFUSED PARTICE LEAVING ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! Set relativistic beta

call convert_pc_to (orbit%p0c * (1 + orbit%vec(6)), orbit%species, beta = orbit%beta)

err_flag = .false.

end subroutine

