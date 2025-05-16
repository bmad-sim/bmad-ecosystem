module time_tracker_mod

use element_at_s_mod

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine odeint_bmad_time (orb, ele, param, t_dir, rf_time, err_flag, track, t_end, dt_step, extra_field)
! 
! Subroutine to do Runge Kutta tracking in time. This routine is adapted from Numerical
! Recipes.  See the NR book for more details.
!
! Tracking is done until the particle is lost or exits the element.
!
! Input: 
!   orb           -- coord_struct: Starting coords: (x, px, y, py, s, ps) [t-based]
!   ele           -- ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field_standard
!   param         -- lat_param_struct: Beam parameters.
!   t_dir         -- real(rp): Direction of time travel = +/-1. Can be negative for patches.
!                     Will be -1 if element has a negative length.
!
!   rf_time       -- real(rp): Time relative to RF clock.
!   track         -- track_struct, optional: Structure holding the track information.
!     %save_track   -- Logical: Set True if track is to be saved.
!   t_end         -- real(rp), optional: If present, maximum time to which the particle will be tracked.
!                     Used for tracking with given time steps. The time orb%t at which tracking stops 
!                     may be less than this if the particle gets to the end of the element
!   extra_field   -- em_field_struct, optional: Static field to be added to the element field. Eg used with space charge.
!
! Output:
!   orb           -- coord_struct: Ending coords: (x, px, y, py, s, ps) [t-based]
!   rf_time       -- real(rp): Updated time.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!   track         -- track_struct, optional: Structure holding the track information.
!   dt_step       -- real(rp), optional: Next RK time step that this tracker would take based on the error tolerance.
!                     Used by track_bunch_time.
!-

subroutine odeint_bmad_time (orb, ele, param, t_dir, rf_time, err_flag, track, t_end, dt_step, extra_field)

use super_recipes_mod, only: super_zbrent

implicit none

type (coord_struct), intent(inout), target :: orb
type (coord_struct), target :: old_orb
type (coord_struct) :: save_orb
type (ele_struct), target :: ele
type (lat_param_struct), target ::  param
type (em_field_struct) :: saved_field
type (track_struct), optional :: track
type (em_field_struct), optional :: extra_field
type (fringe_field_info_struct) fringe_info

real(rp), optional :: t_end, dt_step
real(rp), target :: old_t, dt_tol, s_fringe_edge, t0
real(rp) :: dt, dt_did, dt_next, ds_safe, t_save, dt_save, s_save, dummy, rf_time
real(rp), target  :: dvec_dt(10), vec_err(10), s_target, dt_next_save
real(rp) :: stop_time, s_stop_fwd, s_body_old, ds
real(rp), pointer :: s_body, s_fringe_ptr

integer :: t_dir, n_step, n_pt, old_direction, status

logical :: at_edge_flag, exit_flag, err_flag, err, zbrent_needed, add_ds_safe, has_hit
logical :: edge_kick_applied, track_spin, stop_time_limited

character(*), parameter :: r_name = 'odeint_bmad_time'

! init

s_body => orb%vec(5)
s_fringe_ptr => s_fringe_edge   ! To get around an intel bug: s_fringe_ptr is used in contained routine.

if (ele%key == patch$) then
  s_stop_fwd = 0  ! By convention.
elseif (orb%direction*orb%time_dir*ele%orientation == -1) then
  s_stop_fwd = 0
else
  s_stop_fwd = ele%value(l$)
endif

ds_safe = 0.01_rp * bmad_com%significant_length
dt_next = t_dir * bmad_com%init_ds_adaptive_tracking / c_light  ! Init time step.
if (ele%tracking_method == fixed_step_time_runge_kutta$) then
  if (ele%value(ds_step$) == 0) then
    call out_io (s_error$, r_name, 'FIXED_STEP_TIME_RUNGE_KUTTA TRACKING USED WITHOUT DS_STEP BEING SET!', &
                                   'WILL USE BMAD_COM%INIT_DS_ADAPTIVE_TRACKING AS A FALLBACK FOR ELEMENT: ' // ele%name)
  else
    ds = min(abs(ele%value(ds_step$)), abs(ele%value(l$)))
    dt_next = t_dir * ds / c_light
  endif
endif

stop_time = real_garbage$
if (associated(time_runge_kutta_periodic_kick_hook_ptr)) call time_runge_kutta_periodic_kick_hook_ptr (orb, ele, param, stop_time, true_int$)

call calc_next_fringe_edge (ele, s_fringe_edge, fringe_info, orb, .true., time_tracking = .true.)
old_direction = orb%direction
old_orb = orb

if (present(track)) then
  dt_save = t_dir * track%ds_save/c_light
  t_save = rf_time
endif 

at_edge_flag = .false.
exit_flag = .false.
err_flag = .true.
has_hit = .false. 
edge_kick_applied  = .false.
s_body_old = s_body

! Note: orb%direction * ele%orientation should be equal to sign_of(orb%vec(6)).

do n_step = 1, bmad_com%max_num_runge_kutta_step

  time_runge_kutta_com%num_steps_done = n_step

  ! edge?

  if ((s_body - s_fringe_edge) * sign_of(orb%vec(6)) * t_dir > -ds_safe .and. orb%vec(6) /= 0) then 
    zbrent_needed = .true.
    if (abs(s_body-s_fringe_edge) < ds_safe) zbrent_needed = .false.
    if (n_step == 1) zbrent_needed = .false.

    add_ds_safe = .true.
    if (orb%direction*orb%time_dir*ele%orientation == 1 .and. abs(s_fringe_edge - s_stop_fwd) < ds_safe) then
      if (ele%orientation == 1) then
        orb%location = downstream_end$
      else
        orb%location = upstream_end$
      endif
      add_ds_safe = .false.
      exit_flag = .true.
    elseif (orb%direction*orb%time_dir*ele%orientation == -1 .and. abs(s_fringe_edge) < ds_safe) then
      if (ele%orientation == 1) then
        orb%location = upstream_end$
      else
        orb%location = downstream_end$
      endif
      add_ds_safe = .false.
      exit_flag = .true.
    endif

    ! zbrent

    if (zbrent_needed) then
      dt_tol = ds_safe / (orb%beta * c_light)
      ! The reason for a non-zero t0 is to prevent the situation where a particle starts at the boundary,
      ! gets immediatly turned around, and then super_zbrent chooses the t = 0 solution. 
      ! This happened with an e_gun with particles pushed back to the cathode.
      t0 = 0
      if (old_orb%vec(6)*orb%vec(6) < 0) t0 = min(10.0_rp*dt_tol, 0.1_rp*dt_did)  ! If reflection has happened
      dt = super_zbrent (delta_s_target, t0, dt_did, 1e-15_rp, dt_tol, status)
      dummy = delta_s_target(dt, status) ! Final call to set orb
      dt_did = dt
    endif

    ! Need to apply hard edge kick. 
    ! For super_slaves there may be multipole hard edges at a single s-position.

    edge_kick_applied = .false. 
    do 
      if (.not. fringe_info%has_fringe .or. .not. associated(fringe_info%hard_ele)) exit
      if ((s_body-s_fringe_edge)*orb%time_dir*sign_of(orb%vec(6)) < -ds_safe) exit
      ! Get radius before first edge kick
      if (.not. edge_kick_applied) edge_kick_applied = .true.
      call convert_particle_coordinates_t_to_s(orb, ele, s_save)
      track_spin = (ele%spin_tracking_method == tracking$ .and. ele%field_calc == bmad_standard$)
      call apply_element_edge_kick (orb, fringe_info, ele, param, track_spin, rf_time = rf_time)
      call convert_particle_coordinates_s_to_t(orb, s_save, ele%orientation)
      if (orb%state /= alive$) return
      call calc_next_fringe_edge (ele, s_fringe_edge, fringe_info, orb, time_tracking = .true.)
      ! Trying to take a step through a hard edge can drive Runge-Kutta nuts.
      ! So offset s a very tiny amount to avoid this
      if (add_ds_safe) then
        s_body = s_body + orb%time_dir * sign_of(orb%vec(6)) * ds_safe
        orb%s = orb%s + orb%direction * orb%time_dir * ds_safe
        add_ds_safe = .false.
      endif
    enddo
  endif

  s_body_old = s_body

  ! Check if hit wall.
  ! If so, interpolate position particle at the hit point.
  ! Adapted from runge_kutta_mod's odeint_bmad.

  select case (ele%aperture_at)
  case (continuous$, wall_transition$)
    call check_aperture_limit (orb, ele, in_between$, param, old_orb, check_momentum = .false.)
    if (orb%state /= alive$) then
      ! Cannot do anything if this is the first step
      if (n_step == 1) return
      ! Skip zbrent if the edge kick moved the particle outside the wall
      dt_tol = ds_safe / (orb%beta * c_light)
      if (n_step /= 1  .and. .not. edge_kick_applied) then
        dt = super_zbrent (wall_intersection_func, 0.0_rp, dt_did, 1e-15_rp, dt_tol, status)
        ! Due to the zbrent finite tolerance, the particle may not have crossed the wall boundary.
        ! So step a small amount to make sure that the particle is past the wall.
        dummy = wall_intersection_func(dt+t_dir*ds_safe/c_light, status) ! Final call to set orb
      endif
      orb%location = inside$
      orb%state = lost$
      if (associated(wall_hit_handler_custom_ptr)) then
        ! Convert for wall handler
        call convert_particle_coordinates_t_to_s(orb, ele, s_save)
        call wall_hit_handler_custom_ptr (orb, ele, orb%s)
        call convert_particle_coordinates_s_to_t(orb, s_save, ele%orientation)
      endif
      if (orb%state /= alive$) exit_flag = .true.
    endif
  end select

  ! Check time

  if (present(t_end)) then
    if (abs(orb%t-t_end) < bmad_com%significant_length / c_light) then
      exit_flag = .true.
      orb%location = inside$
    endif
  endif

  ! Save track
  if (present(track) .and. (n_step > 1 .or. orb%state /= alive$)) then
    ! Check if we are past a save time, or if exited
    if (track%ds_save <= 0 .or. (rf_time - t_save) * t_dir >= 0 .or. exit_flag) then
      ! For consistency, convert to s-coordinates for save_a_step
      save_orb = orb
      call convert_particle_coordinates_t_to_s (save_orb, ele, s_body)
      call save_a_step (track, ele, param, .true., save_orb, s_body, .true., rf_time = rf_time)
      ! Set next save time 
      t_save = rf_time + t_dir*dt_save
    endif

    if (dt_did == dt) then
      track%n_ok = track%n_ok + 1
    else
      track%n_bad = track%n_bad + 1
    endif
  endif

  ! Exit when the particle hits an aperture or gets to the end of the element

  if (exit_flag .or. orb%state /= alive$) then
    err_flag = .false. 
    return
  endif

  ! Single Runge-Kutta step. Updates orb% vec(6), s, and t 

  stop_time_limited = .false.
  dt = dt_next
  if (present(t_end)) then
    if (present(dt_step) .and. t_dir*(dt + orb%t) < t_dir*t_end) dt_step = dt
    dt = t_dir * min(t_dir*dt, t_dir*t_end-orb%t)
  endif

  if (stop_time /= real_garbage$ .and. t_dir * dt > t_dir * (stop_time - orb%t)) then
    if (t_dir * stop_time < t_dir * orb%t) then
      call out_io (s_error$, r_name, 'STOP_TIME FROM TIME_RUNGE_KUTTA_PERIODIC_KICK_HOOK IS IN THE PAST!')
      stop_time = real_garbage$
    else
      dt_next_save = dt_next
      dt_next = stop_time - orb%t
      dt = dt_next
      stop_time_limited = .true.
    endif
  endif

  old_orb = orb
  old_t = rf_time

  call rk_adaptive_time_step (ele, param, orb, t_dir, rf_time, dt, dt_did, dt_next, err, extra_field)
  if (err) return
  if (ele%key == patch$) then
    s_stop_fwd = 0  ! By convention.
  elseif (orb%direction*orb%time_dir*ele%orientation == -1) then
    s_stop_fwd = 0
  else
    s_stop_fwd = ele%value(l$)
  endif

  edge_kick_applied = .false.

  if (stop_time_limited) then
    dt_next = dt_next_save
    if (associated(time_runge_kutta_periodic_kick_hook_ptr)) then
      if (abs(orb%t - stop_time) < 0.5_rp * bmad_com%significant_length / c_light) then
        call time_runge_kutta_periodic_kick_hook_ptr (orb, ele, param, stop_time, false_int$)
      endif
    endif
  endif

  if (orb%state /= alive$) then
    err_flag = .false.
    return
  endif

  if (orb%direction /= old_direction) then
    call calc_next_fringe_edge (ele, s_fringe_edge, fringe_info, orb, time_tracking = .true.)
    old_direction = orb%direction
  endif

end do

! Did not get to end.
! If this due to the particle moving too slowly this is not an error.

if (orb%beta < 1e-4) then
  if (time_runge_kutta_com%print_too_many_step_err)  then
    call out_io (s_error$, r_name, 'TOO MANY STEPS WHILE TRACKING THROUGH: ' // ele%name, &
                                   'DUE TO VELOCITY BEING TOO SMALL.', &
                                   'AT (X,Y,Z,T) POSITION FROM ENTRANCE: \4F12.7\ ', &
                                   r_array = [real(rp):: orb%vec(1), orb%vec(3), s_body, orb%t])
  endif
  err_flag = .false.

else
  call out_io (s_error$, r_name, 'STEP SIZE IS TOO SMALL OR TOO MANY STEPS WHILE TRACKING THROUGH: ' // ele%name, &
                                 'AT (X,Y,Z,T) POSITION FROM ENTRANCE: \4F12.7\ ', &
                                 'TYPICALLY THIS IS DUE TO THE FIELD NOT OBEYING MAXWELL''S EQUATIONS.', &
                                 '[OFTEN TIMES THE FIELD IS NOT EVEN BE CONTINUOUS IN THIS CASE!]', &
                                 'THE PARTICLE WILL BE MARKED AS LOST.', &
                                 r_array = [real(rp):: orb%vec(1), orb%vec(3), s_body, orb%t])
endif

orb%location = inside$
orb%state = lost$



!------------------------------------------------------------------------------------------------
contains

! function for zbrent to calculate timestep to exit face surface

function delta_s_target (this_dt, status)

real(rp), intent(in)  :: this_dt
real(rp) :: delta_s_target
integer status
logical err_flag
!
call rk_time_step1 (ele, param, old_t, old_orb, this_dt, orb, vec_err, err_flag = err_flag, extra_field = extra_field)
delta_s_target = s_body - s_fringe_ptr
rf_time = old_t + this_dt
  
end function delta_s_target

!------------------------------------------------------------------------------------------------
! contains

function wall_intersection_func (this_dt, status) result (d_radius)

type (coord_struct) test_orb
real(rp), intent(in) :: this_dt
real(rp) d_radius, s_here
integer status
logical err_flag, no_aperture_here

!
call rk_time_step1 (ele, param, old_t, old_orb, this_dt, orb, vec_err, err_flag = err_flag, extra_field = extra_field)

test_orb = orb
call convert_particle_coordinates_t_to_s(test_orb, ele, s_here)
call offset_particle (ele, unset$, test_orb, s_pos=s_here, set_hvkicks = .false.)
d_radius =  distance_to_aperture (test_orb, in_between$, ele, no_aperture_here)
rf_time = rf_time + this_dt

if (no_aperture_here) then
  call out_io (s_fatal$, r_name, 'CONFUSED APERTURE CALCULATION. PLEASE CONTACT HELP.')
  if (global_com%exit_on_error) call err_exit
endif

end function wall_intersection_func

end subroutine odeint_bmad_time

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine rk_adaptive_time_step (ele, param, orb, t_dir, rf_time, dt_try, dt_did, dt_next, err_flag, extra_field)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, new_orb
type (em_field_struct), optional :: extra_field

real(rp) :: rf_time, dt_try
real(rp) :: dt_did, dt_next

real(rp) :: sqrt_n, err_max, dt, dt_temp, t_new, p2, rel_pc
real(rp) :: r_err(10), r_temp(10), dr_dt(10), r_scaled_err(10), r_scaled_tot(10), abs_scale(10)
real(rp) :: r_scal(10), rel_tol, abs_tol, pc_ref
real(rp), parameter :: safety = 0.9_rp, p_grow = -0.2_rp
real(rp), parameter :: p_shrink = -0.25_rp, err_con = 1.89d-4

integer t_dir

logical err_flag
character(*), parameter :: r_name = 'rk_adaptive_time_step'

! Calc tolerances
! Note that s is in the element frame

call em_field_kick_vector_time (ele, param, rf_time, orb, dr_dt, err_flag, extra_field = extra_field) 
if (err_flag) return

sqrt_N = sqrt(abs(1/(c_light*dt_try)))  ! number of steps we would take to cover 1 meter
rel_tol = bmad_com%rel_tol_adaptive_tracking / sqrt_N
abs_tol = bmad_com%abs_tol_adaptive_tracking / sqrt_N

!

dt = dt_try
new_orb = orb

pc_ref = 0.5_rp * (ele%value(p0c_start$) + ele%value(p0c$))
abs_scale = [1d-2, 1d-6*pc_ref, 1d-2, 1d-6*pc_ref, 1d-2, 1d-2*pc_ref, 1.0_rp, 1.0_rp, 1.0_rp, 1d-4] 

do

  call rk_time_step1 (ele, param, rf_time,  orb, dt, new_orb, r_err, dr_dt, err_flag, .false., extra_field)
  ! Can get errors due to step size too large. For example, for a field map that is only slightly larger than
  ! the aperture, a particle that is outside the aperture and outside of the fieldmap will generate an error.
  ! The solution is to just take a smaller step.
  if (err_flag) then
    if (dt * t_dir < 1d-3/c_light) then
      ! Call rk_time_step1 to generate an error message.
      call rk_time_step1 (ele, param, rf_time,  orb, dt, new_orb, r_err, dr_dt, err_flag, .true., extra_field)
      call out_io (s_fatal$, r_name, 'CANNOT COMPLETE TIME STEP. MARKING PARTICLE AS LOST.')
      orb%state = lost$
      return
    endif
    dt_temp = 0.1_rp * dt
    if (ele%tracking_method == fixed_step_time_runge_kutta$) exit  ! No adaptive step sizing.

  elseif (ele%tracking_method == fixed_step_time_runge_kutta$) then
    exit  ! No adaptive step sizing.

  else
    ! r_scal(7:9) is for spin
    ! Note that cp is in eV, so 1.0_rp is 1 eV
    r_scal(:) = [0.5_rp*(abs(orb%vec) + abs(new_orb%vec)), 1.0_rp, 1.0_rp, 1.0_rp, 0.5_rp*(abs(orb%dt_ref) + abs(new_orb%dt_ref))]
    r_scal(2:6:2) = r_scal(2:6:2) + 1d-6 * (abs(orb%vec(2))+abs(orb%vec(4))+abs(orb%vec(6)))
    r_scaled_tot = r_scal(:) * rel_tol + abs_scale * abs_tol
    r_scaled_err = abs(r_err(:)/r_scaled_tot(:))
    err_max = maxval(r_scaled_err)
    if (err_max <=  1.0_rp) exit
    dt_temp = safety * dt * (err_max**p_shrink)
  endif

  dt = t_dir * max(abs(dt_temp), 0.1_rp*abs(dt))
  t_new = rf_time + dt

  if (t_new == rf_time) then ! Can only happen if dt is very small
    call out_io (s_error$, r_name, 'STEP SIZE UNDERFLOW IN ELEMENT: ' // ele%name, &
                                   'AT (X,Y,Z,T) POSITION FROM ENTRANCE: \4F12.7\ ', &
                                   'TYPICALLY THIS IS DUE TO THE FIELD NOT OBEYING MAXWELL''S EQUATIONS.', &
                                   'OFTEN TIMES THE FIELD IS NOT EVEN CONTINUOUS!', &
                                   'THE PARTICLE WILL BE MARKED AS LOST.', &
                                   r_array = [real(rp):: orb%vec(1), orb%vec(3), orb%vec(5), orb%t])
    err_flag = .true.
    orb%location = inside$
    orb%state = lost$
    return
  endif
end do

! Aperture sanity check

if (orbit_too_large(orb, check_momentum = .false.)) return

! Calculate next step

if (ele%tracking_method /= fixed_step_time_runge_kutta$) then
  if (err_max > err_con) then
    dt_next = safety*dt*(err_max**p_grow)
  else
    dt_next = 5.0_rp * dt
  end if

  if (ele%tracking_method /= fixed_step_time_runge_kutta$) then
    if (abs(dr_dt(5)*dt_next) > ele%value(L$)/4.0_rp) then
      ! Increase step size, limited by an estimated next step ds = L/4
      dt_next = t_dir * abs(ele%value(L$) / (8.0_rp * dr_dt(5)))
    endif
  endif
endif

! finish

dt_did = dt
rf_time = rf_time + dt

orb = new_orb

end subroutine rk_adaptive_time_step

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! Very similar to rk_step1_bmad, except that em_field_kick_vector_time is called
!  and new_orb%s and %t are updated to the global values

subroutine rk_time_step1 (ele, param, rf_time, orb, dt, new_orb, r_err, dr_dt, err_flag, print_err, extra_field)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, new_orb, temp_orb
type (em_field_struct), optional :: extra_field

real(rp), optional, intent(in) :: dr_dt(10)
real(rp), intent(in) :: rf_time, dt
real(rp), intent(out) :: r_err(10)
real(rp) :: dr_dt1(10), dr_dt2(10), dr_dt3(10), dr_dt4(10), dr_dt5(10), dr_dt6(10), r_temp(10), pc
real(rp), parameter :: a2=0.2_rp, a3=0.3_rp, a4=0.6_rp, &
    a5=1.0_rp, a6=0.875_rp, b21=0.2_rp, b31=3.0_rp/40.0_rp, &
    b32=9.0_rp/40.0_rp, b41=0.3_rp, b42=-0.9_rp, b43=1.2_rp, &
    b51=-11.0_rp/54.0_rp, b52=2.5_rp, b53=-70.0_rp/27.0_rp, &
    b54=35.0_rp/27.0_rp, &
    b61=1631.0_rp/55296.0_rp, b62=175.0_rp/512.0_rp, &
    b63=575.0_rp/13824.0_rp, b64=44275.0_rp/110592.0_rp, &
    b65=253.0_rp/4096.0_rp, c1=37.0_rp/378.0_rp, &
    c3=250.0_rp/621.0_rp, c4=125.0_rp/594.0_rp, &
    c6=512.0_rp/1771.0_rp, dc1=c1-2825.0_rp/27648.0_rp, &
    dc3=c3-18575.0_rp/48384.0_rp, dc4=c4-13525.0_rp/55296.0_rp, &
    dc5=-277.0_rp/14336.0_rp, dc6=c6-0.25_rp

real(rp) quat(0:3), dt_now
logical err_flag
logical, optional :: print_err

!

if (present(dr_dt)) then
  dr_dt1 = dr_dt
else
  call em_field_kick_vector_time(ele, param, rf_time, orb, dr_dt1, err_flag, print_err, extra_field)
  if (err_flag) return
endif

temp_orb = orb

dt_now = a2*dt
call transfer_this_orbit (temp_orb, dt_now, orb, b21*dt*dr_dt1)
call em_field_kick_vector_time(ele, param, rf_time + dt_now, temp_orb, dr_dt2, err_flag, print_err, extra_field)
if (err_flag) return

dt_now = a3*dt
call transfer_this_orbit (temp_orb, dt_now, orb, dt*(b31*dr_dt1+b32*dr_dt2))
call em_field_kick_vector_time(ele, param, rf_time + dt_now, temp_orb, dr_dt3, err_flag, print_err, extra_field)
if (err_flag) return

dt_now = a4*dt
call transfer_this_orbit (temp_orb, dt_now, orb, dt*(b41*dr_dt1+b42*dr_dt2+b43*dr_dt3))
call em_field_kick_vector_time(ele, param, rf_time + dt_now, temp_orb, dr_dt4, err_flag, print_err, extra_field)
if (err_flag) return

dt_now = a5*dt
call transfer_this_orbit (temp_orb, dt_now, orb, dt*(b51*dr_dt1+b52*dr_dt2+b53*dr_dt3+b54*dr_dt4))
call em_field_kick_vector_time(ele, param, rf_time + dt_now, temp_orb, dr_dt5, err_flag, print_err, extra_field)
if (err_flag) return

dt_now = a6*dt
call transfer_this_orbit (temp_orb, dt_now, orb, dt*(b61*dr_dt1+b62*dr_dt2+b63*dr_dt3+b64*dr_dt4+b65*dr_dt5))
call em_field_kick_vector_time(ele, param, rf_time + dt_now, temp_orb, dr_dt6, err_flag, print_err, extra_field)
if (err_flag) return

! Output new orb and error vector

call transfer_this_orbit (new_orb, dt, orb, dt*(c1*dr_dt1+c3*dr_dt3+c4*dr_dt4+c6*dr_dt6)) 

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  quat =          omega_to_quat(dt*c1*dr_dt1(7:9))
  quat = quat_mul(omega_to_quat(dt*c3*dr_dt3(7:9)), quat)
  quat = quat_mul(omega_to_quat(dt*c4*dr_dt4(7:9)), quat)
  quat = quat_mul(omega_to_quat(dt*c6*dr_dt6(7:9)), quat)
  new_orb%spin = quat_rotate(quat, orb%spin)
endif

new_orb%t = orb%t + dt
new_orb%s = orb%s + ele%orientation * (new_orb%vec(5) - orb%vec(5))

if (new_orb%vec(6) > 0) then
  new_orb%direction = ele%orientation
else
  new_orb%direction = -ele%orientation
endif

r_err = dt*(dc1*dr_dt1+dc3*dr_dt3+dc4*dr_dt4+dc5*dr_dt5+dc6*dr_dt6)

!------------------------------------------------------------------------------------------------
contains

subroutine transfer_this_orbit (out_orb, dt_now, in_orb, dvec)

type (coord_struct) in_orb, out_orb
real(rp) dvec(10), dt_now

!

out_orb%vec    = in_orb%vec + dvec(1:6)
out_orb%t      = in_orb%t + dt_now
out_orb%s      = in_orb%s + ele%orientation * (out_orb%vec(5) - in_orb%vec(5))
out_orb%dt_ref = in_orb%dt_ref + dt_now - dvec(10)

pc = sqrt(out_orb%vec(2)**2 + out_orb%vec(4)**2 + out_orb%vec(6)**2)
call convert_pc_to (pc, orb%species, beta = out_orb%beta)

end subroutine transfer_this_orbit

end subroutine rk_time_step1

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Subroutine em_field_kick_vector_time (ele, param, rf_time, orbit, dvec_dt, err_flag, print_err, extra_field))
!
! Subroutine to convert particle coordinates from t-based to s-based system. 
!
! Input:
!   ele             -- coord_struct: input particle
!   param           -- real(rp): Reference momentum. The sign indicates direction of p_s. 
!   rf_time         -- real(rp): RF time.
!   orbit           -- coord_struct: in t-based system
!   err_flag        -- logical: Set True if there is an error. False otherwise.
!   print_err       -- logical, optional: Passed to em_field_calc
!   extra_field     -- em_field_struct, optional: Static field to be added to the element field. Eg used with space charge.
!
! Output:
!    dvec_dt(10)  -- real(rp): Derivatives.
!-

subroutine em_field_kick_vector_time (ele, param, rf_time, orbit, dvec_dt, err_flag, print_err, extra_field)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) param
type (em_field_struct), optional :: extra_field
type (em_field_struct) field
type (coord_struct), intent(in) :: orbit

real(rp), intent(out) :: dvec_dt(10)
real(rp) rf_time, s_pos, s_tiny
real(rp) vel(3), force(3), momentum
real(rp) :: pc, e_tot, mc2, gamma, charge, beta, p0, h, beta0, dp_dt, dbeta_dt

integer ie

logical :: err_flag
logical, optional :: print_err

character(28), parameter :: r_name = 'em_field_kick_vector_time'

! Calculate the field. ...
!
! Note that only orbit%vec(1) = x and orbit%vec(3) = y are used in em_field_calc,
!  and they coincide in both coordinate systems, so we can use the 'normal' routine.
!
! Also the hard edge to the field for field_calc = bmad_standard can drive the integrator nuts.
! The solution is to pretend that the particle is inside the element. 
! This will not affect any results since the integrator will always appropriately interpolate to the edge.

s_pos = orbit%vec(5)
s_tiny = 0.01_rp * bmad_com%significant_length

if (ele%key == patch$) then
else
  if (s_pos < min(0.0_rp, ele%value(l$)) + s_tiny) s_pos = min(0.0_rp, ele%value(l$)) + s_tiny
  if (s_pos > max(0.0_rp, ele%value(l$)) - s_tiny) s_pos = max(0.0_rp, ele%value(l$)) - s_tiny
endif

call em_field_calc (ele, param, s_pos, orbit, .true., field, .false., err_flag, rf_time = rf_time, print_err = print_err)

if (err_flag) return

if (present(extra_field)) then
  field%E = field%E + extra_field%E
  field%B = field%B + extra_field%B
endif

! Get e_tot from momentum
! velocities v_x, v_y, v_s:  c*[c*p_x, c*p_y, c*p_s]/e_tot

mc2 = mass_of(orbit%species) ! Note: mc2 is in eV
charge = charge_of(orbit%species) ! Note: charge is in units of |e_charge|

momentum = sqrt(orbit%vec(2)**2 + orbit%vec(4)**2 + orbit%vec(6)**2)
e_tot = sqrt(momentum**2 + mc2**2) 
vel(1:3) = c_light * [orbit%vec(2),  orbit%vec(4),  orbit%vec(6)] / e_tot ! meters/sec

! Computation for dr/dt where r(t) = [x, c*p_x, y, c*p_y, s, c*p_s]
! 
! p_x = m c \beta_x \gamma
! p_y = m c \beta_y \gamma
! p_s = m c h \beta_s \gamma 
!
! Note: v_s = (ds/dt) h, so ds/dt = v_s / h in the equations below
!
! h = 1 + g * x
!
! dx/dt   = v_x 
! dcp_x/dt = c*charge * ( Ex + v_y * Bs - v_s * By ) + cp_s * v_s * g / h
! dy/dt   = v_y
! dcp_y/dt = c*charge * ( Ey + v_s * Bx - v_x * Bs )
! ds/dt = v_s / h 
! dcp_s/dt = -(1/h) * cp_s * v_x * g + c*charge * ( Es + v_x By - v_y Bx )

! Straight coordinate systems have a simple Lorentz force
! If ele%orientation = -1 then +s global direction = -s element frame direction

force = charge * (field%E + cross_product(vel, field%B))

dvec_dt(1) = vel(1)
dvec_dt(2) = c_light*force(1)
dvec_dt(3) = vel(2)
dvec_dt(4) = c_light*force(2)
dvec_dt(5) = vel(3)  ! s_body
dvec_dt(6) = c_light*force(3)

! Curvilinear coordinates have added terms

if (ele%key == sbend$ .or. ele%key == rf_bend$) then
  h = 1 + ele%value(g$) * orbit%vec(1)
  dvec_dt(2) = dvec_dt(2) + orbit%vec(6) * vel(3) * ele%value(g$) / h
  dvec_dt(5) = vel(3) / h
  dvec_dt(6) = dvec_dt(6) - orbit%vec(6) * vel(1) * ele%value(g$) / h
else
  h = 1.0_rp
endif

! The effective reference velocity is different from the velocity of the reference particle for wigglers where the 
! reference particle is not traveling along the reference line and in elements where the reference velocity is not constant.

if (orbit%beta == 0) then
  dvec_dt(10) = 0
else
  dp_dt = dot_product(force, vel) / (orbit%beta * c_light)
  dbeta_dt = mass_of(orbit%species)**2 * dp_dt * c_light / e_tot**3

  ! Note: There is a potential problem with RF and e_gun elements when calculating beta0 when there is slicing and where
  ! the particle is non-relativistic. To avoid z-shifts with slicing, use the lord delta_ref_time.
  ! There is a further potential problem in that a pipe with a lcavity superimposed on it will not have a 
  ! constant reference energy

  if ((ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) .and. ele%key /= pipe$) then
    do ie = 1, ele%n_lord
      ele0 => pointer_to_lord(ele, ie)
      if (ele0%key /= pipe$) exit
    enddo
  else
    ele0 => ele
  endif

  if (ele0%value(delta_ref_time$) == 0 .or. ele0%key == patch$) then
    beta0 = ele0%value(p0c$) / ele0%value(e_tot$) ! Singular case can happen when parsing lattice. 
  else
    beta0 = ele0%value(l$) / (c_light * ele0%value(delta_ref_time$))
  endif

  dvec_dt(10) = orbit%beta * ele%orientation * orbit%vec(6) / (h * beta0 * momentum)  ! dt_ref / dt
endif

! Spin

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  dvec_dt(7:9) = spin_omega (field, orbit, 0, .false.)
  if (ele%key == sbend$ .or. ele%key == rf_bend$) dvec_dt(7:9) = dvec_dt(7:9) + [0.0_rp, ele%value(g$), 0.0_rp] * vel(3) / h
else
  dvec_dt(7:9) = 0
endif

end subroutine em_field_kick_vector_time
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function particle_in_global_frame (orb, in_time_coordinates, in_body_frame, w_mat_out) result (particle) 
!
! Returns the particle in global time coordinates given is coordinates orb in lattice lat.
!   
! Input:
!   orb                 -- Coord_struct: particle in s-coordinates
!   branch              -- branch_struct: branch that contains branch%ele(orb%ix_ele)
!   in_time_coordinates -- Logical (optional): Default is false. If true, orb
!                            will taken as in time coordinates.    
!   in_body_frame        -- Logical (optional): Default is true. If false, ele offsets
!                            will be ignored.
!
! Result:
!   particle            -- Coord_struct: particle in global time coordinates
!-

function particle_in_global_frame (orb, branch, in_time_coordinates, in_body_frame, w_mat_out) result (particle)

implicit none

type (coord_struct) :: orb, particle
type (branch_struct) :: branch
type (floor_position_struct) :: floor_at_particle, global_position
type (ele_struct), pointer :: ele
real(rp) :: w_mat(3,3), s_body
real(rp), optional :: w_mat_out(3,3)
logical, optional :: in_time_coordinates, in_body_frame
character(28), parameter :: r_name = 'particle_in_global_frame'

! Get last tracked element  

ele =>  branch%ele(orb%ix_ele)
particle = orb;

! Convert to time coordinates

if (.not. logic_option(.false., in_time_coordinates)) then
  ! Set s_body to be relative to entrance of ele 
  s_body = particle%s - ele%s_start
  call convert_particle_coordinates_s_to_t (particle, s_body, ele%orientation)
endif

! Set for coords_local_curvilinear_to_floor

floor_at_particle%r = particle%vec(1:5:2)
floor_at_particle%theta = 0.0_rp
floor_at_particle%phi = 0.0_rp
floor_at_particle%psi = 0.0_rp

! Get [X,Y,Z] and w_mat for momenta rotation below
global_position = coords_local_curvilinear_to_floor (floor_at_particle, ele, &
  in_body_frame = logic_option(.true., in_body_frame) , w_mat = w_mat)

!Set x, y, z
particle%vec(1:5:2) = global_position%r

!Rotate momenta 
particle%vec(2:6:2) = matmul(w_mat, particle%vec(2:6:2))

if (present(w_mat_out)) w_mat_out = w_mat

end function particle_in_global_frame

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine drift_orbit_time(orbit, beta0, delta_s, delta_t)
!
! Simple routine to drift a particle orbit in time-based coordinates by a distance delta_s
!   or a time delta_t
!   If the particle has zero longitudinal velocity, then the particle is not drifted
!   and a warning is printed.  
!
! Input:
!   orbit      -- coord_struct: particle orbit in time-based coordinates.
!   beta0      -- real(rp): reference velocity v/c.
!   delta_s    -- real(rp), optional: s-coordinate distance to drift particle.
!   delta_t    -- real(rp), optional: -coordinate distancet to drift particle.
!
! Output:
!   orbit      -- coord_struct: particle orbit in time-based coordinates
!-

subroutine drift_orbit_time (orbit, beta0, delta_s, delta_t)

use bmad_struct
  
implicit none
  
type (coord_struct) :: orbit
real(rp), optional :: delta_s, delta_t
real(rp) :: v_s, e_tot, vel(3), beta0

character(*), parameter :: r_name = 'drift_orbit_time'
  
! Get e_tot from momentum

e_tot = sqrt(orbit%vec(2)**2 + orbit%vec(4)**2 +  orbit%vec(6)**2 + mass_of(orbit%species)**2) 

! velocities v_x, v_y, v_s:  c*[c*p_x, c*p_y, c*p_s]/e_tot

vel(1:3) = c_light * [orbit%vec(2), orbit%vec(4), orbit%vec(6)]/ e_tot 


if (present(delta_s)) then  
  if (vel(3) == 0) then
     ! Do not drift
     call out_io (s_warn$, r_name, 'v_s == 0, will not drift')
     return
  endif 
  delta_t = delta_s / vel(3)
else
  delta_s = vel(3) * delta_t
endif

! Drift x, y, s

orbit%vec(1) = orbit%vec(1) + vel(1)*delta_t  !x
orbit%vec(3) = orbit%vec(3) + vel(2)*delta_t  !y
orbit%vec(5) = orbit%vec(5) + delta_s  !s
orbit%s      = orbit%s + delta_s
orbit%t      = orbit%t + delta_t 
orbit%dt_ref = orbit%dt_ref + delta_t - c_light * delta_s / beta0

end subroutine drift_orbit_time

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_time_particle_distribution  (time_file_unit, bunch, ele, style, branch, format, err)
!
! Subroutine to write a time-based bunch from a standard Bmad bunch
! 
! Note: 'BMAD' style (absolute curvilinear coordinates): 
!       n_particles_alive 
!       x/m  m*c^2 \beta_x*\gamma/eV y/m m*c^2\beta_y*\gamma/eV s/m m*c^2\beta_z*\gamma/eV time/s charge/C
!      
!       'OPAL' style (absolute curvilinear coordinates): 
!       n_particles_alive
!       x/m  \beta_x*\gamma  y/m \beta_y*\gamma s/m \beta_s*\gamma
!
!       'ASTRA' style (global Cartesian coordinates, first line is the reference particle used for z, pz, and t calculation):
!       x/m y/m  z/m  m*c^2 \beta_x*\gamma/eV m*c^2 \beta_y*\gamma/eV m*c^2 \beta_z*\gamma/eV time/ns charge/nC species status
!       
!       'GPT' style (global Cartesian coordinates, with header labeling the columns)
!       x/m y/m z/m \beta_x*\gamma \beta_y*\gamma \beta_z*\gamma t/s elementary_charge/C charge/elementary_charge
!
! Input:
!   time_file_unit -- Integer: unit number to write to, if > 0
!   bunch          -- bunch_struct: bunch to be written.
!                            Particles are drifted to bmad_bunch%t_center for output
!   ele            -- ele_struct: Element being tracked through.
!   style          -- character(16), optional: Style of output file:
!                            'BMAD' (default), 'OPAL', 'ASTRA', 'GPT'
!   branch         -- branch_struct, optional: Required for 'ASTRA' style
!   format         -- character(10): format for numerical output. 
!                                    default: 'es15.7'
!
! Output:          
!   err            -- Logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_time_particle_distribution (time_file_unit, bunch, ele, style, branch, format, err)

implicit none

type (bunch_struct) :: bunch
type (ele_struct) ele
type (branch_struct), optional :: branch
type (coord_struct) :: orb, orb_ref

real(rp) :: dt, pc, gmc, gammabeta(3), charge_alive

character(10) :: rfmt 
integer :: time_file_unit
character(10), optional :: format
integer :: n_alive
integer :: i, i_style, a_species_id, a_status
integer, parameter :: bmad$ = 1, opal$ = 2, astra$ = 3, gpt$ = 4
logical, optional   :: err

character(*), optional  :: style
character(40)  :: r_name = 'write_time_particle_distribution'

!

if (present(style)) then
  call match_word (style, [character(8):: 'BMAD', 'OPAL', 'ASTRA', 'GPT'], i_style)
  if (i_style == 0) then
    call out_io (s_error$, r_name, 'Invalid style: '//trim(style))
  endif
else
  i_style = bmad$
endif

if (present(err)) err = .true.

rfmt = string_option('es15.7', format) ! Format for numbers

n_alive = count(bunch%particle(:)%state == alive$) ! Number of alive particles

! First line
select case (i_style)
  case (bmad$, opal$)
    ! Number of particles
    write(time_file_unit, '(i8)') n_alive !was: size(bunch%particle)  
  case (astra$)
    ! Reference particle is the average of all particles
    if (.not. present(branch)) call out_io (s_error$, r_name, 'Branch must be specified for ASTRA style')
    charge_alive = sum(bunch%particle(:)%charge, mask = (bunch%particle%state == alive$))
    if (charge_alive == 0) then
      call out_io (s_warn$, r_name, 'Zero alive charge in bunch, nothing written to file')
      return
    endif
    
    ! Averages
    do i = 1, 6
      orb_ref%vec(i) = sum( bunch%particle(:)%vec(i) *  bunch%particle(:)%charge, mask = (bunch%particle(:)%state == alive$)) / charge_alive
    enddo  
    
    ! Use the mean time as a reference 
    orb_ref%t = sum( bunch%particle(:)%t *  bunch%particle(:)%charge, mask = (bunch%particle(:)%state == alive$)) / charge_alive
    
    ! For now just use the first particle to copy information:
    orb = bunch%particle(1)
    orb_ref%s = orb%s
    orb_ref%ix_ele = orb%ix_ele
    orb_ref%p0c = orb%p0c
    orb_ref%species = orb%species
    a_species_id = astra_species_id(orb_ref%species)
    orb_ref = particle_in_global_frame (orb_ref,  branch)
    if (orb_ref%p0c == 0) then
      a_status = -1 ! Starting at cathode
    else 
      a_status = 5
    endif
    write(time_file_unit, '(8'//rfmt//', 2i8)') orb_ref%vec(1:5:2), orb_ref%vec(2:6:2), &
                     1e9_rp*orb_ref%t, 1e9_rp*orb_ref%charge, a_species_id, a_status
  case (gpt$)
    write(time_file_unit, '(9a13)') 'x', 'y', 'z', 'GBx', 'GBy', 'GBz', 't', 'q', 'nmacro'
                       
end select

! All particles
do i = 1, size(bunch%particle) 
  orb = bunch%particle(i)
  
  ! Only write live particles
  if (orb%state /= alive$) cycle
  
  ! Get time to track backwards by
  dt = orb%t - orb_ref%t

  ! Get pc before conversion
  pc = (1+orb%vec(6)) * orb%p0c 
  
  ! convert to time coordinates
  call convert_particle_coordinates_s_to_t (orb, orb%s-ele%s_start, ele%orientation)
  
  ! get \gamma m c
  gmc = sqrt(pc**2 + mass_of(orb%species)**2) / c_light
   
  !'track' particles backwards in time and write to file
  ! (x, y, s) - dt mc2 \beta_x \gamma / \gamma m c
  orb%vec(1) = orb%vec(1) - dt*orb%vec(2)/gmc
  orb%vec(3) = orb%vec(3) - dt*orb%vec(4)/gmc
  orb%vec(5) = orb%vec(5) - dt*orb%vec(6)/gmc
  orb%t = orb%t - dt
  
  ! 
  select case (i_style)
  case (bmad$) 
    write(time_file_unit, '(8'//rfmt//')')  orb%vec(1:6), bunch%t_center, bunch%particle(i)%charge 
  
  case (opal$)  
    gammabeta =  orb%vec(2:6:2) / mass_of(orb%species)
    ! OPAL has a problem with zero beta_s
    if ( gammabeta(3) == 0 ) gammabeta(3) = 1d-30 
    write(time_file_unit, '(6'//rfmt//')')  orb%vec(1), gammabeta(1), &
                        orb%vec(3), gammabeta(2), orb%vec(5), gammabeta(3)
  case (astra$)
    orb = particle_in_global_frame (orb,  branch, in_time_coordinates = .true.)
    a_species_id = astra_species_id(orb_ref%species)
    ! The reference particle is used for z, pz, and t
    write(time_file_unit, '(8'//rfmt//', 2i8)')  orb%vec(1), orb%vec(3), orb%vec(5) - orb_ref%vec(5), &
                            orb%vec(2), orb%vec(4), orb%vec(6) - orb_ref%vec(6), &
                            1e9_rp*(orb%t - orb_ref%t), 1e9_rp*orb%charge, a_species_id, a_status
  case (gpt$)
    orb = particle_in_global_frame (orb,  branch, in_time_coordinates = .true.)
    gammabeta =  orb%vec(2:6:2) / mass_of(orb%species) 
    write(time_file_unit, '(9'//rfmt//')')  orb%vec(1), orb%vec(3), orb%vec(5), &
                      gammabeta(1), gammabeta(2), gammabeta(3), &
                      orb%t, charge_of(orb%species)*e_charge, orb%charge/e_charge  
  end select

end do 

if (present(err)) err = .false.

!-----------------------------------

contains

function astra_species_id (species) result (index)
implicit none
integer :: species, index
select case (species)
      case (electron$)
        index = 1
      case (positron$)
        index = 2
      case default
        call out_io (s_warn$, r_name, 'Only electrons or positrons allowed for Astra. Setting index to -1')
        index = -1
end select
end function astra_species_id

end subroutine  write_time_particle_distribution

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! subroutine track_until_dead (start_orb, lat, end_orb, track)
!
! Subroutine to track a particle arbitrarily through a lattice, forwards or backwards,
!   until it is lost or exits the lattice.
!
!   The starting element is located using start_orb%s.
!
! Input:
!   start_orb   -- coord_struct: Starting coords.
!   lat         -- lat_stuct: lattice that contains and element at start_orb%s

! Output:
!   end_orb     -- coord_struct: final coords
!   track       -- track_struct: (optional) 
!
!-

subroutine track_until_dead (start_orb, lat, end_orb, track)

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (coord_struct) :: particle
type (track_struct), optional :: track
type (coord_struct) :: start_orb, start2_orb, end_orb
type (ele_struct), pointer :: ele, ele_now
type (ele_pointer_struct), allocatable :: chain_ele(:)

real(rp) s_rel
  
integer :: ix_ele, point_id, iteration
integer :: ix_pass_now, n_links, ix_pass
integer, parameter :: max_iteration = 100000
  
logical :: wall_flag, exit_flag, err, verbose

character(30), parameter :: r_name = 'track_until_dead'
  
!

start2_orb = start_orb

verbose = bmad_com%debug

branch => lat%branch(0)
ix_ele = element_at_s(branch%lat, start2_orb%s, .true., branch%ix_branch)
wall_flag = .false.
exit_flag = .false.
ele => lat%ele(ix_ele)
 
 
track_loop: do iteration = 1, max_iteration
  
  s_rel = start2_orb%s - ele%s_start
  if ( (s_rel <  -bmad_com%significant_length .or. s_rel > ele%value(L$) + bmad_com%significant_length)) then
  call out_io (s_fatal$, r_name, 'PARTICLE STARTED BEYOND ELEMENT BOUNDS FOR: ' // ele%name)
    print *, 's_rel: ', s_rel
    print *, 'start2_orb vec: ', start2_orb%vec
    print *, 'start2_orb%s: ', start2_orb%s
    print *, 'ele%s_start', ele%s_start
    print *, 'ele%s', ele%s 
    if (global_com%exit_on_error) call err_exit
  endif
 
  if (verbose) print *, "  Tracking element:", ele%name, ele%ix_ele

  !Track particle until it hits something (beginning, end, wall)
  if (present(track)) then
    call track1 (start2_orb, ele, lat%param, end_orb, track = track, init_to_edge = .false.)
  else
    call track1 (start2_orb, ele, lat%param, end_orb, init_to_edge = .false.)
  end if

  !Determine if dead, or next element
  if (end_orb%state /= alive$) wall_flag = .true.

  if (verbose .and. present(track)) print *, "    Number of track points: ", track%n_pt

  if (verbose .and. wall_flag) print *, "  Particle lost in open area or &
    hit wall at element index: ", ele%name, ele%ix_ele

  
  ! Save current element
  ele_now => ele
  
  ! step to next ele
  ele => pointer_to_next_ele (ele, end_orb%direction*end_orb%time_dir)
  if (.not. associated(ele) ) exit
  
  ! Check for wrap around
  if (end_orb%direction*end_orb%time_dir == -1 .and. ele%ix_ele == 0) then
    ! At beginning. Check that the lattice starts with a multipass element
    call multipass_chain (ele_now, ix_pass_now, n_links, chain_ele)
    call multipass_chain (ele, ix_pass, n_links)
    if (ix_pass_now == 1 .and. ix_pass == -1) then
      ele => pointer_to_next_ele(chain_ele(2)%ele, -1) 
    else
      ! No wrap around. Exit.
      exit
    endif
  endif
  
  ! Check for end
  if (ele%ix_ele == branch%n_ele_track) exit_flag = .true. 

  !Exit loop if 
  if (wall_flag .or. exit_flag) exit

  if (verbose) then
    print *, '  entering ele:     = ', ele%name, ele%ix_ele
    print *, '  end_orb%direction = ', end_orb%direction
    print *, '  end_orb%location  = ', end_orb%location
  endif

  !Place end_orb%s correctly 
  if (end_orb%direction*end_orb%time_dir == +1) then
    ! particle arrives at the beginning of the element
    end_orb%s        = ele%s_start
    end_orb%location = upstream_end$
    end_orb%ix_ele   = ele%ix_ele
  else
    ! particle arrives at the beginning of the element
    end_orb%s        = ele%s 
    end_orb%location = downstream_end$
    end_orb%ix_ele   = ele%ix_ele
  endif

  ! Sanity check
  !if (ix_ele /= ele%ix_ele) then
  !  call out_io (s_fatal$, r_name, 'IX_ELE INCONSISTENCY IN ELE: ' // ele%name)
  !   if (global_com%exit_on_error) call err_exit
  ! endif

  !New start coords for ele
  start2_orb = end_orb  

end do track_loop


end subroutine track_until_dead

end module
