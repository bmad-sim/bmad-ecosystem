module time_tracker_mod

use beam_def_struct
use em_field_mod
use wall3d_mod
use geometry_mod
use runge_kutta_mod ! for common struct only
use multipass_mod

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine odeint_bmad_time (orb, ele, param, s1, s2, t_rel, local_ref_frame, err_flag, track)
! 
! Subroutine to do Runge Kutta tracking in time. This routine is adapted from Numerical
! Recipes.  See the NR book for more details.
!
! Input: 
!   orb   -- Coord_struct: Starting coords: (x, px, y, py, s, ps) [t-based]
!   ele     -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field_standard
!   param   -- lat_param_struct: Beam parameters.
!   s1      -- Real: Limit point to stop at.
!   s2      -- Real: Limit point to stop at.
!   t_rel   -- Real: The effective start time.
!   local_ref_frame 
!           -- Logical: If True then take the 
!                input and output coordinates as being with 
!                respect to the frame of referene of the element. 
!
!   track   -- Track_struct: Structure holding the track information.
!     %save_track -- Logical: Set True if track is to be saved.
!
! Output:
!   orb      -- Coord_struct: Ending coords: (x, px, y, py, s, ps) [t-based]
!   t_rel    -- Real: The effective end time.
!   err_flag -- Logical: Set True if there is an error. False otherwise.
!   track    -- Track_struct: Structure holding the track information.
!
!-

subroutine odeint_bmad_time (orb, ele, param, s1, s2, t_rel, local_ref_frame, err_flag, track)

use nr, only: zbrent

implicit none

type (coord_struct), intent(inout), target :: orb
type (coord_struct), target :: orb_old
type (coord_struct) :: orb_save
type (ele_struct), target :: ele
type (lat_param_struct), target ::  param
type (em_field_struct) :: saved_field
type (track_struct), optional :: track
type (fringe_edge_info_struct) fringe_info

real(rp), intent(in) :: s1, s2
real(rp), target :: t_rel, t_old, dt_tol
real(rp) :: dt, dt_did, dt_next, ds_safe, t_save, dt_save, s_save, dummy
real(rp), target  :: dvec_dt(9), vec_err(9), s_target, dt_next_save
real(rp) :: wall_d_radius, old_wall_d_radius = 0
real(rp) :: s_edge_track, ref_time, stop_time

integer, parameter :: max_step = 100000
integer :: n_step, n_pt, old_direction

logical, target :: local_ref_frame
logical :: at_edge_flag, exit_flag, err_flag, err, zbrent_needed, add_ds_safe, has_hit
logical :: edge_kick_applied, track_spin, stop_time_limited

character(30), parameter :: r_name = 'odeint_bmad_time'

! init
ds_safe = bmad_com%significant_length / 10
dt_next = bmad_com%init_ds_adaptive_tracking / c_light  ! Init time step.
call time_runge_kutta_periodic_kick_hook (orb, ele, param, stop_time, true_int$)

! local s coordinates for vec(5)
! Should not need to shift orb%s but, for example, an x_offset in a bend can confuse
! calc_next_fringe_edge.
! WRONG: orb%s = s1 + ele%s + ele%value(z_offset_tot$) - ele%value(l$)

orb%vec(5) = orb%s - (ele%s + ele%value(z_offset_tot$) - ele%value(l$))

call calc_next_fringe_edge (ele, orb%direction, s_edge_track, fringe_info, .true., orb)
old_direction = orb%direction

if ( present(track) ) then
   dt_save = track%ds_save/c_light
   t_save = t_rel
endif 

at_edge_flag = .false.
exit_flag = .false.
err_flag = .true.
has_hit = .false. 
edge_kick_applied  = .false.
  
do n_step = 1, max_step

  ! overstepped edge?
  if ((orb%vec(5) - s_edge_track)*orb%direction > -ds_safe) then
  
    zbrent_needed = .true.
    if ((orb%vec(5)-s_edge_track)*orb%direction < ds_safe) zbrent_needed = .false.

    add_ds_safe = .true.
    if (orb%direction == 1 .and. abs(s_edge_track - ele%value(l$)) < ds_safe) then
      orb%location = downstream_end$
      add_ds_safe = .false.
      exit_flag = .true.
    elseif (orb%direction == -1 .and. abs(s_edge_track) < ds_safe) then
      orb%location = upstream_end$
      add_ds_safe = .false.
      exit_flag = .true.
    endif

    ! zbrent

    dt_tol = ds_safe / (orb%beta * c_light)
    if (zbrent_needed) then
      dt = zbrent (delta_s_target, 0.0_rp, dt_did, dt_tol)
      dummy = delta_s_target(dt) ! Final call to set orb
      dt_did = dt
    endif

    ! Need to apply hard edge kick. 
    ! For super_slaves there may be multipole hard edges at a single s-position.
    edge_kick_applied = .false. 
    do 
      if (.not. associated(fringe_info%hard_ele)) exit
      if ((orb%vec(5)-s_edge_track)*orb%direction < -ds_safe) exit
      ! Get radius before first edge kick
      if (.not. edge_kick_applied) then 
        old_wall_d_radius = ref_frame_wall3d_d_radius (orb)
        edge_kick_applied = .true. 
      endif  
      if (orb%direction == +1) then 
        ref_time = fringe_info%hard_ele%value(ref_time_start$)
      else 
        ref_time = fringe_info%hard_ele%ref_time
      end if
      s_save = orb%vec(5)
      call convert_particle_coordinates_t_to_s(orb, ref_time) 
      track_spin = (ele%spin_tracking_method == tracking$ .and. ele%field_calc == bmad_standard$)
      call apply_element_edge_kick (orb, fringe_info, t_rel, ele, param, track_spin)
      call convert_particle_coordinates_s_to_t(orb)
      orb%vec(5) = s_save
      call calc_next_fringe_edge (ele, orb%direction, s_edge_track, fringe_info)
      ! Trying to take a step through a hard edge can drive Runge-Kutta nuts.
      ! So offset s a very tiny amount to avoid this
      if (add_ds_safe) then
        orb%vec(5) = orb%vec(5) + orb%direction * ds_safe
        orb%s = orb%s + orb%direction * ds_safe
      endif
    enddo      
    
  endif

  ! Wall check
  ! Adapted from runge_kutta_mod's odeint_bmad:
  ! Check if hit wall.
  ! If so, interpolate position particle at the hit point

  if (runge_kutta_com%check_wall_aperture) then
  
    wall_d_radius = ref_frame_wall3d_d_radius (orb)      
    select case (runge_kutta_com%hit_when)
    case (outside_wall$)
      has_hit = (wall_d_radius > 0)
    case (wall_transition$)
      has_hit = (wall_d_radius * old_wall_d_radius < 0 .and. n_step > 1)
      old_wall_d_radius = wall_d_radius
    case default
      call out_io (s_fatal$, r_name, 'BAD RUNGE_KUTTA_COM%HIT_WHEN SWITCH SETTING!')
      if (global_com%exit_on_error) call err_exit
    end select

    ! Cannot do anything if already hit
    if (has_hit .and. n_step == 1) then
      orb%state = lost$
      exit_flag = .true. 
    endif

    if (has_hit) then
      dt_tol = ds_safe / (orb%beta * c_light)
      ! Skip zbrent if the edge kick moved the particle outside the wall
      if (n_step /= 1  .and. .not. (edge_kick_applied .and. old_wall_d_radius < 0)) then
        dt = zbrent (wall_intersection_func, 0.0_rp, dt_did, dt_tol)
        dummy = wall_intersection_func(dt) ! Final call to set orb
      endif
      orb%state = lost$
      ! Convert for wall handler
      call convert_particle_coordinates_t_to_s(orb, ele%ref_time)
      call wall_hit_handler_custom (orb, ele, orb%s, orb%t)
      call convert_particle_coordinates_s_to_t(orb)
      ! Restore vec(5) to relative s 
      orb%vec(5) = orb%s - (ele%s + ele%value(z_offset_tot$) - ele%value(l$))
    endif
  endif

  if (orb%state /= alive$) exit_flag = .true.

  !Save track
  if ( present(track) ) then
    !Check if we are past a save time, or if exited
    if (t_rel >= t_save .or. exit_flag) then
      ! TODO: Set local_ref_frame=.true., and make sure offset_particle does the right thing
      call save_a_step (track, ele, param, .false., orb%vec(5), orb, s_save)
      ! Query the local field to save
      call em_field_calc (ele, param, orb%vec(5), t_rel, orb, local_ref_frame, saved_field, .false., err_flag)
      if (err_flag) return
      track%field(track%n_pt) = saved_field
      ! Set next save time 
      t_save = t_rel + dt_save
    end if
  endif

  ! Exit when the particle hits surface s1 or s2, or hits wall
  if (exit_flag) then
    err_flag = .false. 
    return
  endif

  ! Single Runge-Kutta step. Updates orb% vec(6), s, and t 

  stop_time_limited = .false.
  dt = dt_next

  if (stop_time /= real_garbage$ .and. dt > stop_time - orb%t) then
    if (stop_time < orb%t) then
      call out_io (s_error$, r_name, 'STOP_TIME FROM TIME_RUNGE_KUTTA_PERIODIC_KICK_HOOK IS IN THE PAST!')
      stop_time = real_garbage$
    else
      dt_next_save = dt_next
      dt_next = stop_time - orb%t
      dt = dt_next
      stop_time_limited = .true.
    endif
  endif

  orb_old = orb
  t_old = t_rel

  call rk_adaptive_time_step (ele, param, orb, t_rel, dt, dt_did, dt_next, local_ref_frame, err)

  if (stop_time_limited) then
    dt_next = dt_next_save
    if (abs(orb%t - stop_time) < bmad_com%significant_length / c_light) then
      call time_runge_kutta_periodic_kick_hook (orb, ele, param, stop_time, false_int$)
      if (orb%state /= alive$) return
    endif
  endif

  if (orb%direction /= old_direction) then
    call calc_next_fringe_edge (ele, orb%direction, s_edge_track, fringe_info)
    old_direction = orb%direction
  endif

end do

! Did not get to end

call out_io (s_warn$, r_name, 'STEPS EXCEEDED MAX_STEP FOR ELE: '//ele%name )
orb%location = inside$
orb%state = lost$



!------------------------------------------------------------------------------------------------
contains

! function for zbrent to calculate timestep to exit face surface

function delta_s_target (this_dt)

real(rp), intent(in)  :: this_dt
real(rp) :: delta_s_target
logical err_flag
!
call rk_time_step1 (ele, param, orb_old, t_old, this_dt, orb, vec_err, local_ref_frame, err_flag = err_flag)
delta_s_target = orb%vec(5) - s_edge_track
t_rel = t_old + this_dt
	
end function delta_s_target

!------------------------------------------------------------------------------------------------
! contains

function wall_intersection_func (this_dt) result (d_radius)

real(rp), intent(in) :: this_dt
real(rp) d_radius
logical err_flag
!
call rk_time_step1 (ele, param, orb_old, t_old, this_dt, &
	 				  orb, vec_err, local_ref_frame, err_flag = err_flag)
				  	 				  
d_radius = ref_frame_wall3d_d_radius (orb)
t_rel = t_old + this_dt

end function wall_intersection_func

! Wall sections are in the reference frame, not the local frame, so offset before check
function ref_frame_wall3d_d_radius (orb) result(ref_frame_d_radius)
type(coord_struct) :: orb, test_orb
real(rp) :: ref_frame_d_radius
! Make a copy
test_orb = orb
call offset_particle (ele, param, unset$, test_orb, ds_pos=test_orb%vec(5), set_hvkicks = .false., set_multipoles = .false.)
ref_frame_d_radius = wall3d_d_radius (test_orb%vec, ele)
end function ref_frame_wall3d_d_radius


end subroutine odeint_bmad_time

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine rk_adaptive_time_step (ele, param, orb, t, dt_try, dt_did, dt_next, local_ref_frame, err_flag)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, orb_new

real(rp), intent(inout) :: t
real(rp), intent(in)    :: dt_try
real(rp), intent(out)   :: dt_did, dt_next

real(rp) :: sqrt_n, err_max, dt, dt_temp, t_new, p2, rel_pc
real(rp) :: r_err(9), r_temp(9), dr_dt(9)
real(rp) :: r_scal(9), rel_tol, abs_tol
real(rp), parameter :: safety = 0.9_rp, p_grow = -0.2_rp
real(rp), parameter :: p_shrink = -0.25_rp, err_con = 1.89d-4
real(rp), parameter :: tiny = 1.0e-30_rp

logical local_ref_frame, err_flag
character(24), parameter :: r_name = 'rk_adaptive_time_step'

! Calc tolerances
! Note that s is in the element frame

call em_field_kick_vector_time (ele, param, t, orb, local_ref_frame, dr_dt, err_flag) 
if (err_flag) return

sqrt_N = sqrt(abs(1/(c_light*dt_try)))  ! number of steps we would take to cover 1 meter
rel_tol = bmad_com%rel_tol_adaptive_tracking / sqrt_N
abs_tol = bmad_com%abs_tol_adaptive_tracking / sqrt_N

!

dt = dt_try
orb_new = orb

do

  call rk_time_step1 (ele, param, orb, t, dt, orb_new, r_err, local_ref_frame, dr_dt, err_flag)
  ! Can get errors due to step size too large 
  if (err_flag) then
    if (dt < 1d-3/c_light) then
      call out_io (s_fatal$, r_name, 'CANNOT COMPLETE TIME STEP. ABORTING.')
      if (global_com%exit_on_error) call err_exit
      return
    endif
    dt_temp = dt / 10
  else
    ! r_scal(7:9) is for spin
    r_scal(:) = [abs(orb%vec) + abs(orb_new%vec), c_light, c_light, c_light]  + TINY
    r_scal(1:5:2) = r_scal(1:5:2) + [0.01_rp, 0.01_rp, ele%value(L$)]
    !Note that cp is in eV, so 1.0_rp is 1 eV
    r_scal(2:6:2) = r_scal(2:6:2) + 1.0_rp + 1d-4* (abs(orb%vec(2))+abs(orb%vec(4))+abs(orb%vec(6)))
    err_max = maxval(abs(r_err(:)/(r_scal(:)*rel_tol + abs_tol)))
    if (err_max <=  1.0) exit
    dt_temp = safety * dt * (err_max**p_shrink)
  endif
  dt = sign(max(abs(dt_temp), 0.1_rp*abs(dt)), dt)
  t_new = t + dt

  if (t_new == t) then
    err_flag = .true.
    call out_io (s_fatal$, r_name, 'STEPSIZE UNDERFLOW IN ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

end do

if (err_max > err_con) then
  dt_next = safety*dt*(err_max**p_grow)
else
  dt_next = 5.0_rp * dt
end if

! Increase step size, limited by an estimated next step ds = L/4

if (abs(dr_dt(5)*dt_next) > ele%value(L$)/4.0_rp) then
  dt_next = abs(ele%value(L$)/8.0_rp / dr_dt(5))
endif

! finish

dt_did = dt
t = t+dt

orb = orb_new

end subroutine rk_adaptive_time_step

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! Very similar to rk_step1_bmad, except that em_field_kick_vector_time is called
!  and orb_new%s and %t are updated to the global values

subroutine rk_time_step1 (ele, param, orb, t, dt, orb_new, r_err, local_ref_frame, dr_dt, err_flag)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, orb_new, orb_temp

real(rp), optional, intent(in) :: dr_dt(9)
real(rp), intent(in) :: t, dt
real(rp), intent(out) :: r_err(9)
real(rp) :: dr_dt1(9), dr_dt2(9), dr_dt3(9), dr_dt4(9), dr_dt5(9), dr_dt6(9), r_temp(9), pc
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

logical local_ref_frame, err_flag

!

if (present(dr_dt)) then
  dr_dt1 = dr_dt
else
  call em_field_kick_vector_time(ele, param, t, orb, local_ref_frame, dr_dt1, err_flag)
  if (err_flag) return
endif

orb_temp = orb

call transfer_this_orbit (orb_temp, orb, b21*dt*dr_dt1)
call em_field_kick_vector_time(ele, param, t+a2*dt, orb_temp, local_ref_frame, dr_dt2, err_flag)
if (err_flag) return

call transfer_this_orbit (orb_temp, orb, dt*(b31*dr_dt1+b32*dr_dt2))
call em_field_kick_vector_time(ele, param, t+a3*dt, orb_temp, local_ref_frame, dr_dt3, err_flag) 
if (err_flag) return

call transfer_this_orbit (orb_temp, orb, dt*(b41*dr_dt1+b42*dr_dt2+b43*dr_dt3))
call em_field_kick_vector_time(ele, param, t+a4*dt, orb_temp, local_ref_frame, dr_dt4, err_flag)
if (err_flag) return

call transfer_this_orbit (orb_temp, orb, dt*(b51*dr_dt1+b52*dr_dt2+b53*dr_dt3+b54*dr_dt4))
call em_field_kick_vector_time(ele, param, t+a5*dt, orb_temp, local_ref_frame, dr_dt5, err_flag)
if (err_flag) return

call transfer_this_orbit (orb_temp, orb, dt*(b61*dr_dt1+b62*dr_dt2+b63*dr_dt3+b64*dr_dt4+b65*dr_dt5))
call em_field_kick_vector_time(ele, param, t+a6*dt, orb_temp, local_ref_frame, dr_dt6, err_flag)
if (err_flag) return

! Output new orb and error vector

call transfer_this_orbit (orb_new, orb, dt*(c1*dr_dt1+c3*dr_dt3+c4*dr_dt4+c6*dr_dt6)) 

orb_new%t = orb%t + dt
orb_new%s = orb%s + orb_new%vec(5) - orb%vec(5)

if (orb_new%vec(6) > 0) then
  orb_new%direction = 1
else
  orb_new%direction = -1
endif

pc = sqrt(orb_new%vec(2)**2 +orb_new%vec(4)**2 + orb_new%vec(6)**2)
call convert_pc_to (pc, orb%species, beta = orb_new%beta)

r_err = dt*(dc1*dr_dt1+dc3*dr_dt3+dc4*dr_dt4+dc5*dr_dt5+dc6*dr_dt6)

!------------------------------------------------------------------------------------------------
contains

subroutine transfer_this_orbit (orb_out, orb_in, dvec)

type (coord_struct) orb_in, orb_out
real(rp) dvec(9), a_quat(4), omega(3), angle

!

orb_out%vec = orb_in%vec + dvec(1:6)

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  orb_out%spin = orb_in%spin
  call rotate_spinor(dvec(7:9), orb_out%spin)
endif

end subroutine transfer_this_orbit

end subroutine rk_time_step1

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!+
! Subroutine em_field_kick_vector_time (ele, param, t_rel, orbit, local_ref_frame, dvec_dt, err_flag)
!
! Subroutine to convert particle coordinates from t-based to s-based system. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele             -- coord_struct: input particle
!   param           -- real: Reference momentum. The sign indicates direction of p_s. 
!   t_rel           -- real: element coordinate system: t
!   orbit           -- coord_struct:
!                    %vec(1:6)  in t-based system
!   local_ref_frame --
!   err_flag        -- logical: Set True if there is an error. False otherwise.
! Output:
!    dvec_dt(9)  -- real(rp): Derivatives.
!-

subroutine em_field_kick_vector_time (ele, param, t_rel, orbit, local_ref_frame, dvec_dt, err_flag)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (em_field_struct) field

type (coord_struct), intent(in) :: orbit

real(rp), intent(in) :: t_rel    
real(rp), intent(out) :: dvec_dt(9)

real(rp) f_bend, kappa_x, kappa_y
real(rp) vel(3), force(3)
real(rp) :: pc, e_tot, mc2, gamma, charge, beta, p0, h

logical :: local_ref_frame, err_flag

character(28), parameter :: r_name = 'em_field_kick_vector_time'

! calculate the field. 
! Note that only orbit%vec(1) = x and orbit%vec(3) = y are used in em_field_calc,
!	and they coincide in both coordinate systems, so we can use the 'normal' routine:

call em_field_calc (ele, param, orbit%vec(5), t_rel, orbit, local_ref_frame, field, .false., err_flag)
if (err_flag) return

! Get e_tot from momentum
! velocities v_x, v_y, v_s:  c*[c*p_x, c*p_y, c*p_s]/e_tot

mc2 = mass_of(orbit%species) ! Note: mc2 is in eV
charge = charge_of(orbit%species) ! Note: charge is in units of |e_charge|

e_tot = sqrt( orbit%vec(2)**2 + orbit%vec(4)**2 + orbit%vec(6)**2 + mc2**2) 
vel(1:3) = c_light * [orbit%vec(2),  orbit%vec(4),  orbit%vec(6)] / e_tot 

! Computation for dr/dt where r(t) = [x, c*p_x, y, c*p_y, s, c*p_s]
! 
! p_x = m c \beta_x \gamma
! p_y = m c \beta_y \gamma
! p_s = m c h \beta_s \gamma 
!
! Note: v_s = (ds/dt) h, so ds/dt = v_s / h in the equations below
!
! h = 1 + \kappa_x * x + \kappa_y * y
!
! dx/dt   = v_x 
! dcp_x/dt = cp_s * v_s * \kappa_x / h + c*charge * ( Ex + v_y * Bs - v_s * By )
! dy/dt   = v_y
! dcp_y/dt = cp_s * v_s * \kappa_y / h + c*charge * ( Ey + v_s * Bx - v_x * Bs )
! ds/dt = v_s / h 
! dcp_s/dt = -(1/h) * cp_s * ( v_x * \kappa_x + v_y * \kappa_y ) + c*charge * ( Es + v_x By - v_y Bx )

! Straight coordinate systems have a simple Lorentz force

force = charge * (field%E + cross_product(vel, field%B))
dvec_dt(1) = vel(1)
dvec_dt(2) = c_light*force(1)
dvec_dt(3) = vel(2)
dvec_dt(4) = c_light*force(2)
dvec_dt(5) = vel(3)
dvec_dt(6) = c_light*force(3)

! Curvilinear coordinates have added terms

kappa_x = 0;  kappa_y = 0

if (ele%key == sbend$) then   
  if (ele%value(ref_tilt_tot$) /= 0 .and. .not. local_ref_frame) then
    kappa_x = ele%value(g$) * cos(ele%value(ref_tilt_tot$))
    kappa_y = ele%value(g$) * sin(ele%value(ref_tilt_tot$))
  else
    kappa_x = ele%value(g$)
  endif
  h = 1 + kappa_x * orbit%vec(1) + kappa_y * orbit%vec(3) ! h = 1 + kappa_x * x + kappa_y * y

  dvec_dt(2) = dvec_dt(2) + orbit%vec(6) * vel(3) * kappa_x / h
  dvec_dt(4) = dvec_dt(4) + orbit%vec(6) * vel(3) * kappa_y / h
  dvec_dt(5) = vel(3) / h
  dvec_dt(6) = dvec_dt(6) - orbit%vec(6) * (vel(1)*kappa_x + vel(2)*kappa_y) / h
endif

! Spin

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  dvec_dt(7:9) = spin_omega (field, orbit, .false.) + &
                      ele%orientation * [-kappa_y, kappa_x, 0.0_rp] * vel(3) / (1 + kappa_x * orbit%vec(1) + kappa_y * orbit%vec(3))
else
  dvec_dt(7:9) = 0
endif

end subroutine em_field_kick_vector_time
  
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function particle_in_global_frame (orb, in_time_coordinates, in_ele_frame, w_mat_out) result (particle) 
!
! Returns the particle in global time coordinates given is coordinates orb in lattice lat.
!   
!
! Module needed:
!   geometry_mod
!
! Input:
!   orb                 -- Coord_struct: particle in s-coordinates
!   branch              -- branch_struct: branch that contains branch%ele(orb%ix_ele)
!   in_time_coordinates -- Logical (optional): Default is false. If true, orb
!                            will taken as in time coordinates.    
!   in_ele_frame        -- Logical (optional): Default is true. If false, ele offsets
!                            will be ignored.
!
! Result:
!   particle            -- Coord_struct: particle in global time coordinates
!
!-

function particle_in_global_frame (orb, branch, in_time_coordinates, in_ele_frame, w_mat_out) result (particle)

implicit none

type (coord_struct) :: orb, particle
type (branch_struct) :: branch
type (floor_position_struct) :: floor_at_particle, global_position
type (ele_struct), pointer :: ele
real(rp) :: w_mat(3,3)
real(rp), optional :: w_mat_out(3,3)
logical, optional :: in_time_coordinates, in_ele_frame
character(28), parameter :: r_name = 'particle_in_global_frame'

!Get last tracked element  
ele =>  branch%ele(orb%ix_ele)

!Convert to time coordinates
particle = orb;
if (.not. logic_option( .false., in_time_coordinates)) then
  call convert_particle_coordinates_s_to_t (particle)
  ! Set vec(5) to be relative to entrance of ele 
  particle%vec(5) =  particle%vec(5) - (ele%s - ele%value(L$))
endif

!Set for coords_local_curvilinear_to_floor
floor_at_particle%r = particle%vec(1:5:2)
floor_at_particle%theta = 0.0_rp
floor_at_particle%phi = 0.0_rp
floor_at_particle%psi = 0.0_rp
! Get [X,Y,Z] and w_mat for momenta rotation below
global_position = coords_local_curvilinear_to_floor (floor_at_particle, ele, &
  in_ele_frame = logic_option(.true., in_ele_frame) , w_mat = w_mat)

!Set x, y, z
particle%vec(1:5:2) = global_position%r

!Rotate momenta 
particle%vec(2:6:2) = matmul(w_mat, particle%vec(2:6:2))

if (present(w_mat_out)) w_mat_out = w_mat

end function particle_in_global_frame

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine convert_particle_coordinates_t_to_s (particle, tref)
!
! Subroutine to convert particle coordinates from t-based to s-based system. 
!
! Modules needed:
!   use bmad
!
! Input:
!   particle   -- coord_struct: input particle coordinates
!   tref       -- real: reference time for z coordinate
!
! Output:
!   particle   -- coord_struct: output particle 
!-

subroutine convert_particle_coordinates_t_to_s (particle, tref)

implicit none

type (coord_struct), intent(inout), target ::particle
real(rp) :: p0c
real(rp), intent(in) :: tref
real(rp) :: pctot
real(rp), pointer :: vec(:)

!

vec => particle%vec
p0c = particle%p0c
pctot = sqrt (vec(2)**2 + vec(4)**2 + vec(6)**2)
if (vec(6) >= 0) then
  particle%direction = 1
else
  particle%direction = -1
endif

! Convert t to s. vec(1) and vec(3) are unchanged.

vec(2) = vec(2)/p0c
vec(4) = vec(4)/p0c
vec(5) = -c_light * particle%beta *  (particle%t - tref) 
vec(6) = pctot/p0c - 1.0_rp

end subroutine convert_particle_coordinates_t_to_s

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine convert_particle_coordinates_s_to_t (particle)
!
! Subroutine to convert particle coordinates from s-based to t-based system. 
!
! Note: t coordinates are:            
!     vec(1) = x                              [m]
!     vec(2) = c*p_x = m c^2 \gamma \beta_x   [eV]
!     vec(3) = y                              [m]
!     vec(4) = c*p_y = m c^2 \gamma beta_y    [eV]
!     vec(5) = s                              [m]
!     vec(6) = c*p_s = m c^2 \gamma \beta_s   [eV]
!
! Modules needed:
!   use bmad
!
! Input:
!   particle   -- coord_struct: input particle
!                       %vec(2), %vec(4), %vec(6)
!                       %s, %p0c
! Output:
!    particle   -- coord_struct: output particle 
!-

subroutine convert_particle_coordinates_s_to_t (particle)

implicit none

type (coord_struct), intent(inout), target :: particle
real(rp), pointer :: vec(:)

vec => particle%vec

! Convert s to t
vec(6) = particle%direction * particle%p0c * sqrt( ((1+vec(6)))**2 - vec(2)**2 -vec(4)**2 )
! vec(1) = vec(1) !this is unchanged
vec(2) = vec(2) * particle%p0c
! vec(3) = vec(3) !this is unchanged
vec(4) = vec(4) * particle%p0c
vec(5) = particle%s

end subroutine convert_particle_coordinates_s_to_t

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine drift_orbit_time(orbit, mc2, delta_s)
!
! Simple routine to drift a particle orbit in time-based coordinates by a distance delta_s
!   If the particle has zero longitudinal velocity, then the particle is not drifted
!   and a warning is printed.  
!
! Modules Needed:
!   use bmad_struct
!
! Input:
!   orbit      -- coord_struct: particle orbit in time-based coordinates
!   mc2        -- real(rp): particle mass in eV
!   delta_s    -- real(rp): s-coordinate distance to drift particle
!                  .
!
! Output:
!   orbit      -- coord_struct: particle orbit in time-based coordinates
!                                     
!-
subroutine drift_orbit_time(orbit, mc2, delta_s)
use bmad_struct
  
implicit none
  
type (coord_struct) :: orbit
real(rp) :: mc2, delta_s, delta_t, v_s, e_tot, vel(3)

character(28), parameter :: r_name = 'drift_orbit_time'
  
! Get e_tot from momentum

e_tot = sqrt( orbit%vec(2)**2 + orbit%vec(4)**2 +  orbit%vec(6)**2 + mc2**2) 

! velocities v_x, v_y, v_s:  c*[c*p_x, c*p_y, c*p_s]/e_tot

vel(1:3) = c_light*[  orbit%vec(2), orbit%vec(4), orbit%vec(6) ]/ e_tot 

if( vel(3) == 0 )then
   ! Do not drift
   call out_io (s_warn$, r_name, 'v_s == 0, will not drift')
   return
endif 
  
delta_t = delta_s / vel(3)

! Drift x, y, s
orbit%vec(1) = orbit%vec(1) + vel(1)*delta_t  !x
orbit%vec(3) = orbit%vec(3) + vel(2)*delta_t  !y
orbit%vec(5) = orbit%vec(5) + vel(3)*delta_t  !s
orbit%s =  orbit%s + delta_s
orbit%t =  orbit%t + delta_t 

end subroutine drift_orbit_time





!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine write_time_particle_distribution  (time_file_unit, bunch, style, branch, err)
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
!   style          -- character(16), optional: Style of output file:
!                            'BMAD' (default), 'OPAL', 'ASTRA', 'GPT'
!   branch         -- branch_struct, optional: Required for 'ASTRA' style
!
! Output:          
!   err            -- Logical, optional: Set True if, say a file could not be opened.
!-



subroutine write_time_particle_distribution (time_file_unit, bunch, style, branch, err)

implicit none

integer			    :: time_file_unit
type (bunch_struct) :: bunch
type (branch_struct), optional :: branch



type (coord_struct) :: orb, orb_ref
real(rp)        :: dt, pc, gmc, gammabeta(3), charge_alive

character(10)   ::  rfmt 
integer :: n_alive
integer :: i, i_style, a_species_id, a_status
integer, parameter :: bmad$ = 1, opal$ = 2, astra$ = 3, gpt$ = 4
logical, optional   :: err

character(*), optional  :: style 
character(16) :: style_names(4) = ['BMAD        ', &
								   'OPAL        ', &
								   'ASTRA       ', & 
								   'GPT         ']
character(40)	:: r_name = 'write_time_particle_distribution'

!

if (present(style)) then
  call match_word (style, style_names, i_style)
  if (i_style == 0) then
    call out_io (s_error$, r_name, 'Invalid style: '//trim(style))
  endif
else
  i_style = bmad$
endif

if (present(err)) err = .true.

!Format for numbers
  rfmt = 'es13.5'

! Number of alive particles
n_alive = count(bunch%particle(:)%state == alive$)

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
    do i = 1, 6
      orb_ref%vec(i) = sum( bunch%particle(:)%vec(i) *  bunch%particle(:)%charge, mask = (bunch%particle(:)%state == alive$)) / charge_alive
    enddo  
    ! For now just use the first particle as a reference. 
    orb = bunch%particle(1)
    orb_ref%t = branch%ele(bunch%ix_ele)%ref_time
    orb_ref%ix_ele = bunch%ix_ele
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
  
  !Get time to track backwards by
  dt = orb%t - bunch%t_center
  
  !Get pc before conversion
  pc = (1+orb%vec(6)) * orb%p0c 
  
  !convert to time coordinates
  call convert_particle_coordinates_s_to_t (orb)
  
  !get \gamma m c
  gmc = sqrt(pc**2 + mass_of(orb%species)**2) / c_light
  
  !'track' particles backwards in time and write to file
  ! (x, y, s) - dt mc2 \beta_x \gamma / \gamma m c
  orb%vec(1) = orb%vec(1) - dt*orb%vec(2)/gmc
  orb%vec(3) = orb%vec(3) - dt*orb%vec(2)/gmc
  orb%vec(5) = orb%vec(5) - dt*orb%vec(2)/gmc
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
										    orb%vec(3), gammabeta(2), &
											orb%vec(5), gammabeta(3)
  case (astra$)
    orb = particle_in_global_frame (orb,  branch, in_time_coordinates = .true.)
    a_species_id = astra_species_id(orb_ref%species)
    ! The reference particle is used for z, pz, and t
    write(time_file_unit, '(8'//rfmt//', 2i8)')  orb%vec(1), &
	    	 									 orb%vec(3), &
	    	 									 orb%vec(5) - orb_ref%vec(5), &
	    	 									 orb%vec(2), &
	    	 									 orb%vec(4), &
	    	 									 orb%vec(6) - orb_ref%vec(6), &
                                                 1e9_rp*(orb%t - orb_ref%t), &
                                                 1e9_rp*orb%charge, &
                                                 a_species_id, &
                                                 a_status
  case (gpt$)
    orb = particle_in_global_frame (orb,  branch, in_time_coordinates = .true.)
    gammabeta =  orb%vec(2:6:2) / mass_of(orb%species) 
    write(time_file_unit, '(9'//rfmt//')')  orb%vec(1), orb%vec(3), orb%vec(5), &
                                            gammabeta(1), gammabeta(2), gammabeta(3), &
										    orb%t, &
											charge_of(orb%species)*e_charge, &
											orb%charge/e_charge  
  end select

end do 

if (present(err)) err = .false.

contains

function astra_species_id(species) result (index)
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
end function

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
! Module needed:
!   use track1_mod
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
  
integer :: ix_ele, ix_offset, point_id, iteration
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
  
  s_rel = start2_orb%s - (ele%s - ele%value(L$) )
  if ( (s_rel <  -bmad_com%significant_length .or. s_rel > ele%value(L$) + bmad_com%significant_length)) then
	call out_io (s_fatal$, r_name, 'PARTICLE STARTED BEYOND ELEMENT BOUNDS FOR: ' // ele%name)
    print *, 's_rel: ', s_rel
    print *, 'start2_orb vec: ', start2_orb%vec
    print *, 'start2_orb%s: ', start2_orb%s
    print *, 'ele%s - L ', ele%s -ele%value(L$) 
    print *, 'ele%s', ele%s 
    if (global_com%exit_on_error) call err_exit
  endif
 
  if (verbose) print *, "  Tracking element:", ele%name, ele%ix_ele

  !Track particle until it hits something (beginning, end, wall)
  if (present(track)) then
    call track1( start2_orb, ele, lat%param, end_orb, track = track)
  else
    call track1( start2_orb, ele, lat%param, end_orb)
  end if

  !Determine if dead, or next element
  if (end_orb%state /= alive$) wall_flag = .true.

  if (verbose .and. present(track)) print *, "    Number of track points: ", track%n_pt

  if (verbose .and. wall_flag) print *, "  Particle lost in open area or &
    hit wall at element index: ", ele%name, ele%ix_ele

  
  ! get next index
  ix_offset = end_orb%direction
  
  ! Save current element
  ele_now => ele
  
  ! step to next ele
  ele => pointer_to_next_ele (ele, ix_offset)
  if (.not. associated(ele) ) exit
  
  ! Check for wrap around
  if (ix_offset == -1 .and. ele%ix_ele == 0) then
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
  if (end_orb%direction == +1) then
    ! particle arrives at the beginning of the element
    end_orb%s = ele%s - ele%value(L$) 
  else
    ! particle arrives at the beginning of the element
    end_orb%s = ele%s 
  endif

  ! Sanity check
  !if (ix_ele /= ele%ix_ele) then
  !	call out_io (s_fatal$, r_name, 'IX_ELE INCONSISTENCY IN ELE: ' // ele%name)
  !   if (global_com%exit_on_error) call err_exit
  ! endif
  

  !New start coords for ele
  start2_orb = end_orb  

end do track_loop


end subroutine track_until_dead

end module
