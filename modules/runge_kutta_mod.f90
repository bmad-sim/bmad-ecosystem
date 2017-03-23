module runge_kutta_mod

use track1_mod
use em_field_mod
use fringe_edge_track_mod
use spin_mod

type runge_kutta_common_struct
  integer :: num_steps_done        ! Number of integration steps. Not used by Bmad. For external use.
end type

type (runge_kutta_common_struct), save :: runge_kutta_com

private :: rk_adaptive_step, rk_step1

contains

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine odeint_bmad (orbit, ele, param, s1, s2, local_ref_frame, err_flag, track)
! 
! Subroutine to do Runge Kutta tracking. This routine is adapted from Numerical
! Recipes.  See the NR book for more details.
!
! Notice that this routine has an two tolerances: 
!   bmad_com%rel_tol_adaptive_trackingrel_tol 
!   bmad_com%abs_tol_adaptive_trackingrel_tol 
! %rel_tol is scalled by the step size to to able to relate it to the final accuracy.
!
! Essentually (assuming random errors) one of these conditions holds:
!      %error in tracking < rel_tol
! or
!     absolute error in tracking < abs_tol
!
! Note: For elements where the reference energy is not constant (lcavity, etc.), and 
! with particles where the velocity is energy dependent (ie non ultra-relativistic), 
! the calculation of z is off since the reference velocity is unknown in the body of 
! the element. In this case, the reference velocity is simply taken to be a 
! constant equal to the reference velocity at the exit end. 
! It is up to the calling routine to handle the correction for this.
!
! Modules needed:
!   use bmad
!
! Input: 
!   orbit     -- Coord_struct: Starting coords: (x, px, y, py, z, delta).
!   ele       -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field_standard
!   param     -- lat_param_struct: Beam parameters.
!     %enegy       -- Energy in GeV
!     %particle    -- Particle type [positron$, or electron$]
!   s1        -- Real: Starting point relative to physical entrance.
!   s2        -- Real: Ending point relative physical entrance.
!   local_ref_frame 
!              -- Logical: If True then take the 
!                input and output coordinates as being with 
!                respect to the frame of referene of the element. 
!
! Output:
!   orbit    -- Coord_struct: Ending coords: (x, px, y, py, z, delta).
!   err_flag -- Logical: Set True if there is an error. False otherwise.
!   track    -- Track_struct, optional: Structure holding the track information.
!-

subroutine odeint_bmad (orbit, ele, param, s1, s2, local_ref_frame, err_flag, track)

use nr, only: zbrent

implicit none

type (coord_struct) :: orbit
type (coord_struct) old_orbit
type (ele_struct) ele
type (lat_param_struct) param
type (track_struct), optional :: track
type (fringe_edge_info_struct) fringe_info

real(rp), intent(in) :: s1, s2
real(rp) :: ds, ds_did, ds_next, s, s_last, ds_save
real(rp) :: s_edge_track, position(6)
real(rp) :: old_s, dist_to_wall, ds_tiny

integer :: n_step, s_dir, nr_max, this_state

logical local_ref_frame, err_flag, err, at_hard_edge, track_spin

character(*), parameter :: r_name = 'odeint_bmad'

!

err_flag = .true.
s = s1
s_dir = sign(1.0_rp, s2-s1)
ds_next = bmad_com%init_ds_adaptive_tracking * s_dir
ds_tiny  = bmad_com%significant_length/100
track_spin = (ele%spin_tracking_method == tracking$ .and. ele%field_calc == bmad_standard$)

! Should not need to shift orb%s but, for example, an x_offset in a bend can confuse
! calc_next_fringe_edge.

orbit%s = s1 + ele%s_start + ele%value(z_offset_tot$)

! For elements where the reference energy is changing the reference energy in the body is 
! taken, by convention, to be the reference energy at the exit end.

call reference_energy_correction (ele, orbit, first_track_edge$)

! If the element is using a hard edge model then need to stop at the hard edges
! to apply the appropriate hard edge kick.
! calc_next_fringe_edge assumes that s = 0 is beginning of element which is not true of a patch element.

call calc_next_fringe_edge (ele, s_edge_track, fringe_info, orbit, .true.)
if (ele%key == patch$) s_edge_track = s2

! Save initial point

if (present(track)) then
  s_last = s - 2.0_rp * track%ds_save
  if ((abs(s-s_last) > track%ds_save)) then
    call save_a_step (track, ele, param, local_ref_frame, orbit, s, .true.)
    s_last = s
  endif
endif

! now track

err = .false.

do n_step = 1, bmad_com%max_num_runge_kutta_step

  runge_kutta_com%num_steps_done = n_step
  
  ! Check if we we need to apply a hard edge kick.
  ! For super_slaves there may be multiple hard edges at a single s-position.

  do
    if (.not. associated(fringe_info%hard_ele)) exit
    if ((s-s_edge_track)*s_dir < -ds_tiny) exit
    call apply_element_edge_kick (orbit, fringe_info, ele, param, track_spin)
    call calc_next_fringe_edge (ele, s_edge_track, fringe_info, orbit)
    ! Trying to take a step through a hard edge can drive Runge-Kutta nuts.
    ! So offset s a very tiny amount to avoid this
    s = s + ds_tiny * s_dir
  enddo

  ! Check if we are done.

  if ((s-s2)*s_dir > -ds_tiny) then
    if (present(track)) call save_a_step (track, ele, param, local_ref_frame, orbit, s, .true.)
    call reference_energy_correction (ele, orbit, second_track_edge$)
    err_flag = .false.
    return
  end if

  ! Need to propagate a step. First calc tolerances.

  ds = ds_next
  at_hard_edge = .false.

  if ((s+ds-s_edge_track)*s_dir > 0.0) then
    at_hard_edge = .true.
    ds_save = ds
    ds = s_edge_track - s - ds_tiny*s_dir / 2
  endif

  old_orbit = orbit
  old_s = s
  call rk_adaptive_step (ele, param, orbit, s, ds, abs(s2-s1), ds_did, ds_next, local_ref_frame, err)
  if (err) return

  ! Check x/y limit apertures

  select case (ele%aperture_at)
  case (continuous$, wall_transition$)
    call check_aperture_limit (orbit, ele, in_between$, param, old_orbit)
    if (orbit%state /= alive$) then
      if (n_step == 1) return  ! Cannot do anything if this is the first step
      this_state = orbit%state ! zbrent will overwrite this so save.
      dist_to_wall = zbrent (wall_intersection_func, 0.0_rp, ds_did, ds_tiny)
      orbit%state = this_state
      call wall_hit_handler_custom (orbit, ele, s)
      if (orbit%state /= alive$) return
      if (ele%aperture_at /= wall_transition$) then
        call out_io (s_error$, r_name, 'CUSTOM CODE IS KEEPING A PARTICLE ALIVE ACCROSS A BOUNDARY!', &
                                       'IN THIS CASE, THE APERTURE_AT COMPONENT OF ELEMENT: ' // ele%name, &
                                       'NEEDS TO BE SET TO "WALL_TRANSITION".')
        if (global_com%exit_on_error) call err_exit
        return
      endif
      ! Due to the zbrent finite tolerance, the particle may not have crossed the wall boundary.
      ! So step a small amount to make sure that the particle is past the wall.
      dist_to_wall = wall_intersection_func(ds_did+ds_tiny)
    endif
  end select

  ! Save track

  if (present(track)) then
    if ((abs(s-s_last) > track%ds_save)) then
      call save_a_step (track, ele, param, local_ref_frame, orbit, s, .true.)
      s_last = s
    endif

    if (ds_did == ds) then
      track%n_ok = track%n_ok + 1
    else
      track%n_bad = track%n_bad + 1
    end if
  endif

  ! Calculate next step size. If there was a hard edge then take into account the step that would have
  ! been taken if no hard edge was present.

  if (at_hard_edge .and. abs(ds_next) >= abs(ds)) then
    ds_next = max(abs(ds_save), abs(ds_next)) * s_dir
  endif

  if ((s + ds_next - s2) * s_dir > 0) then
    ds_next = s2 - s
    at_hard_edge = .true.
  endif

  if (abs(ds_next) <  bmad_com%min_ds_adaptive_tracking) ds_next = s_dir * bmad_com%min_ds_adaptive_tracking

  ! Check for step size smaller than minimum. If so we consider the particle lost
 
  if (.not. at_hard_edge .and. abs(ds) < bmad_com%fatal_ds_adaptive_tracking) exit

end do

! Here if step size too small or too many steps.
! One possibility is that the particle is turning around longitudinally.
! A particle is considered to be turning around if the z-velocity is less than 10% of the total velocity.
! Only issue an error message if the particle is *not* turning around since, in this case, there might be an
! error in how the field is calculated and we must warn the user of this.

if (sqrt(orbit%vec(2)**2 + orbit%vec(4)**2) / (1 + orbit%vec(6)) > 0.9) then
  orbit%state = lost_z_aperture$
else
  call out_io (s_error$, r_name, 'STEP SIZE IS TOO SMALL OR TOO MANY STEPS WHILE TRACKING THROUGH: ' // ele%name, &
                                 'AT S-POSITION FROM ENTRANCE: \F10.5\ ', &
                                 'COULD BE DUE TO A DISCONTINUITY IN THE FIELD ', &
                                 r_array = [s])
  orbit%state = lost$
endif

err_flag = .false.

!-----------------------------------------------------------------
contains

function wall_intersection_func (ds) result (d_radius)

real(rp), intent(in) :: ds
real(rp) d_radius
real(rp) t_new, r_err(11), dr_ds(10)

!

call rk_step1 (ele, param, old_orbit, dr_ds, old_s, ds, orbit, r_err, local_ref_frame, err_flag)

s = old_s + ds
orbit%s = s + ele%s_start
call check_aperture_limit (orbit, ele, in_between$, param, old_orbit)

d_radius = param%unstable_factor

end function wall_intersection_func

end subroutine odeint_bmad

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!+
! Subroutine rk_adaptive_step (ele, param, orb, s, ds_try, ds_did, ds_next, local_ref_frame, err_flag)
!
! Private routine used by odeint_bmad.
! Not meant for general use
!-

subroutine rk_adaptive_step (ele, param, orb, s, ds_try, ds12, ds_did, ds_next, local_ref_frame, err_flag)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, orb_new

real(rp), intent(inout) :: s
real(rp), intent(in)    :: ds_try
real(rp), intent(out)   :: ds_did, ds_next

real(rp) :: dr_ds(10), r_scal(10), rel_tol, abs_tol
real(rp) :: err_max, ds, ds_temp, s_new, p2
real(rp) :: r_err(10), r_temp(10), pol
real(rp) :: sqrt_N, rel_pc, t_new, ds12
real(rp), parameter :: safety = 0.9_rp, p_grow = -0.2_rp
real(rp), parameter :: p_shrink = -0.25_rp, err_con = 1.89d-4
real(rp), parameter :: tiny = 1.0e-30_rp

logical local_ref_frame, err_flag
character(20), parameter :: r_name = 'rk_adaptive_step'

! Init. Note:
!   rel_tol: Same as eps for odeint scalled by sqrt(ds/(s2-s1))
!                   where ds is the step size for the current step. rel_tol
!                   sets the %error of the result
!   abs_tol: Sets the absolute error of the result

ds = ds_try
orb_new = orb

do
  call rk_step1 (ele, param, orb, dr_ds, s, ds, orb_new, r_err, local_ref_frame, err_flag)
  ! Can get errors due to step size too large 
  if (err_flag) then
    if (ds < 1d-3) then
      call out_io (s_error$, r_name, 'Problem with field calc. Tracked particle will be marked as dead.')
      orb%state = lost$
     return
    endif
    ds_temp = ds / 10

  else
    sqrt_N = sqrt(abs(ds12/ds))  ! number of steps we would take with this ds
    rel_tol = bmad_com%rel_tol_adaptive_tracking / sqrt_N
    abs_tol = bmad_com%abs_tol_adaptive_tracking / sqrt_N
    pol = 1  ! Spin scale
    r_scal(:) = abs([orb%vec(:), orb%t, pol, pol, pol]) + abs(ds*dr_ds(:)) + TINY
    err_max = maxval(abs(r_err(:)/(r_scal(:)*rel_tol + abs_tol)))
    if (err_max <=  1.0) exit
    ds_temp = safety * ds * (err_max**p_shrink)
  endif
  ds = sign(max(abs(ds_temp), 0.1_rp*abs(ds)), ds)
  s_new = s + ds

  if (s_new == s) then
    err_flag = .true.
    call out_io (s_fatal$, r_name, 'STEPSIZE UNDERFLOW IN ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
end do

!

if (err_max > err_con) then
  ds_next = safety*ds*(err_max**p_grow)
else
  ds_next = 5.0_rp * ds
end if

ds_did = ds
s = s+ds

orb_new%s = orb%s + ds

orb = orb_new
err_flag = .false.

end subroutine rk_adaptive_step

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine rk_step1 (ele, param, orb, dr_ds1, s, ds, orb_new, r_err, local_ref_frame, err)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, orb_new, orb_temp(5)

real(rp) :: dr_ds1(10)
real(rp), intent(in) :: s, ds
real(rp), intent(out) :: r_err(10)
real(rp) :: dr_ds2(10), dr_ds3(10), dr_ds4(10), dr_ds5(10), dr_ds6(10)
real(rp), parameter :: a2=0.2_rp, a3=0.3_rp, a4=0.6_rp, &
    a5=1.0_rp, a6=0.875_rp, b21=0.2_rp, b31=3.0_rp/40.0_rp, &
    b32=9.0_rp/40.0_rp, b41=0.3_rp, b42=-0.9_rp, b43=1.2_rp, &
    b51=-11.0_rp/54.0_rp, b52=2.5_rp, b53=-70.0_rp/27.0_rp, &
    b54=35.0_rp/27.0_rp, b61=1631.0_rp/55296.0_rp, b62=175.0_rp/512.0_rp, &
    b63=575.0_rp/13824.0_rp, b64=44275.0_rp/110592.0_rp, &
    b65=253.0_rp/4096.0_rp, c1=37.0_rp/378.0_rp, &
    c3=250.0_rp/621.0_rp, c4=125.0_rp/594.0_rp, &
    c6=512.0_rp/1771.0_rp, dc1=c1-2825.0_rp/27648.0_rp, &
    dc3=c3-18575.0_rp/48384.0_rp, dc4=c4-13525.0_rp/55296.0_rp, &
    dc5=-277.0_rp/14336.0_rp, dc6=c6-0.25_rp

real(rp) quat(0:3)
logical local_ref_frame, err

!

call kick_vector_calc (ele, param, s, orb, local_ref_frame, dr_ds1, err)
if (err) return

!

call transfer_this_orbit (orb, b21*ds*dr_ds1, orb_temp(1))
call kick_vector_calc(ele, param, s + a2*ds, orb_temp(1), local_ref_frame, dr_ds2, err)
if (err) return

call transfer_this_orbit (orb, ds*(b31*dr_ds1 + b32*dr_ds2), orb_temp(2))
call kick_vector_calc(ele, param, s + a3*ds, orb_temp(2), local_ref_frame, dr_ds3, err)
if (err) return

call transfer_this_orbit (orb, ds*(b41*dr_ds1 + b42*dr_ds2 + b43*dr_ds3), orb_temp(3))
call kick_vector_calc(ele, param, s + a4*ds, orb_temp(3), local_ref_frame, dr_ds4, err)
if (err) return

call transfer_this_orbit (orb, ds*(b51*dr_ds1 + b52*dr_ds2 + b53*dr_ds3 + b54*dr_ds4), orb_temp(4))
call kick_vector_calc(ele, param, s + a5*ds, orb_temp(4), local_ref_frame, dr_ds5, err)
if (err) return

call transfer_this_orbit (orb, ds*(b61*dr_ds1 + b62*dr_ds2 + b63*dr_ds3 + b64*dr_ds4 + b65*dr_ds5), orb_temp(5))
call kick_vector_calc(ele, param, s + a6*ds, orb_temp(5), local_ref_frame, dr_ds6, err)
if (err) return

call transfer_this_orbit (orb, ds*(c1*dr_ds1 + c3*dr_ds3 + c4*dr_ds4 + c6*dr_ds6), orb_new)

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  quat =          omega_to_quat(ds*c1*dr_ds1(8:10))
  quat = quat_mul(omega_to_quat(ds*c3*dr_ds3(8:10)), quat)
  quat = quat_mul(omega_to_quat(ds*c4*dr_ds4(8:10)), quat)
  quat = quat_mul(omega_to_quat(ds*c6*dr_ds6(8:10)), quat)
  orb_new%spin = rotate_vec_given_quat(orb%spin, quat)
endif

r_err=ds*(dc1*dr_ds1 + dc3*dr_ds3 + dc4*dr_ds4 + dc5*dr_ds5 + dc6*dr_ds6)

!----------------------------------------------------------
contains

subroutine transfer_this_orbit (orb_in, dvec, orb_out)

type (coord_struct) orb_in, orb_out
real(rp) dvec(10), t_temp, a_quat(4), omega(3), angle

!

orb_out = orb_in
orb_out%vec = orb_in%vec + dvec(1:6)
orb_out%t = orb_in%t + dvec(7)

if (dvec(6) /= 0) then
  call convert_pc_to (orb_out%p0c * (1 + orb_out%vec(6)), param%particle, beta = orb_out%beta)
endif

end subroutine transfer_this_orbit

end subroutine rk_step1

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine kick_vector_calc (ele, param, s_rel, orbit, local_ref_frame, dr_ds, field, err)
!
! Subroutine to calculate the dr/ds "kick vector" where
!     r = [x, p_x, y, p_y, z, p_z, t, spin1_re, spin1_im, spin2_re, spin2_im] 
!
! Remember: In order to simplify the calculation, in the body of any element, P0 is taken to be 
! the P0 at the exit end of the element.
!
!   dr(1)/ds = dx/dt * dt/ds
!   where:
!     dx/dt = v_x = p_x / (1 + p_z)
!     dt/ds = (1 + g*x) / v_s
!     g = 1/rho, rho = bending radius (nonzero only in a dipole)
!
!   dr(2)/ds = dP_x/dt * dt/ds / P0 + g_x * P_z
!   where:
!     dP_x/dt = EM_Force_x
!     g_x = bending in x-plane.
!
!   dr(3)/ds = dy/dt * dt/ds
!   where:
!     dy/dt = v_x 
! 
!   dr(4)/ds = dP_y/dt * ds/dt / P0 + g_y * P_z
!   where:
!     dP_y/dt = EM_Force_y
!     g_y = bending in y-plane.
!
!   dr(5)/ds = beta * c_light * [dt/ds(ref) - dt/ds] + dbeta/ds * c_light * [t(ref) - t]
!            = beta * c_light * [dt/ds(ref) - dt/ds] + dbeta/ds * vec(5) / beta
!   where:
!     dt/ds(ref) = 1 / beta(ref)
!     Note: dt/ds(ref) formula is inaccurate at low energy in an lcavity.
!
!   dr(6)/ds = (EM_Force dot v_hat) * dt/ds / P0
!   where:
!      v_hat = velocity normalized to 1.
!
!   dr(7)/ds = dt/ds
!
!   dr(8:10)/ds = Spin omega vector
!
! Modules needed:
!   use bmad
!
! Input:
!   ele   -- Ele_struct: Element being tracked thorugh.
!   param -- lat_param_struct: Lattice parameters.
!   s_rel -- Real(rp): Distance from the start of the element to the particle.
!   orbit -- coord_struct: Position of particle.
!   local_ref_frame 
!         -- Logical, If True then take the input coordinates 
!               as being with respect to the frame of referene of the element. 
!
! Output:
!   dr_ds(10)   -- real(rp): Kick vector.
!   field       -- em_field_struct: Local field.
!   err         -- Logical: Set True if there is an error.
!-

subroutine kick_vector_calc (ele, param, s_rel, orbit, local_ref_frame, dr_ds, err)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (em_field_struct) field
type (coord_struct) orbit

real(rp), intent(in) :: s_rel
real(rp), intent(out) :: dr_ds(10)
real(rp) f_bend, gx_bend, gy_bend, dt_ds, dp_ds, dbeta_ds
real(rp) vel(3), E_force(3), B_force(3)
real(rp) e_tot, dt_ds_ref, p0, beta0, v2, pz_p0

integer rel_dir, sign_z_vel

logical :: local_ref_frame, err

character(24), parameter :: r_name = 'kick_vector_calc'

! Init

err = .true.
beta0 = ele%value(p0c$) / ele%value(e_tot$) 
dt_ds_ref = orbit%direction / (beta0 * c_light)
p0 = ele%value(p0c$) / c_light
e_tot = orbit%p0c * (1 + orbit%vec(6)) / orbit%beta
rel_dir = ele%orientation * orbit%direction

! Calculate the field. 
! Important: Field is in frame of element. When ele%orientation = -1 => +z in -s direction.

if (ele%orientation == 1) then
  call em_field_calc (ele, param, s_rel, orbit, local_ref_frame, field, .false., err)
else
  call em_field_calc (ele, param, ele%value(l$)-s_rel, orbit, local_ref_frame, field, .false., err)
endif
if (err) return

! Bend factor

vel(1:2) = [orbit%vec(2), orbit%vec(4)] / (1 + orbit%vec(6))
v2 = vel(1)**2 + vel(2)**2
if (v2 > 0.99999999_rp) return
vel = orbit%beta * c_light * [vel(1), vel(2), sqrt(1 - v2) * rel_dir]
E_force = charge_of(orbit%species) * field%E
B_force = charge_of(orbit%species) * cross_product(vel, field%B)

f_bend = 1 ! Longitudinal distance per unit s-distance. Equal to 1 except when off axis in a bend.
gx_bend = 0; gy_bend = 0  ! bend strengths (1/rho) in x and y planes.

if (ele%key == sbend$) then
  if (ele%value(ref_tilt_tot$) /= 0 .and. .not. local_ref_frame) then
    gx_bend = ele%value(g$) * cos(ele%value(ref_tilt_tot$))
    gy_bend = ele%value(g$) * sin(ele%value(ref_tilt_tot$))
  else
    gx_bend = ele%value(g$)
  endif
  f_bend = 1 + orbit%vec(1) * gx_bend + orbit%vec(3) * gy_bend
endif

dt_ds = orbit%direction * f_bend / abs(vel(3))
dp_ds = dot_product(E_force, vel) * dt_ds / (orbit%beta * c_light)
dbeta_ds = mass_of(orbit%species)**2 * dp_ds * c_light / e_tot**3
pz_p0 = (1 + orbit%vec(6)) * orbit%direction * abs(vel(3)) / (orbit%beta * c_light)  ! Pz / P0

dr_ds(1) = vel(1) * dt_ds
dr_ds(2) = (E_force(1) + B_force(1)) * dt_ds / p0 + gx_bend * pz_p0
dr_ds(3) = vel(2) * dt_ds
dr_ds(4) = (E_force(2) + B_force(2)) * dt_ds / p0 + gy_bend * pz_p0
dr_ds(5) = orbit%beta * c_light * (dt_ds_ref - dt_ds) + dbeta_ds * orbit%vec(5) / orbit%beta
dr_ds(6) = dp_ds / p0
dr_ds(7) = dt_ds

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  ! dr_ds(8:10) = Omega/v_z
  sign_z_vel = orbit%direction * ele%orientation
  dr_ds(8:10) = orbit%direction * (f_bend  * spin_omega (field, orbit, sign_z_vel) + sign_z_vel * [-gy_bend, gx_bend, 0.0_rp])
else
  dr_ds(8:10) = 0
endif

err = .false.

end subroutine kick_vector_calc

end module
