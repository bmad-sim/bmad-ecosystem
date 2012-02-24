module runge_kutta_mod

use em_field_mod
use tracking_integration_mod
use track1_mod

contains

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine odeint_bmad (orb_start, ele, param, orb_end, s1, s2, &
!                            rel_tol, abs_tol, h1, hmin, local_ref_frame, err_flag, track)
! 
! Subroutine to do Runge Kutta tracking. This routine is adapted from Numerical
! Recipes.  See the NR book for more details.
!
! Notice that this routine has an two tolerance arguments rel_tol and abs_tol.
! Odeint only has 1. rel_tol (essentually equivalent to eps in odeint) 
! is scalled by the step size to to able to relate it to the final accuracy.
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
!   orb_start -- Coord_struct: Starting coords: (x, px, y, py, z, delta).
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
!   rel_tol   -- Real: Same as eps for odeint scalled by sqrt(h/(s2-s1))
!                   where h is the step size for the current step. rel_tol
!                   sets the %error of the result
!   abs_tol   -- Real: Sets the absolute error of the result
!   h1        -- Real: Initial guess for a step size.
!   h_min     -- Real: Minimum step size (can be zero).
!   local_ref_frame 
!              -- Logical: If True then take the 
!                input and output coordinates as being with 
!                respect to the frame of referene of the element. 
!
! Output:
!   orb_end  -- Coord_struct: Ending coords: (x, px, y, py, z, delta).
!   err_flag -- Logical: Set True if there is an error. False otherwise.
!   track    -- Track_struct, optional: Structure holding the track information.
!-

subroutine odeint_bmad (orb_start, ele, param, orb_end, s1, s2, &
                    rel_tol, abs_tol, h1, h_min, local_ref_frame, err_flag, track)

implicit none

type (coord_struct), intent(in) :: orb_start
type (coord_struct), intent(out) :: orb_end
type (ele_struct) ele
type (ele_struct), pointer :: hard_ele
type (lat_param_struct) param
type (track_struct), optional :: track

real(rp), intent(in) :: s1, s2, rel_tol, abs_tol, h1, h_min
real(rp), parameter :: tiny = 1.0e-30_rp
real(rp) :: h, h_did, h_next, s, s_sav, rel_tol_eff, abs_tol_eff, sqrt_N, h_save
real(rp) :: dr_ds(7), r_scal(7), t, s_hard_edge

integer, parameter :: max_step = 10000
integer :: n_step, hard_end

logical local_ref_frame, abs_time, err_flag, err

! init

err_flag = .true.
s = s1
h = sign(h1, s2-s1)
orb_end = orb_start
if (s1 == s2) return

orb_end%s = s1 + ele%s + ele%value(s_offset_tot$) - ele%value(l$)

! If the element is using a hard edge model then need to stop at the hard edges
! to apply the appropriate hard edge kick.

nullify (hard_ele)
call calc_next_hard_edge (ele, s_hard_edge, hard_ele, hard_end)

! For elements where the reference energy is changing the reference energy in the body is 
! taken, by convention, to be the reference energy at the exit end.

call lcavity_reference_energy_correction (ele, param, orb_end)

abs_time = .false.
if (associated(ele%lat)) then
  if (ele%lat%absolute_time_tracking) abs_time = .true.
endif

if (abs_time) then
  t = orb_end%t
else
  t = -orb_end%vec(5) / (orb_end%beta * c_light)    ! Time
endif

! Save initial point

if (present(track)) then
  s_sav = s - 2.0_rp * track%ds_save
  call init_saved_orbit (track, 1000)
  if ((abs(s-s_sav) > track%ds_save)) call save_a_step (track, ele, param, local_ref_frame, s, orb_end, s_sav)
endif

! now track

err = .false.

do n_step = 1, max_step

  ! Check if we we need to apply a hard edge kick.
  ! For super_slaves there may be multiple hard edges at a single s-position.

  do
    if (.not. associated(hard_ele)) exit
    if ((s-s_hard_edge)*(s_hard_edge-s1) < 0.0) exit
    call apply_hard_edge_kick (orb_end, t, hard_ele, ele, param, hard_end)
    call calc_next_hard_edge (ele, s_hard_edge, hard_ele, hard_end)
    ! Trying to take a step through a hard edge can drive Runge-Kutta nuts.
    ! So offset s a very tiny amount to avoid this
    s = s + 1d-12
  enddo

  ! Check if we are done.

  if ((s-s2)*(s2-s1) >= 0.0) then
    if (present(track)) call save_a_step (track, ele, param, local_ref_frame, s, orb_end, s_sav)
    err_flag = .false.
    return
  end if

  ! Need to propagate a step. First calc tolerances.

  call kick_vector_calc (ele, param, s, t, orb_end, local_ref_frame, dr_ds)

  sqrt_N = sqrt(abs((s2-s1)/h))  ! number of steps we would take with this h
  rel_tol_eff = rel_tol / sqrt_N
  abs_tol_eff = abs_tol / sqrt_N
  r_scal(:) = abs([orb_end%vec(:), t]) + abs(h*dr_ds(:)) + TINY

  h_save = h
  if ((s+h-s_hard_edge)*(s+h-s1) > 0.0) h = s_hard_edge-s

  call rkqs_bmad (ele, param, orb_end, dr_ds, s, t, h, rel_tol_eff, abs_tol_eff, r_scal, h_did, h_next, local_ref_frame, err)
  if (err) return

  if (present(track)) then
    if ((abs(s-s_sav) > track%ds_save)) call save_a_step (track, ele, param, local_ref_frame, s, orb_end, s_sav)
    if (h_did == h) then
      track%n_ok = track%n_ok + 1
    else
      track%n_bad = track%n_bad + 1
    end if
  endif

  ! Calculate next step size. If there was a hard edge then take into account the step that would have
  ! been taken if no hard edge was present.

  if (h_save == h) then  ! No hard edge case
    h = h_next
  else                   ! Was hard edge case.
    if (h_next < h) then
      h = h_next
    else
      h = max(h_save, h_next)
    endif
  endif

  ! Check for step size smaller than minimum. If so we consider the particle lost
 
  if (abs(h) < h_min) exit

end do

! Here if step size too small or too many steps

param%lost = .true.
param%plane_lost_at = z_plane$

end subroutine odeint_bmad

!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine rkqs_bmad (ele, param, orb, dr_ds, s, t, h_try, rel_tol, abs_tol, &
                                       r_scal, h_did, h_next, local_ref_frame, err_flag)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, orb_new

real(rp), intent(in)    :: dr_ds(7), r_scal(7)
real(rp), intent(inout) :: s, t
real(rp), intent(in)    :: h_try, rel_tol, abs_tol
real(rp), intent(out)   :: h_did, h_next

real(rp) :: err_max, h, h_temp, s_new, p2
real(rp) :: r_err(7), r_temp(7)
real(rp) :: rel_pc, t_new
real(rp), parameter :: safety = 0.9_rp, p_grow = -0.2_rp
real(rp), parameter :: p_shrink = -0.25_rp, err_con = 1.89e-4

logical local_ref_frame, err_flag

!

h = h_try
orb_new = orb

do

  call rkck_bmad (ele, param, orb, dr_ds, s, t, h, orb_new, t_new, r_err, local_ref_frame)
  err_max = maxval(abs(r_err(:)/(r_scal(:)*rel_tol + abs_tol)))
  if (err_max <=  1.0) exit
  h_temp = safety * h * (err_max**p_shrink)
  h = sign(max(abs(h_temp), 0.1_rp*abs(h)), h)
  s_new = s + h

  if (s_new == s) then
    err_flag = .true.
    print *, 'ERROR IN RKQS_BMAD: STEPSIZE UNDERFLOW'
    if (bmad_status%exit_on_error) call err_exit
  endif

end do

if (err_max > err_con) then
  h_next = safety*h*(err_max**p_grow)
else
  h_next = 5.0_rp * h
end if

h_did = h
s = s+h

orb_new%s = orb%s + h
orb_new%t = orb%t + (t_new - t)

orb = orb_new
t = t_new

end subroutine rkqs_bmad

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine rkck_bmad (ele, param, orb, dr_ds, s, t, h, orb_new, t_new, r_err, local_ref_frame)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, orb_new, orb_temp(5)

real(rp), intent(in) :: dr_ds(7)
real(rp), intent(in) :: s, t, h
real(rp), intent(out) :: r_err(7), t_new
real(rp) :: ak2(7), ak3(7), ak4(7), ak5(7), ak6(7), t_temp(5)
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

complex(rp) dspin1(2), dspin3(2), dspin4(2), dspin6(2)

logical local_ref_frame

!

call transfer_this_orbit (orb, b21*h*dr_ds(1:6), orb_temp(1))
t_temp(1) = t + b21*h*dr_ds(7)
call kick_vector_calc(ele, param, s + a2*h, t_temp(1), orb_temp(1), local_ref_frame, ak2)

call transfer_this_orbit (orb, h*(b31*dr_ds(1:6) + b32*ak2(1:6)), orb_temp(2))
t_temp(2) = t + h*(b31*dr_ds(7) + b32*ak2(7))
call kick_vector_calc(ele, param, s + a3*h, t_temp(2), orb_temp(2), local_ref_frame, ak3)

call transfer_this_orbit (orb, h*(b41*dr_ds(1:6) + b42*ak2(1:6) + b43*ak3(1:6)), orb_temp(3))
t_temp(3) = t + h*(b41*dr_ds(7) + b42*ak2(7) + b43*ak3(7))
call kick_vector_calc(ele, param, s + a4*h, t_temp(3), orb_temp(3), local_ref_frame, ak4)

call transfer_this_orbit (orb, h*(b51*dr_ds(1:6) + b52*ak2(1:6) + b53*ak3(1:6) + b54*ak4(1:6)), orb_temp(4))
t_temp(4) = t + h*(b51*dr_ds(7) + b52*ak2(7) + b53*ak3(7) + b54*ak4(7))
call kick_vector_calc(ele, param, s + a5*h, t_temp(4), orb_temp(4), local_ref_frame, ak5)

call transfer_this_orbit (orb, h*(b61*dr_ds(1:6) + b62*ak2(1:6) + b63*ak3(1:6) + b64*ak4(1:6) + b65*ak5(1:6)), orb_temp(5))
t_temp(5) = t + h*(b61*dr_ds(7) + b62*ak2(7) + b63*ak3(7) + b64*ak4(7) + b65*ak5(7))
call kick_vector_calc(ele, param, s + a6*h, t_temp(5), orb_temp(5), local_ref_frame, ak6)

call transfer_this_orbit (orb, h*(c1*dr_ds(1:6) + c3*ak3(1:6) + c4*ak4(1:6) + c6*ak6(1:6)), orb_new)
t_new = t + h*(c1*dr_ds(7) + c3*ak3(7) + c4*ak4(7) + c6*ak6(7))
r_err=h*(dc1*dr_ds + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)

! Spin

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  call dspin_dz (ele, param, s,        t,         orb,         local_ref_frame, dspin1)
  call dspin_dz (ele, param, s + a3*h, t_temp(2), orb_temp(2), local_ref_frame, dspin3)
  call dspin_dz (ele, param, s + a4*h, t_temp(3), orb_temp(3), local_ref_frame, dspin4)
  call dspin_dz (ele, param, s + a6*h, t_temp(5), orb_temp(5), local_ref_frame, dspin6)

  orb_new%spin = orb%spin + h * (c1*dspin1 + c3*dspin3 + c4*dspin4 + c6*dspin6)
endif

!----------------------------------------------------------
contains

subroutine transfer_this_orbit (orb_in, dvec, orb_out)

type (coord_struct) orb_in, orb_out
real(rp) dvec(6)

!

orb_out%vec = orb_in%vec + dvec
orb_out%p0c = orb_in%p0c
if (dvec(6) == 0) then
  orb_out%beta = orb_in%beta
else
  call convert_pc_to (orb_out%p0c * (1 + orb_out%vec(6)), param%particle, beta = orb_out%beta)
endif

end subroutine

!----------------------------------------------------------
! contains

subroutine dspin_dz (ele, param, s, t, orb, local_ref_frame, dspin)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb
type (em_field_struct) :: field

real(rp) s, t
real(rp) :: Omega(3)

complex(rp) :: dspin(2), quaternion(2,2)

logical local_ref_frame

! this uses a modified Omega' = -Omega/v_z

Omega = spin_omega_at (field, orb, ele, param, s)
quaternion = (i_imaginary/2.0_rp)* (pauli(1)%sigma*Omega(1) + pauli(2)%sigma*Omega(2) + pauli(3)%sigma*Omega(3))

!   quaternion = normalized_quaternion (quaternion)

dspin = matmul(quaternion, orb%spin)

end subroutine dspin_dz

end subroutine rkck_bmad

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine kick_vector_calc (ele, param, s_rel, t_rel, orbit, local_ref_frame, dr_ds)
!
! Subroutine to calculate the dr/ds "kick vector" where
!     r = [x, p_x, y, p_y, z, p_z, t] 
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
!   dr(2)/ds = dP_x/dt * dt/ds / P0
!   where:
!     dP_x/dt = EM_Force_x + gamma * mass * vec * (dtheta/dt)^2
!     vec * (dtheta/dt)^2 = g * v_s^2 / (1 + g*x)^2   ! Centrifugal term
!
!   dr(3)/ds = dy/dt * dt/ds
!   where:
!     dy/dt = v_x 
! 
!   dr(4)/ds = dP_y/dt * ds/dt / P0
!   where:
!     dP_y/dt = EM_Force_y
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
! Modules needed:
!   use bmad
!
! Input:
!   ele   -- Ele_struct: Element being tracked thorugh.
!   param -- lat_param_struct: Lattice parameters.
!   s_rel -- Real(rp): Distance from the start of the element to the particle.
!   t_rel  -- Real(rp): Time relative to the reference particle.
!   orbit -- coord_struct: Position of particle.
!   local_ref_frame 
!         -- Logical, If True then take the input coordinates 
!               as being with respect to the frame of referene of the element. 
!
! Output:
!   dr_ds(7) -- Real(rp): Kick vector.
!-

subroutine kick_vector_calc (ele, param, s_rel, t_rel, orbit, local_ref_frame, dr_ds)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (em_field_struct) field
type (coord_struct) orbit

real(rp), intent(in) :: s_rel, t_rel 
real(rp), intent(out) :: dr_ds(7)
real(rp) f_bend, gx_bend, gy_bend, dt_ds, dp_ds, dbeta_ds
real(rp) vel(3), E_force(3), B_force(3)
real(rp) e_tot, dt_ds_ref, p0, beta0

logical :: local_ref_frame

character(24), parameter :: r_name = 'kick_vector_calc'

!

beta0 = ele%value(p0c$) / ele%value(e_tot$) 
dt_ds_ref = 1 / (beta0 * c_light)
p0 = ele%value(p0c$) / c_light
e_tot = orbit%p0c * (1 + orbit%vec(6)) / orbit%beta

! calculate the field

call em_field_calc (ele, param, s_rel, t_rel, orbit, local_ref_frame, field, .false.)

! Bend factor

vel(1:2) = [orbit%vec(2), orbit%vec(4)] / (1 + orbit%vec(6))
vel = orbit%beta * c_light * [vel(1), vel(2), sqrt(1 - vel(1)**2 - vel(2)**2)]
E_force = charge_of(param%particle) * field%E
B_force = charge_of(param%particle) * cross_product(vel, field%B)

f_bend = 1
gx_bend = 0; gy_bend = 0

if (ele%key == sbend$) then
  if (ele%value(tilt_tot$) /= 0 .and. .not. local_ref_frame) then
    gx_bend = ele%value(g$) * cos(ele%value(tilt_tot$))
    gy_bend = ele%value(g$) * sin(ele%value(tilt_tot$))
  else
    gx_bend = ele%value(g$)
  endif
  f_bend = 1 + orbit%vec(1) * gx_bend + orbit%vec(3) * gy_bend
endif

dt_ds = f_bend / vel(3)
dp_ds = dot_product(E_force, vel) * dt_ds / (orbit%beta * c_light)
dbeta_ds = mass_of(param%particle)**2 * dp_ds / e_tot**3

dr_ds(1) = vel(1) * dt_ds
dr_ds(2) = (E_force(1) + B_force(1) + e_tot * gx_bend / (dt_ds * c_light)**2) * dt_ds / p0
dr_ds(3) = vel(2) * dt_ds
dr_ds(4) = (E_force(2) + B_force(2) + e_tot * gy_bend / (dt_ds * c_light)**2) * dt_ds / p0
dr_ds(5) = orbit%beta * c_light * (dt_ds_ref - dt_ds) + dbeta_ds * orbit%vec(5) / orbit%beta
dr_ds(6) = dp_ds / p0
dr_ds(7) = dt_ds

end subroutine kick_vector_calc

end module
