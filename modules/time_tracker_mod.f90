module time_tracker_mod

use bmad_struct
use beam_def_struct
use em_field_mod

type (ele_struct), save, private, pointer :: ele_com
type (lat_param_struct), save, private, pointer :: param_com
type (coord_struct), save, private, pointer ::  orb_old_com
type (coord_struct), save, private, pointer :: orb_com
real(rp), save, private, pointer :: vec_err_com(:)
real(rp), save, private, pointer :: s_target_com, t_old_com
real(rp), save, private, pointer :: t_rel_com
logical, save, private, pointer :: local_ref_frame_com

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine odeint_bmad_time (orb, ele, param, s1, s2, t_rel, &
!                             dt1, local_ref_frame, err_flag, track)
! 
! Subroutine to do Runge Kutta tracking in time. This routine is adapted from Numerical
! Recipes.  See the NR book for more details.
!
!
! Modules needed:
!   use bmad
!
! Input: 
!   orb   -- Coord_struct: Starting coords: (x, px, y, py, s, ps) [t-based]
!   ele     -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field_standard
!   param   -- lat_param_struct: Beam parameters.
!     %enegy       -- Energy in GeV
!     %particle    -- Particle type [positron$, or electron$]
!   s1      -- Real: Starting point.
!   s2      -- Real: Ending point.
!   t_rel   -- Real: time relative to entering reference time
!   dt1      -- Real: Initial guess for a time step size.
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
!   err_flag -- Logical: Set True if there is an error. False otherwise.
!   track    -- Track_struct: Structure holding the track information.
!
!-

subroutine odeint_bmad_time (orb, ele, param, s1, s2, t_rel, &
                                dt1, local_ref_frame, err_flag, track)

use nr, only: zbrent

implicit none

type (coord_struct), intent(inout), target :: orb
type (coord_struct), target :: orb_old
type (ele_struct), target :: ele
type (lat_param_struct), target ::  param
type (em_field_struct) :: saved_field
type (track_struct), optional :: track

real(rp), intent(in) :: s1, s2, dt1
real(rp), target :: t_rel, t_old, dt_tol
real(rp) :: dt, dt_did, dt_next, ds_safe, t_save, dt_save
real(rp), target  :: dvec_dt(6), vec_err(6), s_target

integer, parameter :: max_step = 100000
integer :: n_step, n_pt

logical, target :: local_ref_frame
logical :: exit_flag, err_flag, err, zbrent_needed, add_ds_safe

character(30), parameter :: r_name = 'odeint_bmad_time'

! init
ds_safe = bmad_com%significant_length / 10
dt_next = dt1

! local s coordinates for vec(5)
orb%vec(5) = orb%s - (ele%s - ele%value(l$))

! Allocate track arrays

!n_pt = max_step
if ( present(track) ) then
   dt_save = track%ds_save/c_light
   t_save = t_rel
endif 

! Now Track

exit_flag = .false.
err_flag = .true.

do n_step = 1, max_step

  ! Single Runge-Kutta step. Updates orb% vec(6), s, and t 

  dt = dt_next
  orb_old = orb
  t_old = t_rel
  call rk_adaptive_time_step (ele, param, orb, t_rel, dt, dt_did, dt_next, local_ref_frame, err)

  ! Check entrance and exit faces
  if (orb%vec(5) > s2 - ds_safe .or. orb%vec(5) < s1 + ds_safe) then
    zbrent_needed = .true.
    add_ds_safe = .true.

    if (orb%vec(5) < s1 + ds_safe) then 
      s_target = s1
      if (orb%vec(5) > s1 - ds_safe) zbrent_needed = .false.
      if (s1 == 0) add_ds_safe = .false.
    else
      s_target = s2
      if (orb%vec(5) < s2 + ds_safe) zbrent_needed = .false.
      if (abs(s2 - ele%value(l$)) < ds_safe) add_ds_safe = .false.
    endif

    exit_flag = .true.
    orb%status = outside$

    !Set common structures for zbrent's internal functions 
    ele_com => ele
    param_com => param
    orb_old_com => orb_old
    orb_com => orb
    vec_err_com => vec_err
    local_ref_frame_com => local_ref_frame
    s_target_com => s_target
    t_rel_com => t_rel
    t_old_com => t_old

    !---
    dt_tol = ds_safe / (orb%beta * c_light)
    if (zbrent_needed) dt = zbrent (delta_s_target, 0.0_rp, dt, dt_tol)
    ! Trying to take a step through a hard edge can drive Runge-Kutta nuts.
    ! So offset s a very tiny amount to avoid this
    orb%vec(5) = s_target 
    if (add_ds_safe) orb%vec(5) = orb%vec(5) + sign(ds_safe, orb%vec(6))
    orb%s = orb%vec(5) + ele%s - ele%value(l$)
  endif

  ! Check wall or aperture at every intermediate step and flag for exit if wall is hit

  call  particle_hit_wall_check_time(orb_old, orb, param, ele)
  if (orb%status == dead$) exit_flag = .true.
  
  !Save track
  if ( present(track) ) then
    !Check if we are past a save time, or if exited
    if (t_rel >= t_save .or. exit_flag) then
      track%n_pt = track%n_pt + 1
      n_pt = track%n_pt
      track%orb(n_pt) = orb
      !Query the local field to save
      call em_field_calc (ele, param, orb%vec(5), t_rel, orb, local_ref_frame, saved_field, .false., err_flag)
      if (err_flag) return
      track%field(n_pt) = saved_field
       !Set next save time 
       t_save = t_rel + dt_save
    end if
  endif

  ! Exit when the particle hits surface s1 or s2, or hits wall
  if (exit_flag) then
    err_flag = .false. 
    return
  endif

end do

if (bmad_status%type_out) then
  call out_io (s_warn$, r_name, 'STEPS EXCEEDED MAX_STEP FOR ELE: '//ele%name )
  !print *, '  Skipping particle; coordinates will not be saved'
  orb%status = inside$
  return
end if

if (bmad_status%exit_on_error) call err_exit

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
real(rp) :: r_err(6), r_temp(6), dr_dt(6)
real(rp) :: r_scal(6), rel_tol, abs_tol
real(rp), parameter :: safety = 0.9_rp, p_grow = -0.2_rp
real(rp), parameter :: p_shrink = -0.25_rp, err_con = 1.89e-4
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
  if (err_flag) return
  r_scal(:) = abs(orb%vec) + abs(orb_new%vec) + TINY
  r_scal(5) = ele%value(l$)
  err_max = maxval(abs(r_err(:)/(r_scal(:)*rel_tol + abs_tol)))
  if (err_max <=  1.0) exit
  dt_temp = safety * dt * (err_max**p_shrink)
  dt = sign(max(abs(dt_temp), 0.1_rp*abs(dt)), dt)
  t_new = t + dt

  if (t_new == t) then
    err_flag = .true.
    call out_io (s_fatal$, r_name, 'STEPSIZE UNDERFLOW IN ELEMENT: ' // ele%name)
    if (bmad_status%exit_on_error) call err_exit
    return
  endif

end do

if (err_max > err_con) then
  dt_next = safety*dt*(err_max**p_grow)
else
  dt_next = 5.0_rp * dt
end if

dt_did = dt
t = t+dt

orb = orb_new

end subroutine rk_adaptive_time_step

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine rk_time_step1 (ele, param, orb, t, dt, orb_new, r_err, local_ref_frame, dr_dt, err_flag)

!Very similar to rk_step1_bmad, except that em_field_kick_vector_time is called
!  and orb_new%s and %t are updated to the global values

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, orb_new, orb_temp

real(rp), optional, intent(in) :: dr_dt(6)
real(rp), intent(in) :: t, dt
real(rp), intent(out) :: r_err(6)
real(rp) :: dr_dt1(6), dr_dt2(6), dr_dt3(6), dr_dt4(6), dr_dt5(6), dr_dt6(6), r_temp(6), pc
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

orb_temp%vec = orb%vec + b21*dt*dr_dt1
call em_field_kick_vector_time(ele, param, t+a2*dt, orb_temp, local_ref_frame, dr_dt2, err_flag)
if (err_flag) return

orb_temp%vec = orb%vec + dt*(b31*dr_dt1+b32*dr_dt2)
call em_field_kick_vector_time(ele, param, t+a3*dt, orb_temp, local_ref_frame, dr_dt3, err_flag) 
if (err_flag) return

orb_temp%vec = orb%vec + dt*(b41*dr_dt1+b42*dr_dt2+b43*dr_dt3)
call em_field_kick_vector_time(ele, param, t+a4*dt, orb_temp, local_ref_frame, dr_dt4, err_flag)
if (err_flag) return

orb_temp%vec = orb%vec + dt*(b51*dr_dt1+b52*dr_dt2+b53*dr_dt3+b54*dr_dt4)
call em_field_kick_vector_time(ele, param, t+a5*dt, orb_temp, local_ref_frame, dr_dt5, err_flag)
if (err_flag) return

orb_temp%vec = orb%vec + dt*(b61*dr_dt1+b62*dr_dt2+b63*dr_dt3+b64*dr_dt4+b65*dr_dt5)
call em_field_kick_vector_time(ele, param, t+a6*dt, orb_temp, local_ref_frame, dr_dt6, err_flag)
if (err_flag) return

!Output new orb and error vector

orb_new%vec = orb%vec +dt*(c1*dr_dt1+c3*dr_dt3+c4*dr_dt4+c6*dr_dt6)
orb_new%t = orb%t + dt
orb_new%s = orb%s + orb_new%vec(5) - orb%vec(5)
pc = sqrt(orb_new%vec(2)**2 +orb_new%vec(4)**2 + orb_new%vec(6)**2)
call convert_pc_to (pc, param%particle, beta = orb_new%beta)

r_err = dt*(dc1*dr_dt1+dc3*dr_dt3+dc4*dr_dt4+dc5*dr_dt5+dc6*dr_dt6)

end subroutine rk_time_step1

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! function for zbrent to calculate timestep to exit face surface

function delta_s_target (this_dt)
  real(rp), intent(in)  :: this_dt
  real(rp) :: delta_s_target
  logical err_flag
  !
  call rk_time_step1 (ele_com, param_com, orb_old_com, t_old_com, this_dt, &
                          orb_com, vec_err_com, local_ref_frame_com, err_flag = err_flag)
  delta_s_target = orb_com%vec(5) - s_target_com
  t_rel_com = t_old_com + this_dt

end function delta_s_target
  
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine particle_hit_wall_check_time
!
! Subroutine to check whether particle has collided with element walls,
! and to calculate location of impact if it has
!
! Modules needed:
!   use bmad
!   use capillary_mod
!
! Input
!   orb     -- coord_struct: Previous particle coordinates  [t-based]
!   orb_new -- coord_struct: Current particle coordinates [t-based]
!   param   -- lat_param_struct: Lattice parameters
!    %particle -- integer: Type of particle
!   ele     -- ele_struct: Lattice element
!
! Output
!   param   -- lat_param_struct: Lattice parameters
!    %lost -- logical: True if orbit hit wall
!   orb_new -- coord_struct: Location of hit
!    %phase_x -- real(rp): Used to store hit angle
!-

subroutine  particle_hit_wall_check_time(orb, orb_new, param, ele)
 
use capillary_mod

implicit none

type (coord_struct) :: orb, orb_new
type (coord_struct), pointer :: old_orb, now_orb
type (lat_param_struct) :: param
type (ele_struct) :: ele
type (photon_track_struct), target :: particle
integer :: section_ix
real(rp) :: norm, perp(3), dummy_real
real(rp) :: edge_tol = 1e-8

!-----------------------------------------------
!Do nothing if there is no wall
if ( .not. associated(ele%wall3d%section) ) then
   return
end if

   
!Prepare coordinate structures for capillary_photon_d_radius
old_orb => particle%old%orb
now_orb => particle%now%orb
old_orb = orb
now_orb = orb_new


!If now_orb is before element, change it
!We can do this because it can't hit the element wall outside of the
!element, and these do not affect tracks
if (now_orb%vec(5) < 0) then
   now_orb%vec(5) = 0
end if


!If now_orb is too close to the wall, move it edge_tol away
if (abs(capillary_photon_d_radius(particle%now, ele)) < edge_tol) then
   now_orb%vec(1) = now_orb%vec(1) + sign(edge_tol, now_orb%vec(1))
   now_orb%vec(3) = now_orb%vec(3) + sign(edge_tol, now_orb%vec(3))
end if


!Change from particle coordinates to photon coordinates
! (coord_struct to photon_coord_struct)

!Get e_tot from momentum, calculate beta_i = c*p_i / p_tot, pretending that these are traveling at v=c
!p_tot = sqrt(orb%vec(2)**2 + orb%vec(4)**2 + orb%vec(6)**2)
!old_orb%vec(2) = orb%vec(2) / p_tot
!old_orb%vec(4) = orb%vec(4) / p_tot
!old_orb%vec(6) = orb%vec(6) / p_tot

!p_tot = sqrt(orb_new%vec(2)**2 + orb_new%vec(4)**2 + orb_new%vec(6)**2 )
!now_orb%vec(2) = orb_new%vec(2) / p_tot
!now_orb%vec(4) = orb_new%vec(4) / p_tot
!now_orb%vec(6) = orb_new%vec(6) / p_tot


!More coordinate changes
!Equations taken from track_a_capillary in capillary_mod
!particle%old%energy = ele%value(e_tot$) * (1 + orb%vec(6))

!particle%old%ix_section = 1

!particle%now%energy = ele%value(e_tot$) * (1 + orb_new%vec(6))

!Pretend that now_orb and old_orb are photons to calculate the wall intersection
particle%old%track_len = 0
particle%now%track_len = sqrt( &
     (now_orb%vec(1) - old_orb%vec(1))**2 + &
     (now_orb%vec(3) - old_orb%vec(3))**2 + &
     (now_orb%vec(5) - old_orb%vec(5))**2)

old_orb%vec(2) = (now_orb%vec(1) - old_orb%vec(1)) /particle%now%track_len
old_orb%vec(4) = (now_orb%vec(3) - old_orb%vec(3)) /particle%now%track_len
old_orb%vec(6) = (now_orb%vec(5) - old_orb%vec(5)) /particle%now%track_len
now_orb%vec(2) = old_orb%vec(2) 
now_orb%vec(4) = old_orb%vec(4) 
now_orb%vec(6) = old_orb%vec(6) 

!If particle hit wall, find out where
if (capillary_photon_d_radius(particle%now, ele) > 0) then

   call capillary_photon_hit_spot_calc (particle, ele)

   orb_new = now_orb

   !Calculate perpendicular to get angle of impact
   dummy_real = capillary_photon_d_radius(particle%now, ele, perp)

   !Calculate angle of impact; cos(hit_angle) = norm_photon_vec \dot perp
   !****
   !Notice we store this in orb_new%phase_x; this is so that we don't have
   !to add another argument- yeah, it's a hack.
   !****
   orb_new%phase_x = acos( &
        ((now_orb%vec(1) - old_orb%vec(1)) * perp(1) + &
        (now_orb%vec(3) - old_orb%vec(3)) * perp(2) + &
        (now_orb%vec(5) - old_orb%vec(5)) * perp(3)) / particle%now%track_len)

   !Restore momenta from original orb 
   orb_new%vec(2) = orb%vec(2) 
   orb_new%vec(4) = orb%vec(4) 
   orb_new%vec(6) = orb%vec(6) 

   !Set orb%s
   orb_new%s = orb_new%vec(5) + ele%s - ele%value(l$)

   !Note that the time is not set!

   param%lost = .true.
   orb_new%status = dead$
endif

end subroutine  particle_hit_wall_check_time

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
!    dvec_dt(6)  -- real(rp): Derivatives.
!-

subroutine em_field_kick_vector_time (ele, param, t_rel, orbit, local_ref_frame, dvec_dt, err_flag)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (em_field_struct) field

type (coord_struct), intent(in) :: orbit

real(rp), intent(in) :: t_rel    
real(rp), intent(out) :: dvec_dt(6)

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

mc2 = mass_of(param%particle) ! Note: mc2 is in eV
charge = charge_of(param%particle) ! Note: charge is in units of |e_charge|

e_tot = sqrt( orbit%vec(2)**2 +  orbit%vec(4)**2 +  orbit%vec(6) **2 + mc2**2) 
vel(1:3) = c_light*[  orbit%vec(2),  orbit%vec(4),  orbit%vec(6) ]/ e_tot 

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
! dcp_y/dt = cp_s * v_s * \kappa_y / h + c*charge * ( Ey + * Bx - v_x * Bs )
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

if (ele%key == sbend$) then   
  if (ele%value(tilt_tot$) /= 0 .and. .not. local_ref_frame) then
    kappa_x = ele%value(g$) * cos(ele%value(tilt_tot$))
    kappa_y = ele%value(g$) * sin(ele%value(tilt_tot$))
  else
    kappa_x = ele%value(g$)
    kappa_y = 0
  endif
  h = 1 + kappa_x *  orbit%vec(1) + kappa_y *  orbit%vec(3) ! h = 1 + \kappa_x * x + \kappa_y * y

  dvec_dt(2) = dvec_dt(2) + orbit%vec(6) * vel(3) * kappa_x / h
  dvec_dt(4) = dvec_dt(4) + orbit%vec(6) * vel(3) * kappa_y / h
  dvec_dt(5) = vel(3) / h
  dvec_dt(6) = dvec_dt(6) - orbit%vec(6) * (vel(1)*kappa_x + vel(2)*kappa_y) / h
endif

end subroutine em_field_kick_vector_time

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine convert_particle_coordinates_t_to_s (particle, mc2, tref)
!
! Subroutine to convert particle coordinates from t-based to s-based system. 
!
! Modules needed:
!   use bmad
!
! Input:
!   particle   -- coord_struct: input particle coordinates
!                    %vec(:)
!                    %t 
!                    %p0c
!   mc2        -- real: particle rest mass in eV
!   tref       -- real: reference time for z coordinate
! Output:
!    particle   -- coord_struct: output particle 
!-

subroutine convert_particle_coordinates_t_to_s (particle, mc2, tref)

!use bmad_struct

implicit none

type (coord_struct), intent(inout), target ::particle
real(rp) :: p0c
real(rp), intent(in) :: mc2
real(rp), intent(in) :: tref

real(rp) :: pctot

real(rp), pointer :: vec(:)
vec => particle%vec
p0c=abs(particle%p0c)

! Convert t to s
pctot = sqrt (vec(2)**2 + vec(4)**2 + vec(6)**2)
! vec(1) = vec(1)   !this is unchanged
vec(2) = vec(2)/p0c
! vec(3) = vec(3)   !this is unchanged
vec(4) = vec(4)/p0c
! z \equiv -c \beta(s)  (t(s) - t_0(s)) 
vec(5) = -c_light * (pctot/sqrt(pctot**2 +mc2**2)) *  (particle%t - tref) 
vec(6) = pctot/p0c - 1.0_rp

end subroutine convert_particle_coordinates_t_to_s

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine convert_particle_coordinates_s_to_t (particle)
!
! Subroutine to convert particle coordinates from s-based to t-based system. 
!     The sign of particle%p0c indicates the direction of p_s
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
!   p0c        -- real: Reference momentum. The sign indicates direction of p_s 
! Output:
!    particle   -- coord_struct: output particle 
!-

subroutine convert_particle_coordinates_s_to_t (particle)

implicit none

type (coord_struct), intent(inout), target :: particle
real(rp), pointer :: vec(:)

vec => particle%vec

! Convert s to t
vec(6) = particle%p0c * sqrt( ((1+vec(6)))**2 - vec(2)**2 -vec(4)**2 )
! vec(1) = vec(1) !this is unchanged
vec(2) = vec(2)*abs(particle%p0c)
! vec(3) = vec(3) !this is unchanged
vec(4) = vec(4)*abs(particle%p0c)
vec(5) = particle%s

end subroutine convert_particle_coordinates_s_to_t

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine drift_orbit_time(orbit, mc2, delta_s)
!
! Simple routine to drift a particle orbit in time-based coordinates by a distance delta_s
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
  
  ! Get e_tot from momentum

  e_tot = sqrt( orbit%vec(2)**2 + orbit%vec(4)**2 +  orbit%vec(6)**2 + mc2**2) 

  ! velocities v_x, v_y, v_s:  c*[c*p_x, c*p_y, c*p_s]/e_tot

  vel(1:3) = c_light*[  orbit%vec(2), orbit%vec(4), orbit%vec(6) ]/ e_tot 

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
! Subroutine write_time_particle_distribution  (time_file_unit, bunch, mc2, err)
!
! Subroutine to write an time-based bunch from a standard Bmad bunch
! 
! Note: The time-based file format is
!       n_particles
!       x/m  m*c^2 \beta_x*\gamma/eV  y/m m*c^2\beta_y*\gamma/eV s/m m*c^2\beta_z*\gamma/eV time/s charge/C
!       . . .
!       all at the same time. 
!       This is very similar to subroutine write_opal_particle_distribution
!
! Input:
!   time_file_unit -- Integer: unit number to write to, if > 0
!   bunch          -- bunch_struct: bunch to be written.
!                            Particles are drifted to bmad_bunch%t_center for output
!   mc2            -- real(rp): particle mass in eV
!
! Output:          
!   err            -- Logical, optional: Set True if, say a file could not be opened.
!-



subroutine write_time_particle_distribution (time_file_unit, bunch, mc2, err)

implicit none

integer			    :: time_file_unit
type (bunch_struct) :: bunch
real(rp)            :: mc2
logical, optional   :: err

type (coord_struct) :: orb
real(rp)        :: dt, pc, gmc
character(40)	:: r_name = 'write_time_particle_distribution'
character(10)   ::  rfmt 
integer n_particle, i


!
if (present(err)) err = .true.

n_particle = size(bunch%particle)

!Format for numbers
  rfmt = 'es13.5'

!Write number of particles to first line
write(time_file_unit, '(i8)') n_particle

!\gamma m c

!Write out all particles to file
do i = 1, n_particle
  orb = bunch%particle(i)
  
  !Get time to track backwards by
  dt = orb%t - bunch%t_center
  
  !Get pc before conversion
  pc = (1+orb%vec(6))*orb%p0c 
  
  !convert to time coordinates
  call convert_particle_coordinates_s_to_t (orb)
  
  !get \gamma m c
  gmc = sqrt(pc**2 + mc2**2) / c_light
  
  !'track' particles backwards in time and write to file
  ! (x, y, s) - dt mc2 \beta_x \gamma / \gamma m c
  write(time_file_unit, '(8'//rfmt//')')  orb%vec(1) - dt*orb%vec(2)/gmc, orb%vec(2), &
                                          orb%vec(3) - dt*orb%vec(4)/gmc, orb%vec(4), &
                                          orb%vec(5) - dt*orb%vec(6)/gmc, orb%vec(6), &
                    										  bunch%t_center, bunch%particle(i)%charge 
end do 

end subroutine  write_time_particle_distribution


end module
