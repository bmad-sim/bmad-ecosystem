module runge_kutta_mod

use bmad_interface

type runge_kutta_common_struct
  integer :: num_steps_done = -1              ! Number of integration steps. Not used by Bmad. For external use.
  logical :: calc_field_derivatives = .false. ! Experimental: For use with extensions to the spin BMT equation that involve field derivatives.
end type

type (runge_kutta_common_struct), save :: runge_kutta_com

private :: rk_adaptive_step, rk_step1

contains

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine odeint_bmad (orbit, ele, param, s1_body, s2_body, err_flag, track, mat6, make_matrix)
! 
! Subroutine to do Runge Kutta tracking. This routine is adapted from Numerical
! Recipes.  See the NR book for more details.
!
! Notice that this routine has an two tolerances: 
!   bmad_com%rel_tol_adaptive_tracking
!   bmad_com%abs_tol_adaptive_tracking
!
! Note: For elements where the reference energy is not constant (lcavity, etc.), and 
! with elements where the reference particle does not follow the reference trajectory (wigglers for example),
! the calculation of z is "off" while the particle is inside the element. At the ends there is no problem.
!
! Input: 
!   orbit        -- coord_struct: Starting coords: (x, px, y, py, z, delta) in element body coords.
!   ele          -- Ele_struct: Element to track through.
!   param        -- lat_param_struct: Lattice parameters.
!   s1_body      -- real: Starting point relative to physical entrance.
!   s2_body      -- real: Ending point relative physical entrance.
!   mat6(6,6)    -- real(rp), optional: Transfer matrix before the element.
!   make_matrix  -- logical, optional: If True then make the 6x6 transfer matrix.
!
! Output:
!   orbit       -- Coord_struct: Ending coords: (x, px, y, py, z, delta) in element body coords.
!   err_flag    -- Logical: Set True if there is an error. False otherwise. Note: a particle getting
!                    lost, for example hitting an aperture, is *not* an error.
!   track       -- Track_struct, optional: Structure holding the track information.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix propagated through the element.
!-

recursive subroutine odeint_bmad (orbit, ele, param, s1_body, s2_body, err_flag, track, mat6, make_matrix)

use super_recipes_mod, only: super_zbrent

implicit none

type (coord_struct) :: orbit
type (coord_struct) old_orbit
type (ele_struct) ele
type (lat_param_struct) param
type (track_struct), optional :: track
type (fringe_field_info_struct) fringe_info

real(rp), intent(in) :: s1_body, s2_body
real(rp), optional :: mat6(6,6)
real(rp), target :: ds, ds_did, ds_next, s_body, s_last_save, ds_saved, s_edge_body
real(rp) :: old_s, ds_zbrent, dist_to_wall, ds_tiny, mat6_old(6,6)
real(rp), pointer :: s_body_ptr

integer :: n_step, s_dir, nr_max, n_step_max, status

logical, optional :: make_matrix
logical err_flag, err, at_hard_edge, track_spin, too_large

character(*), parameter :: r_name = 'odeint_bmad'

!

err_flag = .true.
s_body = s1_body
s_last_save = s_body
s_body_ptr => s_body     ! s_body_ptr used to get around an ifort bug.
s_dir = sign(1.0_rp, s2_body-s1_body)
ds_next = bmad_com%init_ds_adaptive_tracking * s_dir
if (ele%tracking_method == fixed_step_runge_kutta$) ds_next = ele%value(ds_step$) * s_dir
ds_tiny  = bmad_com%significant_length/100
track_spin = (ele%spin_tracking_method == tracking$ .and. &
                                (ele%field_calc == bmad_standard$ .or. ele%field_calc == fieldmap$))

if (ele%orientation == 1) then
  orbit%s = s1_body + ele%s_start
else
  orbit%s = ele%s - s1_body
endif

! For elements where the reference energy is changing the reference energy in the body is 
! taken, by convention, to be the reference energy at the exit end.

call reference_energy_correction (ele, orbit, first_track_edge$, mat6, make_matrix)

! If the element is using a hard edge model then need to stop at the hard edges
! to apply the appropriate hard edge kick.
! calc_next_fringe_edge assumes that s = 0 is beginning of element which is not true of a patch element.

call calc_next_fringe_edge (ele, s_edge_body, fringe_info, orbit, .true.)
if (ele%key == patch$) s_edge_body = s2_body

! Save initial point

if (present(track)) then
  call save_a_step (track, ele, param, .true., orbit, s_body, .true., mat6, make_matrix)
  s_last_save = s_body
endif

! now track

err = .false.

n_step_max = bmad_com%max_num_runge_kutta_step
if (ele%tracking_method == fixed_step_runge_kutta$) n_step_max = max(n_step_max, 2*nint(ele%value(num_steps$)))

do n_step = 1, n_step_max

  runge_kutta_com%num_steps_done = n_step
  
  ! Check if we we need to apply a hard edge kick.
  ! For super_slaves there may be multiple hard edges at a single s-position.

  do
    if (.not. fringe_info%has_fringe .or. .not. associated(fringe_info%hard_ele)) exit
    if ((s_body-s_edge_body)*s_dir < -ds_tiny) exit

    if (present(track) .and. s_body /= s_last_save) call save_a_step (track, ele, param, .true., orbit, s_body, .true., mat6, make_matrix)

    call apply_element_edge_kick (orbit, fringe_info, ele, param, track_spin, mat6, make_matrix)

    if (present(track)) then
      call save_a_step (track, ele, param, .true., orbit, s_body, .true., mat6, make_matrix)
      s_last_save = s_body
    endif

    call calc_next_fringe_edge (ele, s_edge_body, fringe_info, orbit)
    ! Trying to take a step through a hard edge can drive Runge-Kutta nuts.
    ! So offset s a very tiny amount to avoid this
    if (abs(s_body - s2_body) > ds_tiny) then
      s_body = s_body + ds_tiny * s_dir
      orbit%s = orbit%s + ds_tiny * s_dir * ele%orientation
    endif
  enddo

  ! Check if we are done.
  ! Use s = s_body - ds_tiny for save_a_step to make sure s-position is not outside the element (don't want
  ! to have the field = 0 due to being outside).

  if ((s_body-s2_body)*s_dir > -ds_tiny) then
    if (present(track) .and. abs(s_body-s_last_save) > 2*ds_tiny) call save_a_step (track, ele, param, .true., orbit, s_body-ds_tiny*s_dir, .true., mat6, make_matrix)
    call reference_energy_correction (ele, orbit, second_track_edge$, mat6, make_matrix)
    err_flag = .false.
    return
  end if

  ! Need to propagate a step. First calc tolerances.

  ds = ds_next
  at_hard_edge = .false.

  if ((s_body+ds-s_edge_body)*s_dir > 0.0) then
    at_hard_edge = .true.
    ds_saved = ds
    ds = s_edge_body - s_body - 0.5_rp*ds_tiny*s_dir
  endif

  old_orbit = orbit
  old_s = s_body
  call rk_adaptive_step (ele, param, orbit, s_body, ds, abs(s2_body-s1_body), ds_did, ds_next, err, mat6, make_matrix)
  if (err) return

  ! Check x/y limit apertures

  select case (ele%aperture_at)
  case (continuous$, wall_transition$)
    call check_aperture_limit (orbit, ele, in_between$, param, old_orbit, check_momentum = .false.)
    if (orbit%state /= alive$) then
      if (n_step == 1) return  ! Cannot do anything if this is the first step
      ! Due to the zbrent finite tolerance, the particle may not have crossed the wall boundary.
      ! So step a small amount to make sure that the particle is past the wall.
      ds_zbrent = super_zbrent (wall_intersection_func, 0.0_rp, ds_did, 1e-15_rp, ds_tiny, status)
      dist_to_wall = wall_intersection_func(ds_zbrent+ds_tiny, status)

      if (associated(wall_hit_handler_custom_ptr)) call wall_hit_handler_custom_ptr (orbit, ele, s_body)
      if (orbit%state /= alive$) return
      if (ele%aperture_at /= wall_transition$) then
        call out_io (s_error$, r_name, 'CUSTOM CODE IS KEEPING A PARTICLE ALIVE ACCROSS A BOUNDARY!', &
                                       'IN THIS CASE, THE APERTURE_AT COMPONENT OF ELEMENT: ' // ele%name, &
                                       'NEEDS TO BE SET TO "WALL_TRANSITION".')
        if (global_com%exit_on_error) call err_exit
        return
      endif
    endif
  end select

  ! Check for a crazy orbit

  if (orbit%state == alive$) too_large = orbit_too_large (orbit, param)

  ! Save track

  if (present(track)) then
    if (track%ds_save <= 0 .or. (abs(s_body-s_last_save) > track%ds_save)) then
      call save_a_step (track, ele, param, .true., orbit, s_body, .true., mat6, make_matrix)
      s_last_save = s_body
    endif

    if (ds_did == ds) then
      track%n_ok = track%n_ok + 1
    else
      track%n_bad = track%n_bad + 1
    endif
  endif

  ! Exit?

  if (orbit%state /= alive$) then
    err_flag = .false.
    return
  endif

  ! Calculate next step size. If there was a hard edge then take into account the step that would have
  ! been taken if no hard edge was present.

  if (at_hard_edge .and. abs(ds_next) >= abs(ds)) then
    ds_next = max(abs(ds_saved), abs(ds_next)) * s_dir
  endif

  if ((s_body + ds_next - s2_body) * s_dir > 0) then
    ds_next = s2_body - s_body
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

if (sqrt(orbit%vec(2)**2 + orbit%vec(4)**2) > 0.9_rp * (1.0_rp + orbit%vec(6)) .or. orbit%vec(6) < -0.99_rp) then
  orbit%state = lost_pz$
else
  call out_io (s_error$, r_name, 'STEP SIZE IS TOO SMALL OR TOO MANY STEPS WHILE TRACKING THROUGH: ' // ele%name, &
                                 'AT (X,Y,Z) POSITION FROM ENTRANCE: \3F14.7\ ', &
                                 'TYPICALLY THIS IS DUE TO THE FIELD NOT OBEYING MAXWELL''S EQUATIONS.', &
                                 '[OFTEN TIMES THE FIELD IS NOT EVEN BE CONTINUOUS IN THIS CASE!]', &
                                 'THE PARTICLE WILL BE MARKED AS LOST.', &
                                 r_array = [orbit%vec(1), orbit%vec(3), s_body])
  orbit%state = lost$
endif

err_flag = .false.

!-----------------------------------------------------------------
contains

function wall_intersection_func (ds, status) result (d_radius)

real(rp), intent(in) :: ds
real(rp) d_radius
real(rp) t_new, r_err(11), dr_ds(11)

integer :: status
logical no_aperture_here

!

call rk_step1 (ele, param, old_orbit, dr_ds, old_s, ds, orbit, r_err, err_flag, .true.)

s_body_ptr = old_s + ds
orbit%s = s_body_ptr + ele%s_start
d_radius = distance_to_aperture (orbit, in_between$, ele, no_aperture_here)

if (no_aperture_here) then
  call out_io (s_fatal$, r_name, 'CONFUSED APERTURE CALCULATION. PLEASE CONTACT HELP.')
  if (global_com%exit_on_error) call err_exit
  status = 1
endif

end function wall_intersection_func

end subroutine odeint_bmad

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!+
! Subroutine rk_adaptive_step (ele, param, orb, s_body, ds_try, ds_did, ds_next, err_flag, mat6, make_matrix)
!
! Private routine used by odeint_bmad.
! Not meant for general use
!-

recursive subroutine rk_adaptive_step (ele, param, orb, s_body, ds_try, ds12, ds_did, ds_next, err_flag, mat6, make_matrix)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, orb_new, end1, end2

real(rp), intent(inout) :: s_body
real(rp), intent(in)    :: ds_try
real(rp), intent(out)   :: ds_did, ds_next
real(rp), optional :: mat6(6,6)

real(rp) :: dr_ds(11), r_scal(10), rel_tol, abs_tol
real(rp) :: err_max, ds, ds_temp, p2, dmat(6,6)
real(rp) :: r_err(11), pol
real(rp) :: sqrt_N, rel_pc, t_new, ds12
real(rp), parameter :: safety = 0.9_rp, p_grow = -0.2_rp
real(rp), parameter :: p_shrink = -0.25_rp, err_con = 1.89e-4_rp
real(rp), parameter :: tiny = 1.0e-30_rp

integer ii

logical err_flag
logical, optional :: make_matrix
character(20), parameter :: r_name = 'rk_adaptive_step'

! Init. Note:
!   rel_tol: Same as eps for odeint scalled by sqrt(ds/(s2_body-s1_body))
!                   where ds is the step size for the current step. rel_tol
!                   sets the %error of the result
!   abs_tol: Sets the absolute error of the result

ds = ds_try
orb_new = orb

do
  call rk_step1 (ele, param, orb, dr_ds, s_body, ds, orb_new, r_err, err_flag, .false.)
  ! Can get errors due to step size too large. For example, for a field map that is only slightly larger than
  ! the aperture, a particle that is outside the aperture and outside of the fieldmap will generate an error.
  ! The solution is to just take a smaller step.
  if (err_flag) then
    if (ds < 1d-3) then
      ! call rk_step1 to generate an error message.
      call rk_step1 (ele, param, orb, dr_ds, s_body, ds, orb_new, r_err, err_flag, .true.)
      call out_io (s_error$, r_name, 'Problem with field calc. Tracked particle will be marked as dead.')
      orb%state = lost$
     return
    endif
    ds_temp = ds / 10
    if (ele%tracking_method == fixed_step_runge_kutta$) exit  

  elseif (ele%tracking_method == fixed_step_runge_kutta$) then
    exit  ! No adaptive step sizing.

  else
    sqrt_N = sqrt(abs(ds12/ds))  ! number of steps we would take with this ds
    rel_tol = bmad_com%rel_tol_adaptive_tracking / sqrt_N
    abs_tol = bmad_com%abs_tol_adaptive_tracking / sqrt_N
    pol = 1.0_rp  ! Spin scale
    r_scal(:) = abs([real(rp):: orb%vec(:), 1e-9_rp, pol, pol, pol]) + abs(ds*dr_ds(1:10)) + TINY
    err_max = maxval(abs(r_err(1:10)/(r_scal(:)*rel_tol + abs_tol)))
    if (err_max <=  1.0) exit
    ds_temp = safety * ds * (err_max**p_shrink)
  endif

  ds = sign(max(abs(ds_temp), 0.1_rp*abs(ds)), ds)

  if (s_body + ds == s_body) then
    err_flag = .true.
    call out_io (s_error$, r_name, 'STEPSIZE UNDERFLOW IN ELEMENT: ' // ele%name, &
                                   'AT (X,Y,Z) POSITION FROM ENTRANCE: \3F12.7\ ', &
                                   'TYPICALLY THIS IS DUE TO THE FIELD NOT OBEYING MAXWELL''S EQUATIONS.', &
                                   'OFTEN TIMES THE FIELD IS NOT EVEN CONTINUOUS!', &
                                   'THE PARTICLE WILL BE MARKED AS LOST.', &
                                   r_array = [orb%vec(1), orb%vec(3), s_body])
    orb%state = lost$
    return
  endif
end do

!

if (logic_option(.false., make_matrix)) then
  do ii = 1, 6
    end2 = orb
    end2%vec(6) = end2%vec(6)
    end2%vec(ii) = end2%vec(ii) + bmad_com%d_orb(ii)
    call adjust_start(orb, end2)
    call odeint_bmad (end2, ele, param, s_body, s_body+ds, err_flag)
    if (end2%state /= alive$ .or. err_flag) then
      call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING (+). MATRIX NOT CALCULATED FOR ELEMENT: ' // ele%name)
      return
    endif

    end1 = orb
    end1%vec(6) = end1%vec(6)
    end1%vec(ii) = end1%vec(ii) - bmad_com%d_orb(ii)
    call adjust_start(orb, end1)
    call odeint_bmad (end1, ele, param, s_body, s_body+ds, err_flag)
    if (end1%state /= alive$ .or. err_flag) then
      call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING (-). MATRIX NOT CALCULATED FOR ELEMENT: ' // ele%name)
      return
    endif

    dmat(1:6, ii) = (end2%vec - end1%vec) / (2.0_rp * bmad_com%d_orb(ii))
  enddo
  mat6 = matmul(dmat, mat6)
endif

!

if (ele%tracking_method == fixed_step_runge_kutta$) then
  ds_next = ele%value(ds_step$) * sign_of(ds)
elseif (err_max > err_con) then
  ds_next = safety*ds*(err_max**p_grow)
else
  ds_next = 5.0_rp * ds
end if

ds_did = ds
s_body = s_body + ds

orb_new%s = orb%s + ds * ele%orientation

orb = orb_new

!--------------------------------------------
contains

subroutine adjust_start (start0, start)

type (coord_struct) start0, start

call convert_pc_to (start%p0c * (1.0_rp + start%vec(6)), start%species, beta = start%beta)

if (start0%beta == 0) then
  start%t = start0%t + start%vec(5) / (c_light * start%beta)
else
  start%t = start0%t + start%vec(5) / (c_light * start%beta) - start0%vec(5) / (c_light * start0%beta)
endif

end subroutine adjust_start

end subroutine rk_adaptive_step

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine rk_step1 (ele, param, orb, dr_ds1, s_body, ds, orb_new, r_err, err, print_err)
!
! Fifth-order Cash–Karp Runge-Kutta step. Adapted from the routine rkqs from NR. Not meant for general use.
!
! Input:
!   ele         -- ele_struct: Elemnet being tracked through.
!   param       -- lat_param_struct: Some lattice parameters.
!   orb         -- coord_struct: Initial position.
!   s_body      -- real(rp): Initial longitudinal position.
!   ds          -- real(rp): Step to do.
!
! Output:
!   dr_ds1(11)  -- real(rp): Derivative vector at the initial position.
!                   See kick_vector_calc for details.
!   orb_new     -- coord_struct: Ending coords.
!   r_err(11)   -- real(rp): Error estimates.
!   err         -- logical: Set True if there is an error.
!   print_err   -- logical, optional: Print message if there is an error? Default is True.
!-

recursive subroutine rk_step1 (ele, param, orb, dr_ds1, s_body, ds, orb_new, r_err, err, print_err)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orb, orb_new, orb_temp(5)

real(rp) :: dr_ds1(11)
real(rp), intent(in) :: s_body, ds
real(rp), intent(out) :: r_err(11)
real(rp) :: dr_ds2(11), dr_ds3(11), dr_ds4(11), dr_ds5(11), dr_ds6(11)
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

real(rp) quat(0:3), s1
logical err
logical, optional :: print_err

!

call kick_vector_calc (ele, param, s_body, orb, dr_ds1, err, print_err)
if (err) return

!

s1 = s_body + a2*ds
call transfer_this_orbit (orb, s1, b21*ds*dr_ds1, orb_temp(1), .true.)
call kick_vector_calc(ele, param, s1, orb_temp(1), dr_ds2, err, print_err)
if (err) return

s1 = s_body + a3*ds
call transfer_this_orbit (orb, s1, ds*(b31*dr_ds1 + b32*dr_ds2), orb_temp(2), .true.)
call kick_vector_calc(ele, param, s1, orb_temp(2), dr_ds3, err, print_err)
if (err) return

s1 = s_body + a4*ds
call transfer_this_orbit (orb, s1, ds*(b41*dr_ds1 + b42*dr_ds2 + b43*dr_ds3), orb_temp(3), .true.)
call kick_vector_calc(ele, param, s1, orb_temp(3), dr_ds4, err, print_err)
if (err) return

s1 = s_body + a5*ds
call transfer_this_orbit (orb, s1, ds*(b51*dr_ds1 + b52*dr_ds2 + b53*dr_ds3 + b54*dr_ds4), orb_temp(4), .true.)
call kick_vector_calc(ele, param, s1, orb_temp(4), dr_ds5, err, print_err)
if (err) return

s1 = s_body + a6*ds
call transfer_this_orbit (orb, s1, ds*(b61*dr_ds1 + b62*dr_ds2 + b63*dr_ds3 + b64*dr_ds4 + b65*dr_ds5), orb_temp(5), .true.)
call kick_vector_calc(ele, param, s1, orb_temp(5), dr_ds6, err, print_err)
if (err) return

call transfer_this_orbit (orb, s1, ds*(c1*dr_ds1 + c3*dr_ds3 + c4*dr_ds4 + c6*dr_ds6), orb_new, .false.)

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  quat =          omega_to_quat(ds*c1*dr_ds1(8:10))
  quat = quat_mul(omega_to_quat(ds*c3*dr_ds3(8:10)), quat)
  quat = quat_mul(omega_to_quat(ds*c4*dr_ds4(8:10)), quat)
  quat = quat_mul(omega_to_quat(ds*c6*dr_ds6(8:10)), quat)
  orb_new%spin = quat_rotate(quat, orb%spin)
endif

r_err=ds*(dc1*dr_ds1 + dc3*dr_ds3 + dc4*dr_ds4 + dc5*dr_ds5 + dc6*dr_ds6)

!----------------------------------------------------------
contains

subroutine transfer_this_orbit (orb_in, s1, dvec, orb_out, all_transfer)

type (coord_struct) orb_in, orb_out
real(rp) dvec(11), t_temp, s1
logical all_transfer

! Ignore dvec(5) and compute the change in z from the change in t to keep z and t consistent.

if (all_transfer) orb_out = orb_in
orb_out%vec = orb_in%vec + dvec(1:6)
orb_out%t = orb_in%t + dvec(7)
orb_out%s = s1 + ele%s_start

if (dvec(6) /= 0) then
  call convert_pc_to (orb_out%p0c * (1.0_rp + orb_out%vec(6)), param%particle, beta = orb_out%beta)
endif

orb_out%vec(5) = orb_in%vec(5) * orb_out%beta / orb_in%beta + c_light * orb_out%beta * (dvec(11) - dvec(7))

end subroutine transfer_this_orbit

end subroutine rk_step1

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine kick_vector_calc (ele, param, s_rel, orbit, dr_ds, field, err, print_err)
!
! Subroutine to calculate the dr/ds "kick vector" where
!     r = [x, p_x, y, p_y, z, p_z, t, spin_x,y,z]
!
! Remember: In order to simplify the calculation, in the body of any element, P0 is taken to be 
! the P0 at the exit end of the element.
!
!   dr(1)/ds = dx/ds = dx/dt * dt/ds
!   where:
!     dx/dt = v_x = p_x / (1 + p_z)
!     dt/ds = (1 + g*x) / v_s
!     g = 1/rho, rho = bending radius (nonzero only in a dipole)
!
!   dr(2)/ds = dp_x/ds = dP_x/dt * dt/ds / P0 + g_x * P_z
!   where:
!     dP_x/dt = EM_Force_x
!     g_x = bending in x-plane.
!
!   dr(3)/ds = dy/ds = dy/dt * dt/ds
!   where:
!     dy/dt = v_x 
! 
!   dr(4)/ds = dp_y/ds = dP_y/dt * ds/dt / P0 + g_y * P_z
!   where:
!     dP_y/dt = EM_Force_y
!     g_y = bending in y-plane.
!
!   NOTE: dr(5)/ds IS IGNORED WHEN CALCULATING Z. SEE TRANSFER_THIS_ORBIT ABOVE.
!   dr(5)/ds = dz/ds = beta * c_light * [dt/ds(ref) - dt/ds] + dbeta/ds * c_light * [t(ref) - t]
!                    = beta * c_light * [dt/ds(ref) - dt/ds] + dbeta/ds * vec(5) / beta
!   where:
!     dt/ds(ref) = 1 / beta(ref)
!
!   dr(6)/ds = dp_z/ds = d(EM_Force dot v_hat) * dt/ds / P0
!   where:
!      v_hat = velocity normalized to 1.
!
!   dr(7)/ds = dt/ds
!
!   dr(8:10)/ds = Spin omega vector
!
!   dr(11)/ds = dt_ref/ds
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
!   dr_ds(11)   -- real(rp): Kick vector.
!   field       -- em_field_struct: Local field.
!   err         -- Logical: Set True if there is an error.
!-

recursive subroutine kick_vector_calc (ele, param, s_body, orbit, dr_ds, err, print_err)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
type (lat_param_struct) param
type (em_field_struct) field
type (coord_struct) orbit, orbit2

real(rp), intent(in) :: s_body
real(rp), intent(out) :: dr_ds(11)
real(rp) g_bend, dt_ds, dp_ds, dbeta_ds, s2
real(rp) vel(3), E_force(3), B_force(3), w_mat(3,3), ww(3,3), ww_inv(3,3), r_vec(3)
real(rp) e_tot, dt_ds_ref, p0, beta0, v2, pz_p0, mc2, delta, dh_bend, rel_pc

integer rel_dir, i

logical :: err, calc_dfield
logical, optional :: print_err

character(*), parameter :: r_name = 'kick_vector_calc'

! Init.
! The effective reference velocity is different from the velocity of the reference particle for wigglers where the 
! reference particle is not traveling along the reference line and in elements where the reference velocity is not constant.

! Note: There is a potential problem with RF and e_gun elements when calculating beta0 when there is slicing and where
! the particle is non-relativistic. To avoid z-shifts with slicing, use the lord delta_ref_time.
! Exception: pipe element.

if ((ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) .and. ele%key /= pipe$) then
  ele0 => pointer_to_super_lord(ele)
else
  ele0 => ele
endif

if (ele0%value(delta_ref_time$) == 0 .or. ele0%key == patch$) then
  beta0 = ele0%value(p0c$) / ele0%value(e_tot$) ! Singular case can happen when parsing lattice. 
else
  beta0 = ele0%value(l$) / (c_light * ele0%value(delta_ref_time$))
endif

err = .true.
rel_dir = ele%orientation * orbit%direction
dt_ds_ref = ele%orientation / (beta0 * c_light)
p0 = ele%value(p0c$) / c_light
e_tot = orbit%p0c * (1.0_rp + orbit%vec(6)) / orbit%beta

! Calculate the field. 
! Note: Field is in frame of element. When ele%orientation = -1 => +z in -s direction.

calc_dfield = runge_kutta_com%calc_field_derivatives

if (ele%key == patch$) then
  ! Particle going towards an end uses the coordinate of that end.
  ! But the field is specified in the exit end coords so must transform if particle is traveling towards the entrance end.
  if (rel_dir == 1) then
    call em_field_calc (ele, param, s_body, orbit, .true., field, calc_dfield, err, print_err = print_err)
  else
    call floor_angles_to_w_mat(ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), w_mat = ww, w_mat_inv = ww_inv)
    orbit2 = orbit
    r_vec = matmul(ww_inv, [orbit%vec(1), orbit%vec(3), s_body] - [ele%value(x_offset$), ele%value(y_offset$), ele%value(z_offset$)])
    orbit2%vec(1:3:2) = r_vec(1:2)
    call em_field_calc (ele, param, r_vec(3)+ele%value(l$), orbit2, .true., field, calc_dfield, err, print_err = print_err)
    field%B = matmul(ww, field%B)
    field%E = matmul(ww, field%E)
    if (calc_dfield) then
      do i = 1, 3
        field%dB(:,i) = matmul(ww, field%dB(:,i))
        field%dE(:,i) = matmul(ww, field%dE(:,i))
      enddo
    endif
  endif

else
  call em_field_calc (ele, param, s_body, orbit, .true., field, calc_dfield, err, print_err = print_err)
endif

if (err) return

! Bend factor

delta = orbit%vec(6)
rel_pc = 1.0_rp + delta
vel(1:2) = [orbit%vec(2), orbit%vec(4)] / rel_pc
v2 = vel(1)**2 + vel(2)**2
if (v2 > 0.99999999_rp) return
vel = orbit%beta * c_light * [vel(1), vel(2), sqrt(1.0_rp - v2) * rel_dir]
E_force = charge_of(orbit%species) * field%E
B_force = charge_of(orbit%species) * cross_product(vel, field%B)


if (ele%key == sbend$ .or. ele%key == rf_bend$) then
  g_bend = ele%value(g$)
  dh_bend = orbit%vec(1) * g_bend
else
  g_bend = 0
  dh_bend = 0 ! Longitudinal distance deviation per unit s-distance. Equal to 0 except when off axis in a bend.
endif

mc2 = mass_of(orbit%species)
dt_ds = rel_dir * (1.0_rp + dh_bend) / abs(vel(3))
dp_ds = dot_product(E_force, vel) * dt_ds / (orbit%beta * c_light)
dbeta_ds = mc2**2 * dp_ds * c_light / e_tot**3
pz_p0 = rel_pc * rel_dir * abs(vel(3)) / (orbit%beta * c_light)  ! Pz / P0

dr_ds(1) = vel(1) * dt_ds
dr_ds(2) = (E_force(1) + B_force(1)) * dt_ds / p0 + g_bend * pz_p0
dr_ds(3) = vel(2) * dt_ds
dr_ds(4) = (E_force(2) + B_force(2)) * dt_ds / p0
! Modified dr_ds(5) to reduce roundoff. The old formula is:
! dr_ds(5) = orbit%beta * c_light * (dt_ds_ref - dt_ds) + dbeta_ds * orbit%vec(5) / orbit%beta
! Note that beta0 in a wiggler is not the velocity one would calculated based upon the ref energy.
dr_ds(5) = rel_dir * (orbit%beta / beta0 - 1) + &
           rel_dir * (sqrt_one(-v2) - dh_bend) / sqrt(1-v2) + dbeta_ds * orbit%vec(5) / orbit%beta
dr_ds(6) = dp_ds / p0
dr_ds(7) = dt_ds
dr_ds(11) = dt_ds_ref

if (bmad_com%spin_tracking_on .and. ele%spin_tracking_method == tracking$) then
  ! dr_ds(8:10) = Omega/v_z
  dr_ds(8:10) = rel_dir * (1.0_rp + dh_bend) * spin_omega (field, orbit, rel_dir)
  if (ele%key == sbend$ .or. ele%key == rf_bend$) dr_ds(8:10) = dr_ds(8:10) + [0.0_rp, g_bend, 0.0_rp]
else
  dr_ds(8:10) = 0
endif

err = .false.

end subroutine kick_vector_calc

end module
