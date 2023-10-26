!+                                
! Subroutine ele_to_fibre (ele, ptc_fibre, param, use_offsets, err_flag, integ_order, steps, for_layout, ref_in)
!
! Routine to convert a Bmad element to a PTC fibre element.
!
! Note: You need to call set_ptc before using this routine.
!
! Note: If ele contains a gen_grad_map, this routine may not be called in between calls to 
!   FPP alloc and kill since the setting up of the PTC pancake uses FPP.
!
! Input:
!   ele             -- Ele_struct: Bmad element.
!   param           -- lat_param_struct: 
!   use_offsets     -- Logical: Does ptc_fibre include element offsets, pitches and tilt?
!   integ_order     -- Integer, optional: Order for the 
!                        sympletic integrator. Possibilities are: 2, 4, or 6
!                        Overrides ele%value(integrator_order$).
!                        default = 2 (if not set with set_ptc).
!   steps           -- Integer, optional: Number of integration steps.
!                        Overrides ele%value(ds_step$).
!   for_layout      -- Logical, optional: If True then fibre will be put in the PTC layout.
!                        Default is False.
!   ref_in          -- coord_struct, optional: Particle to be tracked. ref_particle$, electron$, etc.
!                        This argument should only be present when the fibre is not to be put in a layout.
!
! Output:
!   err_flag        -- logical: Set True if setup OK. False otherwise.
!   ptc_fibre       -- Fibre: PTC fibre element.
!+

subroutine ele_to_fibre (ele, ptc_fibre, param, use_offsets, err_flag, integ_order, steps, for_layout, ref_in)

use ptc_interface_mod, dummy => ele_to_fibre
use madx_ptc_module, pi_dum => pi, pi2_dum => twopi
use dabnew_pancake

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param
type (coord_struct), optional :: ref_in
type (ele_struct), pointer :: field_ele, ele2
type (cartesian_map_term1_struct), pointer :: wt
type (cartesian_map_struct), pointer :: cm
type (cylindrical_map_struct), pointer :: cy
type (gen_grad_map_struct), pointer :: gg_map
type (gen_grad1_struct), pointer :: gg
type (em_taylor_struct), target :: em_taylor(3)
type (em_taylor_term_struct), pointer :: tm
type (fibre), pointer :: ptc_fibre
type (keywords) ptc_key, key2
type (work) energy_work
type (tree_element), pointer :: arbre(:)
type (c_damap) ptc_c_damap
type (real_8) ptc_re8(6)
type (taylor) ptc_taylor

real(rp) leng, hk, vk, s_rel, z_patch, phi_tot, norm, rel_charge, kl(0:n_pole_maxx), t(0:n_pole_maxx)
real(rp) dx, dy, cos_t, sin_t, coef, coef_e, coef_b, kick_magnitude, ap_lim(2), ap_dxy(2), e1, e2
real(rp) beta0, beta1, ref0(6), ref1(6), fh, dz_offset, ff, z0, z, dz_step, vec(6)
real(rp), pointer :: val(:)
real(rp), target :: value0(num_ele_attrib$)
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) ld, hd, lc, hc, angc, xc, dc
real(rp), parameter :: dz_pan7(0:6) = [0.0_rp, 1.0_rp/9.0_rp, 1.0_rp/6.0_rp, 1.0_rp/3.0_rp, 1.0_rp/2.0_rp, 2.0_rp/3.0_rp, 5.0_rp/6.0_rp]

complex(rp) k_0

integer, optional :: integ_order, steps
integer i, ii, j, k, m, n, key, n_term, exception, ix, met, net, ap_type, ap_pos, ns, n_map
integer np, max_order, ix_pole_max, nn, n_period, icoef, n_step, n_pan, field(8)
integer, allocatable :: pancake_field(:,:)

logical use_offsets, err_flag, kill_spin_fringe, onemap, found, is_planar_wiggler, use_taylor, done_it, change
logical, optional :: for_layout

character(24) pancake_name
character(*), parameter :: r_name = 'ele_to_fibre'

!

err_flag = .true.
val => ele%value
key = ele%key

if (.not. ele%is_on) then
  val => value0
  val = ele%value

  select case (ele%key)
  case (sbend$)
    val = 0
    val(l$)            = ele%value(l$)
    val(g$)            = ele%value(g$)
    val(dg$)           = -ele%value(g$)
    val(angle$)        = ele%value(angle$)
    val(rho$)          = ele%value(rho$)
    val(ref_tilt_tot$) = val(ref_tilt_tot$)

  case (lcavity$, rfcavity$)
    val(voltage$)  = 0
    val(gradient$) = 0

  case (taylor$)
    key     = drift$
    val(l$) = 0

  case default
    key             = drift$
    val             = 0
    val(l$)         = ele%value(l$)
    val(ds_step$)   = val(l$)
    val(num_steps$) = 1
  end select
endif

!

ele2 => pointer_to_field_ele(ele, 1, s_rel)

n_map = 0
if (associated(ele2%cylindrical_map)) n_map = n_map + 1
if (associated(ele2%cartesian_map)) n_map = n_map + 1
if (associated(ele2%gen_grad_map)) n_map = n_map + 1

if (n_map > 1) then
  call out_io (s_fatal$, r_name, 'PTC TRACKING IS ONLY ABLE TO HANDLE A SINGLE FIELD MAP IN AN ELEMENT.', &
                                 'ELEMENT HAS MULTIPLE FIELD MAPS: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

use_taylor = (n_map == 0 .and. (key == wiggler$ .or. key == undulator$) .and. ele%field_calc == helical_model$)
if (use_taylor) key = match$

! 

call zero_key(ptc_key)  ! init key

select case (ele%ptc_integration_type)
case (drift_kick$);  ptc_key%model = 'DRIFT_KICK'
case (matrix_kick$); ptc_key%model = 'MATRIX_KICK'
case (ripken_kick$); ptc_key%model = 'DELTA_MATRIX_KICK'
end select

if (key == sbend$ .and. val(angle$) == 0 .and. ptc_key%model /= 'DRIFT_KICK') then
  ptc_key%model = 'DRIFT_KICK'
  ! Only need to issue a warning if K1 is nonzero.
  !if (val(k1$) /= 0) call out_io (s_warn$, r_name, &
  !          'BEND WITH ZERO BENDING ANGLE WILL USE PTC_INTEGRATION_TYPE OF DRIFT_KICK: ' // ele%name)
endif

if (key == sbend$ .and. val(g$) + val(dg$) == 0 .and. ptc_key%model /= 'DRIFT_KICK') then
  ptc_key%model = 'DRIFT_KICK'
  ! Only need to issue a warning if K1 is nonzero.
  !if (val(k1$) /= 0) call out_io (s_warn$, r_name, &
  !          'BEND WITH ZERO NET BENDING FIELD WILL USE PTC_INTEGRATION_TYPE OF DRIFT_KICK: ' // ele%name)
endif

if (present(ref_in)) then
  rel_charge = charge_of(ref_in%species) / charge_of(param%particle)
else
  rel_charge = charge_of(default_tracking_species(param)) / charge_of(param%particle)
endif

leng = val(l$)
if (use_offsets .and. key == sbend$) then
  ptc_key%tiltd = val(ref_tilt_tot$)
else
  ptc_key%tiltd = 0
endif

ptc_key%method = nint(val(integrator_order$))
if (ptc_key%method == 0) ptc_key%method = bmad_com%default_integ_order 
if (present(integ_order)) ptc_key%method = integ_order

if (present(steps)) then
  ptc_key%nstep = steps
elseif (leng == 0) then
  ptc_key%nstep = 1
elseif (key == taylor$ .or. key == match$ .or. key == multipole$ .or. &
                                key == ab_multipole$ .or. key == patch$) then
  ptc_key%nstep = 1
  leng = 0  ! Problem is that PTC will not ignore the length in tracking which is different from the Bmad convention.
else
  if (val(ds_step$) == 0) then
    call out_io (s_fatal$, r_name, 'DS_STEP IS ZERO FOR ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
  endif

  ptc_key%nstep = nint(abs(leng) / val(ds_step$))
  if (ptc_key%nstep == 0) ptc_key%nstep = 1
endif

ptc_key%list%name = ele%name
ptc_key%list%l    = leng

! Fringes

ix = both_ends$
if (attribute_index(ele, 'FRINGE_AT') > 0) ix = nint(ele%value(fringe_at$))
kill_spin_fringe = is_false(ele%value(spin_fringe_on$))

ptc_key%list%kill_ent_fringe = (ix == exit_end$ .or. ix == no_end$)
ptc_key%list%kill_exi_fringe = (ix == entrance_end$ .or. ix == no_end$)
ptc_key%list%kill_ent_spin = (ix == exit_end$ .or. ix == no_end$ .or. kill_spin_fringe)
ptc_key%list%kill_exi_spin = (ix == entrance_end$ .or. ix == no_end$ .or. kill_spin_fringe)

!

if (key == sbend$ .and. val(l$) == 0) key = kicker$
if (ele2%field_calc == fieldmap$ .and. ele2%tracking_method /= bmad_standard$) key = wiggler$

select case (key)

!------------------------------
case (ac_kicker$)
  if (.not. allocated(ele%ac_kick%frequency)) then
    call out_io (s_error$, r_name, 'AC_KICKER CANNOT BE TRANSLATED TO PTC IF THE ELEMENT USES (TIME, AMP)', &
                                   'TO SPECIFY THE FIELD TIME DEPENDENCE: ' // ele%name)
    return
  endif

  if (size(ele%ac_kick%frequency) > 1) then
    call out_io (s_error$, r_name, 'AC_KICKER CANNOT BE TRANSLATED TO PTC IF MORE THAN FREQUENCY IS USED. ' // ele%name)
    return
  endif

  ptc_key%magnet = 'rfcavity'
  ptc_key%list%n_bessel = 0
  ptc_key%list%permfringe = 0
  ptc_key%list%cavity_totalpath = 0
  ptc_key%list%freq0 = -ele%ac_kick%frequency(1)%f
  ptc_key%list%delta_e = 0     ! For radiation calc.
  ptc_key%list%lag = twopi * (val(phi0_multipass$) + ele%ac_kick%frequency(1)%phi + val(t_offset$)*ptc_key%list%freq0)

!------------------------------
case (crab_cavity$)
  if (val(rf_frequency$) == 0) then
    call out_io (s_fatal$, r_name, 'RF FREQUENCY IS ZERO FOR: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  if (ele%value(crab_x1$) /= 0 .or. ele%value(crab_x2$) /= 0) then
    call out_io (s_warn$, r_name, 'CRAB_X1, CRAB_X2, etc. DO NOT HAVE A PTC EQUIVALENT!')
  endif

  ptc_key%magnet = 'rfcavity'
  ptc_key%list%n_bessel = 0
  !!ptc_key%list%volt = 1d-6 * e_accel_field(ele, voltage$)
  ptc_key%list%permfringe = 0
  ptc_key%list%cavity_totalpath = 0
  ptc_key%list%freq0 = val(rf_frequency$)
  phi_tot = val(phi0$) + val(phi0_multipass$) + 0.25_rp 
  ptc_key%list%delta_e = 0     ! For radiation calc.
  ptc_key%list%lag = twopi * phi_tot

!------------------------------
case (drift$, rcollimator$, ecollimator$, monitor$, instrument$, pipe$) 
  if (val(hkick$) == 0 .and. val(vkick$) == 0) then
    ptc_key%magnet = 'drift'
  else
    ptc_key%magnet = 'quadrupole'
  endif

!------------------------------
case (quadrupole$) 
  ptc_key%magnet = 'quadrupole'
  ptc_key%list%usethin = .false.  ! So zero length element is not treated as a multipole

!------------------------------
case (sad_mult$)
  if (val(l$) == 0) then
    ptc_key%magnet = 'multipole'  ! No Bz field if zero length.
  else
    ptc_key%magnet = 'solenoid'
    ptc_key%list%bsol = val(ks$)
    ! PTC tracking uses matrix-kick where the "matrix" step uses the solenoid component and the Kick step
    ! uses the quadrupole (and higher order) components. SAD and Bmad combine the quadrupole component
    ! in the Matrix step. Thus PTC may need a smaller step size.
    if (.not. present(integ_order) .and. .not. present(steps)) then
      call multipole_ele_to_kt (ele, .false., ix_pole_max, kl, t, magnetic$)
      ptc_key%nstep = max(ptc_key%nstep, nint(ele%value(l$) * abs(kl(1)) / (ele%value(eps_step_scale$) * ptc_com%cut_factor)))
      ptc_key%method = 2
      if (ptc_key%nstep > 18) then
        ptc_key%nstep = nint(ptc_key%nstep / 7.0)
        ptc_key%method = 6
      elseif (ptc_key%nstep > 4) then
        ptc_key%nstep = nint(ptc_key%nstep / 3.0)
        ptc_key%method = 4
      endif
      switch_to_drift_kick = .false.
      call change_method_in_create_fibre(ptc_key, 1, change)
    endif
  endif

!------------------------------
case (sbend$) 
  ! PTC does not consider a finite e1/e2 part of the fringe so must zero e1/e2 if needed.
  ix = nint(ele%value(fringe_type$))

  e1 = ele%value(e1$)
  if (ptc_key%list%kill_ent_fringe .or. ix == none$) e1 = 0

  e2 = ele%value(e2$)
  if (ptc_key%list%kill_exi_fringe .or. ix == none$) e2 = 0

  ptc_key%magnet = 'sbend'
  ptc_key%list%t1   = e1
  ptc_key%list%t2   = e2

  ptc_key%list%b0   = ele%value(g$) * leng
  if (abs(ptc_key%list%b0) < 1e-30_rp) ptc_key%list%b0 = 0  ! To get around PTC bug.
  ptc_key%list%hgap = ele%value(hgap$)
  ptc_key%list%fint = ele%value(fint$)
  ptc_key%list%hgap2 = ele%value(hgapx$)
  ptc_key%list%fint2 = ele%value(fintx$)

  if (nint(ele%value(ptc_field_geometry$)) == straight$) then
    ptc_key%magnet = 'wedgrbend'
    ptc_key%list%t1   = e1 - ele%value(angle$)/2
    ptc_key%list%t2   = e2 - ele%value(angle$)/2
  endif

!------------------------------
case (sextupole$)
  ptc_key%magnet = 'sextupole'
  ptc_key%list%usethin = .false.  ! So zero length element is not treated as a multipole

!------------------------------
case (octupole$)
  ptc_key%magnet = 'octupole'
  ptc_key%list%usethin = .false.  ! So zero length element is not treated as a multipole

!------------------------------
case (solenoid$)
  ptc_key%magnet = 'solenoid'
  ptc_key%list%bsol = val(ks$)
  ptc_key%list%usethin = .false.  ! So zero length element is not treated as a multipole

!------------------------------
case (sol_quad$)
  ptc_key%magnet = 'solenoid'
  ptc_key%list%bsol = val(ks$)
  ptc_key%list%usethin = .false.  ! So zero length element is not treated as a multipole

!------------------------------
case (marker$, detector$, fork$, photon_fork$, beginning_ele$, patch$, floor_shift$, fiducial$, taylor$, match$)
  ptc_key%magnet = 'marker'
  ptc_key%nstep = 1

!------------------------------
case (kicker$, hkicker$, vkicker$)
  ptc_key%magnet = 'kicker'

case (gkicker$)
  ptc_key%magnet = 'marker'
  ptc_key%nstep = 1
  

!------------------------------
case (rfcavity$, lcavity$)
  if (ele%value(rf_frequency$) == 0) then
    call out_io (s_fatal$, r_name, 'RF FREQUENCY IS ZERO FOR: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  ix = nint(ele%value(cavity_type$))
  if (ele%tracking_method == bmad_standard$) ix = ptc_standard$

  select case (ix)
  case (traveling_wave$)
    ptc_key%magnet = 'twcavity'
    ptc_key%list%volt = 1d-6 * e_accel_field(ele, voltage$)
    ptc_key%list%cavity_totalpath = 1  ! 
  case (standing_wave$)
    ptc_key%magnet = 'rfcavity'
    ptc_key%list%cavity_totalpath = 1  ! 
    if (ptc_key%nstep == 1) ptc_key%nstep = 5  ! Avoid bug with nstep = 1.
    ptc_key%list%n_bessel = -1   ! Pillbox cavity.
    ptc_key%list%volt = 2d-6 * e_accel_field(ele, voltage$)
    ptc_key%list%l = val(l_active$)
    ptc_key%list%h1 = (val(l$) - val(l_active$)) / 2
    ptc_key%list%h2 = (val(l$) - val(l_active$)) / 2
  case (ptc_standard$)
    ptc_key%magnet = 'rfcavity'
    ptc_key%list%volt = 1d-6 * e_accel_field(ele, voltage$)
    ptc_key%list%n_bessel = 0
    ptc_key%list%permfringe = 0
    ptc_key%list%cavity_totalpath = 0
  end select

  ptc_key%list%freq0 = val(rf_frequency$)
  phi_tot = val(phi0$) + val(phi0_multipass$) + val(phi0_err$) + val(phi0_autoscale$)

  if (key == lcavity$) then
    ptc_key%list%lag = pi / 2 - twopi * phi_tot
  else
    ptc_key%list%lag = twopi * phi_tot
  endif

  ptc_key%list%delta_e = 0     ! For radiation calc.

!------------------------------
! ptc elsep cannot do spin tracking so use general electrostatic element instead.
case (elseparator$)
!  ptc_key%magnet = 'elseparator'
!  hk = val(hkick$) / leng
!  vk = val(vkick$) / leng
!  if (hk == 0 .and. vk == 0) then
!    ptc_key%tiltd = 0
!  else
!    if (param%particle < 0) then
!      hk = -hk
!      vk = -vk
!    endif
!    ptc_key%tiltd = -atan2 (hk, vk) + ele%value(tilt_tot$)
!  endif
!  ptc_key%list%volt = 1d-6 * ele%value(e_tot$) * sqrt(hk**2 + vk**2)

!------------------------------
case (ab_multipole$, multipole$)
  ptc_key%magnet = 'multipole'

!------------------------------
! beambeam element in PTC is a special drift that must be setup after the integration 
! node array of the fibre is created.

case (beambeam$)
  ptc_key%magnet = 'drift'

!------------------------------
case (wiggler$, undulator$)
  ptc_key%magnet = 'wiggler'

!------------------------------
case default
  call out_io (s_fatal$, r_name, 'CONVERSION TO PTC NOT IMPLEMENTED FOR ELEMENTS OF TYPE ' // trim(key_name(ele%key)), &
                                 'FOR ELEMENT: ' // trim(ele%name))
  if (global_com%exit_on_error) call err_exit
end select

!------------------------------
! Fringe

if (key == sbend$ .and. ele%value(l$) /= 0) then

  ix = nint(ele%value(ptc_fringe_geometry$))
  ptc_key%list%bend_fringe = (ix == x_invariant$)

  ix = nint(ele%value(fringe_type$))
  select case (ix)
  case (none$)
    ptc_key%list%permfringe = 0
  case (basic_bend$, linear_edge$)
    ptc_key%list%permfringe = 0
  case (full$)
    ptc_key%list%permfringe = 1
  case (hard_edge_only$)
    ptc_key%list%permfringe = 1
  case (soft_edge_only$)
    ptc_key%list%permfringe = 2
  case (sad_full$)
    ptc_key%list%permfringe = 3
  end select

elseif (attribute_index(ele, 'FRINGE_TYPE') > 0) then  ! If fringe_type is a valid attribute

  ptc_key%list%bend_fringe = .false.

  ix = nint(ele%value(fringe_type$))
  select case (ix)
  case (none$)
    ptc_key%list%permfringe = 0
  case (hard_edge_only$)
    ptc_key%list%permfringe = 1
  case (soft_edge_only$)
    ptc_key%list%permfringe = 2
  case (full$)
    ptc_key%list%permfringe = 3
  end select

  if (key == sad_mult$ .and. ele%value(l$) == 0) ptc_key%list%permfringe = 0
endif

! Electric fields present? Everything is an sbend!

if (key /= multipole$ .and. (associated(ele%a_pole_elec) .or. key == elseparator$)) then
  ptc_key%magnet = 'sbend'
  ptc_key%model = 'DRIFT_KICK'   ! PTC demands this.
  ptc_key%exact = .true.  ! PTC does not implement a non-exact model when there are electric fields.
  SOLVE_ELECTRIC = .true.
  if (leng == 0) then
    call out_io (s_fatal$, r_name, 'ZERO LENGTH ELEMENT WITH AN ELECTRIC FIELD NOT ALLOWED IN PTC: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

! Magnetic multipole components

call ele_to_ptc_magnetic_an_bn (ele, ptc_key%list%k, ptc_key%list%ks, ptc_key%list%nmul)

! Cylindrical_map

if (associated(ele2%cylindrical_map) .and. ele2%field_calc == fieldmap$) then
  PUT_A_ABELL = 1
  ptc_key%magnet = 'abell_dragt'
  n_abell = 0
  do i = 1, size(ele%cylindrical_map)
    n_abell = max(2, n_abell, size(ele%cylindrical_map(i)%ptr%term))
  enddo
  m_abell = maxval(ele%cylindrical_map%m)
endif

! Gen_grad_map

if (associated(ele2%gen_grad_map) .and. ele2%field_calc == fieldmap$) then

  if (nint(ele2%value(integrator_order$)) /= 6) ele2%value(integrator_order$) = 4
  ptc_key%method = nint(ele2%value(integrator_order$))

  if (nint(ele2%value(integrator_order$)) == 4) then
    n_step = max(nint(ele%value(l$) / (4.0_rp*ele2%value(ds_step$))), 1)
    n_pan = n_step * 4 + 1
  else
    n_step = max(nint(ele%value(l$) / (7.0_rp*ele2%value(ds_step$))), 1)
    n_pan = n_step * 7 + 1
    dz_step = ele%value(l$) / n_step
  endif

  gg_map => ele2%gen_grad_map(1)

  if (key == sbend$ .and. ele%value(g$) /= 0) then
    ld = ele%value(l$)
    hd = ele%value(g$)
    lc = ele%value(l$)   ! Integration length

    if (gg_map%curved_ref_frame) then
      hc = ele%value(g$) ! pancake curvature = ele curvature
      angc = (ele%value(angle$) - (n_pan-1) * ele2%value(ds_step$) * ele%value(g$)) / 2
      xc = gg_map%r0(1) * cos(ele%value(angle$)/2) - ele%value(rho$) * (1 - cos(angc))
      dc = gg_map%r0(1) * sin(ele%value(angle$)/2) + ele%value(rho$) * sin(angc)
    else
      hc = 0.d0     ! pancake curvature
      angc = ele%value(angle$) / 2
      xc = gg_map%r0(1) + ele%value(rho$) * (1 - cos(ele%value(angle$)/2))
      dc = (ele%value(l$) - (n_pan-1) * ele2%value(ds_step$)) / 2
    endif

  else
    ld = ele%value(l$)
    lc = ele%value(l$)

    hd = 0
    hc = 0                ! pancake curvature

    angc = 0
    xc = gg_map%r0(1)
    dc = 0
  endif

  if (gg_map%field_type /= magnetic$) then
    call out_io (s_fatal$, r_name, 'FIELD TYPE MUST BE MAGNETIC FOR A TAYLOR_FIELD IF USED WITH PTC.', &
                                   'TAYLOR_FIELD USED IN ELEMENT: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  n = max(len_trim(gg_map%file)-24, str_last_in_set(gg_map%file, '/'))
  m = min(n+24, len(gg_map%file))
  pancake_name = gg_map%file(n+1:m)
  call str_substitute(pancake_name, ':', '-')  ! Make valid file name

  ! Note: True => Not canonical tracking.
  call set_pancake_constants(n_pan, angc, xc, dc, gg_map%r0(2), hc, lc, hd, ld, .true., pancake_name)

  max_order = 0
  do i = 1, size(gg_map%gg)
    max_order = max(max_order, gg_map%gg(i)%m+gg_map%gg(i)%n_deriv_max)
  enddo

  ptc_key%magnet = 'PANCAKEBMADZERO'

else
  ! This needed since corresponding create_fibre_append dummy arg is an assumed shape array.
  allocate (pancake_field(0,0))  
endif


! Check number of slices and integrator order.
! Wigglers are special since they do not have a uniform field so ignore these.

select case (ele%key)
case (wiggler$, undulator$, gkicker$, rfcavity$, lcavity$)
case default
  switch_to_drift_kick = .false.
  key2 = ptc_key  ! Don't want to change ptc_key
  call change_method_in_create_fibre(key2, 1, change)
  if (change .and. ptc_com%print_step_warning) then
    call out_io (s_warn$, r_name, 'Number of integration steps looks large for element: ' // ele%name, &
              ' Present num_steps: ' // int_str(ptc_key%nstep) // ' at integrator_order: ' // int_str(ptc_key%method), &
              ' Suggested num_steps: ' // int_str(key2%nstep) // ' at integrator_order: ' // int_str(key2%method), &
              ' [Note: To turn off this message set ptc_com%print_step_warning = False.]')
  endif
end select

!--------------------------------------------
! Create ptc_fibre
! The EXCEPTION argument is an error_flag. Set to 1 if there is an error. Never reset.

! The integration node array of the fibre is not setup if for_layout = True.
! [In this case the setup is done on the entire layout later.]
! Therefore only do the beambeam setup (which needs this array to exist) if for_layout = False.

call set_ptc_quiet(12, set$, n)

if (logic_option(.false., for_layout)) then
  call create_fibre_append (.true., m_u%end, ptc_key, EXCEPTION, bri = pancake_field)   ! ptc routine
  ptc_fibre => m_u%end%end

  if (present (ref_in)) then
    call out_io (s_fatal$, r_name, 'REF_IN ARGUMENT SHOULD NOT BE PRESENT WHEN FOR_LAYOUT IS TRUE!')
    if (global_com%exit_on_error) call err_exit
    return
  endif

  ptc_fibre%dir = ele%orientation
  if (present(ref_in)) ptc_fibre%dir = ptc_fibre%dir * ref_in%direction

! Non-layout case

else
  call set_madx (energy = ele%value(e_tot$), method = ptc_key%method , step = ptc_key%nstep)

  if (associated(bmadl%start)) then
    call kill (bmadl)
    call set_up(bmadl)
  endif
  call create_fibre_append (.false., bmadl, ptc_key, EXCEPTION, bri = pancake_field)   ! ptc routine

  ptc_fibre => bmadl%end

  ! NB: Set of pointers only needed if doing stuff other than tracking (like calculating misalignments).

  ptc_fibre%mag%p%dir     => ptc_fibre%dir
  ptc_fibre%mag%p%beta0   => ptc_fibre%beta0
  ptc_fibre%mag%p%gamma0i => ptc_fibre%gamma0i
  ptc_fibre%mag%p%gambet  => ptc_fibre%gambet
  ptc_fibre%mag%p%mass    => ptc_fibre%mass
  ptc_fibre%mag%p%charge  => ptc_fibre%charge

  ptc_fibre%magp%p%dir     => ptc_fibre%dir
  ptc_fibre%magp%p%beta0   => ptc_fibre%beta0
  ptc_fibre%magp%p%gamma0i => ptc_fibre%gamma0i
  ptc_fibre%magp%p%gambet  => ptc_fibre%gambet
  ptc_fibre%magp%p%mass    => ptc_fibre%mass
  ptc_fibre%magp%p%charge  => ptc_fibre%charge

  if (key == beambeam$) call beambeam_fibre_setup(ele, ptc_fibre)

  ptc_fibre%dir = ele%orientation
  if (present(ref_in)) ptc_fibre%dir = ptc_fibre%dir * ref_in%direction

  call survey (ptc_fibre, ent = global_frame, a = global_origin)
  ! Note: Misalignments/patch setup for the layout is handled after the layout is instantiated.
  if (ele%key == gkicker$) then
    vec = ele%value(x_kick$:pz_kick$)
    beta0 = ele%value(p0c$) / ele%value(E_tot$)
    call convert_bmad_to_ptc(vec, beta0, .true.)
    ptc_fibre%patch%b_ang = vec(2:6:2)
    ptc_fibre%patch%b_d   = vec(1:5:2)
    ptc_fibre%patch%patch = 20

  else
    call misalign_ptc_fibre (ele, use_offsets, ptc_fibre, .false.)
  endif

endif

!----------------------------------------------

if (ptc_key%magnet == 'PANCAKEBMADZERO') then
  call init_pancake (max_order+1, 2)

  field = 0
  call alloc_pancake(field)
  call daall0_pancake(icoef)

  do j = 1, n_pan
    call kill(ptc_fibre%mag%pa%b(j))
    call kill(ptc_fibre%magp%pa%b(j))
  enddo

  z0 = ele%s_start - ele2%s_start

  do i = 0, n_pan-1
    if (nint(ele2%value(integrator_order$)) == 4) then
      z = z0 + i * ele%value(l$) / (n_pan - 1)
    else
      k = i / 7
      z = z0 + (k + dz_pan7(i-7*k)) * dz_step
    endif

    call gen_grad_at_s_to_em_taylor(ele2, gg_map, z, em_taylor)

    do j = 1, 3
      call dacon_pancake(field(j), 0.0_rp)
      do k = 1, size(em_taylor(j)%term)
        tm => em_taylor(j)%term(k)
        coef = tm%coef * c_light / ele2%value(p0c$)
        call dacon_pancake(icoef, 0.0_rp)
        call dapok_pancake(icoef, tm%expn, coef)
        call daadd_pancake(field(j), icoef, field(j))
      enddo
      !! call dapri_pancake(field(j), 6)  ! Taylor print
    enddo

    call set_tree_g_pancake(ptc_fibre%mag%pa%b(i+1), field)
    call set_tree_g_pancake(ptc_fibre%magp%pa%b(i+1), field)
  enddo

  call kill_pancake(field)
  call dadal1_pancake(icoef)
endif


if (ele%key == floor_shift$) Then
  ptc_fibre%patch%track = .false.
endif

!

call set_ptc_quiet(12, unset$, n)

! Cylindrical field

if (associated(ele2%cylindrical_map) .and. ele2%field_calc == fieldmap$) then
  do m = 0, m_abell
    found = .false.
    do j = 1, size(ele%cylindrical_map)
      cy => ele%cylindrical_map(j)
      if (m /= cy%m) cycle
      found = .true.
      exit
    enddo

    if (found) then
      if (cy%r0(1) /= 0 .or. cy%r0(2) /= 0) then
        call out_io (s_error$, r_name, 'TRANSVERSE R0 OFFSETS NOT YET IMPLEMENTED FOR PTC TRACKING FOR: ' // ele%name)
        if (global_com%exit_on_error) call err_exit
        return
      endif
      coef_e = cy%field_scale * master_parameter_value(cy%master_parameter, ele)
      if (key == lcavity$ .or. key == rfcavity$) coef_e = coef_e * ele%value(field_autoscale$)
      coef_b = coef_e * c_light / ele%value(p0c$)
      coef_e = -coef_e * 1d-6  ! Notice negative sign.

      ptc_fibre%mag%ab%t(m)  = cy%theta0_azimuth  ! Magnetic theta0
      ptc_fibre%mag%ab%te(m) = cy%theta0_azimuth  ! Electric theta0
      ptc_fibre%mag%ab%dz(m) = cy%dz

      n = size(cy%ptr%term%b_coef)
      select case (cy%ele_anchor_pt)
      case (anchor_beginning$)
        k_0 = -I_imaginary * twopi * cy%r0(3) / (n * cy%dz)
      case (anchor_center$)
        k_0 = -I_imaginary * twopi * (cy%r0(3) + ele%value(l$)/2) / (n * cy%dz)
      case (anchor_end$)
        k_0 = -I_imaginary * twopi * (cy%r0(3) + ele%value(l$)) / (n * cy%dz)
      end select

      nn = max(n, 2)  ! n = 1 is a singular case so treat it as n = 2
      do ii = 1, nn
        if (ii == 2 .and. n == 1) then
          k = ii - 1 - nn  ! = -1
          ptc_fibre%mag%ab%b(m,k) = 0
          ptc_fibre%mag%ab%e(m,k) = 0
        elseif (ii <= nn/2) then
          k = ii - 1 
          ptc_fibre%mag%ab%b(m,k) = coef_b * exp(k_0*k) * cy%ptr%term(ii)%b_coef
          ptc_fibre%mag%ab%e(m,k) = coef_e * exp(k_0*k) * cy%ptr%term(ii)%e_coef
        else
          k = ii - 1 - nn
          ptc_fibre%mag%ab%b(m,k) = coef_b * exp(k_0*k) * cy%ptr%term(ii)%b_coef
          ptc_fibre%mag%ab%e(m,k) = coef_e * exp(k_0*k) * cy%ptr%term(ii)%e_coef
        endif
      enddo
      !! ptc_fibre%mag%ab%b(m,:) = [cy%ptr%term(n/2+1:n)%b_coef, cy%ptr%term(1:n/2)%b_coef] * coef
    else
      ptc_fibre%mag%ab%t(m)   = 0
      ptc_fibre%mag%ab%dz(m)  = 1
      ptc_fibre%mag%ab%b(m,:) = 0
      ptc_fibre%mag%ab%e(m,:) = 0
    endif

    ptc_fibre%magp%ab%t(m)   = ptc_fibre%mag%ab%t(m)
    ptc_fibre%magp%ab%te(m)  = ptc_fibre%mag%ab%te(m)
    ptc_fibre%magp%ab%dz(m)  = ptc_fibre%mag%ab%dz(m)
    do ii = lbound(ptc_fibre%magp%ab%b, 2), ubound(ptc_fibre%magp%ab%b, 2)
      ptc_fibre%magp%ab%b(m,ii) = ptc_fibre%mag%ab%b(m,ii)
      ptc_fibre%magp%ab%e(m,ii) = ptc_fibre%mag%ab%e(m,ii)
    enddo

    ptc_fibre%mag%ab%xprime = xprime_abell   ! is_false(ele%value(ptc_canonical_coords$))
    ptc_fibre%magp%ab%xprime = ptc_fibre%mag%ab%xprime
  enddo

endif

!-----------------------------
! The E-field units that PTC wants on input are MV/m (MAD convention). 
! Note: PTC convert MV/m to GV/m internally.

if (key /= multipole$ .and. (associated(ele%a_pole_elec) .or. key == elseparator$)) then
  if (key /= sbend$) then
    ptc_fibre%mag%p%bend_fringe = .false.
    ptc_fibre%magp%p%bend_fringe = .false.
  endif
  fh = 1d-9 * sign_of(charge_of(param%particle)) / VOLT_C
  if (key == sbend$ .and. nint(ele%value(exact_multipoles$)) == vertically_pure$ .and. ele%value(g$) /= 0) then
    call multipole_ele_to_ab(ele, .false., ix_pole_max, a_pole, b_pole, electric$, include_kicks$)
    ! Notice that a_pole and b_pole are reversed for electric fields.
    call convert_bend_exact_multipole(ele%value(g$), horizontally_pure$, b_pole, a_pole)
  else
    call multipole_ele_to_ab(ele, .false., ix_pole_max, a_pole, b_pole, electric$)
  endif

  if (key == elseparator$) then
    a_pole(0) = a_pole(0) + val(vkick$) * ele%value(p0c$) / leng
    b_pole(0) = b_pole(0) + val(hkick$) * ele%value(p0c$) / leng
  endif

  do i = 0, n_pole_maxx
    if (a_pole(i) /= 0) call add (ptc_fibre, -(i+1), 0, fh*a_pole(i), electric = .true.)
    if (b_pole(i) /= 0) call add (ptc_fibre,  (i+1), 0, fh*b_pole(i), electric = .true.)
  enddo
  SOLVE_ELECTRIC = .false.
endif

! Aperture

if (ele%aperture_at /= no_aperture$ .and. (ele%value(x1_limit$) /= 0 .or. ele%value(x2_limit$) /= 0) .and. &
            (ele%value(y1_limit$) /= 0 .or. ele%value(y2_limit$) /= 0) .and. &
            (ele%aperture_type == rectangular$ .or. ele%aperture_type == elliptical$)) then

  select case (ele%aperture_type)
  case (rectangular$);    ap_type = 2
  case (elliptical$);     ap_type = 1
  end select

  select case (ele%aperture_at)
  case (both_ends$);      ap_pos =  0
  case (entrance_end$);   ap_pos = -1
  case (exit_end$);       ap_pos =  1
  end select

  ap_lim = 0.5_rp * [ele%value(x1_limit$) + ele%value(x2_limit$), ele%value(y1_limit$) + ele%value(y2_limit$)]
  ap_dxy = 0.5_rp * [ele%value(x1_limit$) - ele%value(x2_limit$), ele%value(y1_limit$) - ele%value(y2_limit$)]
  call assign_aperture (ptc_fibre, ap_type, ap_lim, ap_lim(1), ap_lim(2), ap_dxy(1), ap_dxy(2), pos = ap_pos)
endif

! sad_mult & quadrupole
! Following the SAD convention: A zero length sad_mult has no fringe.

if ((key == sad_mult$ .and. ele%value(l$) /= 0) .or. key == quadrupole$) then
  ptc_fibre%mag%va  = -sign(sqrt(24.0_rp * abs(ele%value(fq1$))), ele%value(fq1$))
  ptc_fibre%magp%va = -sign(sqrt(24.0_rp * abs(ele%value(fq1$))), ele%value(fq1$))
  ptc_fibre%mag%vs  = ele%value(fq2$)
  ptc_fibre%magp%vs = ele%value(fq2$)
endif

if (key == sad_mult$) then
  cos_t = cos(ele%value(tilt_tot$))
  sin_t = sin(ele%value(tilt_tot$))
  dx =  ele%value(x_offset_mult$) * cos_t + ele%value(y_offset_mult$) * sin_t 
  dy = -ele%value(x_offset_mult$) * sin_t + ele%value(y_offset_mult$) * cos_t

  if (ptc_fibre%mag%kind == kind5) then  ! Non-zero length.
    ptc_fibre%mag%s5%dx = dx
    ptc_fibre%mag%s5%dy = dy

    ptc_fibre%magp%s5%dx = dx
    ptc_fibre%magp%s5%dy = dy
 elseif (ptc_fibre%mag%kind /= kind3) then
    call out_io (s_fatal$, r_name, 'INTERNAL ERROR SETTING MULT OFFSET. PLEASE CONTACT DAVID SAGAN.')
  endif
endif

! Set reference energy to the exit reference energy.

energy_work = 0
call find_energy (energy_work, p0c = 1d-9 * ele%value(p0c$))
ptc_fibre = energy_work

! FieldMap cartesian_map element. 
! Include all wiggler elements even planar_model with field_calc = bmad_standard$

if (.not. associated(ele2%gen_grad_map) .and. .not. associated(ele2%cylindrical_map) .and. (key == wiggler$ .or. &
                          key == undulator$ .or. (associated(ele2%cartesian_map) .and. ele2%field_calc == fieldmap$))) then

  
  is_planar_wiggler = ((key == wiggler$ .or. key == undulator$) .and. ele2%field_calc == planar_model$) 

  if (associated(ele2%grid_field)) then
    call out_io (s_error$, r_name, 'PTC TRACKING IS NOT ABLE TO USE GRID_FIELDS. FOR ELEMENT: ' // ele%name)
    return
  endif

  if (n_map == 0) then
    if (is_planar_wiggler) then
      allocate(cm)
      call create_wiggler_cartesian_map(ele2, cm)
    else
      call out_io (s_error$, r_name, 'NOT ABLE TO DO PTC TRACKING FOR NON-PLANAR WIGGLER WITHOUT A CARTESIAN (OR OTHER TYPE OF) MAP.', &
                                     'FOR ELEMENT: ' // ele%name)
      return
    endif
  else
    cm => ele2%cartesian_map(1)
  endif

  if (any(abs(cm%ptr%term%kz) < 1e-100_rp)) then
    call out_io (s_error$, r_name, 'CARTESIAN MAP KZ VALUES MUST BE NON-ZERO (MAKE > 1E-100).', &
                                   'FOR ELEMENT: ' // ele%name)
    return
  endif

  if (cm%field_type == electric$) then
    call out_io (s_fatal$, r_name, 'PTC IS NOT ABLE TO HANDLE CARTESIAN_MAP WITH ELECTRIC FIELDS. FOR ELEMENT: ' // ele%name)
    return
  endif  

  select case (cm%ele_anchor_pt)
  case (anchor_beginning$); s_rel = s_rel - cm%r0(3)
  case (anchor_center$);    s_rel = s_rel - cm%r0(3) + ele%value(l$) / 2
  case (anchor_end$);       s_rel = s_rel - cm%r0(3) + ele%value(l$)
  end select

  n_term = size(cm%ptr%term)
  call POINTERS_W (ptc_fibre%mag%wi%w, n_term, 0)  ! n_term_electric needed   

  ptc_fibre%mag%wi%w%k(1,1:n_term)   = cm%ptr%term%kx
  ptc_fibre%mag%wi%w%k(2,1:n_term)   = cm%ptr%term%ky
  ptc_fibre%mag%wi%w%k(3,1:n_term)   = cm%ptr%term%kz
  ptc_fibre%mag%wi%w%f(1:n_term)     = cm%ptr%term%phi_z + s_rel * cm%ptr%term%kz
  ptc_fibre%mag%wi%w%x0(1:n_term)    = cm%ptr%term%x0
  ptc_fibre%mag%wi%w%y0(1:n_term)    = cm%ptr%term%y0
  ptc_fibre%mag%wi%w%form(1:n_term)  = 3*(cm%ptr%term%family - 1) + cm%ptr%term%form

  if (ele%is_on) then
    ff = c_light * cm%field_scale * charge_of(ele2%ref_species) / ele%value(p0c$)
    if (has_attribute(ele2, 'POLARITY')) ff = ff * ele2%value(polarity$)
    do i = 1, size(ptc_fibre%mag%wi%w%a(1:n_term))
      wt => cm%ptr%term(i)
      ptc_fibre%mag%wi%w%a(i) = ff * wt%coef
    enddo
  else
    ptc_fibre%mag%wi%w%a(1:n_term) = 0
  endif

  ! Correct z-position. This is only used when totalpath is not used.
  ! This is important with the ele_compute_ref_energy_and_time routine.

  if (ptc_private%base_state%totalpath == 0) then
    z_patch = ele%value(delta_ref_time$) * c_light * ele%value(p0c$) / ele%value(e_tot$) - ele%value(l$)
    ptc_fibre%mag%wi%internal(6) = z_patch
  endif

  call copy (ptc_fibre%mag, ptc_fibre%magp)

  if (n_map == 0 .and. is_planar_wiggler) then
    deallocate(cm%ptr)
    deallocate(cm)
  endif
endif

! Set charge

ptc_fibre%charge = rel_charge

! Taylor maps
! In theory can put in a taylor map for any element but for now only setup Bmad taylor and match elements.

if (key == taylor$ .or. key == match$) then
  ! The map can be split into pieces by taking the log of the map.
  ! onemap = T means do not split.
  ! At this point in time there is no splitting allowed.

  if (.not. associated(ele%taylor(1)%term)) call ele_to_taylor (ele, param)

  onemap = .true.

  call alloc(ptc_c_damap)

  ! 

  beta0 = ele%value(p0c_start$)/ele%value(e_tot_start$)
  beta1 = ele%value(p0c$)/ele%value(e_tot$)  

  call alloc (ptc_re8)
  call alloc (ptc_taylor)

  call taylor_to_real_8 (ele%taylor, beta0, beta1, ptc_re8, ref0, ref1)
  ptc_c_damap = ptc_re8

  if (associated(ele%spin_taylor(1)%term)) then
    do j = 0, 3
      ptc_taylor = ele%spin_taylor(j)
      ptc_c_damap%q%x(j) = ptc_taylor
    enddo
  endif

  call kill (ptc_taylor)
  call kill (ptc_re8);

  ! 

  if (.not. onemap) then
    call nth_root(ptc_c_damap, ptc_c_damap, ptc_fibre%mag%p%nst)
  endif

  if (ptc_fibre%dir == 1) then
    if (.not.associated(ptc_fibre%mag%forward)) then 
      allocate(ptc_fibre%mag%forward(3))
    else
      call KILL(ptc_fibre%mag%forward)  ! The kill only zeros the Berz-part, it stays associated.
    endif

    call SET_TREE_G_complex(ptc_fibre%mag%forward, ptc_c_damap)
    ptc_fibre%mag%do1mapf = onemap
    ptc_fibre%mag%usef = .true.
    arbre => ptc_fibre%mag%forward
  else
    if (.not. associated(ptc_fibre%mag%backward)) then 
      allocate(ptc_fibre%mag%backward(3))
    else
      call KILL(ptc_fibre%mag%backward)  ! The kill only zeros the Berz-part, it stays associated.
    endif
    call SET_TREE_G_complex(ptc_fibre%mag%backward, ptc_c_damap)
    ptc_fibre%mag%do1mapb = onemap
    ptc_fibre%mag%useb = .true.
    arbre => ptc_fibre%mag%backward
  endif

  !! arbre => ptc_fibre%mag%forward
  call mat_make_unit(arbre(1)%rad)    ! Radiation damping matrix. Unit matrix  => radiation off
  arbre(1)%fix0(1:6) = ref0
  arbre(1)%fixr(1:6) = ref1
  arbre(1)%fix(1:6) = ref1

  if (onemap) then
    arbre(1)%ds = ptc_fibre%mag%p%ld
  else
    arbre(1)%ds = ptc_fibre%mag%p%ld/ptc_fibre%mag%p%nst 
  endif

  arbre(1)%beta0 = ptc_fibre%beta0

  if (ptc_fibre%dir == 1) then
    if (.not.associated(ptc_fibre%magp%forward)) then 
      allocate(ptc_fibre%magp%forward(3))
      allocate(ptc_fibre%magp%usef)
    else
      call KILL(ptc_fibre%magp%forward)
    endif
    !call SET_TREE_G_complex(ptc_fibre%magp%forward,m)
    do i = 1, 3
      call alloc_tree(ptc_fibre%magp%forward(i), ptc_fibre%mag%forward(i)%n, ptc_fibre%mag%forward(i)%np)
      call copy_tree(ptc_fibre%mag%forward(i), ptc_fibre%magp%forward(i))
    enddo
    ptc_fibre%magp%do1mapf = onemap
    ptc_fibre%magp%usef = .true.             ! Use the Taylor map by default 
    arbre => ptc_fibre%magp%forward

  else
    if (.not.associated(ptc_fibre%magp%backward)) then 
      allocate(ptc_fibre%magp%backward(3))
      allocate(ptc_fibre%magp%useb)
    else
      call KILL(ptc_fibre%magp%backward)
    endif
    ! call SET_TREE_G_complex(ptc_fibre%magp%backward,m)
    do i = 1, 3
      call alloc_tree(ptc_fibre%magp%backward(i), ptc_fibre%mag%backward(i)%n, ptc_fibre%mag%backward(i)%np)
      call copy_tree(ptc_fibre%mag%backward(i), ptc_fibre%magp%backward(i))
    enddo
    ptc_fibre%magp%do1mapb = onemap
    ptc_fibre%magp%useb = .true.
    arbre => ptc_fibre%magp%backward
  endif

  ! File names in case a flat file is to be created.
  write (ptc_fibre%mag%filef, '(2a, i0, a)') trim(downcase(ele%name)), '_', ele%ix_ele, '.txf'
  !!write (ptc_fibre%mag%fileb, '(2a, i0, a)') trim(downcase(ele%name)), '_', ele%ix_ele, '.txb'

  !

  call mat_make_unit(arbre(1)%rad)    ! Radiation damping matrix. Unit matrix  => radiation off
  arbre(1)%fix0(1:6) = ref0
  arbre(1)%fixr(1:6) = ref1
  arbre(1)%fix(1:6) = ref1

  if (onemap) then
    arbre(1)%ds = ptc_fibre%mag%p%ld    ! Length of magnet (arc of circle if a bend)
  else
    arbre(1)%ds = ptc_fibre%mag%p%ld/ptc_fibre%mag%p%nst
  endif

  arbre(1)%beta0 = ptc_fibre%beta0
   
  call kill(ptc_c_damap)
endif

! Customization if wanted

err_flag = .false.

call ele_to_fibre_hook (ele, ptc_fibre, param, use_offsets, err_flag)

end subroutine ele_to_fibre
