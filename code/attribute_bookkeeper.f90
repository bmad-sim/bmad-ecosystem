
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine attribute_bookkeeper (ele, param, force_bookkeeping)
!
! Routine to recalculate the dependent attributes of an element.
! If the attributes have changed then any Taylor Maps will be killed.
!
! Note: This routine does not do any other bookkeeping. Consider using
! control_bookkeeper or lattice_bookkeeper instead.
! 
! Note: The following is an old, incomplete list.
!
! BEAMBEAM:   
!     bbi_const$ = param%n_part * charge$ * classical_radius_factor /
!                           (2 * pi * p0c$ * (sig_x$ + sig_y$)
! CRYSTAL:
!     bragg_angle_in$
!     bragg_angle_out$
!     tilt_corr$
!
! ELSEPARATOR:
!     e_field$ = sqrt(hkick$**2 + vkick$**2) * p0c$ / l$
!     voltage$ = e_field$ * gap$ 
!
! LCAVITY:    
!     voltage$ = gradient$ * L$ 
! 
! RFCAVITY:   
!     harmon$  = rf_frequency$ / T0
!
! SBEND:      
!     angle$   = L$ * G$
!     l_chord$ = 2 * sin(Angle$/2) / G$
!     rho$     = 1 / G$
!
! WIGGLER (map_type):
!     B_MAX$    
!     k1$  = -0.5 * (c_light * b_max$ / p0c$)**2
!     rho$ = p0c$ / (c_light * b_max$)
!
! WIGGLER (periodic_type):
!     k1$  = -0.5 * (c_light * b_max$ / p0c$)**2
!     rho$ = p0c$ / (c_light * b_max$)
!     n_pole$ = L$ / l_pole$
!
! Modules needed:
!   use bmad
!
! Input:
!   ele            -- Ele_struct: Element with attributes 
!   param          -- lat_param_struct: 
!   force_bookkeeping 
!                  -- Logical, optional: If present and True then force
!                       attribute bookkeeping to be done independent of
!                       the state of ele%bookkeeping_stat%attributes.
! Output:
!   ele            -- Ele_struct: Element with self-consistant attributes.
!
! Programming Note: If the dependent attributes are changed then 
!       the attribute_free routine must be modified.
!-

subroutine attribute_bookkeeper (ele, param, force_bookkeeping)

use s_fitting, only: check_bend
use bookkeeper_mod, except_dummy => attribute_bookkeeper

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord, slave, slave2
type (lat_param_struct) param
type (coord_struct) start, end
type (em_field_struct) field
type (branch_struct), pointer :: branch
type (cartesian_map_term1_struct), pointer :: term

real(rp) factor, gc, f2, phase, E_tot, polarity, dval(num_ele_attrib$), time
real(rp) w_inv(3,3), len_old, f
real(rp), pointer :: val(:), tt
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx), eps6
real(rp) kick_magnitude, bend_factor, quad_factor, radius0, step_info(7), dz_dl_max_err
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer i, n, n_div, ixm, ix_pole_max

character(20) ::  r_name = 'attribute_bookkeeper'

logical, optional :: force_bookkeeping
logical err_flag, set_l
logical non_offset_changed, offset_changed, offset_nonzero, is_on
logical :: v_mask(num_ele_attrib$), vv_mask(num_ele_attrib$), offset_mask(num_ele_attrib$)
logical :: dval_change(num_ele_attrib$)

! Note: If the dependent attributes are changed then attribute_free must be modified.

! Some init

val => ele%value

! Overlay and group and hybrid elements do not have any dependent attributes

select case (ele%key)
case (overlay$, group$, hybrid$)
  ele%bookkeeping_state%attributes = ok$
  return
end select

! 

if (bmad_com%auto_bookkeeper) then
  call attributes_need_bookkeeping(ele, dval)
  if (ele%bookkeeping_state%attributes /= stale$ .and. .not. logic_option(.false., force_bookkeeping)) return

else
  if (ele%bookkeeping_state%attributes /= stale$ .and. .not. logic_option(.false., force_bookkeeping)) return

  if (ele%lord_status /= not_a_lord$) then
    call set_ele_status_stale (ele, control_group$)
  endif

  if (ele%old_value(l$) /= val(l$)) then
    call set_ele_status_stale (ele, s_and_floor_position_group$)
  endif

  if (ele%lord_status /= multipass_lord$) then
    call set_ele_status_stale (ele, mat6_group$)
  endif

  dval = abs(val - ele%old_value)
endif

ele%bookkeeping_state%attributes = ok$
ele%bookkeeping_state%rad_int = stale$
ele%bookkeeping_state%ptc     = stale$

! Transfer tilt to tilt_tot, etc.

if (.not. on_a_girder(ele) .and. has_orientation_attributes(ele) .and. ele%slave_status /= multipass_slave$) then
  select case (ele%key)
  case (sbend$)
    val(roll_tot$)     = val(roll$)
    val(ref_tilt_tot$) = val(ref_tilt$)
  case (crystal$, mirror$, multilayer_mirror$)
    val(tilt_tot$)     = val(tilt$)
    val(ref_tilt_tot$) = val(ref_tilt$)
  case default
    val(tilt_tot$)     = val(tilt$)
  end select

  val(x_offset_tot$) = val(x_offset$)
  val(y_offset_tot$) = val(y_offset$)
  val(z_offset_tot$) = val(z_offset$)
  val(x_pitch_tot$)  = val(x_pitch$)
  val(y_pitch_tot$)  = val(y_pitch$)
endif

! Super_lord length change is put in last slave

if (ele%lord_status == super_lord$ .and. dval(l$) /= 0) then
  slave => pointer_to_slave(ele, ele%n_slave)
  len_old = slave%value(l$)
  slave%value(l$) = ele%value(l$) + ele%value(lord_pad1$) + ele%value(lord_pad2$)
  do i = 1, ele%n_slave - 1
    slave2 => pointer_to_slave(ele, i)
    slave%value(l$) = slave%value(l$) - slave2%value(l$)
  enddo
  ! Now adjust the lengths of all of the lords of the slave whose length was changed.
  do i = 1, slave%n_lord
    lord => pointer_to_lord(slave, i)
    if (lord%ix_ele == ele%ix_ele) cycle   ! Do not adjust ele length
    lord%value(l$) = lord%value(l$) + (slave%value(l$) - len_old)
  enddo
endif

! Taylor elements need no more bookkeeping

if (ele%key == taylor$) then
  ele%old_value = val
  return
endif

! Field_master...

if (ele%field_master) then

  if (val(p0c$) == 0) then
    factor = 0
  else
    factor = charge_of(param%particle) * c_light / val(p0c$)
  endif

  select case (ele%key)
  case (quadrupole$)
    val(k1$) = factor * val(B1_gradient$)
  case (sextupole$)
    val(k2$) = factor * val(B2_gradient$)
  case (octupole$)
    val(k3$) = factor * val(B3_gradient$)
  case (solenoid$)
    val(ks$) = factor * val(Bs_field$)
  case (sad_mult$)
    val(ks$) = factor * val(Bs_field$)
  case (sol_quad$)
    val(ks$) = factor * val(Bs_field$)
    val(k1$) = factor * val(B1_gradient$)
  case (sbend$)
    val(g$)     = factor * val(B_field$)
    val(g_err$) = factor * val(B_field_err$)
    val(k1$)    = factor * val(B1_gradient$)
    val(k2$)    = factor * val(B2_gradient$)
  case (hkicker$)
    val(kick$) = factor * val(BL_kick$)
  case (vkicker$)
    val(kick$) = factor * val(BL_kick$)
  case (bend_sol_quad$)
    val(g$)  = factor * val(B_field$)
    val(k1$) = factor * val(B1_gradient$)
    val(ks$) = factor * val(Bs_field$)
  end select

  if (has_hkick_attributes(ele%key) .and. ele%key /= elseparator$) then
    val(hkick$) = factor * val(BL_hkick$)
    val(vkick$) = factor * val(BL_vkick$)
  endif

else

  if (charge_of(param%particle) == 0) then
    factor = 0
  else
    factor = val(p0c$) / (charge_of(param%particle) * c_light)
  endif

  select case (ele%key)
  case (quadrupole$)
    val(B1_gradient$) = factor * val(k1$)
  case (sextupole$)
    val(B2_gradient$) = factor * val(k2$)
  case (octupole$)
    val(B3_gradient$) = factor * val(k3$)
  case (solenoid$)
    val(Bs_field$)    = factor * val(ks$)
  case (sad_mult$)
    val(Bs_field$)    = factor * val(ks$)
  case (sol_quad$)
    val(Bs_field$)    = factor * val(ks$)
    val(B1_gradient$) = factor * val(k1$)
  case (sbend$)
    val(B_field$)     = factor * val(g$)
    val(B_field_err$) = factor * val(g_err$)
    val(B1_gradient$) = factor * val(k1$)
    val(B2_gradient$) = factor * val(k2$)
  case (hkicker$)
    val(BL_kick$) = factor * val(kick$)
  case (vkicker$) 
    val(BL_kick$) = factor * val(kick$)
  case (bend_sol_quad$)
    val(B_field$)     = factor * val(g$)
    val(B1_gradient$) = factor * val(k1$)
    val(Bs_field$)    = factor * val(ks$)
  end select

  if (has_hkick_attributes(ele%key) .and. ele%key /= elseparator$) then
    val(BL_hkick$) = factor * val(hkick$)
    val(BL_vkick$) = factor * val(vkick$)
  endif

endif

! If the ref energy has not been set then, for example, k1 for a bend with field_master = T has not been set.
! So have to wait until the energy is set to figure out ds_step.

if (attribute_index(ele, 'DS_STEP') > 0 .and. val(p0c$) > 0) then  ! If this is an attribute for this element...

  dz_dl_max_err = 1d-10

  ! Set ds_step and/or num_steps if not already set.

  if (val(ds_step$) == 0 .and. val(num_steps$) == 0) then
    select case (ele%key)
    case (wiggler$, undulator$) 
      if (val(l_pole$) /= 0) val(ds_step$) = val(l_pole$) / 10

    case (sbend$, quadrupole$, sextupole$)
      if (val(integrator_order$) == 0 .and. val(l$) /= 0) then
        bend_factor = sqrt(val(hkick$)**2 + val(vkick$)**2) / val(l$)
        radius0 = ele%value(r0_mag$)
        if (radius0 == 0) radius0 = 0.01   ! Use a 1 cm scale default

        select case (ele%key)
        case (sbend$)
          bend_factor = bend_factor + abs(val(g$)) + abs(val(g_err$)) 
          quad_factor = val(l$) * (abs(val(k1$)) + abs(val(k2$)) * radius0 + bend_factor**2)
        case (quadrupole$)
          quad_factor = val(l$) * (abs(val(k1$)) + bend_factor**2)
        case (sextupole$)
          quad_factor = val(l$) * (abs(val(k2$)) * radius0 + bend_factor**2)
        end select

        if (associated(ele%a_pole)) then
          call multipole_ele_to_ab (ele, .false., ix_pole_max, a_pole, b_pole)
          quad_factor = quad_factor + abs(a_pole(1)) + abs(b_pole(1)) + radius0 * (abs(a_pole(2)) + abs(b_pole(2)))
        endif

        if (associated(ele%a_pole_elec)) then
          radius0 = ele%value(r0_elec$)
          if (radius0 == 0) radius0 = 0.01   ! Use a 1 cm scale default
          call multipole_ele_to_ab (ele, .false., ix_pole_max, a_pole, b_pole, electric$)
          bend_factor = bend_factor + (abs(a_pole(0)) + abs(b_pole(0))) / ele%value(p0c$)
          quad_factor = quad_factor + (abs(a_pole(1)) + abs(b_pole(1)) + radius0 * (abs(a_pole(2)) + abs(b_pole(2)))) / ele%value(p0c$)
        endif

        ! check_bend is a PTC routine
        call check_bend (val(l$), quad_factor, bend_factor, dz_dl_max_err, step_info, ixm)
        val(integrator_order$) = ixm
        val(num_steps$) = max(nint(step_info(ixm)), 1)
        val(ds_step$) = abs(val(l$)) / val(num_steps$)
      endif

    case (kicker$, hkicker$, vkicker$)
      if (val(l$) /= 0) then
        if (ele%key == kicker$) then
          kick_magnitude = sqrt(val(hkick$)**2 + val(vkick$)**2) / val(l$)
        else
          kick_magnitude = val(kick$) / val(l$)
        endif
        call check_bend (val(l$), 0.0_rp, kick_magnitude, dz_dl_max_err, step_info, ixm)
        val(integrator_order$) = ixm
        val(num_steps$) = max(nint(step_info(ixm+1)), 1)
        val(ds_step$) = abs(val(l$)) / val(num_steps$)
      endif

    case (lcavity$, rfcavity$)
      if (val(l$) /= 0) then
        val(num_steps$) = 10
        val(ds_step$) = abs(val(l$)) / val(num_steps$) 
      endif
    end select
  endif

  if (val(ds_step$) <= 0) then
    if (val(num_steps$) <= 0 .or. abs(val(l$)) == 0) then
      val(ds_step$) = bmad_com%default_ds_step
    else
      val(ds_step$) = abs(val(l$)) / val(num_steps$)
    endif
  endif
   
  val(num_steps$) = max(1, nint(abs(val(l$) / val(ds_step$))))

endif

!----------------------------------
! General bookkeeping...

if (attribute_index(ele, 'L_HARD_EDGE') /= 0) val(l_hard_edge$) = val(l$)

select case (ele%key)

! Bend_sol_quad

case (bend_sol_quad$)

! BeamBeam

case (beambeam$)

  if (val(n_slice$) == 0) val(n_slice$) = 1.0 ! revert to default

  if (val(charge$) == 0 .or. param%n_part == 0) then
    val(bbi_const$) = 0

  else

    if (val(sig_x$) == 0 .or. val(sig_y$) == 0) then
      call out_io(s_abort$, r_name, 'ZERO SIGMA IN BEAMBEAM ELEMENT!')
      call type_ele(ele, .true., 0, .false., 0, .false.)
      if (global_com%exit_on_error) call err_exit
    endif

    if (val(p0c$) /= 0) then  ! Can happen when parsing lattice file.
      val(bbi_const$) = -param%n_part * val(charge$) * classical_radius_factor /  &
                             (2 * pi * val(p0c$) * (val(sig_x$) + val(sig_y$)))
    endif

  endif

! Crystal

case (crystal$, multilayer_mirror$)

  if (ele%key == crystal$) then
    call crystal_type_to_crystal_params (ele, err_flag)
    call crystal_attribute_bookkeeper (ele)
  else
    call multilayer_type_to_multilayer_params (ele, err_flag)
  endif

  ele%photon%surface%has_curvature = (any(ele%photon%surface%curvature_xy /= 0))

! E_Gun

case (e_gun$)
  if (ele%lord_status /= multipass_lord$) then
    if (val(gradient$) /= ele%old_value(gradient$) .or. val(l$) /= ele%old_value(l$)) then
      call set_ele_status_stale (ele, ref_energy_group$)
      val(voltage$) = val(gradient$) * val(l$)
      val(voltage_err$) = val(gradient_err$) * val(l$)
    endif
  endif

! Elseparator

case (elseparator$)

  if (ele%field_master) then
    if (val(p0c$) /= 0) then
      val(hkick$) = 0
      val(vkick$) = val(l$) * val(e_field$) / val(p0c$)
      val(voltage$) = val(e_field$) * val(gap$)
    endif

  else
    if (val(l$) == 0) then
      val(e_field$) = 0
      val(voltage$) = 0
    else
      val(e_field$) = sqrt(val(hkick$)**2 + val(vkick$)**2) * val(p0c$) / val(l$)
      val(voltage$) = val(e_field$) * val(gap$) 
    endif
  endif

! Lcavity

case (lcavity$)
  if (ele%lord_status /= multipass_lord$) then
    if (val(phi0$) /= ele%old_value(phi0$) .or. val(phi0_multipass$) /= ele%old_value(phi0_multipass$) .or. &
        val(gradient$) /= ele%old_value(gradient$) .or. val(e_loss$) /= ele%old_value(e_loss$) .or. &
        val(l$) /= ele%old_value(l$)) then
      call set_ele_status_stale (ele, ref_energy_group$)
    endif
  endif

  if (val(l$) /= 0) then
    val(voltage$) = val(gradient$) * val(l$)
    val(voltage_err$) = val(gradient_err$) * val(l$)
  endif

  if (val(rf_frequency$) /= 0 .and. ele%field_calc == bmad_standard$) then
    val(l_hard_edge$) = c_light * nint(val(n_cell$)) / (2 * val(rf_frequency$))
  endif

! Patch

case (patch$) 

  call floor_angles_to_w_mat (val(x_pitch$), val(y_pitch$), val(tilt$), w_mat_inv = w_inv)
  val(l$) = w_inv(3,1) * val(x_offset$) + w_inv(3,2) * val(y_offset$) + w_inv(3,3) * val(z_offset$)
  val(ds_step$) = val(l$)
  val(num_steps$) = 1

! Quadrupole

case (quadrupole$)

! RFcavity

case (rfcavity$)
  if (param%geometry == closed$ .and. associated(ele%branch) .and. val(p0c$) /= 0) then
    branch => ele%branch
    time = branch%ele(branch%n_ele_track)%ref_time
    if (time /= 0) then
      if (ele%field_master) then
        val(rf_frequency$) = val(harmon$) / time
      else
        val(harmon$) = val(rf_frequency$) * time
      endif
    endif
  endif

  if (val(rf_frequency$) /= 0 .and. ele%field_calc == bmad_standard$) then
    val(l_hard_edge$) = c_light * nint(val(n_cell$)) / (2 * val(rf_frequency$))
  endif

  if (val(voltage$) == 0) then
    val(gradient$) = 0
  elseif (val(l$) == 0) then
    val(gradient$) = 1d30    ! Something large
  else
    val(gradient$) = val(voltage$) / val(l$)
  endif

! sad_mult

case (sad_mult$)

  if (ele%value(eps_step_scale$) > 0) then
    call multipole_ele_to_kt (ele, .true., ix_pole_max, knl, tilt)
    eps6 = 6 * ele%value(eps_step_scale$) * bmad_com%sad_eps_scale
    n_div = 1
    ! This is the same algorithm as in SAD to determine the step size.
    do n = 2, ix_pole_max
      if (knl(n) == 0) cycle  
      n_div = max(n_div, int(sqrt(abs(knl(n)) * ele%value(l$) * bmad_com%sad_amp_max**(n-1) / (eps6 * factorial(n-1)))) + 1)
    enddo

    ele%value(num_steps$) = min(n_div, bmad_com%sad_n_div_max)
    ele%value(ds_step$) = ele%value(l$) / ele%value(num_steps$)
  endif


! Sbend

case (sbend$)

  val(angle$) = val(l$) * val(g$)

  if (val(g$) == 0) then
    val(rho$) = 0
    val(l_chord$) = val(l$)
    val(l_sagitta$) = 0
  else
    val(rho$) = 1 / val(g$)
    val(l_chord$) = 2 * val(rho$) * sin(val(angle$)/2)
    val(l_sagitta$) = val(rho$) * (1 - cos(val(angle$)/2))
  endif

  if (ele_value_has_changed(ele, [g$], [1e-10_rp], .false.)) then
    call set_ele_status_stale (ele, floor_position_group$)
  endif

! Sol_quad

case (sol_quad$)

! Solenoid

case (solenoid$)

! Wiggler

case (wiggler$, undulator$) 

  ! Calculate b_max for map_type wigglers. 

  if (ele%sub_key == map_type$ .and. val(b_max$) == 0) then
    is_on = ele%is_on  ! Save
    polarity = val(polarity$)
    ele%is_on = .true.
    val(polarity$) = 1
    start%vec = 0
    val(b_max$) = 0
    n = nint(val(num_steps$))
    do i = 0, n
      call em_field_calc (ele, param, i * val(l$) / n, start, .true., field, rf_time = 0.0_rp)
      val(b_max$) = max(val(b_max$), sqrt(sum(field%b**2)))
    enddo
    ele%is_on = is_on
    val(polarity$) = polarity
  endif

  if (val(p0c$) == 0) then
    val(k1$) = 0
  else
    val(k1$) = -0.5 * (c_light * val(b_max$) / val(p0c$))**2
  endif

  if (val(b_max$) == 0) then
    val(rho$) = 0
  else
    val(rho$) = val(p0c$) / (c_light * val(b_max$))
  endif

  if (val(l_pole$) == 0) then
    val(n_pole$) = 0
  else
    val(n_pole$) = val(l$) / val(l_pole$)
  endif

  ! Periodic_type wigglers have a single term %cylindrical_map(1)%ptr%term(1) for use with tracking, etc.
  ! The phase of this term is set so that tracking with a particle starting
  ! on-axis ends on-axis. For this to be true, there must be an integer number
  ! of poles.

  ! For super_slave and sliced elements, the phi_z is set by the position with respect to the lord in
  ! the routine makeup_super_slave1 and so should not be touched here.

  if (ele%sub_key == periodic_type$ .and. ele%slave_status /= super_slave$ .and. &
      ele%slave_status /= multipass_slave$ .and. ele%slave_status /= slice_slave$) then
    if (.not. associated(ele%cartesian_map)) then
      allocate (ele%cartesian_map(1))
      allocate (ele%cartesian_map(1)%ptr)
      allocate (ele%cartesian_map(1)%ptr%term(1))
      ele%cartesian_map(1)%master_parameter = polarity$
      write (ele%cartesian_map(1)%ptr%file, '(a, i0, a, i0)') 'attribute_bookkeeper:', &
                                                       ele%ix_branch, '>>', ele%ix_ele ! Unique name
    endif

    term => ele%cartesian_map(1)%ptr%term(1)
    if (val(l$) == 0) then
      term%ky = 0
    else
      if (val(n_pole$) == 0) then
        call out_io (s_error$, r_name, 'NUMBER OF POLES NOT SET FOR WIGGLER: ' // trim(ele%name))
        term%ky = pi / val(l$)
      else
        term%ky = pi * val(n_pole$) / val(l$)
      endif
    endif
    term%coef   = val(b_max$)
    term%kx     = 0
    term%kz     = term%ky
    term%x0     = 0
    term%y0     = 0
    term%phi_z  = -term%kz * val(l$) / 2 
    term%type   = hyper_y_family_y$
  endif

end select

! It might be true that value == old_value here even though this was not true at the start
! of this routine. For example: A super_slave, via makeup_super_slave (which will be
! called by control_bookkeeper1 before calling this routine), will inherit its lords 
! num_steps value but then this routine will reset num_steps to be the correct value.

! So stop here if nothing has truely changed.

if (bmad_com%auto_bookkeeper) then
  if (all(val == ele%old_value)) return
endif

! Since things have changed we need to kill the Taylor Map and ptc_genfield.
! The factor of 1d-15 is to avoid negligible changes which can be caused if the digested 
! file was created on a different machine from the machine where the code is run.

v_mask = .true.
v_mask([x_offset$, y_offset$, z_offset$, &
            tilt$, x_pitch$, y_pitch$, x_offset_tot$, y_offset_tot$, z_offset_tot$, &
            tilt_tot$, x_pitch_tot$, y_pitch_tot$]) = .false.
offset_mask = .not. v_mask
v_mask( [x1_limit$, x2_limit$, y1_limit$, y2_limit$] ) = .false.

! With runge_kutta tracking, must use a less stringent tolerance for what is a significant change.

dval = abs(val - ele%old_value)
dval_change = (dval > small_rel_change$ * abs(val))
dval_change(p0c_start$)   = (dval(p0c_start$)   > (1d-3 + bmad_com%rel_tol_adaptive_tracking * val(p0c_start$)))
dval_change(E_tot_start$) = (dval(E_tot_start$) > (1d-3 + bmad_com%rel_tol_adaptive_tracking * val(E_tot_start$)))
dval_change(p0c$)         = (dval(p0c$)         > (1d-3 + bmad_com%rel_tol_adaptive_tracking * val(p0c$)))
dval_change(E_tot$)       = (dval(E_tot$)       > (1d-3 + bmad_com%rel_tol_adaptive_tracking * val(E_tot$)))

! delta_ref_time can have relatively large changes since this is computed 
! as an absolute time difference. Also it is a dependent attribute.
dval_change(delta_ref_time$) = .false.  

if (has_orientation_attributes(ele)) then
  ! non_offset_changed is used to determine if maps should be killed.
  ! Since maps do not depend upon custom attributes, veto them in the mask.
  vv_mask = v_mask
  vv_mask(custom_attribute1$:custom_attribute_max$) = .false.
  non_offset_changed = (any(dval_change .and. v_mask))
  offset_changed =  (any(dval_change .and. offset_mask))
  offset_nonzero = (any(val /= 0 .and. offset_mask))
else
  non_offset_changed = (any(dval_change))
  offset_changed = .false.
  offset_nonzero = .false.
endif

! If an element has just been offset and bmad_com%conserve_taylor_map = T then 
! conserve the taylor map.

if (associated(ele%taylor(1)%term) .and. ele%taylor_map_includes_offsets .and. &
        offset_nonzero .and. offset_changed .and. .not. non_offset_changed .and. &
        bmad_com%conserve_taylor_maps .and. ele%key /= patch$) then
  ele%taylor_map_includes_offsets = .false.
  if (associated(ele%branch) .and. ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$) then
    do i = 1, ele%n_lord
      lord => pointer_to_lord(ele, i)
      if (lord%slave_status == multipass_slave$) lord => pointer_to_lord(lord, 1)
      lord%taylor_map_includes_offsets = .false.
    enddo
  endif
  if (any(ele%old_value /= 0 .and. offset_mask)) non_offset_changed = .true.  ! To trigger kill_taylor below
  call out_io (s_info$, r_name, &
      'Note: bmad_com%conserve_taylor_maps = True (this is the default)', &
      'But: Element has just been offset: ' // ele%name, &
      "To conserve the element's Taylor map, I will set ele%taylor_map_includes_offsets = False.")
endif

! Kill the taylor map and ptc_genfield if necessary.

if (non_offset_changed .or. (offset_changed .and. ele%taylor_map_includes_offsets)) then
  if (associated(ele%taylor(1)%term)) call kill_taylor(ele%taylor, ele%spin_taylor)
  if (associated(ele%ptc_genfield%field)) call kill_ptc_genfield(ele%ptc_genfield%field)
endif

! Make stale ele%rad_int_cache if allocated

if (associated(ele%rad_int_cache)) ele%rad_int_cache%stale = .true.  ! Forces recalc

! Set old_value = value

do i = 1, num_ele_attrib$
  if (.not. dval_change(i) .or. .not. v_mask(i)) cycle
  if (i == check_sum$) then
    if (ele%key /= sad_mult$)   call set_flags_for_changed_attribute (ele, ele%tracking_method, .false.)
  else
    call set_flags_for_changed_attribute (ele, ele%value(i), .false.)
  endif
enddo
 
ele%bookkeeping_state%attributes = ok$
ele%old_value = val

end subroutine attribute_bookkeeper

