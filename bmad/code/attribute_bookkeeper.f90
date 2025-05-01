!+
! Subroutine attribute_bookkeeper (ele, force_bookkeeping)
!
! Routine to recalculate the dependent attributes of an element.
! If the attributes have changed then any Taylor Maps will be killed.
!
! Note: This routine does not do any other bookkeeping. Consider using
! control_bookkeeper or lattice_bookkeeper instead.
!
! Input:
!   ele            -- Ele_struct: Element with attributes 
!   force_bookkeeping 
!                  -- Logical, optional: If present and True then force
!                       attribute bookkeeping to be done independent of
!                       the state of ele%bookkeeping_stat%attributes.
!                       This will also cause attribute_bookkeeper to assume intelligent bookkeeping.
! Output:
!   ele            -- Ele_struct: Element with self-consistant attributes.
!
! Programming Note: If the dependent attributes are changed then 
!       the attribute_free routine must be modified.
!-

subroutine attribute_bookkeeper (ele, force_bookkeeping)

use s_fitting, only: check_bend
use bookkeeper_mod, except_dummy => attribute_bookkeeper
use xraylib_interface, except_dummy2 => attribute_bookkeeper
use xraylib, dummy => r_e
use super_recipes_mod, only: super_brent
use ptc_layout_mod, only: update_ele_from_fibre
use particle_species_mod
use bmad_parser_struct, only: bp_com
use bmad_parser_mod, only: settable_dep_var_bookkeeping

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord, slave, slave2
type (coord_struct) start, end
type (em_field_struct) field
type (branch_struct), pointer :: branch
type (photon_element_struct), pointer :: ph
type (wake_lr_mode_struct), pointer :: lr
type (converter_prob_pc_r_struct), pointer :: ppcr
type (molecular_component_struct), allocatable :: component(:)
type (material_struct), pointer :: material
type (material_struct) :: materi

real(rp) factor, e_factor, gc, f2, phase, E_tot, polarity, dval(num_ele_attrib$), time, beta
real(rp) w_inv(3,3), len_old, f, dl, b_max, zmin, ky, kz, tot_mass
real(rp), pointer :: val(:), tt
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx), eps6
real(rp) kick_magnitude, bend_factor, quad_factor, radius0, step_info(7), dz_dl_max_err
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)

integer i, j, ix, ig, n, n_div, ixm, ix_pole_max, particle, geometry, i_max, status, z_material

character(20) ::  r_name = 'attribute_bookkeeper'

logical, optional :: force_bookkeeping
logical err_flag, set_l, field_bookkeeping_doable
logical non_offset_changed, offset_changed, offset_nonzero, is_on
logical :: v_mask(num_ele_attrib$), vv_mask(num_ele_attrib$), offset_mask(num_ele_attrib$)
logical :: dval_change(num_ele_attrib$)

! Note: If the dependent attributes are changed then attribute_free must be modified.

! Some init

val => ele%value
branch => pointer_to_branch(ele)

if (associated(branch)) then
  particle = branch%param%particle
  geometry = branch%param%geometry
else
  particle = positron$  ! Just to keep the code happy
  geometry = open$
endif
  
! Overlay and group and hybrid elements do not have any dependent attributes

select case (ele%key)
case (overlay$, group$, hybrid$)
  ele%bookkeeping_state%attributes = ok$
  return
end select

! 

if (bmad_com%auto_bookkeeper .and. .not. logic_option(.false., force_bookkeeping)) then
  call attributes_need_bookkeeping(ele, dval)
  if (ele%bookkeeping_state%attributes /= stale$) return

  if (.false. .and. bp_com%parser_name == '') then   ! If not parsing should not be here
    call out_io (s_warn$, r_name, &
      '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', &
      '!!!!! Using intelligent bookkeeping is now mandated for all Bmad based programs.               !!!!!', &
      '!!!!! See the "Intelligent Bookkeeping" section in the Bmad manual.                            !!!!!', &
      '!!!!! This program will run now but if this program modifies any lattice parameters, the       !!!!!', &
      '!!!!! correctness of the results is questionable.                                              !!!!!', &
      '!!!!! Contact the maintainer of this program with this information.                            !!!!!', &
      '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
  endif
else
  call attributes_need_bookkeeping(ele)
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

! lr wake

if (associated(ele%wake)) then
  do i = 1, size(ele%wake%lr%mode)
    lr => ele%wake%lr%mode(i)
    if (lr%freq_in < 0) lr%freq = ele%value(rf_frequency$)

    ! Old style lattice files set Q and not damp.
    if (lr%q /= real_garbage$) then  
      if (lr%q == 0) then
        call out_io (s_error$, r_name, 'Q factor for LR wake mode is zero which does not make sense!', &
                                       'For element: ' // ele%name)
      else
        lr%damp = pi * lr%freq / lr%q
      endif
    endif

    if (lr%damp == 0) then
      lr%q = 1d100
    else
      lr%q = pi * lr%freq / lr%damp
    endif

  enddo
endif

! Transfer tilt to tilt_tot, etc.

if (.not. associated(pointer_to_girder(ele)) .and. has_orientation_attributes(ele) &
                                                  .and. ele%slave_status /= multipass_slave$) then
  select case (ele%key)
  case (sbend$, rf_bend$)
    val(roll_tot$)     = val(roll$)
    val(ref_tilt_tot$) = val(ref_tilt$)
    ele%bookkeeping_state%has_misalign = (val(ref_tilt_tot$) /= 0 .or. val(roll_tot$) /= 0)

  case (crystal$, mirror$, multilayer_mirror$)
    val(tilt_tot$)     = val(tilt$)
    val(ref_tilt_tot$) = val(ref_tilt$)
    ele%bookkeeping_state%has_misalign = (val(ref_tilt_tot$) /= 0 .or. val(tilt_tot$) /= 0)

  case default
    val(tilt_tot$)     = val(tilt$)
    ele%bookkeeping_state%has_misalign = (val(tilt_tot$) /= 0)
  end select

  val(x_offset_tot$) = val(x_offset$)
  val(y_offset_tot$) = val(y_offset$)
  val(z_offset_tot$) = val(z_offset$)
  val(x_pitch_tot$)  = val(x_pitch$)
  val(y_pitch_tot$)  = val(y_pitch$)

  ele%bookkeeping_state%has_misalign = (ele%bookkeeping_state%has_misalign .or. val(x_offset_tot$) /= 0 .or. &
                                        val(y_offset_tot$) /= 0 .or. val(z_offset_tot$) /= 0 .or. &
                                        val(x_pitch_tot$) /= 0 .or. val(y_pitch_tot$) /= 0)
  if (ele%key == sad_mult$) ele%bookkeeping_state%has_misalign = (ele%bookkeeping_state%has_misalign .or. &
                                                         val(x_offset_mult$) /= 0 .or. val(y_offset_mult$) /= 0)
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

field_bookkeeping_doable = (val(p0c$) /= 0 .and. particle /= not_set$)

if (field_bookkeeping_doable) then
  if (ele%field_master) then

    factor = charge_of(particle) * c_light / val(p0c$)

    select case (ele%key)
    case (beambeam$)
      val(ks$) = factor * val(Bs_field$)
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
    case (rf_bend$)
      if (is_true(ele%value(init_needed$))) call settable_dep_var_bookkeeping(ele)
      val(g$)     = factor * val(B_field$)
    case (sbend$)
      if (is_true(ele%value(init_needed$))) call settable_dep_var_bookkeeping(ele)
      val(g$)     = factor * val(B_field$)
      val(dg$)    = factor * val(dB_field$)
      val(k1$)    = factor * val(B1_gradient$)
      val(k2$)    = factor * val(B2_gradient$)
    case (hkicker$, vkicker$)
      val(kick$) = factor * val(BL_kick$)
    end select

    if (has_hkick_attributes(ele%key) .and. ele%key /= elseparator$) then
      val(hkick$) = factor * val(BL_hkick$)
      val(vkick$) = factor * val(BL_vkick$)
    endif

  else

    if (charge_of(particle, 0) == 0) then
      factor = 0
    else
      factor = val(p0c$) / (charge_of(particle) * c_light)
    endif

    select case (ele%key)
    case (beambeam$)
      val(Bs_field$)    = factor * val(ks$)
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
    case (rf_bend$)
      if (is_true(ele%value(init_needed$))) call settable_dep_var_bookkeeping(ele)
      val(B_field$)     = factor * val(g$)
    case (sbend$)
      if (is_true(ele%value(init_needed$))) call settable_dep_var_bookkeeping(ele)
      val(B_field$)     = factor * val(g$)
      val(dB_field$)    = factor * val(dg$)
      val(B1_gradient$) = factor * val(k1$)
      val(B2_gradient$) = factor * val(k2$)
    case (hkicker$)
      val(BL_kick$) = factor * val(kick$)
    case (vkicker$) 
      val(BL_kick$) = factor * val(kick$)
    end select

    if (has_hkick_attributes(ele%key) .and. ele%key /= elseparator$) then
      val(BL_hkick$) = factor * val(hkick$)
      val(BL_vkick$) = factor * val(vkick$)
    endif

  endif
endif

! If the ref energy has not been set then, for example, k1 for a bend with field_master = T has not been set.
! So have to wait until the energy is set to figure out ds_step.

if (attribute_index(ele, 'DS_STEP') > 0 .and. val(p0c$) > 0) then  ! If this is an attribute for this element...

  dz_dl_max_err = 1d-8

  ! Set ds_step and/or num_steps if not already set.

  if (val(ds_step$) == 0 .and. val(num_steps$) == 0) then
    if (ele%field_calc == fieldmap$ .and. ele%tracking_method /= bmad_standard$ .and. associated(ele%gen_grad_map)) then
      n = ele%gen_grad_map(1)%iz1 + 1 - ele%gen_grad_map(1)%iz0
      if (nint(val(integrator_order$)) /= 6) val(integrator_order$) = 4
      if (nint(val(integrator_order$)) == 4) then
        val(num_steps$) = max(1, nint((n-1)/4.0_rp))
      else 
        val(num_steps$) = max(1, nint((n-1)/7.0_rp))
      endif

      val(ds_step$) = val(l$) / val(num_steps$)

    else
      select case (ele%key)
      case (drift$, pipe$, ecollimator$, rcollimator$, instrument$, monitor$)
        val(num_steps$) = 1
        val(ds_step$) = val(l$)

      case (wiggler$, undulator$) 
        if (val(l_period$) /= 0) val(ds_step$) = val(l_period$) / 20

      case (sbend$, quadrupole$, sextupole$)
        if (val(l$) /= 0) then
          call multipole_ele_to_kt (ele, .false., ix, knl, tilt, magnetic$, include_kicks$)
          knl = abs(knl)
          bend_factor = knl(0) / val(l$)
          radius0 = ele%value(r0_mag$)
          if (radius0 == 0) radius0 = 0.01   ! Use a 1 cm scale default

         select case (ele%key)
         case (sbend$)
           bend_factor = bend_factor + max(abs(val(g$)), abs(val(g$) + val(dg$)))
           quad_factor = knl(1) + knl(2) * radius0 + ele%value(l$) * bend_factor**2
         case (quadrupole$)
           ! The factor of 50 here is empirical based upon simulations with the canonical
           ! CESR bmad_L9a18A000-_MOVEREC lattice.
           quad_factor = 50 * (knl(1) + val(l$) * bend_factor**2)
         case (sextupole$)
           quad_factor = knl(2) * radius0 + val(l$) * bend_factor**2
         end select

         if (associated(ele%a_pole_elec)) then
           radius0 = ele%value(r0_elec$)
           if (radius0 == 0) radius0 = 0.01   ! Use a 1 cm scale default
           call multipole_ele_to_ab (ele, .false., ix_pole_max, a_pole, b_pole, electric$)
           bend_factor = bend_factor + (abs(a_pole(0)) + abs(b_pole(0))) / ele%value(p0c$)
           quad_factor = quad_factor + (abs(a_pole(1)) + abs(b_pole(1)) + radius0 * (abs(a_pole(2)) + abs(b_pole(2)))) / ele%value(p0c$)
         endif

          ! check_bend is a PTC routine
          ix = nint(val(integrator_order$))
          if (ix /= 2 .and. ix /= 4 .and. ix /= 6) val(integrator_order$) = 0  ! Reset if current value is not valid
         call check_bend (val(l$), quad_factor, bend_factor, dz_dl_max_err, step_info, ixm)
          if (val(integrator_order$) == 0) then
            ! Since num_steps is used by Bmad routines, do not use order 6 which can give two few steps for Bmad.
            ixm = min(ixm, 4)
            val(integrator_order$) = ixm
          else
            ixm = val(integrator_order$)
          endif
          val(num_steps$) = max(nint(step_info(ixm+1)), 1)
          val(ds_step$) = abs(val(l$)) / val(num_steps$)
        endif

      case (ac_kicker$, kicker$, hkicker$, vkicker$)
        if (val(l$) /= 0) then
          if (ele%key == ac_kicker$ .or. ele%key == kicker$) then
            kick_magnitude = sqrt(val(hkick$)**2 + val(vkick$)**2) / val(l$)
          else
            kick_magnitude = val(kick$) / val(l$)
          endif

          ix = nint(val(integrator_order$))
          if (ix /= 2 .and. ix /= 4 .and. ix /= 6) val(integrator_order$) = 0  ! Reset if current value is not valid
          call check_bend (val(l$), 0.0_rp, kick_magnitude, dz_dl_max_err, step_info, ixm)
          if (val(integrator_order$) == 0) then
            val(integrator_order$) = ixm
          else
            ixm = val(integrator_order$)
          endif
          val(num_steps$) = max(nint(step_info(ixm+1)), 1)
          val(ds_step$) = abs(val(l$)) / val(num_steps$)
        endif

      end select
    endif
  endif

  if (val(ds_step$) <= 0) then
    if (val(num_steps$) <= 0 .or. abs(val(l$)) == 0) then
      val(ds_step$) = bmad_com%default_ds_step
      if (ele%key == wiggler$ .or. ele%key == undulator$) call out_io(s_warn$, r_name, &
                    'Warning! Element: ' // trim(ele%name) // ' which is a ' // key_name(ele%key), &
                    'is using the bmad_com%default_ds_step since ds_step is not set in the element.', &
                    'this may be very inaccurate.')
    else
      val(ds_step$) = abs(val(l$)) / val(num_steps$)
    endif
  endif

  val(num_steps$) = max(1, nint(abs(val(l$) / val(ds_step$))))
endif

!----------------------------------
! General bookkeeping...

select case (ele%key)

! BeamBeam

case (beambeam$)

  if (val(n_slice$) == 0) val(n_slice$) = 1.0 ! revert to default

  if (strong_beam_strength(ele) == 0) then
    val(bbi_const$) = 0

  else
    if (val(sig_x$) == 0 .or. val(sig_y$) == 0) then
      call out_io(s_error$, r_name, 'ZERO SIGMA IN BEAMBEAM ELEMENT!' // ele%name)
      call type_ele(ele, .true., 0, .false., 0)
      return
    endif

    if (val(p0c$) /= 0) then  ! Can happen when parsing lattice file.
      val(bbi_const$) = -strong_beam_strength(ele) * classical_radius_factor /  &
                                  (2 * pi * val(p0c$) * (val(sig_x$) + val(sig_y$)))
    endif
  endif

! Converter. Note: Reference energy bookkeeping handled in ele_compute_ref_energy_and_time.

case (converter$)
  if (allocated(ele%converter%dist) .and. (dval(angle_out_max$) /= 0 .or. &
                   dval(thickness$) /= 0 .or. dval(pc_out_min$) /= 0 .or. dval(pc_out_max$) /= 0)) then
    do i = 1, size(ele%converter%dist)
      do j = 1, size(ele%converter%dist(i)%sub_dist)
        ppcr => ele%converter%dist(i)%sub_dist(j)%prob_pc_r
        if (allocated(ppcr%p_norm)) deallocate (ppcr%p_norm, ppcr%integ_pc_out, ppcr%integ_r)
      enddo
    enddo
  endif

! Crab_Cavity

case (crab_cavity$)

  if (val(voltage$) == 0) then
    val(gradient$) = 0
  elseif (val(l$) == 0) then
    val(gradient$) = 1d30    ! Something large
  else
    val(gradient$) = val(voltage$) / val(l$)
  endif

  if (val(e_tot$) /= 0) then
    beta = ele%value(p0c$) / ele%value(e_tot$)
    time = branch%param%total_length / (c_light * beta)
    if (time /= 0) then
      if (is_true(val(harmon_master$))) then
        val(rf_frequency$) = val(harmon$) / time
      else
        val(harmon$) = val(rf_frequency$) * time
      endif
    endif
  endif

  if (val(rf_frequency$) /= 0) then
    val(rf_wavelength$) = c_light / val(rf_frequency$)
  else
    val(rf_wavelength$) = 0
  endif

! Foil

case (foil$)

  call molecular_components(ele%component_name, component)
  n = size(component)
  if (.not. allocated(ele%foil%material)) allocate(ele%foil%material(n))
  materi = ele%foil%material(1)

  if (n > 1 .and. size(ele%foil%material) == 1) then
    deallocate (ele%foil%material)
    allocate(ele%foil%material(n))
    ele%foil%material(1) = materi
  endif

  if (n /= size(ele%foil%material)) then
    call out_io(s_error$, r_name, 'NUMBER OF COMPONENTS IN: ' // quote(ele%component_name) // ' (' // int_str(n) // &
                    ') IS NOT THE SAME AS OTHER PARAMETERS IN ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif

  tot_mass = 0
  do ix = 1, n
    material => ele%foil%material(ix)
    material%species = species_id(component(ix)%atom)
    material%number = component(ix)%number
    tot_mass = tot_mass + mass_of(material%species) * material%number
  enddo

  if (n > 1) then
    if (materi%density /= real_garbage$ .and. ele%foil%material(2)%density == real_garbage$) then
      do ix = 1, n
        material => ele%foil%material(ix)
        material%density = materi%density * mass_of(material%species) * material%number / tot_mass
      enddo  
    endif

    if (materi%area_density /= real_garbage$ .and. ele%foil%material(2)%area_density == real_garbage$) then
      do ix = 1, n
        material => ele%foil%material(ix)
        material%area_density = materi%area_density * mass_of(material%species) * material%number / tot_mass
      enddo  
    endif
  endif

  do ix = 1, n
    material => ele%foil%material(ix)
    z_material = atomic_number(material%species)

    if (material%radiation_length == real_garbage$) then
      material%radiation_length_used = x0_radiation_length(material%species)
    else
      material%radiation_length_used = material%radiation_length
    endif

    if (material%density /= real_garbage$ .and. material%area_density /= real_garbage$) then
      call out_io(s_error$, r_name, 'SETTING BOTH DENSITY AND AREA_DENSITY IS NOT PERMITTED FOR: ' // ele%name)
      return
    endif

    if (material%density == real_garbage$) then
      if (material%area_density /= real_garbage$) then
        if (ele%value(thickness$) == 0) then
          material%density_used = real_garbage$
        else
          material%density_used = material%area_density / ele%value(thickness$)
        endif
      else
        material%density_used = ElementDensity(z_material) * 1e3_rp  ! From xraylib. Convert to kg/m^3
      endif
    else
      material%density_used = material%density
    endif

    if (material%area_density == real_garbage$) then
      material%area_density_used = material%density_used * ele%value(thickness$)
    else
      material%area_density_used = material%area_density
    endif

    if (material%density == real_garbage$ .and. material%area_density == real_garbage$ .and. n > 1) then
      call out_io(s_warn$, r_name, 'Neither foil DENSITY(s) nor AREA_DENSITY(s) set for compound material for: ' // ele%name, &
                                   'This will produce HIGHLY inaccurate results!')
    endif
  enddo

! Crystal

case (crystal$, multilayer_mirror$, mirror$, detector$, sample$, diffraction_plate$)

  if (ele%key == crystal$) then
    call crystal_type_to_crystal_params (ele, err_flag)
    call crystal_attribute_bookkeeper (ele)
  elseif (ele%key == multilayer_mirror$) then
    call multilayer_type_to_multilayer_params (ele, err_flag)
  endif

  ph => ele%photon
  ph%curvature%has_curvature = (any(ph%curvature%xy /= 0) .or. ph%curvature%spherical /= 0 .or. any(ph%curvature%elliptical /= 0))

! E_Gun

case (e_gun$)
  if (val(rf_frequency$) /= 0) then
    val(rf_wavelength$) = c_light / val(rf_frequency$)
  else
    val(rf_wavelength$) = 0
  endif

  if (ele%lord_status /= multipass_lord$) then
    if (val(gradient$) /= ele%old_value(gradient$) .or. val(l$) /= ele%old_value(l$)) then
      call set_ele_status_stale (ele, ref_energy_group$)
      val(voltage$) = val(gradient$) * val(l$)
      val(voltage_err$) = val(gradient_err$) * val(l$)
    endif
  endif

  val(voltage_tot$)  = val(voltage$)  + val(voltage_err$)
  val(gradient_tot$) = val(gradient$) + val(gradient_err$)

case (multipole$)
  if (associated(ele%a_pole)) then
    if (ele%a_pole(0) /= 0) then
      call out_io(s_error$, r_name, 'MULTIPOLE: ' // ele_full_name(ele, '@N (&#)'), &
                                    'CANNOT HAVE A FINITE K0L VALUE. WILL SET TO ZERO. SEE THE BMAD MANUAL FOR DETAILS.')
      ele%a_pole(0) = 0
    endif
  endif

! Elseparator

case (elseparator$)
  if (val(E_tot$) == 0) then ! Can happen during lattice parsing.
    e_factor = 1
  else
    e_factor = val(p0c$)**2 / val(E_tot$) ! = p0c * (v/c)
  endif

  if (ele%field_master) then
    if (val(p0c$) /= 0) then
      val(hkick$) = 0
      val(vkick$) = val(l$) * val(e_field$) / e_factor
      val(voltage$) = val(e_field$) * val(gap$)
    endif

  else
    if (val(l$) == 0) then
      val(e_field$) = 0
      val(voltage$) = 0
    else
      val(e_field$) = sqrt(val(hkick$)**2 + val(vkick$)**2) * e_factor / val(l$)
      val(voltage$) = val(e_field$) * val(gap$) 
    endif
  endif

! EM_field

case (em_field$)
  if (val(rf_frequency$) /= 0) then
    val(rf_wavelength$) = c_light / val(rf_frequency$)
  else
    val(rf_wavelength$) = 0
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

  if (val(rf_frequency$) /= 0) then
    val(rf_wavelength$) = c_light / val(rf_frequency$)
  else
    val(rf_wavelength$) = 0
  endif

  ! multipass_slaves will inherit from lord
  if (ele%slave_status /= multipass_slave$) then
    if (val(rf_frequency$) /= 0 .and. ele%field_calc == bmad_standard$ .and. nint(ele%value(cavity_type$)) == standing_wave$) then
      val(l_active$) = 0.5_rp * val(rf_wavelength$) * nint(val(n_cell$))
    else
      val(l_active$) = val(l$)
    endif
  endif

  val(voltage_tot$)  = val(voltage$)  + val(voltage_err$)
  val(gradient_tot$) = val(gradient$) + val(gradient_err$)

! Patch

case (patch$)
  val(l$) = patch_length (ele)

! Quadrupole

case (quadrupole$)

! RFcavity

case (rfcavity$)
  if (associated(branch) .and. val(e_tot$) /= 0) then
    beta = ele%value(p0c$) / ele%value(e_tot$)
    time = branch%param%total_length / (c_light * beta)
    if (time /= 0) then
      if (is_true(val(harmon_master$))) then
        val(rf_frequency$) = val(harmon$) / time
      else
        val(harmon$) = val(rf_frequency$) * time
      endif
    endif
  endif

  if (val(rf_frequency$) /= 0) then
    val(rf_wavelength$) = c_light / val(rf_frequency$)
  else
    val(rf_wavelength$) = 0
  endif

  ! multipass_slaves will inherit from lord
  if (ele%slave_status /= multipass_slave$) then
    if (val(rf_frequency$) /= 0 .and. ele%field_calc == bmad_standard$ .and. nint(ele%value(cavity_type$)) == standing_wave$) then
      val(l_active$) = 0.5_rp * val(rf_wavelength$) * nint(val(n_cell$))
    else
      val(l_active$) = val(l$)
    endif
  endif

  if (val(voltage$) == 0) then
    val(gradient$) = 0
  elseif (val(l$) == 0) then
    val(gradient$) = 1d30    ! Something large
  else
    val(gradient$) = val(voltage$) / val(l$)
  endif

! rf_bend

case (rf_bend$)

  if (associated(branch) .and. val(e_tot$) /= 0) then
    beta = ele%value(p0c$) / ele%value(e_tot$)
    time = branch%param%total_length / (c_light * beta)
    if (time /= 0) then
      if (is_true(val(harmon_master$))) then
        val(rf_frequency$) = val(harmon$) / time
      else
        val(harmon$) = val(rf_frequency$) * time
      endif
    endif
  endif

  if (val(rf_frequency$) /= 0) then
    val(rf_wavelength$) = c_light / val(rf_frequency$)
  else
    val(rf_wavelength$) = 0
  endif

  ! multipass_slaves will inherit from lord
  if (ele%slave_status /= multipass_slave$) then
    if (val(rf_frequency$) /= 0 .and. ele%field_calc == bmad_standard$ .and. nint(ele%value(cavity_type$)) == standing_wave$) then
      val(l_active$) = 0.5_rp * val(rf_wavelength$) * nint(val(n_cell$))
    else
      val(l_active$) = val(l$)
    endif
  endif

  if (field_bookkeeping_doable) then
    val(angle$) = val(l$) * val(g$)

    if (val(g$) == 0) then
      val(rho$) = 0
      val(l_chord$) = val(l$)
      val(l_sagitta$) = 0
    else
      val(rho$) = 1 / val(g$)
      val(l_chord$) = 2 * val(rho$) * sin(val(angle$)/2)
      val(l_sagitta$) = -val(rho$) * cos_one(val(angle$)/2)
    endif

    select case (nint(val(fiducial_pt$)))
    case (none_pt$, center_pt$)
      val(l_rectangle$) = val(l_chord$)
    case default
      val(l_rectangle$) = sinc(val(angle$)) * val(l$)
    end select

    if (ele_value_has_changed(ele, [g$, rho$], [1e-10_rp, 1e-10_rp], .false.)) then
      call set_ele_status_stale (ele, floor_position_group$)
    endif
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

  if (field_bookkeeping_doable) then
    val(angle$) = val(l$) * val(g$)

    if (val(g$) == 0) then
      val(rho$) = 0
      val(l_chord$) = val(l$)
      val(l_sagitta$) = 0
    else
      val(rho$) = 1 / val(g$)
      val(l_chord$) = 2 * val(rho$) * sin(val(angle$)/2)
      val(l_sagitta$) = -val(rho$) * cos_one(val(angle$)/2)
    endif

    select case (nint(val(fiducial_pt$)))
    case (none_pt$, center_pt$)
      val(l_rectangle$) = val(l_chord$)
    case default
      val(l_rectangle$) = sinc(val(angle$)) * val(l$)
    end select

    if (ele_value_has_changed(ele, [g$, rho$], [1e-10_rp, 1e-10_rp], .false.)) then
      call set_ele_status_stale (ele, floor_position_group$)
    endif

    val(g_tot$)       = val(g$)       + val(dg$)
    val(b_field_tot$) = val(b_field$) + val(db_field$)
  endif

! Sol_quad

case (sol_quad$)

! Solenoid

case (solenoid$)

! Wiggler

case (wiggler$, undulator$) 

  ! Calculate b_max for fieldmap wigglers. 

  if (ele%field_calc == fieldmap$ .and. val(b_max$) == 0 .and. val(p0c$) > 0) then
    is_on = ele%is_on  ! Save
    polarity = val(polarity$)
    ele%is_on = .true.
    val(polarity$) = 1
    start%vec = 0
    val(b_max$) = 0
    n = nint(val(num_steps$))
    dl = val(l$) / n
    do i = 0, n
      b_max = -wiggler_field(i * dl)
      if (val(b_max$) > b_max) cycle
      val(b_max$) = b_max
      i_max = i      
    enddo
    if (i_max /= 0 .and. i_max /= n) then
      val(b_max$) = -super_brent((i_max-1)*dl, i_max*dl, (i_max+1)*dl, wiggler_field, 1.0e-8_rp, 0.0_rp, zmin, status)
    endif
    ele%is_on = is_on
    val(polarity$) = polarity
  endif

  val(g_max$) = 0
  val(k1x$) = 0
  val(k1y$) = 0

  if (val(p0c$) /= 0 .and. val(l_period$) /= 0) then
    val(g_max$) = c_light * val(b_max$) / val(p0c$)
    if (ele%field_calc == helical_model$) then
      val(k1x$) =  0.5 * val(g_max$)**2
      val(k1y$) =  0.5 * val(g_max$)**2
    elseif (ele%field_calc == planar_model$) then
      kz = twopi / val(l_period$)
      
      val(k1x$) =  0.5 * (val(g_max$) * val(kx$) / kz)**2
      val(k1y$) = -0.5 * val(g_max$)**2 * (kz**2 + val(kx$)**2) / kz**2
    endif
  endif

  if (val(l_period$) == 0) then
    val(n_period$) = 0
  else
    val(n_period$) = val(l$) / val(l_period$)
  endif

  val(osc_amplitude$) = val(g_max$) * (val(l_period$) / twopi)**2

end select

! It might be true that value == old_value here even though this was not true at the start
! of this routine. For example: A super_slave, via makeup_super_slave (which will be
! called by control_bookkeeper1 before calling this routine), will inherit its lords 
! num_steps value but then this routine will reset num_steps to be the correct value.

! So stop here if nothing has truely changed.

if (bmad_com%auto_bookkeeper) then
  if (all(val == ele%old_value)) return
endif

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

if (ele%tracking_method == taylor$ .and. ele%taylor_map_includes_offsets .and. &
        offset_nonzero .and. offset_changed .and. .not. non_offset_changed .and. &
        bmad_com%conserve_taylor_maps .and. ele%key /= patch$) then
  ele%taylor_map_includes_offsets = .false.
  if (associated(branch) .and. ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$) then
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

! Kill the taylor map if necessary.

if (non_offset_changed .or. (offset_changed .and. ele%taylor_map_includes_offsets)) then
  if (associated(ele%taylor(1)%term)) then
    call kill_taylor(ele%taylor)
  endif
  if (associated(ele%spin_taylor(0)%term)) then
    call kill_taylor(ele%spin_taylor)
    ele%spin_q(0,0) = real_garbage$
  endif
endif

! Make stale ele%rad_map if allocated

if (associated(ele%rad_map)) ele%rad_map%stale = .true.  ! Forces recalc

if (allocated(ele%multipole_cache)) then
  ele%multipole_cache%mag_valid = .false.
  ele%multipole_cache%elec_valid = .false.
endif

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

!----------------------------------------------------
contains

function wiggler_field (z, status) result (b_field)

real(rp), intent(in) :: z
real(rp) b_field
integer, optional :: status


call em_field_calc (ele, branch%param, z, start, .true., field, rf_time = 0.0_rp)
b_field = -norm2(field%b)

end function wiggler_field

end subroutine attribute_bookkeeper

