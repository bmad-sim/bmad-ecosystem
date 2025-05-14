!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)
!
! Routine to mark an element as modified for use with "intelligent" bookkeeping.
!
! This routine is overloaded by set_flags_for_changed_attribute. 
! See set_flags_for_changed_attribute for more details.
!-

subroutine set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)

use bmad_routine_interface, dummy => set_flags_for_changed_all_attribute

implicit none

type (ele_struct), target :: ele
type (all_pointer_struct) all_attrib
logical, optional :: set_dependent

!

if (associated(all_attrib%r)) call set_flags_for_changed_real_attribute(ele, all_attrib%r, set_dependent)
if (associated(all_attrib%i)) call set_flags_for_changed_integer_attribute(ele, all_attrib%i, set_dependent)
if (associated(all_attrib%l)) call set_flags_for_changed_logical_attribute(ele, all_attrib%l, set_dependent)

end subroutine set_flags_for_changed_all_attribute

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_integer_attribute (ele, attrib, set_dependent)
!
! Routine to mark an element as modified for use with "intelligent" bookkeeping.
!
! This routine is overloaded by set_flags_for_changed_attribute. 
! See set_flags_for_changed_attribute for more details.
!-

subroutine set_flags_for_changed_integer_attribute (ele, attrib, set_dependent)

use bmad_routine_interface, dummy => set_flags_for_changed_integer_attribute

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: slave, lord

integer, target :: attrib
integer, pointer :: a_ptr
integer i, ix_pass

real(rp) real_dummy

logical, optional :: set_dependent

! This will set some generic flags

call set_flags_for_changed_real_attribute (ele, real_dummy, set_dependent)

!

a_ptr => attrib

if (ele%key /= taylor$ .or. ele%tracking_method /= taylor$) then
  if (associated(a_ptr, ele%spin_tracking_method)) then
    call kill_taylor (ele%spin_taylor)
  endif
endif

!-------------------------------------------------------------------

select case (ele%key)
case (rfcavity$, lcavity$, e_gun$)
  if (associated(a_ptr, ele%tracking_method) .or. associated(a_ptr, ele%field_calc)) then
    call set_ele_status_stale (ele, ref_energy_group$)
  endif
end select

! Set independent stuff in multipass lord

if (ele%lord_status == multipass_lord$) then
  do i = 1, ele%n_slave
    slave => pointer_to_slave(ele, i)
  
    if (associated(a_ptr, ele%aperture_at)) then
      slave%aperture_at = a_ptr
    elseif (associated(a_ptr, ele%ptc_integration_type)) then
      slave%ptc_integration_type = a_ptr
    elseif (associated(a_ptr, ele%aperture_type)) then
      slave%aperture_type = a_ptr
    elseif (associated(a_ptr, ele%mat6_calc_method)) then
      slave%mat6_calc_method = a_ptr
    elseif (associated(a_ptr, ele%tracking_method)) then
      slave%tracking_method = a_ptr
    elseif (associated(a_ptr, ele%spin_tracking_method)) then
      slave%spin_tracking_method = a_ptr
    elseif (associated(a_ptr, ele%field_calc)) then
      slave%field_calc = a_ptr
    elseif (associated(a_ptr, ele%csr_method)) then
      slave%csr_method = a_ptr
    elseif (associated(a_ptr, ele%space_charge_method)) then
      slave%space_charge_method = a_ptr
    else
      exit
    endif
  enddo
endif

!-------------------------------------------------------------------

if (ele%slave_status == multipass_slave$) then
  lord => pointer_to_multipass_lord(ele, ix_pass)
  if (ix_pass == 1) then
    if (associated(a_ptr, ele%space_charge_method)) then
      lord%space_charge_method = ele%space_charge_method
    endif
  endif
endif

end subroutine set_flags_for_changed_integer_attribute

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_logical_attribute (ele, attrib, set_dependent)
!
! Routine to mark an element as modified for use with "intelligent" bookkeeping.
!
! This routine is overloaded by set_flags_for_changed_attribute. 
! See set_flags_for_changed_attribute for more details.
!-

subroutine set_flags_for_changed_logical_attribute (ele, attrib, set_dependent)

use bmad_routine_interface, dummy => set_flags_for_changed_logical_attribute

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: slave

integer i

real(rp) real_dummy, f

logical, target :: attrib
logical, pointer :: a_ptr
logical, optional :: set_dependent

! Call to set_flags_for_changed_real_attribute will set some generic flags

a_ptr => attrib

if (ele%key /= taylor$ .or. ele%tracking_method /= taylor$) then
  if (associated(a_ptr, ptc_com%exact_model) .or. associated(a_ptr, ptc_com%exact_misalign)) then
    call kill_taylor (ele%spin_taylor)
    call kill_taylor (ele%taylor)
  endif
endif

!-------------------------------------------------------------------

if (associated(a_ptr, ele%field_master)) then
  f = ele%value(p0c$) / (charge_of(ele%ref_species) * c_light)
  if (ele%key == multipole$) then
    if (attrib) then
      ele%a_pole = ele%a_pole * f
    else
      ele%a_pole = ele%a_pole / f
    endif

  elseif (ele%key == ab_multipole$) then
    if (attrib) then
      ele%a_pole = ele%a_pole * f
      ele%b_pole = ele%b_pole * f
    else
      ele%a_pole = ele%a_pole / f
      ele%b_pole = ele%b_pole / f
    endif
  endif

else
  call set_flags_for_changed_real_attribute (ele, real_dummy, set_dependent)
endif

! Set independent stuff in multipass lord

if (ele%lord_status == multipass_lord$) then 

  do i = 1, ele%n_slave
    slave => pointer_to_slave(ele, i)
  
    if (associated(a_ptr, ele%offset_moves_aperture)) then
      slave%offset_moves_aperture = attrib
    else
      exit
    endif
  enddo

endif

end subroutine set_flags_for_changed_logical_attribute

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_lat_attribute (lat, set_dependent)
!
! Routine to mark a lattice as modified for use with "intelligent" bookkeeping.
!
! This routine is overloaded by set_flags_for_changed_attribute. 
! See set_flags_for_changed_attribute for more details.
!-

subroutine set_flags_for_changed_lat_attribute (lat, set_dependent)

use bmad_routine_interface, dummy => set_flags_for_changed_lat_attribute

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch

integer i, j
logical, optional :: set_dependent

!

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  call set_status_flags (branch%param%bookkeeping_state, stale$)
  do j = 0, ubound(branch%ele, 1)
    call set_status_flags (branch%ele(j)%bookkeeping_state, stale$)
  enddo
enddo

end subroutine set_flags_for_changed_lat_attribute

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_real_attribute (ele, attrib, set_dependent)
!
! Routine to mark an element as modified for use with "intelligent" bookkeeping.
!
! This routine is overloaded by set_flags_for_changed_attribute. 
! See set_flags_for_changed_attribute for more details.
!-

subroutine set_flags_for_changed_real_attribute (ele, attrib, set_dependent)

use bmad_routine_interface, dummy => set_flags_for_changed_real_attribute

implicit none

type (ele_struct), target :: ele
type (branch_struct), pointer :: branch
type (grid_field_struct), pointer :: g_field
type (cylindrical_map_struct), pointer :: cl_map

real(rp), optional, target :: attrib
real(rp), pointer :: a_ptr
real(rp) p0c_factor, ff, rel_p1
real(rp), target :: unknown_attrib, dangle

integer i

logical coupling_change, found, dep_set, dep2_set
logical, optional :: set_dependent

!-------------------
! For a particular elemement...

branch => pointer_to_branch(ele)
dep_set = logic_option(.true., set_dependent)
dep2_set = (dep_set .and. ele%value(p0c$) /= 0 .and. charge_of(branch%param%particle, 0) /= 0)
p0c_factor = 0
if (dep2_set) p0c_factor = ele%value(p0c$) / (c_light * charge_of(branch%param%particle))

! If a lord then set the control flag stale

if (ele%lord_status /= not_a_lord$) call set_ele_status_stale (ele, control_group$)

! Groups and overlays do not have any dependent attributes. 
! For all others set the attributes flag stale.

if (ele%key /= group$ .and. ele%key /= overlay$ .and. dep_set) then
  call set_ele_status_stale (ele, attribute_group$)
endif

! Transfer matrix calc needs to be flagged

if (ele%key /= overlay$ .and. ele%key /= group$ .and. ele%lord_status /= multipass_lord$) then
  call set_ele_status_stale (ele, mat6_group$)
endif

! If attrib is not present then point to a dummy location which will not match when 
! the associated() function is used below.

if (.not. present(attrib)) then
  call set_ele_status_stale (ele, all_groups$)
endif

! Use a_ptr with the associated function to see which attribute has been changed.

a_ptr => attrib

! A limit change does not need any bookkeeping

if (associated(a_ptr, ele%value(x1_limit$)) .or. associated(a_ptr, ele%value(x2_limit$)) .or. &
    associated(a_ptr, ele%value(y1_limit$)) .or. associated(a_ptr, ele%value(y2_limit$))) return

! Associated taylor map

if (ele%key /= taylor$ .and. associated(ele%taylor(1)%term)) then
  call kill_taylor (ele%taylor)
  call kill_taylor (ele%spin_taylor)
endif

! delta_ref_time change

if (associated(a_ptr, ele%value(delta_ref_time$))) then
  call set_ele_status_stale (ele, ref_energy_group$)  ! Energy & time
  return
endif

! A length change involves changes in the floor position.

if (associated(a_ptr, ele%value(l$))) then
  if (ele%key /= overlay$ .and. ele%key /= group$) then
    call set_ele_status_stale (ele, s_and_floor_position_group$)
    call set_ele_status_stale (ele, floor_position_group$)
    call set_ele_status_stale (ele, ref_energy_group$)
  endif
  if (ele%value(p0c$) /= ele%value(p0c_start$)) call set_ele_status_stale (ele, ref_energy_group$)
  return
endif

! E_tot and p0c can be varied in an init_ele or a multipass lord with multipass_ref_energy = user_set$.
! In addition, for an init_ele, must also set e_tot_start and p0c_start. This is important
! for lattices with an e_gun element

if (associated(a_ptr, ele%value(e_tot$)) .and. associated(branch)) then
  select case (ele%key)
  case (lcavity$, rfcavity$, e_gun$, em_field$)
    call set_ele_status_stale (ele, ref_energy_group$, .true.)
  case default
    ! Lord energy is set from slave. Not other way around.
    call set_ele_status_stale (ele, ref_energy_group$, .false.)
  end select

  if (dep2_set) then
    call convert_total_energy_to (ele%value(e_tot$), branch%param%particle, pc = ele%value(p0c$))
    ! If there is an e_gun then actually want to vary the start energy and e_tot/p0c will be
    ! dependent parameters dependent upon the gun voltage and starting energy.
    if (ele%key /= lcavity$ .and. ele%key /= e_gun$) then
      ele%value(e_tot_start$) = ele%value(e_tot$)
      ele%value(p0c_start$) = ele%value(p0c$)
    endif
  endif
  return
endif

if (associated(a_ptr, ele%value(p0c$)) .and. associated(branch)) then
  select case (ele%key)
  case (lcavity$, rfcavity$, e_gun$, em_field$)
    call set_ele_status_stale (ele, ref_energy_group$, .true.)
  case default
    ! Lord energy is set from slave. Not other way around.
    call set_ele_status_stale (ele, ref_energy_group$, .false.)
  end select

  if (dep2_set) then
    call convert_pc_to (ele%value(p0c$), branch%param%particle, e_tot = ele%value(e_tot$))
    ! If there is an e_gun then actually want to vary the start energy and e_tot/p0c will be
    ! dependent parameters dependent upon the gun voltage and starting energy.
    if (ele%key /= lcavity$ .and. ele%key /= e_gun$) then
      ele%value(e_tot_start$) = ele%value(e_tot$)
      ele%value(p0c_start$) = ele%value(p0c$)
    endif
  endif
  return
endif

if (associated(a_ptr, ele%value(e_tot_start$)) .and. associated(branch)) then
  select case (ele%key)
  case (lcavity$, rfcavity$, e_gun$, em_field$)
    call set_ele_status_stale (ele, ref_energy_group$, .true.)
  case default
    ! Lord energy is set from slave. Not other way around.
    call set_ele_status_stale (ele, ref_energy_group$, .false.)
  end select
  if (dep2_set) then
    call convert_total_energy_to (ele%value(e_tot_start$), branch%param%particle, pc = ele%value(p0c_start$))
  endif
  return
endif

if (associated(a_ptr, ele%value(p0c_start$)) .and. associated(branch)) then
  select case (ele%key)
  case (lcavity$, rfcavity$, e_gun$, em_field$)
    call set_ele_status_stale (ele, ref_energy_group$, .true.)
  case default
    ! Lord energy is set from slave. Not other way around.
    call set_ele_status_stale (ele, ref_energy_group$, .false.)
  end select
  if (dep2_set) then
    call convert_pc_to (ele%value(p0c_start$), branch%param%particle, e_tot = ele%value(e_tot_start$))
  endif
  return
endif

!

if (associated(a_ptr, ele%value(ds_step$)) .or. associated(a_ptr, ele%value(num_steps$))) then
  if (dep_set .and. ele%value(num_steps$) /= 0 .and. associated(a_ptr, ele%value(num_steps$))) then
     ele%value(ds_step$) = abs(ele%value(l$)) / ele%value(num_steps$)
  endif

  if (ele%key == e_gun$ .or. ele%key == lcavity$ .or. &
                      (ele%key == em_field$ .and. is_false(ele%value(constant_ref_energy$)))) then
    call set_ele_status_stale (ele, ref_energy_group$)
  endif

  if (ele%lord_status == multipass_lord$) then
    call this_multipass_slave_xfer(ele, num_steps$)
    call this_multipass_slave_xfer(ele, ds_step$)
  endif
  return
endif

if (dep2_set .and. has_hkick_attributes(ele%key) .and. ele%key /= elseparator$) then
  if (associated(a_ptr, ele%value(bl_hkick$))) then
    ele%value(hkick$) = ele%value(bl_hkick$) / p0c_factor
    return
  elseif (associated(a_ptr, ele%value(bl_vkick$))) then
    ele%value(vkick$) = ele%value(bl_vkick$) / p0c_factor
    return
  elseif (associated(a_ptr, ele%value(hkick$))) then
    ele%value(bl_hkick$) = ele%value(hkick$) * p0c_factor
    return
  elseif (associated(a_ptr, ele%value(vkick$))) then
    ele%value(bl_vkick$) = ele%value(vkick$) * p0c_factor
    return
  endif
endif

!------------------------------------------------
! By element type

select case (ele%key)

! BeamBeam

case (beambeam$)
  if (dep2_set) then
    if (associated(a_ptr, ele%value(ks$))) then
      ele%value(bs_field$) = ele%value(ks$) * p0c_factor
    elseif (associated(a_ptr, ele%value(bs_field$))) then
      ele%value(ks$) = ele%value(bs_field$) / p0c_factor
    endif

    if (associated(a_ptr, ele%value(species_strong$))) then
      if (ele%value(pc_strong$) >= 0) then
        call convert_pc_to(ele%value(pc_strong$), nint(ele%value(species_strong$)), E_tot = ele%value(E_tot_strong$))
      elseif (ele%value(E_tot_strong$) > 0) then
        call convert_total_energy_to(ele%value(E_tot_strong$), nint(ele%value(species_strong$)), pc = ele%value(pc_strong$))
      endif

    elseif (ele%value(species_strong$) /= real_garbage$) then
      if (associated(a_ptr, ele%value(E_tot_strong$))) then
        call convert_total_energy_to(ele%value(E_tot_strong$), nint(ele%value(species_strong$)), pc = ele%value(pc_strong$))
      elseif (associated(a_ptr, ele%value(pc_strong$))) then
        call convert_pc_to(ele%value(pc_strong$), nint(ele%value(species_strong$)), E_tot = ele%value(E_tot_strong$))
      endif
    endif
  endif

! Beginning_ele

case (beginning_ele$) 
  coupling_change = .false.

  if (associated(a_ptr, ele%a%beta) .or. associated(a_ptr, ele%a%alpha)) then
    if (dep_set) then
      if (ele%a%beta /= 0) ele%a%gamma = (1 + ele%a%alpha**2) / ele%a%beta
    endif
    return
  endif

  if (associated(a_ptr, ele%b%beta) .or. associated(a_ptr, ele%b%alpha)) then
    if (dep_set) then
      if (ele%b%beta /= 0) ele%b%gamma = (1 + ele%b%alpha**2) / ele%b%beta
    endif
    return
  endif

  if (associated(a_ptr, ele%c_mat(1,1)) .or. associated(a_ptr, ele%c_mat(1,2)) .or. & 
      associated(a_ptr, ele%c_mat(2,1)) .or. associated(a_ptr, ele%c_mat(2,2))) then
    if (dep_set) then
      ele%gamma_c = sqrt(1 - ele%c_mat(1,1)*ele%c_mat(2,2) + ele%c_mat(1,2)*ele%c_mat(2,1))
      coupling_change = .true.
    endif
  endif

  if (associated(a_ptr, ele%x%deta_ds) .or. associated(a_ptr, ele%y%deta_ds)) then
    ele%value(deta_ds_master$) = true$
    rel_p1 = 1 + ele%map_ref_orb_out%vec(6)
    ele%x%etap = ele%x%deta_ds * rel_p1 + ele%map_ref_orb_out%vec(2) / rel_p1
    ele%y%etap = ele%y%deta_ds * rel_p1 + ele%map_ref_orb_out%vec(4) / rel_p1
    coupling_change = .true.
  endif

  if (associated(a_ptr, ele%x%etap) .or. associated(a_ptr, ele%y%etap)) then
    ele%value(deta_ds_master$) = false$
    rel_p1 = 1 + ele%map_ref_orb_out%vec(6)
    ele%x%deta_ds = ele%x%etap / rel_p1 - ele%map_ref_orb_out%vec(2) / rel_p1**2
    ele%y%deta_ds = ele%y%etap / rel_p1 - ele%map_ref_orb_out%vec(4) / rel_p1**2
    coupling_change = .true.
  endif

  if (associated(a_ptr, ele%x%eta) .or. associated(a_ptr, ele%y%eta) .or. coupling_change) then
    if (dep_set) then
      call normal_mode_dispersion(ele)
    endif
    return
  endif

  if (associated(a_ptr, ele%a%eta) .or. associated(a_ptr, ele%a%etap) .or. &
      associated(a_ptr, ele%b%eta) .or. associated(a_ptr, ele%b%etap)) then 
    if (dep_set) then
      call normal_mode_dispersion(ele, .true.)
    endif
    return
  endif

  if (associated(a_ptr, ele%floor%r(1)) .or. associated(a_ptr, ele%floor%r(2)) .or. &
      associated(a_ptr, ele%floor%r(3)) .or. associated(a_ptr, ele%floor%theta) .or. &
      associated(a_ptr, ele%floor%phi) .or. associated(a_ptr, ele%floor%psi) .or. &
      associated(a_ptr, ele%s)) then
    call set_ele_status_stale (ele, s_and_floor_position_group$)
    call floor_angles_to_w_mat (ele%floor%theta, ele%floor%phi, ele%floor%psi, ele%floor%w)
    ele%s_start = ele%s
    return
  endif

  if (associated(a_ptr, ele%ref_time)) then
    call set_ele_status_stale (ele, ref_energy_group$)
  endif

! Converter

case (converter$)

! Crystal

case (crystal$)
  if (associated(a_ptr, ele%value(graze_angle_in$)) .or. associated(a_ptr, ele%value(graze_angle_out$))) then
    call set_ele_status_stale (ele, floor_position_group$)
    return
  endif

! Mirror, multilayer_mirror

case (mirror$, multilayer_mirror$)
  if (associated(a_ptr, ele%value(graze_angle$))) then
    call set_ele_status_stale (ele, floor_position_group$)
    return
  endif

! fork, photon_fork

case (fork$, photon_fork$)

! hkicker, vkicker

case (hkicker$, vkicker$)
  if (dep2_set) then
    if (associated(a_ptr, ele%value(kick$))) then
      ele%value(bl_kick$) = ele%value(kick$) * p0c_factor
    elseif (associated(a_ptr, ele%value(bl_kick$))) then
      ele%value(kick$) = ele%value(bl_kick$) / p0c_factor
    endif
  endif

! rfcavity, crab_cavity

case (rfcavity$, crab_cavity$)
  if (dep_set) then
    if (associated(a_ptr, ele%value(harmon$))) then
      if (ele%value(e_tot$) /= 0 .and. branch%param%total_length /= 0) then
        ff = c_light * ele%value(p0c$) / (branch%param%total_length * ele%value(e_tot$))
        ele%value(rf_frequency$) = ele%value(harmon$) * ff
      else
        ele%value(rf_frequency$) = 0
      endif
    elseif (associated(a_ptr, ele%value(rf_frequency$))) then
      if (ele%value(p0c$) /= 0) then
        ff = branch%param%total_length * ele%value(e_tot$) / (c_light * ele%value(p0c$))
        ele%value(harmon$) = ele%value(rf_frequency$) * ff
      else
        ele%value(harmon$) = 0
      endif
    elseif (ele%value(l$) /= 0) then
      if (associated(a_ptr, ele%value(voltage$))) then
        ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)
      elseif (associated(a_ptr, ele%value(gradient$))) then
        ele%value(voltage$) = ele%value(gradient$) * ele%value(l$)
      endif
    endif
  endif

! lcavity, e_gun

case (lcavity$, e_gun$)
  if (associated(a_ptr, ele%value(gradient$)) .or. associated(a_ptr, ele%value(phi0$)) .or. &
      associated(a_ptr, ele%value(voltage$)) .or. associated(a_ptr, ele%value(rf_frequency$)) .or. &
      associated(a_ptr, ele%value(phi0_autoscale$)) .or. associated(a_ptr, ele%value(field_autoscale$))) then
    call set_ele_status_stale (ele, ref_energy_group$)
  endif

  if (dep_set .and. ele%value(l$) /= 0) then
    if (associated(a_ptr, ele%value(voltage$))) then
      ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)
      ele%value(voltage_tot$) = ele%value(voltage$) + ele%value(voltage_err$)
      ele%value(gradient_tot$) = ele%value(gradient$) + ele%value(gradient_err$)
    elseif (associated(a_ptr, ele%value(voltage_err$))) then
      ele%value(gradient_err$) = ele%value(voltage_err$) / ele%value(l$)
      ele%value(voltage_tot$) = ele%value(voltage$) + ele%value(voltage_err$)
      ele%value(gradient_tot$) = ele%value(gradient$) + ele%value(gradient_err$)
    elseif (associated(a_ptr, ele%value(gradient$))) then
      ele%value(voltage$) = ele%value(gradient$) * ele%value(l$)
      ele%value(voltage_tot$) = ele%value(voltage$) + ele%value(voltage_err$)
      ele%value(gradient_tot$) = ele%value(gradient$) + ele%value(gradient_err$)
    elseif (associated(a_ptr, ele%value(gradient_err$))) then
      ele%value(voltage_err$) = ele%value(gradient_err$) * ele%value(l$)
      ele%value(voltage_tot$) = ele%value(voltage$) + ele%value(voltage_err$)
      ele%value(gradient_tot$) = ele%value(gradient$) + ele%value(gradient_err$)
    endif
  endif

  if (ele%key == lcavity$) then 
    if (associated(a_ptr, ele%value(phi0_multipass$)) .or. associated(a_ptr, ele%value(e_loss$))) then
       call set_ele_status_stale (ele, ref_energy_group$)
    endif
  endif

  found = .false.

  if (associated(ele%cylindrical_map)) then
    do i = 1, size(ele%cylindrical_map)
      cl_map => ele%cylindrical_map(i)
      if (associated(a_ptr, cl_map%phi0_fieldmap)) found = .true.
      if (associated(a_ptr, cl_map%field_scale)) found = .true.
      if (cl_map%master_parameter > 0) found = (found .or. associated(a_ptr, ele%value(cl_map%master_parameter)))
      if (associated(a_ptr, cl_map%phi0_fieldmap)) found = .true.
    enddo
  endif

  if (associated(ele%grid_field)) then
    do i = 1, size(ele%grid_field)
      g_field => ele%grid_field(i)
      if (associated(a_ptr, g_field%phi0_fieldmap)) found = .true.
      if (associated(a_ptr, g_field%field_scale)) found = .true.
      if (g_field%master_parameter > 0) found = (found .or. associated(a_ptr, ele%value(g_field%master_parameter)))
      if (associated(a_ptr, g_field%phi0_fieldmap)) found = .true.
    enddo
  endif

  if (found) call set_ele_status_stale (ele, ref_energy_group$)

! Patch

case (patch$)
  ! Any attribute change will shift the reference time.
  call set_ele_status_stale (ele, ref_energy_group$)
  call set_ele_status_stale (ele, floor_position_group$)

! Quadrupole

case (quadrupole$)
  if (dep2_set) then
    if (associated(a_ptr, ele%value(k1$))) then
      ele%value(b1_gradient$) = ele%value(k1$) * p0c_factor
    elseif (associated(a_ptr, ele%value(b1_gradient$))) then
      ele%value(k1$) = ele%value(b1_gradient$) / p0c_factor
    endif
  endif

! Floor_shift, fiducial

case (floor_shift$, fiducial$)
  call set_ele_status_stale (ele, floor_position_group$)

! Octupole

case (octupole$)
  if (dep2_set) then
    if (associated(a_ptr, ele%value(k3$))) then
      ele%value(b3_gradient$) = ele%value(k3$) * p0c_factor
    elseif (associated(a_ptr, ele%value(b3_gradient$))) then
      ele%value(k3$) = ele%value(b3_gradient$) / p0c_factor
    endif
  endif

! Sad_mult

case (sad_mult$)
  if (dep2_set) then
    if (associated(a_ptr, ele%value(ks$))) then
      ele%value(bs_field$) = ele%value(ks$) * p0c_factor
    elseif (associated(a_ptr, ele%value(bs_field$))) then
      ele%value(ks$) = ele%value(bs_field$) / p0c_factor
    endif
  endif

! Sbend

case (sbend$, rf_bend$)
  if (associated(a_ptr, ele%value(angle$)) .or. associated(a_ptr, ele%value(g$)) .or. &
      associated(a_ptr, ele%value(rho$)) .or. associated(a_ptr, ele%value(b_field$))) then
    call set_ele_status_stale (ele, floor_position_group$)
  endif

  ! Attribute_bookkeeper takes care of some stuff so just have to make sure the final attribute state
  ! is independent of the setting for %field_master.
  if (dep_set) then
    if (associated(a_ptr, ele%value(angle$))) then
      call bend_angle_fiducial(ele, dep2_set, p0c_factor)

    elseif (associated(a_ptr, ele%value(l_rectangle$))) then
      select case (nint(ele%value(fiducial_pt$)))
      case (none_pt$)
        ele%value(l$) = ele%value(l_rectangle$) * asinc(0.5_rp * ele%value(l_rectangle$) * ele%value(g$))
      case (center_pt$)
        ele%value(l$) = ele%value(l_rectangle$) * asinc(0.5_rp * ele%value(l_rectangle$) * ele%value(g$))
        dangle = ele%value(l$) * ele%value(g$) - ele%value(angle$)
        ele%value(e1$) = ele%value(e1$) + 0.5_rp * dangle
        ele%value(e2$) = ele%value(e2$) + 0.5_rp * dangle
      case (entrance_end$)
        ele%value(l$) = ele%value(l_rectangle$) * asinc(ele%value(l_rectangle$) * ele%value(g$))
        dangle = ele%value(l$) * ele%value(g$) - ele%value(angle$)
        ele%value(e2$) = ele%value(e2$) + dangle
      case (exit_end$)
        ele%value(l$) = ele%value(l_rectangle$) * asinc(ele%value(l_rectangle$) * ele%value(g$))
        dangle = ele%value(l$) * ele%value(g$) - ele%value(angle$)
        ele%value(e1$) = ele%value(e1$) + dangle
      end select

    elseif (associated(a_ptr, ele%value(rho$))) then
      if (ele%value(rho$) /= 0) then
        ele%value(g$) = 1.0_rp / ele%value(rho$)
        call bend_g_fiducial(ele, dep2_set, p0c_factor)
      endif

    elseif (associated(a_ptr, ele%value(g$))) then
      ele%value(b_field$) = ele%value(g$) * p0c_factor 
      call bend_g_fiducial(ele, dep2_set, p0c_factor)

    elseif (dep2_set) then
      if (associated(a_ptr, ele%value(b_field$))) then
        ele%value(g$) = ele%value(b_field$) / p0c_factor
        call bend_g_fiducial(ele, dep2_set, p0c_factor)

      elseif (associated(a_ptr, ele%value(db_field$))) then
        ele%value(dg$) = ele%value(db_field$) / p0c_factor

      elseif (associated(a_ptr, ele%value(dg$))) then
        ele%value(db_field$) = ele%value(dg$) * p0c_factor

      elseif (associated(a_ptr, ele%value(k1$))) then
        ele%value(b1_gradient$) = ele%value(k1$) * p0c_factor

      elseif (associated(a_ptr, ele%value(b1_gradient$))) then
        ele%value(k1$) = ele%value(b1_gradient$) / p0c_factor

      elseif (associated(a_ptr, ele%value(k2$))) then
        ele%value(b2_gradient$) = ele%value(k2$) * p0c_factor

      elseif (associated(a_ptr, ele%value(b2_gradient$))) then
        ele%value(k2$) = ele%value(b2_gradient$) / p0c_factor
      endif
    endif

    ele%value(g_tot$) = ele%value(g$) + ele%value(dg$)
    ele%value(b_field_tot$) = ele%value(b_field$) + ele%value(db_field$)
  endif

! Sextupole

case (sextupole$)
  if (dep2_set) then
    if (associated(a_ptr, ele%value(k2$))) then
      ele%value(b2_gradient$) = ele%value(k2$) * p0c_factor
    elseif (associated(a_ptr, ele%value(b2_gradient$))) then
      ele%value(k2$) = ele%value(b2_gradient$) / p0c_factor
    endif
  endif

! Sol_Quad

case (sol_quad$)
  if (dep2_set) then
    if (associated(a_ptr, ele%value(ks$))) then
      ele%value(bs_field$) = ele%value(ks$) * p0c_factor
    elseif (associated(a_ptr, ele%value(bs_field$))) then
      ele%value(ks$) = ele%value(bs_field$) / p0c_factor
    elseif (associated(a_ptr, ele%value(k1$))) then
      ele%value(b1_gradient$) = ele%value(k1$) * p0c_factor
    elseif (associated(a_ptr, ele%value(b1_gradient$))) then
      ele%value(k1$) = ele%value(b1_gradient$) / p0c_factor
    endif
  endif

! Solenoid

case (solenoid$)
  if (dep2_set) then
    if (associated(a_ptr, ele%value(ks$))) then
      ele%value(bs_field$) = ele%value(ks$) * p0c_factor
    elseif (associated(a_ptr, ele%value(bs_field$))) then
      ele%value(ks$) = ele%value(bs_field$) / p0c_factor
    endif
  endif

end select

!--------------------------------------------------------------
contains

subroutine this_multipass_slave_xfer(lord, ix_attrib)

type (ele_struct) lord
type (ele_struct), pointer :: slave

integer ix_attrib, i

!

do i = 1, lord%n_slave
  slave => pointer_to_slave(lord, i)
  slave%value(ix_attrib) = lord%value(ix_attrib)
enddo

end subroutine this_multipass_slave_xfer

!--------------------------------------------------------------
! contains

subroutine bend_g_fiducial(ele, dep2_set, p0c_factor)

type (ele_struct) ele
real(rp) p0c_factor, len1, len2
logical dep2_set

! There has been a change in g (and not any other parameters) so need to calculate shifts in L, e1, and e2

select case (nint(ele%value(fiducial_pt$)))
case (none_pt$)

case (center_pt$)
  call bend_g_fiducial2(ele, 0.5_rp, ele%value(e1$), len1)
  call bend_g_fiducial2(ele, 0.5_rp, ele%value(e2$), len2)
  ele%value(l$) = len1 + len2


case (entrance_end$)
  call bend_g_fiducial2(ele, 1.0_rp, ele%value(e2$), ele%value(l$))

case (exit_end$)
  call bend_g_fiducial2(ele, 1.0_rp, ele%value(e1$), ele%value(l$))
end select

if (dep2_set) ele%value(b_field$) = ele%value(g$) * p0c_factor 

end subroutine bend_g_fiducial

!--------------------------------------------------------------
! contains

subroutine bend_g_fiducial2(ele, factor, e_edge, len_out, angle0_in)

type (ele_struct) ele
real(rp) factor, e_edge, len_out, g0, angle0, rp2x, r1(2), theta1, l_rect, length, a, b, c
real(rp), optional :: angle0_in

!

g0 = ele%value(l$) * ele%value(angle$)  ! Old g
if (present(angle0_in)) then
  angle0 = angle0_in                    ! Old angle
else
  angle0 = factor * ele%value(angle$)   ! Old angle
endif
theta1 = angle0 - e_edge

r1 = factor * [-ele%value(l_rectangle$), -cos_one(angle0) * ele%value(l$)]
r1 = rot_2d(r1, theta1)

a = ele%value(g$)
b = 2.0_rp * cos(theta1)
c = ele%value(g$) * r1(1)**2 + 2.0_rp * r1(1) * sin(theta1)
rp2x = r1(2) - 2.0_rp * c / (b + sqrt(b*b - 4.0_rp*a*c))

l_rect = factor * ele%value(l_rectangle$) + rp2x * sin(theta1)
length = asinc(ele%value(g$) * l_rect) * l_rect

e_edge = e_edge + length * ele%value(g$) - angle0
len_out = length

end subroutine bend_g_fiducial2

!--------------------------------------------------------------
! contains

subroutine bend_angle_fiducial(ele, dep2_set, p0c_factor)

type (ele_struct) ele
real(rp) p0c_factor, angle0, len1, len2
logical dep2_set

! There has been a change in angle (and not any other parameters) so need to calculate shifts in L, e1, and e2

if (ele%value(l$) == 0) return

select case (nint(ele%value(fiducial_pt$)))
case (none_pt$)
  ele%value(g$) = ele%value(angle$) / ele%value(l$)

case (center_pt$)
  call bend_angle_fiducial2(ele, 0.5_rp, ele%value(e2$), len1)
  call bend_angle_fiducial2(ele, 0.5_rp, ele%value(e2$), len2)
  ele%value(l$) = len1 + len2
  ele%value(g$) = ele%value(angle$) / ele%value(l$)

case (entrance_end$)
  call bend_angle_fiducial2(ele, 1.0_rp, ele%value(e2$), ele%value(l$))

case (exit_end$)
  call bend_angle_fiducial2(ele, 1.0_rp, ele%value(e1$), ele%value(l$))
end select

if (dep2_set) ele%value(b_field$) = ele%value(g$) * p0c_factor 

end subroutine bend_angle_fiducial

!--------------------------------------------------------------
! contains

subroutine bend_angle_fiducial2(ele, factor, e_edge, len_out)

type (ele_struct) ele
real(rp) factor, e_edge, angle0, theta1, r1(2), len_out

! fiducial_pt = center is not considered here.

angle0 = factor * ele%value(g$) * ele%value(l$)   ! Old angle
theta1 = angle0 - e_edge

r1 = factor * [-ele%value(l_rectangle$), -cos_one(angle0) * ele%value(l$)]
r1 = rot_2d(r1, theta1)

ele%value(g$) = -(sin(theta1) + sin(ele%value(angle$) - theta1)) / r1(1)

call bend_g_fiducial2(ele, factor, e_edge, len_out, angle0)

end subroutine bend_angle_fiducial2

end subroutine set_flags_for_changed_real_attribute
