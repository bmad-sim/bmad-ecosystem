module bookkeeper_mod

use wall3d_mod
use equality_mod
use expression_mod
use attribute_mod

implicit none

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_group_lord (lat, lord, err_flag)
!
! Subroutine to calculate the attributes of group slave elements.
! This routine is private to bookkeeper_mod.
!-

Subroutine makeup_group_lord (lat, lord, err_flag)  

type (lat_struct), target :: lat
type (ele_struct) :: lord
type (ele_struct), pointer :: slave, slave2
type (control_struct), pointer :: control

real(rp), pointer :: r_ptr

integer i, j

logical moved, err_flag

character(40) :: attrib_name
character(20) :: r_name = 'makeup_group_lord'

!

err_flag = .false.
moved = .false.   ! have we longitudinally moved an element?
if (.not. lord%is_on) return

do i = 1, lord%n_slave
  slave => pointer_to_slave(lord, i, control)
  attrib_name = control%attribute
  if (attrib_name == 'L') moved = .true.

  select case (attrib_name)

  !---------
  ! Edge: Varying lengths takes special code.

  case ('START_EDGE', 'END_EDGE', 'S_POSITION', 'ACCORDION_EDGE', 'L')

    if (slave%lord_status == multipass_lord$) then
      do j = 1, slave%n_slave
        slave2 => pointer_to_slave (slave, j)
        call change_this_edge (slave2, control);  if (err_flag) return
      enddo
    else
      call change_this_edge (slave, control);  if (err_flag) return
    endif

  !---------
  ! x_limit, y_limit, aperture

  case ('X_LIMIT')
    call group_change_this (slave, 'X1_LIMIT', control, 1);  if (err_flag) return
    call group_change_this (slave, 'X2_LIMIT', control, 1);  if (err_flag) return

  case ('Y_LIMIT')
    call group_change_this (slave, 'Y1_LIMIT', control, 1);  if (err_flag) return
    call group_change_this (slave, 'Y2_LIMIT', control, 1);  if (err_flag) return

  case ('APERTURE') 
    call group_change_this (slave, 'X1_LIMIT', control, 1);  if (err_flag) return
    call group_change_this (slave, 'X2_LIMIT', control, 1);  if (err_flag) return
    call group_change_this (slave, 'Y1_LIMIT', control, 1);  if (err_flag) return
    call group_change_this (slave, 'Y2_LIMIT', control, 1);  if (err_flag) return

  !---------
  ! All else

  case default

    call group_change_this (slave, attrib_name, control, 1);  if (err_flag) return

  end select

enddo

!---------
! End stuff

if (moved) then
  call s_calc (lat)       ! recompute s distances
  call lat_geometry (lat)
endif

do i = 1, size(lord%control%var)
  lord%control%var(i)%old_value = lord%control%var(i)%value ! update old
enddo

lord%bookkeeping_state%control = ok$

!---------------------------------------------------------------------------------
contains

subroutine change_this_edge (this_slave, ctl)

type (ele_struct) this_slave
type (ele_struct), pointer :: this_slave2
type (branch_struct), pointer :: branch
type (control_struct) ctl
integer ix_min, ix_max, ix1, ix2

!

if (this_slave%lord_status == super_lord$) then
  this_slave2 => pointer_to_slave (this_slave, 1)
  ix_min = this_slave2%ix_ele
  this_slave2 => pointer_to_slave (this_slave, this_slave%n_slave)
  ix_max = this_slave2%ix_ele
  branch => lat%branch(this_slave2%ix_branch)
elseif (this_slave%ix_ele < this_slave%branch%n_ele_track) then
  ix_min = this_slave%ix_ele
  ix_max = this_slave%ix_ele
  branch => lat%branch(this_slave%ix_branch)
else
  call out_io (s_error$, r_name, &
                'A GROUP IS NOT ALLOWED TO CONTROL', &
                'A ' // control_name(slave%slave_status), &
                'YOU TRIED TO CONTROL: ' // slave%name)
  err_flag = .true.
  return
endif

! now that we have the ends we find the elements to either side whose length
! the group can adjust

if (attrib_name /= 'END_EDGE' .and. attrib_name /= 'L') then
  ix1 = ix_min - 1
  do
    if (attribute_name(branch%ele(ix1), l$) == 'L') exit  ! If has length attribute
    ix1 = ix1 - 1
    if (ix1 < 0) then
      call out_io (s_error$, r_name, &
                    'START_EDGE OF CONTROLED', &
                    'ELEMENT IS AT BEGINNING OF LAT AND CANNOT BE', &
                    'VARIED FOR GROUP: ' // lord%name)
      err_flag = .true.
      return
    endif
  enddo
endif

if (attrib_name /= 'START_EDGE' .and. attrib_name /= 'L') then
  ix2 = ix_max + 1 
  do
    if (attribute_name(branch%ele(ix2), l$) == 'L') exit  ! If has length attribute
    ix2 = ix2 + 1
    if (ix2 > branch%n_ele_track) then
      call out_io (s_error$, r_name, &
                    'END_EDGE OF CONTROLED', &
                    'ELEMENT IS AT END OF LAT AND CANNOT BE', &
                    'VARIED FOR GROUP: ' // lord%name)
      err_flag = .true.
      return
    endif
  enddo
endif

! put in changes

select case (attrib_name)

case ('L')
  call group_change_this (branch%ele(ix_max), 'L', ctl, 1);  if (err_flag) return

case ('START_EDGE')
  call group_change_this (branch%ele(ix_min), 'L', ctl, -1);  if (err_flag) return
  call group_change_this (branch%ele(ix1), 'L', ctl, 1);  if (err_flag) return

case ('END_EDGE')
  call group_change_this (branch%ele(ix_max), 'L', ctl, 1);  if (err_flag) return
  call group_change_this (branch%ele(ix2), 'L', ctl, -1);  if (err_flag) return

case ('ACCORDION_EDGE')
  call group_change_this (branch%ele(ix_min), 'L', ctl, 1);  if (err_flag) return
  call group_change_this (branch%ele(ix1), 'L', ctl, -1);  if (err_flag) return

  call group_change_this (branch%ele(ix_max), 'L', ctl, 1);  if (err_flag) return
  call group_change_this (branch%ele(ix2), 'L', ctl, -1);  if (err_flag) return

case ('S_POSITION')
  call group_change_this (branch%ele(ix1), 'L', ctl, 1);  if (err_flag) return
  call group_change_this (branch%ele(ix2), 'L', ctl, -1);  if (err_flag) return
case ('LORD_PAD1')
  call group_change_this (branch%ele(ix1), 'L', ctl, 1, this_slave, 'LORD_PAD1');  if (err_flag) return

case ('LORD_PAD2')
  call group_change_this (branch%ele(ix2), 'L', ctl, 1, this_slave, 'LORD_PAD1');  if (err_flag) return

end select

end subroutine change_this_edge

!---------------------------------------------------------------------------------
! contains 
!+
! Note: It is assumed that edge, super_lord, and varied_length optional args are all present if
! any one of them is present.
!-

recursive subroutine group_change_this (ele, attrib_name, ctl, dir, this_lord, this_pad)

type (ele_struct) ele
type (ele_struct), optional :: this_lord
type (ele_struct), pointer :: my_lord
type (control_struct) ctl
type (all_pointer_struct) :: a_ptr

integer dir, il, ix_slave

real(rp) coef, val_old, delta

logical ok

character(*) attrib_name
character(*), optional :: this_pad
character(100) err_str

!

call pointer_to_attribute (ele, attrib_name, .false., a_ptr, err_flag, do_unlink = .true.)
if (err_flag) then
  if (global_com%exit_on_error) call err_exit
  return
endif

! Evaluate value and old value.

if (allocated(ctl%stack)) then
  ctl%value = expression_stack_value (ctl%stack, err_flag, err_str, lord%control%var, .false.)
  val_old   = expression_stack_value (ctl%stack, err_flag, err_str, lord%control%var, .true.)
  if (err_flag) then
    call out_io (s_error$, r_name, err_str, 'FOR SLAVE: ' // slave%name, 'OF LORD: ' // lord%name)
    return
  endif

else
  val_old = knot_interpolate(lord%control%x_knot, ctl%y_knot, lord%control%var(1)%old_value, &
                                                                  nint(lord%value(interpolation$)), err_flag)
  ctl%value = knot_interpolate(lord%control%x_knot, ctl%y_knot, lord%control%var(1)%value, &
                                                                  nint(lord%value(interpolation$)), err_flag)
  if (err_flag) then
    call out_io (s_error$, r_name, 'EVALUATION PROBLEM FOR GROUP ELEMENT: ' // lord%name, &
                    'WHILE CALCULATING VALUE FOR: ' // trim(ele%name) // '[' // trim(attrib_name) // ']')
    return
  endif
endif

!

delta = ctl%value - val_old
a_ptr%r = a_ptr%r + delta * dir

call set_flags_for_changed_attribute (ele, a_ptr%r)
! super_slave length can be varied by a group so don't check this.
if ((ele%slave_status /= super_slave$ .and. ele%slave_status /= multipass_slave$) .or. attrib_name /= 'L') then
  err_flag = .not. attribute_free (ele, attrib_name, .true., .false., .true.)
  if (err_flag) then
    call out_io (s_blank$, r_name, 'GROUP_LORD TRYING TO CONTROL THIS ATTRIBUTE IS:' // lord%name)
    return
  endif
endif

! Pad check

if (ele%lord_status == super_lord$ .and. a_ptr%r < 0) then
  if (attrib_name == 'LORD_PAD1') then
    call out_io (s_error$, r_name, 'GROUP ELEMENT: ' // lord%name, &
                                   'CONTROLS SUPER_LORD: ' // ele%name, &
                                   'AND LORD_PAD1 IS NOW NEGATIVE: \f8.3\ ', r_array = [a_ptr%r])
    err_flag = .true.
    return
  elseif (attrib_name == 'LORD_PAD2') then
    call out_io (s_error$, r_name, 'GROUP ELEMENT: ' // lord%name, &
                                   'CONTROLS SUPER_LORD: ' // ele%name, &
                                   'AND LORD_PAD2 IS NOW NEGATIVE: \f8.3\ ', r_array = [a_ptr%r])
    err_flag = .true.
    return
  endif
endif

! ele is a super_slave...

if (ele%slave_status == super_slave$) then
  if (attrib_name /= 'L') then
    call out_io (s_error$, r_name, &
                  'CONFUSED GROUP IS TRYING TO VARY SUPER_SLAVE ATTRIBUTE: ' // attrib_name)
    if (global_com%exit_on_error) call err_exit
    err_flag = .true.
    return
  endif

  do il = 1, ele%n_lord
    my_lord => pointer_to_lord(ele, il)
    if (my_lord%lord_status /= super_lord$) cycle
    call group_change_this (my_lord, attrib_name, ctl, dir);  if (err_flag) return
  enddo

  if (present(this_lord)) then
    call group_change_this (my_lord, attrib_name, ctl, -dir);  if (err_flag) return  ! Take out length change.
    call group_change_this (my_lord, this_pad, ctl, dir);  if (err_flag) return    ! And change pad length instead.
  endif

endif

! ele is a multipass_slave...
! In the loop over all multipass_slaves, only modify the multipass_lord once

if (ele%slave_status == multipass_slave$) then
  my_lord => pointer_to_lord(ele, 1, ix_slave_back = ix_slave)
  if (ix_slave == 1) call group_change_this (my_lord, attrib_name, ctl, 1);  if (err_flag) return
endif

end subroutine group_change_this

end subroutine makeup_group_lord

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_multipass_slave (lat, slave, err_flag)
!
! Subroutine to calcualte the attributes of multipass slave elements.
! This routine is not meant for guse.
!-

subroutine makeup_multipass_slave (lat, slave, err_flag)

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord
type (branch_struct), pointer :: branch
type (floor_position_struct), pointer :: f0, f1
type (coord_struct) start, end

real(rp) s, slave_val(num_ele_attrib$), arg
real(rp) d, e, r_lord, r_slave, cos_lord, cos_e, sin_lord, sin_lorde
real(rp) ang_slave, ang_lord, ang_slave_old, d1, d2
real(rp) cos_e2, d_theta, ang_dlord, cos_lorde1, cos_dlord
real(rp) w0_mat(3,3), w1_mat(3,3), w1_inv_mat(3,3), offset(3), dw_mat(3,3)
real(rp) theta, phi, psi, w0_inv_mat(3,3)

integer i, j, ix_slave, ic, ix_s0, n_pass
character(40) :: r_name = 'makeup_multipass_slave'
logical err_flag

!

err_flag = .false.
branch => lat%branch(slave%ix_branch)
call set_ele_status_stale (slave, attribute_group$)
slave%bookkeeping_state%control = ok$

ix_slave = slave%ix_ele
j =  lat%ic(slave%ic1_lord)
lord => lat%ele(lat%control(j)%lord%ix_ele)
n_pass = j - lord%ix1_slave + 1  ! pass number for slave

slave_val = slave%value  ! save

slave%value = lord%value
if (lord%key == lcavity$ .or. lord%key == rfcavity$) then
  slave%value(phi0_multipass$) = slave_val(phi0_multipass$)
  slave%value(l_active$)       = lord%value(l_active$)
endif

! A slave's field_master = T irregardless of the lord's setting.
! This is to make attribute_bookkeeper compute the correct normalized field strength.

slave%field_calc             = refer_to_lords$
slave%value(E_tot_start$)    = slave_val(E_tot_start$)
slave%value(p0c_start$)      = slave_val(p0c_start$)
slave%value(e_tot$)          = slave_val(e_tot$)
slave%value(p0c$)            = slave_val(p0c$)
slave%value(delta_ref_time$) = slave_val(delta_ref_time$)
slave%value(ref_time_start$) = slave_val(ref_time_start$)

if (associated(slave%a_pole)) deallocate(slave%a_pole, slave%b_pole)
if (associated(slave%a_pole_elec)) deallocate(slave%a_pole_elec, slave%b_pole_elec)
if (allocated(slave%multipole_cache)) deallocate(slave%multipole_cache)

! A match element with recalc = True: Restore initial Twiss parameters (which
! are calculated in twiss_propagate1).

if (lord%key == match$) then
  if (is_true(lord%value(recalc$))) then
    if (nint(lord%value(matrix$)) == match_twiss$) then
      slave%value(beta_a0$)    = slave_val(beta_a0$)
      slave%value(beta_b0$)    = slave_val(beta_b0$)
      slave%value(alpha_a0$)   = slave_val(alpha_a0$)
      slave%value(alpha_b0$)   = slave_val(alpha_b0$)
      slave%value(eta_x0$)     = slave_val(eta_x0$)
      slave%value(eta_y0$)     = slave_val(eta_y0$)
      slave%value(etap_x0$)    = slave_val(etap_x0$)
      slave%value(etap_y0$)    = slave_val(etap_y0$)
      slave%value(c11_mat0$:mode_flip1$) = slave_val(c11_mat0$:mode_flip1$)
    endif

    if (nint(lord%value(kick0$)) == match_orbit$) then
      slave%value(x0$)  = slave_val(x0$)
      slave%value(px0$) = slave_val(px0$)
      slave%value(y0$)  = slave_val(y0$)
      slave%value(py0$) = slave_val(py0$)
      slave%value(z0$)  = slave_val(z0$)
      slave%value(pz0$) = slave_val(pz0$)
    endif
  endif
endif

! Sbend field: The design bending strength is same for slave as lord. 
! So the error field must be adjusted so that total_field = design_field + err_field is the same. 
! Note: The lord's energy may not yet be set if bmad_parser is active. So only do calc if p0c is set.

if (slave%key == sbend$ .and. lord%value(p0c$) /= 0) then
  slave%value(b_field$) = lord%value(b_field$) * slave%value(p0c$) / lord%value(p0c$)
  slave%value(db_field$) = (lord%value(b_field$) + lord%value(db_field$)) - slave%value(b_field$)
  slave%value(dg$) = (lord%value(g$) + lord%value(dg$)) - slave%value(g$)
endif

! Multipoles. Note: p0c = 0 Can happen if not finished parsing lattice file.

if (associated (slave%a_pole) .and. slave%value(p0c$) /= 0) then
  if (lord%key == multipole$) then
    slave%a_pole = lord%a_pole * lord%value(p0c$) / slave%value(p0c$)

  elseif (lord%field_master) then  
    slave%a_pole = lord%a_pole
    slave%b_pole = lord%b_pole

  else
    slave%a_pole = lord%a_pole * lord%value(p0c$) / slave%value(p0c$)
    slave%b_pole = lord%b_pole * lord%value(p0c$) / slave%value(p0c$)
  endif
  slave%multipoles_on    = lord%multipoles_on
  slave%scale_multipoles = lord%scale_multipoles
endif

! Electric Multipoles

if (associated (slave%a_pole_elec)) then
  slave%a_pole_elec      = lord%a_pole_elec
  slave%b_pole_elec      = lord%b_pole_elec
  slave%multipoles_on    = lord%multipoles_on
endif

! Custom attributes

if (associated(slave%custom)) then
  slave%custom = lord%custom
endif

! RF wakes

call transfer_wake (lord%wake, slave%wake)

if (associated (slave%wake)) then
  slave%wake%lr%t_ref = lord%wake%lr%t_ref - slave%ref_time
endif

! Methods

if (attribute_index(slave, 'FIELD_MASTER') /= 0) slave%field_master = .true.
slave%taylor_map_includes_offsets = lord%taylor_map_includes_offsets
slave%symplectify                 = lord%symplectify
slave%is_on                       = lord%is_on

! The following is handled by set_flags_for_changed_attribute

!! slave%ptc_integration_type = lord%ptc_integration_type
!! slave%csr_method           = lord%csr_method
!! slave%space_charge_method  = lord%space_charge_method
!! slave%aperture_at          = lord%aperture_at
!! slave%aperture_type        = lord%aperture_type
!! slave%mat6_calc_method     = lord%mat6_calc_method
!! slave%tracking_method      = lord%tracking_method
!! slave%spin_tracking_method = lord%spin_tracking_method

! A multipass_slave is allowed to have num_steps and ds_step set different from the lord.

slave%value(ds_step$)   = slave_val(ds_step$)
slave%value(num_steps$) = slave_val(num_steps$)

end subroutine makeup_multipass_slave

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_super_slave (lat, slave, err_flag)
!
! Subroutine to calcualte the attributes of superposition slave elements.
! This routine is not meant for general use.
!-
       
subroutine makeup_super_slave (lat, slave, err_flag)

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord, slave0, lord1
type (ele_struct) :: sol_quad
type (branch_struct), pointer :: branch

integer i, j, ix, ix_lord, ix_order, ix_slave, n_major_lords, at

real(rp) tilt, k_x, k_y, x_kick, y_kick, ks, k1, coef
real(rp) x_o, y_o, x_p, y_p, s_slave, s_del, k2, k3, c, s
real(rp) sin_n, cos_n, a(0:n_pole_maxx), b(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), t(0:n_pole_maxx), value(num_ele_attrib$), old_value(num_ele_attrib$)
real(rp) sum_1, sum_2, sum_3, sum_4, ks_sum, ks_xp_sum, ks_xo_sum
real(rp) ks_yp_sum, ks_yo_sum, l_slave, r_off(4), leng, offset
real(rp) t_1(4), t_2(4), T_end(4,4), mat4(4,4), mat4_inv(4,4), beta(4)
real(rp) T_tot(4,4), x_o_sol, x_p_sol, y_o_sol, y_p_sol

logical is_first, is_last, err_flag, fixed_step

character(20) :: r_name = 'makeup_super_slave'

! Super_slave:

branch => lat%branch(slave%ix_branch)
ix_slave = slave%ix_ele
err_flag = .false.

if (slave%slave_status /= super_slave$) then
  call out_io(s_abort$, r_name, "ELEMENT IS NOT A SUPER SLAVE: " // slave%name)
  if (global_com%exit_on_error) call err_exit
  err_flag = .true.
  return
endif

slave%field_calc    = refer_to_lords$
slave%aperture_at   = lord_defined$
slave%aperture_type = lord_defined$

if (associated(slave%a_pole)) deallocate(slave%a_pole, slave%b_pole)
if (associated(slave%a_pole_elec)) deallocate(slave%a_pole_elec, slave%b_pole_elec)

! Bookkeeping for EM_Field slave is mostly independent of the lords.
! Exception: If only one lord then treat em_field slave same as other slaves.

if (slave%key == em_field$ .and. slave%n_lord > 1) then
  value = slave%value
  slave%value = 0
  slave%value(l$)                     = value(l$)
  slave%value(delta_ref_time$)        = value(delta_ref_time$)
  slave%value(p0c_start$)             = value(p0c_start$)
  slave%value(e_tot_start$)           = value(e_tot_start$)
  slave%value(p0c$)                   = value(p0c$)
  slave%value(e_tot$)                 = value(e_tot$)
  slave%value(ref_time_start$)        = value(ref_time_start$)
  slave%value(check_sum$)             = value(check_sum$)
  slave%value(multipass_ref_energy$)  = first_pass$
  slave%value(fringe_type$)           = none$
  slave%value(fringe_at$)             = no_end$
  slave%value(spin_fringe_on$)        = false$
  slave%value(autoscale_phase$)       = false$
  slave%value(autoscale_amplitude$)   = false$
  slave%mat6_calc_method              = tracking$
  slave%spin_tracking_method          = tracking$
  slave%tracking_method               = runge_kutta$
  ! Use time_runge_kutta over runge_kutta since runge_kutta is not able to handle a particle starting from rest in an e_gun.
  ! Used fixed_step over non-fixed step since fixed_step is basically only used for testing.
  fixed_step = .false.
  do i = 1, slave%n_lord
    lord => pointer_to_lord(slave, i)
    if (slave%value(ds_step$) == 0 .or. lord%value(ds_step$) < slave%value(ds_step$)) slave%value(ds_step$) = lord%value(ds_step$)
    if (lord%tracking_method == fixed_step_runge_kutta$ .or. lord%tracking_method == fixed_step_time_runge_kutta$) fixed_step = .true.
    if (lord%tracking_method == time_runge_kutta$) slave%tracking_method = time_runge_kutta$
    if (lord%tracking_method == fixed_step_time_runge_kutta$) slave%tracking_method = fixed_step_time_runge_kutta$
  enddo
  if (fixed_step .and. slave%tracking_method == runge_kutta$) slave%tracking_method = fixed_step_runge_kutta$
  if (fixed_step .and. slave%tracking_method == time_runge_kutta$) slave%tracking_method = fixed_step_time_runge_kutta$

  slave%value(constant_ref_energy$) = true$
  do i = 1, slave%n_lord
    lord => pointer_to_lord(slave, i)
    if (lord%lord_status /= super_lord$) cycle
    lord%tracking_method = slave%tracking_method  ! Important for lcavity element to make sure 
                                                  !   tracking along entire element is self-consistent.
    lord%mat6_calc_method = tracking$
    lord%spin_tracking_method = tracking$
    select case (lord%key)
    case (lcavity$)
      slave%value(constant_ref_energy$) = false$
    case (em_field$)
      if (is_false(lord%value(constant_ref_energy$))) slave%value(constant_ref_energy$) = false$
    end select
  enddo

  return
endif

!-----------------------------------------------------------------------
! A "major" super_lord is something other than a pipe.
! If only one "major" super_lord for this super_slave: just transfer attributes except length

ix_lord = 1  ! Index of major lord or first lord 
n_major_lords = 0
do j = 1, slave%n_lord
  lord => pointer_to_lord (slave, j)
  select case (lord%key)
  case (pipe$, instrument$, monitor$, ecollimator$, rcollimator$)
  case default
    n_major_lords = n_major_lords + 1
    ix_lord = j
  end select
enddo

if (n_major_lords < 2) then
  lord => pointer_to_lord (slave, ix_lord, ix_slave_back = ix_order)
  old_value = slave%value

  is_first = (ix_order == 1)
  is_last  = (ix_order == lord%n_slave)

  ! If this is not the first slave: Transfer reference orbit from previous slave

  if (.not. is_first) then
    if (.not. all(slave%map_ref_orb_in%vec == branch%ele(ix_order-1)%map_ref_orb_out%vec)) then
      slave0 => pointer_to_slave(lord, ix_order-1)
      slave%map_ref_orb_in = slave0%map_ref_orb_out
      if (allocated(slave%multipole_cache)) then
        slave%multipole_cache%mag_valid = .false.
        slave%multipole_cache%elec_valid = .false.
      endif
      if (associated(slave%rad_map)) slave%rad_map%stale = .true. ! Forces recalc
    endif
  endif

  ! Find the offset from the longitudinal start of the lord to the start of the slave

  offset = 0 ! length of all slaves before this one
  do i = 1, ix_order-1
    slave0 => pointer_to_slave(lord, i)
    offset = offset + slave0%value(l$)
  enddo

  call makeup_super_slave1 (slave, lord, offset, branch%param, is_first, is_last, err_flag)

  ! A pipe may have a kick so add that in.
  ! Note: check_aperture_limits knows to look at the lords for apertures so we do not need to fiddle with apertures.

  do j = 1, slave%n_lord
    if (j == ix_lord) cycle  ! Do not double count
    lord => pointer_to_lord (slave, j, ix_slave_back = ix_order)
    slave%value(hkick$) = slave%value(hkick$) + lord%value(hkick$)
    slave%value(vkick$) = slave%value(vkick$) + lord%value(vkick$)
    slave%value(bl_hkick$) = slave%value(bl_hkick$) + lord%value(bl_hkick$)
    slave%value(bl_vkick$) = slave%value(bl_vkick$) + lord%value(bl_vkick$)
  enddo

  if (any(slave%value /= old_value)) call attribute_bookkeeper (slave, .true.)

  return
endif

!-----------------------------------------------------------------------
! Multiple super_lords for this super_slave: 
! combine the lord elements.

if (allocated(slave%multipole_cache)) deallocate(slave%multipole_cache)

k_x = 0
k_y = 0
x_kick = 0
y_kick = 0
sum_1 = 0
sum_2 = 0
sum_3 = 0
sum_4 = 0
ks_sum = 0
ks_xp_sum = 0
ks_xo_sum = 0
ks_yp_sum = 0
ks_yo_sum = 0

!

old_value = slave%value

value = 0
value(l$)                = slave%value(l$)
value(E_tot_start$)      = slave%value(E_tot_start$)
value(p0c_start$)        = slave%value(p0c_start$)
value(E_tot$)            = slave%value(E_tot$)
value(p0c$)              = slave%value(p0c$)
value(delta_ref_time$)   = slave%value(delta_ref_time$)
value(ref_time_start$)   = slave%value(ref_time_start$)
value(fringe_at$)        = no_end$
value(fringe_type$)      = none$
value(integrator_order$) = 0

slave%value(x1_limit$:y2_limit$) = 0  ! check_aperture_limits knows to look at the lords for apertures.
slave%is_on                      = .false.

s_slave = slave%s - value(l$)/2  ! center of slave

! A "major" element is something other than a pipe, monitor, etc.
! n_major_lords counts how many major lords there are.

n_major_lords = 0

! sum over all lords...

do j = 1, slave%n_lord

  lord => pointer_to_lord(slave, j, ix_slave_back = ix_order)
  if (lord%lord_status /= super_lord$) cycle

  is_first = (ix_order == 1)
  is_last  = (ix_order == lord%n_slave)

  ! Do some error checking.

  if (associated(lord%wake)) then
    call out_io (s_abort$, r_name, &
            'SUPERPOSITION OF ELEMENTS WITH WAKES NOT YET IMPLEMENTED!', &
            'SUPER_LORD: ' // lord%name)
    err_flag = .true.
    return
  endif

  ! Physically, the lord length cannot be less than the slave length.
  ! In case we are dealing with a non-physical situation, arbitrarily set coef = 1.

  if (abs(slave%value(l$)) >= abs(lord%value(l$))) then
    coef = 1
  else
    coef = slave%value(l$) / lord%value(l$) 
  endif

  ! If this is not the first slave: Transfer reference orbit from previous slave

  if (.not. is_first) then
    if (.not. all(slave%map_ref_orb_in%vec == branch%ele(ix_slave-1)%map_ref_orb_out%vec)) then
      slave%map_ref_orb_in = branch%ele(ix_slave-1)%map_ref_orb_out
      if (allocated(slave%multipole_cache)) then
        slave%multipole_cache%mag_valid = .false.
        slave%multipole_cache%elec_valid = .false.
      endif
      if (associated(slave%rad_map)) slave%rad_map%stale = .true. ! Forces recalc
    endif
  endif

  ! Choose the largest integrator order

  value(integrator_order$) = max(value(integrator_order$), lord%value(integrator_order$))

  ! Choose the smallest ds_step of all the lords.

  if (value(ds_step$) == 0 .or. lord%value(ds_step$) < value(ds_step$)) value(ds_step$) = lord%value(ds_step$)

  ! Methods...
  ! A "major" element is something other than a pipe.
  ! n_major_lords counts how many major lords there are.

  if (n_major_lords == 0) then
    slave%mat6_calc_method            = lord%mat6_calc_method
    slave%tracking_method             = lord%tracking_method
    slave%taylor_map_includes_offsets = lord%taylor_map_includes_offsets
    slave%csr_method                  = lord%csr_method
    slave%space_charge_method         = lord%space_charge_method
  endif

  if (has_attribute (lord, 'FRINGE_TYPE')) then
    if (is_first .and. at_this_ele_end(entrance_end$, nint(lord%value(fringe_at$)))) then
      call set_fringe_on_off(value(fringe_at$), entrance_end$, on$)
      value(fringe_type$) = lord%value(fringe_type$)
    endif

    if (is_last .and. at_this_ele_end(exit_end$, nint(lord%value(fringe_at$)))) then
      call set_fringe_on_off(value(fringe_at$), exit_end$, on$)
      value(fringe_type$) = lord%value(fringe_type$)
    endif
  endif

  select case (lord%key)
  case (pipe$, instrument$, monitor$, ecollimator$, rcollimator$)
  case default
    if (n_major_lords > 0) then
 
     if (slave%mat6_calc_method /= lord%mat6_calc_method) then
        call out_io(s_error$, r_name, &
              'MAT6_CALC_METHOD DOES NOT AGREE FOR DIFFERENT SUPERPOSITION LORDS FOR SLAVE: ' // slave%name, &
              'Conflicting methods are: ' // trim(mat6_calc_method_name(lord%mat6_calc_method)) // ',  ' // & 
              mat6_calc_method_name(slave%mat6_calc_method))
      endif

      if (slave%tracking_method /= lord%tracking_method) then
        if (slave%tracking_method == fixed_step_runge_kutta$ .and. lord%tracking_method == runge_kutta$) then
          ! Do nothing
        elseif (slave%tracking_method == runge_kutta$ .and. lord%tracking_method == fixed_step_runge_kutta$) then
          slave%tracking_method = fixed_step_runge_kutta$
        elseif (slave%tracking_method == fixed_step_time_runge_kutta$ .and. lord%tracking_method == time_runge_kutta$) then
          ! Do nothing
        elseif (slave%tracking_method == time_runge_kutta$ .and. lord%tracking_method == fixed_step_time_runge_kutta$) then
          slave%tracking_method = fixed_step_time_runge_kutta$
        else
          call out_io(s_error$, r_name, &
             'TRACKING_METHOD DOES NOT AGREE FOR DIFFERENT SUPERPOSITION LORDS FOR SLAVE: ' // slave%name, &
             'Conflicting methods are: ' // trim(tracking_method_name(lord%tracking_method)) // ',  ' // & 
             tracking_method_name(slave%tracking_method))
        endif
      endif

      if (slave%csr_method == off$) slave%csr_method = lord%csr_method
      if (slave%space_charge_method == off$) slave%space_charge_method = lord%space_charge_method

      if (slave%taylor_map_includes_offsets .neqv. lord%taylor_map_includes_offsets) then
        call out_io(s_error$, r_name, &
            'TAYLOR_MAP_INCLUDES_OFFSETS DOES NOT AGREE FOR DIFFERENT SUPERPOSITION LORDS FOR SLAVE: ' // slave%name)
      endif

      if ((is_first .or. is_last) .and. has_attribute (lord, 'FRINGE_TYPE')) then
       if (value(fringe_type$) /= lord%value(fringe_type$)) then
         call out_io(s_error$, r_name, &
            'FRINGE_TYPE DOES NOT AGREE FOR DIFFERENT SUPERPOSITION LORDS FOR SLAVE: ' // slave%name)
       endif

     endif
    endif

    n_major_lords = n_major_lords + 1
  end select

  ! descriptive strings.

  if (associated(lord%descrip)) then
    if (.not. associated(slave%descrip)) allocate (slave%descrip)
    slave%descrip = lord%descrip
  endif

  if (lord%type /= '') slave%type = lord%type
  if (lord%alias /= '') slave%alias = lord%alias

  !----------------------------------------------------
  ! kicks, etc.

  if (.not. lord%is_on) cycle
  slave%is_on = .true.  ! Slave is on if at least one lord is on

  if (slave%key == em_field$) cycle  ! Field info is stored in the lord elements.

  tilt = lord%value(tilt_tot$)

  if (lord%key == hkicker$) then
    x_kick = x_kick + lord%value(kick$) * cos(tilt) * coef
    y_kick = y_kick + lord%value(kick$) * sin(tilt) * coef
  elseif (lord%key == vkicker$) then
    x_kick = x_kick - lord%value(kick$) * sin(tilt) * coef
    y_kick = y_kick + lord%value(kick$) * cos(tilt) * coef
  elseif (lord%key == ac_kicker$ .or. lord%key == kicker$) then
    c = cos(tilt) * coef
    s = sin(tilt) * coef
    x_kick = x_kick + c * lord%value(hkick$) - s * lord%value(vkick$)
    y_kick = y_kick + s * lord%value(hkick$) + c * lord%value(vkick$)
  else
    x_kick = x_kick + lord%value(hkick$) * coef
    y_kick = y_kick + lord%value(vkick$) * coef
  endif

  !------

  select case (slave%key)

  ! sextupole

  case (sextupole$) 

    cos_n = lord%value(k2$) * cos(3 * tilt)
    sin_n = lord%value(k2$) * sin(3 * tilt)            
    
    k_x = k_x + cos_n
    k_y = k_y + sin_n

  ! octupole

  case (octupole$)

    cos_n = lord%value(k3$) * cos(4 * tilt)
    sin_n = lord%value(k3$) * sin(4 * tilt)        
    
    k_x = k_x + cos_n
    k_y = k_y + sin_n

  ! solenoid/quadrupole combo.

  case (solenoid$, sol_quad$, quadrupole$)

    x_p = lord%value(x_pitch_tot$);  x_o = lord%value(x_offset_tot$)
    y_p = lord%value(y_pitch_tot$);  y_o = lord%value(y_offset_tot$)

    s_del = s_slave - (lord%s_start + lord%value(z_offset_tot$) + lord%value(l$)/2)
    s_del = modulo2 (s_del, branch%param%total_length/2)

    ks = lord%value(ks$)

    ks_sum = ks_sum + ks

    ks_xp_sum = ks_xp_sum + ks * x_p
    ks_yp_sum = ks_yp_sum + ks * y_p

    ks_xo_sum = ks_xo_sum + ks * (x_o + x_p * s_del)
    ks_yo_sum = ks_yo_sum + ks * (y_o + y_p * s_del)

    cos_n = lord%value(k1$) * cos(2 * tilt)
    sin_n = lord%value(k1$) * sin(2 * tilt)

    k_x = k_x + cos_n
    k_y = k_y + sin_n

    sum_1 = sum_1 + cos_n * x_p + sin_n * y_p
    sum_2 = sum_2 + sin_n * x_p - cos_n * y_p

    sum_3 = sum_3 + cos_n * (x_o + x_p * s_del) + sin_n * (y_o + y_p * s_del)
    sum_4 = sum_4 + sin_n * (x_o + x_p * s_del) - cos_n * (y_o + y_p * s_del)

  ! Everything else

  case default
    ! Everything else has no special needs. 

  end select

enddo

if (slave%tracking_method == bmad_standard$ .and. slave%key == em_field$) slave%tracking_method = runge_kutta$
if (slave%mat6_calc_method == bmad_standard$ .and. slave%key == em_field$) slave%mat6_calc_method = tracking$

!-------------------------------------------------------------------------------
! stuff sums into slave element

if (slave%key == em_field$) then
  slave%value = value
  goto 8000  ! Field info is stored in the lord elements.
endif

slave%value = value

! Kick values

if (x_kick == 0 .and. y_kick == 0) then
  if (slave%key == hkicker$ .or. slave%key == vkicker$) then
    slave%value(kick$) = 0
  else
    slave%value(hkick$) = 0
    slave%value(vkick$) = 0
  endif
elseif (slave%key == hkicker$) then
  slave%value(kick$) = sqrt(x_kick**2 + y_kick**2)
  slave%value(tilt$) = atan2(y_kick, x_kick)
elseif (slave%key == vkicker$) then
  slave%value(kick$) = sqrt(x_kick**2 + y_kick**2)
  slave%value(tilt$) = atan2(-x_kick, y_kick)
elseif (slave%key == ac_kicker$ .or. slave%key == kicker$) then
  slave%value(tilt$) = 0
  slave%value(hkick$) = x_kick
  slave%value(vkick$) = y_kick
else
  slave%value(hkick$) = x_kick
  slave%value(vkick$) = y_kick
endif

!-----------------------------

select case (slave%key)

case (sextupole$) 

  if (k_x == 0 .and. k_y == 0) goto 8000

  k2 = sqrt(k_x**2 + k_y**2)
  tilt = atan2(k_y, k_x) / 3

  if (tilt > pi/6) then
    k2 = -k2
    tilt = tilt - pi/3
  elseif (tilt < -pi/6) then
    k2 = -k2
    tilt = tilt + pi/3
  endif

  slave%value(k2$) = k2
  slave%value(tilt$) = tilt

! octupole

case (octupole$)

  if (k_x == 0 .and. k_y == 0 .and. ks == 0) goto 8000

  k3 = sqrt(k_x**2 + k_y**2)
  tilt = atan2(k_y, k_x) / 4

  if (tilt > pi/8) then
    k3 = -k3
    tilt = tilt - pi/4
  elseif (tilt < -pi/8) then
    k3 = -k3
    tilt = tilt + pi/4
  endif

  slave%value(k3$) = k3
  slave%value(tilt$) = tilt

! sol_quad, etc.

case (solenoid$, sol_quad$, quadrupole$)

  ks = ks_sum
  slave%value(ks$) = ks

  if (k_x == 0 .and. k_y == 0 .and. ks == 0) goto 8000

  if (ks /= 0) then
    x_o_sol = ks_xo_sum / ks
    x_p_sol = ks_xp_sum / ks
    y_o_sol = ks_yo_sum / ks
    y_p_sol = ks_yp_sum / ks
  endif

  if (k_x == 0 .and. k_y == 0) then  ! pure solenoid
    slave%value(k1$) = 0
    slave%value(tilt$) = 0
    slave%value(x_offset$) = x_o_sol
    slave%value(y_offset$) = y_o_sol
    slave%value(x_pitch$)  = x_p_sol
    slave%value(y_pitch$)  = y_p_sol
  endif   

  ! here if have quadrupole component

  if (k_x /= 0 .or. k_y /= 0) then
    k1 = sqrt(k_x**2 + k_y**2)
    tilt = atan2(k_y, k_x) / 2

    if (tilt > pi/4) then
      k1 = -k1
      tilt = tilt - pi/2
    elseif (tilt < -pi/4) then
      k1 = -k1
      tilt = tilt + pi/2
    endif

    slave%value(k1$) = k1
    slave%value(tilt$) = tilt

    cos_n = k_x / (k_x**2 + k_y**2)
    sin_n = k_y / (k_x**2 + k_y**2)

    slave%value(x_pitch$)  = cos_n * sum_1 + sin_n * sum_2
    slave%value(y_pitch$)  = sin_n * sum_1 - cos_n * sum_2
    slave%value(x_offset$) = cos_n * sum_3 + sin_n * sum_4
    slave%value(y_offset$) = sin_n * sum_3 - cos_n * sum_4
  endif

  ! if ks /= 0 then we have to recalculate the offsets and pitches.

  if (ks /= 0 .and. (k_x /= 0 .or. k_y /= 0)) then

    x_p = slave%value(x_pitch$) - x_p_sol; x_o = slave%value(x_offset$) - x_o_sol
    y_p = slave%value(y_pitch$) - y_p_sol; y_o = slave%value(y_offset$) - y_o_sol

    if (x_p == 0 .and. x_o == 0 .and. y_p == 0 .and. y_o == 0) goto 8000

    t_2 = [x_o, x_p, y_o, y_p]
    call tilt_coords (tilt, t_2)

    l_slave = slave%value(l$)

    t_1 = [t_2(2), 0.0_rp, t_2(4), 0.0_rp]
    t_2(1) = t_2(1) + ks * t_2(4) / k1 
    t_2(3) = t_2(3) + ks * t_2(2) / k1
             
    call mat_make_unit (T_end)
    T_end(4,1) =  ks / 2
    T_end(2,3) = -ks / 2

    call transfer_ele (slave, sol_quad)
    sol_quad%key = sol_quad$
    sol_quad%value(ks$) = ks
    sol_quad%value(k1$) = k1
    sol_quad%value(l$)  = l_slave
    call set_flags_for_changed_attribute(sol_quad, sol_quad%value(ks$))
    call set_flags_for_changed_attribute(sol_quad, sol_quad%value(k1$))
    call set_flags_for_changed_attribute(sol_quad, sol_quad%value(l$))
    call make_mat6 (sol_quad, branch%param)
    T_tot = sol_quad%mat6(1:4,1:4)

    r_off = matmul (T_end, l_slave * t_1 / 2 - t_2) 
    r_off = matmul (T_tot, r_off) + matmul (T_end, l_slave * t_1 / 2 + t_2)

    call mat_make_unit (mat4)
    mat4(:,2) = mat4(:,2) + l_slave * T_tot(:,1) / 2
    mat4(:,4) = mat4(:,4) + l_slave * T_tot(:,3) / 2
    mat4(1,2) = mat4(1,2) + l_slave / 2
    mat4(3,4) = mat4(3,4) + l_slave / 2
    mat4 = mat4 - T_tot

    call mat_inverse (mat4, mat4_inv)
    beta = matmul (mat4_inv, r_off)

    call tilt_coords (-tilt, beta)

    slave%value(x_offset$) = beta(1) + x_o_sol
    slave%value(x_pitch$)  = beta(2) + x_p_sol
    slave%value(y_offset$) = beta(3) + y_o_sol
    slave%value(y_pitch$)  = beta(4) + y_p_sol
  endif

end select

! If the slave has %field_master = T then we need to convert k1, etc values to field quantities.

8000 continue

! Coupler and aperture calc.

if (slave%key == lcavity$ .or. slave%key == rfcavity$) call compute_slave_coupler (slave)
if (all(slave%value == old_value)) return

call set_ele_status_stale (slave, attribute_group$)

if (slave%field_master) then
  slave%field_master = .false.   ! So attribute_bookkeeper will do the right thing.
  call attribute_bookkeeper (slave, .true.)
  slave%field_master = .true.
else
  call attribute_bookkeeper (slave, .true.)
endif

end subroutine makeup_super_slave

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_super_slave1 (slave, lord, offset, param, include_upstream_end, include_downstream_end)
!
! Routine to construct a super_slave from a super_lord when the slave has only one lord.
! Note: Reference energy and times are not computed in this routine.
!
! Input:
!   slave                   -- Ele_struct: Slave element.
!   lord                    -- Ele_struct: Lord element.
!   offset                  -- Real(rp): offset of entrance end of slave from entrance end of the lord.
!   param                   -- Lat_param_struct: lattice paramters.
!   include_upstream_end    -- Logical: Slave contains the lord's entrance end?
!   include_downstream_end  -- Logical: Slave contains the lord's exit end?
!
! Output:
!   slave                   -- Ele_struct: Slave element with appropriate values set.
!   err_flag                -- Logical: Set true if there is an error. False otherwise.
!-

subroutine makeup_super_slave1 (slave, lord, offset, param, include_upstream_end, include_downstream_end, err_flag)

type (ele_struct), target :: slave, lord
type (ele_struct), pointer :: slave2, lord2
type (lat_param_struct) param
type (floor_position_struct) from_pos, to_pos

real(rp) offset, s_del, coef, lord_ang, slave_ang, angle
real(rp) value(num_ele_attrib$), cos_a, sin_a, dr(3), w_mat(3,3)
real(rp) off(3), rot(3), cos_t, sin_t, m_trans(3,3)
real(rp) xp, yp, roll, r_roll, tilt, dx, dy, len_slave, len_lord
real(rp) w_inv(3,3), dl, vel

integer i, ifr, ixs
logical include_upstream_end, include_downstream_end, err_flag, include_entrance, include_exit, err2_flag
logical has_fringe
character(24) :: r_name = 'makeup_super_slave1'

! Physically, the lord length cannot be less than the slave length.
! In case we are dealing with a non-physical situation, arbitrarily set coef = 1.

err_flag = .true.
len_slave = slave%value(l$)
len_lord = lord%value(l$)

if (abs(len_slave) >= abs(len_lord)) then
  coef = 1
else
  coef = len_slave / len_lord 
endif

if (lord%orientation == 1) then
  include_entrance = include_upstream_end
  include_exit     = include_downstream_end
else
  include_entrance = include_downstream_end
  include_exit     = include_upstream_end
endif

! Reference energy and time computed in ele_compute_ref_energy_and_time.

value = lord%value

value(l$)              = len_slave                ! do not change slave length, etc.
value(delta_ref_time$) = slave%value(delta_ref_time$)
value(ref_time_start$) = slave%value(ref_time_start$)
value(E_tot_start$)    = slave%value(E_tot_start$)
value(p0c_start$)      = slave%value(p0c_start$)
value(E_tot$)          = slave%value(E_tot$)
value(p0c$)            = slave%value(p0c$)
value(num_steps$)      = slave%value(num_steps$)
if (attribute_name(slave%key, split_id$, .true.) == 'SPLIT_ID') value(split_id$) = slave%value(split_id$)

! Taylor element has a zero length map. Rule: The map gets applied at the entrance end.
! There is no reason why the entrance end was chosen over the exit end.
! Note: multipole, and ab_multipole have zero length maps but in this case multipole_ele_to_kt
! will scale the multipole strength proportional to the slave length.

has_fringe = .true.

if (lord%key == taylor$) then
  if (.not. include_entrance) slave%is_on = .false.
  has_fringe = .false.
elseif (lord%key == multipole$ .or. lord%key == ab_multipole$) then
  has_fringe = .false.
endif

! Ref energy shift for e_gun only happens at start of element.

if (lord%key == e_gun$) then
  if (.not. include_upstream_end) then
    value(e_tot_ref_init$) = lord%value(e_tot_start$)
    value(p0c_ref_init$)   = lord%value(p0c_start$)
  endif
endif

! fringe fields 

if (has_fringe) then
  ifr = nint(value(fringe_at$))
  if (include_entrance .and. include_exit) then
    ! Inherit from lord
  elseif (include_entrance) then
    if (ifr == entrance_end$ .or. ifr == both_ends$) then
      value(fringe_at$) = entrance_end$
    else
      value(fringe_at$) = no_end$
    endif
  elseif (include_exit) then
    if (ifr == exit_end$ .or. ifr == both_ends$) then
      value(fringe_at$) = exit_end$
    else
      value(fringe_at$) = no_end$
    endif
  else
    value(fringe_at$) = no_end$
  endif
endif

!

if (lord%key == hkicker$ .or. lord%key == vkicker$) then
  value(kick$)    = lord%value(kick$) * coef
  value(bl_kick$) = lord%value(bl_kick$) * coef
elseif (has_hkick_attributes(lord%key)) then
  value(hkick$)    = lord%value(hkick$) * coef
  value(vkick$)    = lord%value(vkick$) * coef
  value(bl_hkick$) = lord%value(bl_hkick$) * coef
  value(bl_vkick$) = lord%value(bl_vkick$) * coef
endif

select case (slave%key)
case (crab_cavity$)
  value(voltage$)     = lord%value(voltage$) * coef
  value(voltage_tot$) = lord%value(voltage_tot$) * coef
case (lcavity$, rfcavity$, e_gun$)
  value(voltage$)     = lord%value(voltage$) * coef
  value(voltage_tot$) = lord%value(voltage_tot$) * coef
  value(voltage_err$) = lord%value(voltage_err$) * coef
end select

!

if (allocated(lord%multipole_cache)) then
  slave%multipole_cache = lord%multipole_cache
  if (allocated(slave%multipole_cache%a_pole_mag)) then
    slave%multipole_cache%a_pole_mag = slave%multipole_cache%a_pole_mag * coef
    slave%multipole_cache%b_pole_mag = slave%multipole_cache%b_pole_mag * coef
  endif

  if (allocated(slave%multipole_cache%a_kick_mag)) then
    slave%multipole_cache%a_kick_mag = slave%multipole_cache%a_kick_mag * coef
    slave%multipole_cache%b_kick_mag = slave%multipole_cache%b_kick_mag * coef
  endif

  if (allocated(slave%multipole_cache%a_pole_elec)) then
    slave%multipole_cache%a_pole_elec = slave%multipole_cache%a_pole_elec
    slave%multipole_cache%b_pole_elec = slave%multipole_cache%b_pole_elec
  endif

  if (allocated(slave%multipole_cache%a_kick_elec)) then
    slave%multipole_cache%a_kick_elec = slave%multipole_cache%a_kick_elec
    slave%multipole_cache%b_kick_elec = slave%multipole_cache%b_kick_elec
  endif

else
  if (allocated(slave%multipole_cache)) deallocate(slave%multipole_cache)
endif

! s_del is the distance between lord and slave centers

if (has_orientation_attributes(slave)) then

  if ((slave%key == sbend$ .or. slave%key == rf_bend$) .and. value(g$) /= 0) then

    roll = value(roll_tot$);     tilt = value(ref_tilt_tot$)
    off = [value(x_offset_tot$), value(y_offset_tot$), value(z_offset_tot$)]
    xp  = value(x_pitch_tot$);   yp = value(y_pitch_tot$)

    value(ref_tilt$) = tilt

    if (any(off /= 0) .or. xp /= 0 .or. yp /= 0 .or. roll /= 0) then
      from_pos = floor_position_struct([0,0,0], mat3_unit$, 0, 0, 0)
      from_pos%r(3) = offset + len_slave/2
      to_pos = coords_body_to_rel_exit (from_pos, lord)
      to_pos = bend_shift (to_pos, lord%value(g$), offset + len_slave/2 - len_lord, ref_tilt = tilt)

      w_mat = to_pos%w
      if (tilt /= 0) call rotate_mat (w_mat, z_axis$, -tilt, right_multiply = .true.)
      call floor_w_mat_to_angles (w_mat, value(x_pitch$), value(y_pitch$), value(roll$))

      off = to_pos%r

      if (roll /= 0) then
        rot = [-value(rho$) * cos_one(value(g$)*len_slave/2), 0.0_rp, 0.0_rp]
        call rotate_vec(rot, z_axis$, tilt)
        off = off + rot
        call rotate_vec(rot, z_axis$, roll)
        off = off - rot
      endif

      value(x_offset$) = off(1)
      value(y_offset$) = off(2)
      value(z_offset$) = off(3)
    endif

  ! Not an sbend

  else
    s_del = offset + len_slave/2 - len_lord/2
    value(tilt$)     = value(tilt_tot$)
    value(x_pitch$)  = value(x_pitch_tot$)
    value(y_pitch$)  = value(y_pitch_tot$)
    
    value(x_offset$) = value(x_offset_tot$) + s_del * sin(value(x_pitch_tot$)) * cos(value(y_pitch_tot$))
    value(y_offset$) = value(y_offset_tot$) + s_del * sin(value(y_pitch_tot$))
    value(z_offset$) = value(z_offset_tot$) + s_del * (cos(value(x_pitch_tot$)) * cos(value(y_pitch_tot$)) - 1)
  endif
endif

! Patch
! The rotation part of the patch is applied at the entrance end of the patch.
! Excluding the rotation, a patch is just a drift.

if (lord%key == patch$) then
  if ((include_upstream_end .and. lord%orientation == 1) .or. (include_downstream_end .and. lord%orientation == -1)) then
    call floor_angles_to_w_mat (lord%value(x_pitch$), lord%value(y_pitch$), lord%value(tilt$), w_mat_inv = w_inv)
    dl = len_lord - len_slave
    value(x_offset$)     = lord%value(x_offset$) - dl * w_inv(3,1)
    value(y_offset$)     = lord%value(y_offset$) - dl * w_inv(3,2)
    value(z_offset$)     = lord%value(z_offset$) - dl * w_inv(3,3)
    value(t_offset$)     = lord%value(t_offset$)
    value(e_tot_offset$) = lord%value(e_tot_offset$)
    value(e_tot_set$)    = lord%value(e_tot_set$)
    value(p0c_set$)      = lord%value(p0c_set$)
  else
    value(x_pitch$)      = 0
    value(y_pitch$)      = 0
    value(tilt$)         = 0
    value(x_offset$)     = 0
    value(y_offset$)     = 0
    value(z_offset$)     = len_slave ! L is set by create_element_slice
    value(t_offset$)     = 0
    value(e_tot_offset$) = 0
    value(e_tot_set$)    = 0
    value(p0c_set$)      = 0
  endif

  value(x_pitch_tot$)     = value(x_pitch$)
  value(y_pitch_tot$)     = value(y_pitch$)
  value(tilt_tot$)        = value(tilt$)
  value(x_offset_tot$)    = value(x_offset$)
  value(y_offset_tot$)    = value(y_offset$)
  value(z_offset_tot$)    = value(z_offset$)
  value(flexible$) = false$  ! Flexible calc must be handled by the lord.

  ! During parsing the reference energy may not be set.
  ! In this case, do not try to compute things since will get a divide by zero.

  if (lord%value(p0c$) /= 0) then
    vel = c_light * value(p0c$) / value(E_tot$)
    value(ref_time_start$) = lord%value(ref_time_start$) + offset / vel
    value(delta_ref_time$) = value(t_offset$) + value(l$) / vel
    slave%ref_time = value(ref_time_start$) + value(delta_ref_time$)
  endif
endif

!

slave%value = value
slave%mat6_calc_method            = lord%mat6_calc_method
slave%tracking_method             = lord%tracking_method
slave%field_master                = lord%field_master
slave%ptc_integration_type        = lord%ptc_integration_type
slave%taylor_map_includes_offsets = lord%taylor_map_includes_offsets
slave%symplectify                 = lord%symplectify
slave%multipoles_on               = lord%multipoles_on
slave%scale_multipoles            = lord%scale_multipoles
slave%is_on                       = lord%is_on
slave%csr_method                  = lord%csr_method
slave%space_charge_method         = lord%space_charge_method

if (slave%tracking_method == bmad_standard$ .and. slave%key == em_field$) slave%tracking_method = runge_kutta$
if (slave%mat6_calc_method == bmad_standard$ .and. slave%key == em_field$) slave%mat6_calc_method = tracking$

! The slave can have more than one lord here if the other lords are pipes.

if (slave%n_lord == 1) then
  slave%aperture_type               = lord%aperture_type
  slave%aperture_at                 = no_aperture$

  select case (lord%aperture_at)
  case (continuous$, wall_transition$, surface$)
    slave%aperture_at = lord%aperture_at
  case (entrance_end$)
    if (include_entrance) slave%aperture_at = entrance_end$
  case (exit_end$)
    if (include_exit) slave%aperture_at = exit_end$
  case (both_ends$)
    if (include_entrance) slave%aperture_at = entrance_end$
    if (include_exit) slave%aperture_at = exit_end$
  end select
endif

! wiggler fields and electro-magnetic fields

if (slave%key == wiggler$ .or. slave%key == undulator$) slave%value(n_period$) = lord%value(n_period$) * coef

! If an sbend:
!     1) renormalize the angles
!     2) zero the face angles next to the split

if (slave%key == sbend$) then
  ifr = nint(value(fringe_at$))
  if (ifr == no_end$ .or. ifr == exit_end$) then
    slave%value(e1$)    = 0
    slave%value(h1$)    = 0
    slave%value(fint$)  = 0
    slave%value(hgap$)  = 0
  endif

  if (ifr == no_end$ .or. ifr == entrance_end$) then
    slave%value(e2$)    = 0
    slave%value(h2$)    = 0
    slave%value(fintx$) = 0
    slave%value(hgapx$) = 0
  endif
endif                       

! If there are long range wakes they must be scaled.

if (associated (slave%wake)) then
  slave%wake%lr%mode%freq_in   = lord%wake%lr%mode%freq_in
  slave%wake%lr%mode%freq      = lord%wake%lr%mode%freq
  slave%wake%lr%mode%damp      = lord%wake%lr%mode%damp
  slave%wake%lr%mode%phi       = lord%wake%lr%mode%phi
  slave%wake%lr%mode%angle     = lord%wake%lr%mode%angle
  slave%wake%lr%mode%m         = lord%wake%lr%mode%m
  slave%wake%lr%mode%polarized = lord%wake%lr%mode%polarized
  slave%wake%lr%mode%r_over_q  = lord%wake%lr%mode%r_over_q * coef
endif

!

if (slave%key == lcavity$ .or. slave%key == rfcavity$) call compute_slave_coupler (slave)

if (slave%key == lcavity$) then
  slave%value(coupler_at$) = no_end$
  slave%value(e_loss$) = lord%value(e_loss$) * coef
endif

call attribute_bookkeeper (slave, .true.)

err_flag = .false.

end subroutine makeup_super_slave1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine compute_slave_coupler (slave)
!
! This routine is not meant for general use.
!-

subroutine compute_slave_coupler (slave)

type (ele_struct) slave
type (ele_struct), pointer :: lord
real(rp) value(num_ele_attrib$)
integer i
logical entrance_end, exit_end

!

do i = 1, slave%n_lord
  lord => pointer_to_lord (slave, i)
  if (lord%key /= rfcavity$ .and. lord%key /= lcavity$) cycle
  entrance_end = lord_edge_aligned (slave, entrance_end$, lord)
  exit_end = lord_edge_aligned (slave, exit_end$, lord)
  if (entrance_end .and. exit_end) then
    slave%value(coupler_at$) = both_ends$
  elseif (entrance_end) then
    slave%value(coupler_at$) = entrance_end$
  elseif (exit_end) then
    slave%value(coupler_at$) = exit_end$
  else
    slave%value(coupler_at$) = no_end$
  endif
enddo

end subroutine compute_slave_coupler

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_control_slave (lat, slave, err_flag)
!
! This routine is not meant for general use.
!-

subroutine makeup_control_slave (lat, slave, err_flag)

type (lat_struct), target :: lat
type (ele_struct), target :: slave
type (ele_struct), pointer :: lord, slave0, my_lord, my_slave
type (branch_struct), pointer :: branch
type (floor_position_struct) slave_floor
type (all_pointer_struct) ptr_attrib(20), a_ptr
type (control_struct), pointer :: control

real(rp) ds, s_slave, val_attrib(20)
real(rp) t, x_off, y_off, x_pitch, y_pitch, l_gs(3), l_g_off(3), l_slave_off_tot(3)
real(rp) w_slave_inv(3,3), w_gird(3,3), w_gs(3,3), w_gird_mis_tot(3,3)
real(rp) w_slave_mis_tot(3,3), w_slave_mis(3,3), dr, length
real(rp), pointer :: v(:), vs(:), tt

integer i, j, ix, iv, ix_slave, icom, l_stat, n_attrib
logical err_flag, on_an_offset_girder

character(*), parameter :: r_name = 'makeup_control_slave'

!

err_flag = .false.
branch => slave%branch
n_attrib = 0
val_attrib = 0

l_stat = slave%lord_status
ix_slave = slave%ix_ele

on_an_offset_girder = .false.

do i = 1, slave%n_lord
  lord => pointer_to_lord(slave, i, control)

  if (lord%lord_status == control_lord$) cycle
  if (lord%lord_status == multipass_lord$) cycle
  if (lord%key == group$) cycle

  if (lord%key == girder$) then
    if (.not. has_orientation_attributes(slave)) cycle   ! Example: Match element does not have orientation.
    v => lord%value
    vs => slave%value

    if (v(x_offset_tot$) == 0 .and. v(y_offset_tot$) == 0 .and. v(z_offset_tot$) == 0 .and. &
        v(x_pitch_tot$) == 0 .and. v(y_pitch_tot$) == 0 .and. v(tilt_tot$) == 0) cycle
    ! Transformation to get the total misalignment:
    !   T_slave_mis_tot = G_slave^-1 . G_gird . T_gird_mis_tot . G_gird^-1 . G_slave . T_slave_mis
    ! where G = transformation wrt Global coordinate system.
    select case (slave%key)
    case (crystal$, mirror$, multilayer_mirror$)
      slave0 => pointer_to_next_ele (slave, -1) 
      slave_floor = slave0%floor
    case default
      call ele_geometry (slave%floor, slave, slave_floor, -0.5_rp)
    end select
    call floor_angles_to_w_mat (slave_floor%theta, slave_floor%phi, slave_floor%psi, w_mat_inv = w_slave_inv)
    call floor_angles_to_w_mat (lord%floor%theta, lord%floor%phi, lord%floor%psi, w_gird)
    w_gs = matmul(w_slave_inv, w_gird)
    l_gs = matmul(w_slave_inv, (lord%floor%r - slave_floor%r))

    call floor_angles_to_w_mat (v(x_pitch_tot$), v(y_pitch_tot$), v(tilt_tot$), w_gird_mis_tot)
    l_slave_off_tot = matmul(w_gs, [v(x_offset_tot$), v(y_offset_tot$), v(z_offset_tot$)]) + l_gs
    w_slave_mis_tot = matmul(w_gs, w_gird_mis_tot)

    w_slave_mis_tot = matmul(w_slave_mis_tot, transpose(w_gs))     ! Transpose = inverse
    l_slave_off_tot = matmul(w_slave_mis_tot, -l_gs) + l_slave_off_tot

    call floor_angles_to_w_mat (vs(x_pitch$), vs(y_pitch$), vs(tilt$), w_slave_mis)
    l_slave_off_tot = matmul(w_slave_mis_tot, [vs(x_offset$), vs(y_offset$), vs(z_offset$)]) + l_slave_off_tot
    w_slave_mis_tot = matmul(w_slave_mis_tot, w_slave_mis)

    ! If slave is an sbend then correct offsets since roll axis is displaced from the bend center.
    if ((slave%key == sbend$ .or. slave%key == rf_bend$) .and. vs(g$) /= 0) then
      call floor_w_mat_to_angles (w_slave_mis_tot, vs(x_pitch_tot$), vs(y_pitch_tot$), vs(roll_tot$))
      dr = (1 - cos(vs(angle$)/2)) / vs(g$)
      vs(x_offset_tot$) = l_slave_off_tot(1) + dr * (1 - cos(vs(roll_tot$)))
      vs(y_offset_tot$) = l_slave_off_tot(2) - dr * sin(vs(roll_tot$))  
      vs(z_offset_tot$) = l_slave_off_tot(3) 
    else
      call floor_w_mat_to_angles (w_slave_mis_tot, vs(x_pitch_tot$), vs(y_pitch_tot$), vs(tilt_tot$))
      vs(x_offset_tot$) = l_slave_off_tot(1)
      vs(y_offset_tot$) = l_slave_off_tot(2)
      vs(z_offset_tot$) = l_slave_off_tot(3)
    endif

    slave%bookkeeping_state%has_misalign = .true.
    on_an_offset_girder = .true.

    cycle
  endif

  if (lord%key /= overlay$) then
    call out_io (s_abort$, r_name, 'THE LORD IS NOT AN OVERLAY: ', lord%name)
    err_flag = .true.
    return
  endif

  ! overlay lord

  select case (control%attribute)
  case ('X_LIMIT')
    call overlay_change_this(lord, slave%value(x1_limit$), control, val_attrib, ptr_attrib);  if (err_flag) return
    call overlay_change_this(lord, slave%value(x2_limit$), control, val_attrib, ptr_attrib);  if (err_flag) return
  case ('Y_LIMIT')
    call overlay_change_this(lord, slave%value(y1_limit$), control, val_attrib, ptr_attrib);  if (err_flag) return
    call overlay_change_this(lord, slave%value(y2_limit$), control, val_attrib, ptr_attrib);  if (err_flag) return
  case ('APERTURE')
    call overlay_change_this(lord, slave%value(x1_limit$), control, val_attrib, ptr_attrib);  if (err_flag) return
    call overlay_change_this(lord, slave%value(x2_limit$), control, val_attrib, ptr_attrib);  if (err_flag) return
    call overlay_change_this(lord, slave%value(y1_limit$), control, val_attrib, ptr_attrib);  if (err_flag) return
    call overlay_change_this(lord, slave%value(y2_limit$), control, val_attrib, ptr_attrib);  if (err_flag) return
  case default
    err_flag = .not. attribute_free (slave, control%attribute, .true., .true., .true.)
    if (err_flag) then
      call out_io (s_abort$, r_name, 'OVERLAY LORD: ' // lord%name, &
           'IS TRYING TO VARY NON-FREE ATTRIBUTE: ' // trim(slave%name) // '[' // trim(control%attribute) // ']')
      return
    endif

    call pointer_to_attribute (slave, control%attribute, .true., a_ptr, err_flag, do_unlink = .true.)
    if (err_flag) then
      if (global_com%exit_on_error) call err_exit
      return
    endif
    call overlay_change_this(lord, a_ptr%r, control, val_attrib, ptr_attrib)
  end select

enddo

! Transfer values from val_attrib to slave elements

do iv = 1, n_attrib

  a_ptr%r => ptr_attrib(iv)%r
  ! If there is no significant change in the attribute then do not set bookkeeping flags stale.
  if (abs(a_ptr%r - val_attrib(iv)) <= small_rel_change$ * (abs(a_ptr%r) + abs(val_attrib(iv)))) cycle
  a_ptr%r = val_attrib(iv)
  call set_ele_status_stale (slave, attribute_group$)
  call set_flags_for_changed_attribute (slave, a_ptr%r)

  ! If varying length then must update any associated super_lords and super_slaves

  if (associated(a_ptr%r, slave%value(l$))) then

    ! If varying a  super_lord length then adjust last super_slave length to match.
    if (slave%lord_status == super_lord$) then
      length = 0
      do i = 1, slave%n_slave-1
        my_slave => pointer_to_slave(slave, i)
        length = length + my_slave%value(l$)
      enddo
      my_slave => pointer_to_slave(slave, slave%n_slave)
      my_slave%value(l$) = a_ptr%r + slave%value(lord_pad1$) + slave%value(lord_pad2$) - length
      if (my_slave%value(l$) < 0 .and. all(my_slave%key /= [drift$, pipe$, patch$])) then
        call out_io (s_error$, r_name, 'APPLICATION OF OVERLAY LORD: ' // lord%name, &
                'IS MAKING THE LENGTH ELEMENT: ' // my_slave%name, 'LESS THAN ZERO')
      endif
      call set_flags_for_changed_attribute (my_slave, my_slave%value(l$))
    else
      my_slave => slave
    endif

    ! If varying a super_slave length then vary all associated super_lord lengths to match.
    if (my_slave%slave_status == super_slave$) then
      do i = 1, my_slave%n_lord
        my_lord => pointer_to_lord(my_slave, i)
        if (my_lord%lord_status /= super_lord$) cycle
        length = 0
        do j = 1, my_lord%n_slave
          slave0 => pointer_to_slave(my_lord, j)
          length = length + slave0%value(l$)
        enddo
        my_lord%value(l$) = length - my_lord%value(lord_pad1$) - my_lord%value(lord_pad2$)
        call set_flags_for_changed_attribute (my_lord, my_lord%value(l$))
      enddo
    endif
  endif

  call s_calc (lat)

enddo

! If no girder then simply transfer tilt to tilt_tot, etc.

if (.not. on_an_offset_girder .and. has_orientation_attributes(slave)) then
  select case (slave%key)
  case (sbend$, rf_bend$)
    slave%value(roll_tot$)     = slave%value(roll$)
    slave%value(ref_tilt_tot$) = slave%value(ref_tilt$)
  case (crystal$, mirror$, multilayer_mirror$)
    slave%value(tilt_tot$)     = slave%value(tilt$)
    slave%value(ref_tilt_tot$) = slave%value(ref_tilt$)
  case default
    slave%value(tilt_tot$)     = slave%value(tilt$)
  end select

  slave%value(x_offset_tot$) = slave%value(x_offset$)
  slave%value(y_offset_tot$) = slave%value(y_offset$)
  slave%value(z_offset_tot$) = slave%value(z_offset$)
  slave%value(x_pitch_tot$)  = slave%value(x_pitch$)
  slave%value(y_pitch_tot$)  = slave%value(y_pitch$)
endif

slave%bookkeeping_state%control = ok$

!-------------------------------------------------------------------------------
contains

! Add to slave attribute the contribution from a lord overlay.

subroutine overlay_change_this (lord, r_attrib, c, val_attrib, ptr_attrib)

type (ele_struct), pointer :: lord
type (control_struct) c
type (all_pointer_struct) ptr_attrib(:)

real(rp), target :: r_attrib
real(rp) val_attrib(:)
integer iv
logical err_flag, ok

character(100) err_str

! First evaluate the contribution from the overlay lord

if (.not. lord%is_on) return

if (allocated(c%stack)) then
  c%value = expression_stack_value(c%stack, err_flag, err_str, lord%control%var, .false.)
  if (err_flag) then
    call out_io (s_error$, r_name, err_str, 'FOR SLAVE: ' // slave%name, 'OF LORD: ' // lord%name)
    err_flag = .true.
    return
  endif

else
  c%value = knot_interpolate (lord%control%x_knot, c%y_knot, lord%control%var(1)%value, &
                                                                      nint(lord%value(interpolation$)), err_flag)
  if (err_flag) then
    call out_io (s_error$, r_name, 'VARIABLE VALUE OUTSIDE OF SPLINE KNOT RANGE.')
    return
  endif
endif

! If the contribution (c%value) contributes to a slave attribute that is on the val_attrib list then
! just add c%value to the slave attribute

do iv = 1, n_attrib
  if (.not. associated(ptr_attrib(iv)%r, r_attrib)) cycle
  val_attrib(iv) = val_attrib(iv) + c%value
  return
enddo

! Must be a slave attribute that is not in the val_attrib list
! So add this slave attribute to the list and set the value to c%value

n_attrib = n_attrib + 1
ptr_attrib(n_attrib)%r => r_attrib
val_attrib(n_attrib) = c%value

end subroutine overlay_change_this

end subroutine makeup_control_slave 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine aperture_bookkeeper (ele)
!
! Routine to calculate aperture limits when ele%attribute_type is set to auto_aperture$
!
! Input:
!   ele   -- ele_struct: Element with aperture.
!
! Output:
!   ele   -- ele_struct: Element with apertures set.
!-

subroutine aperture_bookkeeper (ele)

type (ele_struct), target :: ele
type (surface_displacement_struct), pointer :: displacement
type (surface_h_misalign_struct), pointer :: h_misalign
type (surface_segmented_struct), pointer :: segmented
type (pixel_detec_struct), pointer :: pixel
type (wall3d_section_struct), pointer :: sec

real(rp) angle, r_wall, dr_dtheta, x, y

integer i, j
logical is_set

character(*), parameter :: r_name = 'aperture_bookkeeper'

!

select case (ele%key)
case (diffraction_plate$, mask$) 

  ele%value(x1_limit$) = -100
  ele%value(y1_limit$) = -100
  ele%value(x2_limit$) = -100
  ele%value(y2_limit$) = -100

  do i = 1, size(ele%wall3d(1)%section)
    sec => ele%wall3d(1)%section(i)
    if (sec%type == opaque$) cycle
    do j = 1, 100
      angle = twopi * j / 100.0_rp
      call calc_wall_radius (sec%v, cos(angle), sin(angle), r_wall, dr_dtheta)
      x = sec%r0(1) + 1.001 * r_wall * cos(angle)
      y = sec%r0(2) + 1.001 * r_wall * sin(angle)
      ele%value(x1_limit$) = max(ele%value(x1_limit$), -x)
      ele%value(y1_limit$) = max(ele%value(y1_limit$), -y)
      ele%value(x2_limit$) = max(ele%value(x1_limit$), x)
      ele%value(y2_limit$) = max(ele%value(y1_limit$), y)
    enddo
  enddo


! Element not a diffraction_plate nor a mask.

case default
  if (.not. associated (ele%photon)) then
    call out_io (s_error$, r_name, 'ELEMENT APERTURE TYPE SET TO "SURFACE" BUT', &
                                   'THERE IS NOT A SURFACE ASSOCIATED WITH ELEMENT: ' // ele%name)
    return
  endif

  is_set = .false.
  ele%value(x1_limit$) = -1e30
  ele%value(y1_limit$) = -1e30
  ele%value(x2_limit$) = -1d30
  ele%value(y2_limit$) = -1d30

  if (allocated(ele%photon%pixel%pt)) then
    pixel => ele%photon%pixel
    ele%value(x1_limit$) = max(ele%value(x1_limit$), -(pixel%r0(1) + (lbound(pixel%pt, 1) - 0.5) * pixel%dr(1)))
    ele%value(y1_limit$) = max(ele%value(y1_limit$), -(pixel%r0(2) + (lbound(pixel%pt, 2) - 0.5) * pixel%dr(2)))
    ele%value(x2_limit$) = max(ele%value(x2_limit$),  (pixel%r0(1) + (ubound(pixel%pt, 1) + 0.5) * pixel%dr(1)))
    ele%value(y2_limit$) = max(ele%value(y2_limit$),  (pixel%r0(2) + (ubound(pixel%pt, 2) + 0.5) * pixel%dr(2)))
    is_set = .true.
  endif

  if (allocated(ele%photon%segmented%pt)) then
    segmented => ele%photon%segmented
    ele%value(x1_limit$) = max(ele%value(x1_limit$), -(segmented%r0(1) + (lbound(segmented%pt, 1) - 0.5) * segmented%dr(1)))
    ele%value(y1_limit$) = max(ele%value(y1_limit$), -(segmented%r0(2) + (lbound(segmented%pt, 2) - 0.5) * segmented%dr(2)))
    ele%value(x2_limit$) = max(ele%value(x2_limit$),  (segmented%r0(1) + (ubound(segmented%pt, 1) + 0.5) * segmented%dr(1)))
    ele%value(y2_limit$) = max(ele%value(y2_limit$),  (segmented%r0(2) + (ubound(segmented%pt, 2) + 0.5) * segmented%dr(2)))
    is_set = .true.
  endif

  if (allocated(ele%photon%h_misalign%pt)) then
    h_misalign => ele%photon%h_misalign
    ele%value(x1_limit$) = max(ele%value(x1_limit$), -(h_misalign%r0(1) + (lbound(h_misalign%pt, 1) - 0.5) * h_misalign%dr(1)))
    ele%value(y1_limit$) = max(ele%value(y1_limit$), -(h_misalign%r0(2) + (lbound(h_misalign%pt, 2) - 0.5) * h_misalign%dr(2)))
    ele%value(x2_limit$) = max(ele%value(x2_limit$),  (h_misalign%r0(1) + (ubound(h_misalign%pt, 1) + 0.5) * h_misalign%dr(1)))
    ele%value(y2_limit$) = max(ele%value(y2_limit$),  (h_misalign%r0(2) + (ubound(h_misalign%pt, 2) + 0.5) * h_misalign%dr(2)))
    is_set = .true.
  endif

  if (allocated(ele%photon%displacement%pt)) then
    displacement => ele%photon%displacement
    ele%value(x1_limit$) = max(ele%value(x1_limit$), -(displacement%r0(1) + (lbound(displacement%pt, 1) - 0.5) * displacement%dr(1)))
    ele%value(y1_limit$) = max(ele%value(y1_limit$), -(displacement%r0(2) + (lbound(displacement%pt, 2) - 0.5) * displacement%dr(2)))
    ele%value(x2_limit$) = max(ele%value(x2_limit$),  (displacement%r0(1) + (ubound(displacement%pt, 1) + 0.5) * displacement%dr(1)))
    ele%value(y2_limit$) = max(ele%value(y2_limit$),  (displacement%r0(2) + (ubound(displacement%pt, 2) + 0.5) * displacement%dr(2)))
    is_set = .true.
  endif

  if (.not. is_set) then
    call out_io (s_error$, r_name, 'ELEMENT APERTURE TYPE SET TO "SURFACE" BUT', &
                                   'NO GRID IS DEFINED: ' // ele%name)
    return
  endif

end select

end subroutine aperture_bookkeeper

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine attributes_need_bookkeeping (ele, dval) result (is_needed)
!
! Routine, when bmad_com%auto_bookkeeper = True, to decide if attribute
! bookkeeping needs to be done for an element.
!
! Input:
!   ele         -- ele_struct: Element under consideration.
!
! Output:
!   dval(:)     -- real(rp), optional: Array of differences between old and new ele%value(:) values.
!   ele%bookkeeping_state%attributes 
!                    -- Set ok$ if not needed, stale$ otherwise
!-

subroutine attributes_need_bookkeeping (ele, dval)

type (ele_struct), target :: ele
type (ele_struct), pointer :: ele0
real(rp), optional :: dval(:)
real(rp) dv(num_ele_attrib$)
integer i

!

select case (ele%key)
case (overlay$, group$, hybrid$)
  ele%bookkeeping_state%attributes = ok$
  if (present(dval)) dval = 0
  return
end select

! Check_sum is a hash number that is used to see if a value has been changed.
! This is used implicitly in attribute_bookkeeper.

ele0 => ele
if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) ele0 => pointer_to_super_lord(ele0)
ele%value(check_sum$) = 0

if (associated(ele0%a_pole)) then
  do i = 0, ubound(ele0%a_pole, 1)
    ele%value(check_sum$) = ele%value(check_sum$) + fraction(ele0%a_pole(i)) + fraction(ele0%b_pole(i))
    ele%value(check_sum$) = ele%value(check_sum$) + (exponent(ele0%a_pole(i)) + exponent(fraction(ele0%b_pole(i)))) / 10
  enddo
endif


if (associated(ele0%a_pole_elec)) then
  do i = 0, ubound(ele0%a_pole_elec, 1)
    ele%value(check_sum$) = ele%value(check_sum$) + fraction(ele0%a_pole_elec(i)) + fraction(ele0%b_pole_elec(i))
    ele%value(check_sum$) = ele%value(check_sum$) + (exponent(ele0%a_pole_elec(i)) + exponent(fraction(ele0%b_pole_elec(i)))) / 10
  enddo
endif

if (associated(ele%cartesian_map)) then
  do i = 1, size(ele%cartesian_map)
    ele%value(check_sum$) = ele%value(check_sum$) + ele%cartesian_map(i)%field_scale
  enddo
endif

if (associated(ele%cylindrical_map)) then
  do i = 1, size(ele%cylindrical_map)
    ele%value(check_sum$) = ele%value(check_sum$) + ele%cylindrical_map(i)%field_scale
  enddo
endif

if (associated(ele%gen_grad_map)) then
  do i = 1, size(ele%gen_grad_map)
    ele%value(check_sum$) = ele%value(check_sum$) + ele%gen_grad_map(i)%field_scale
  enddo
endif

if (associated(ele%grid_field)) then
  do i = 1, size(ele%grid_field)
    ele%value(check_sum$) = ele%value(check_sum$) + ele%grid_field(i)%field_scale
  enddo
endif

!

if (present(dval)) then
  dv = abs(ele%value - ele%old_value)
  dv(x1_limit$:y2_limit$) = 0  ! Limit changes do not need bookkeeping
  dval = dv

  if (all(dv == 0) .and. ele%key /= capillary$) then
    ele%bookkeeping_state%attributes = ok$
  else
    ele%bookkeeping_state%attributes = stale$
  endif
endif

end subroutine attributes_need_bookkeeping

end module
