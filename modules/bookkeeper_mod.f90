module bookkeeper_mod

use wall3d_mod
use bmad_utils_mod
use multipole_mod
use geometry_mod
use equality_mod
use em_field_mod
use xraylib_interface
use expression_mod

implicit none

integer, parameter :: save_state$ = 3, restore_state$ = 4, off_and_save$ = 5

private control_bookkeeper1, makeup_control_slave
private makeup_group_lord, makeup_super_slave1, makeup_super_slave

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine set_flags_for_changed_attribute (...)
!
! Routine to mark an element or lattice as modified for use with "intelligent" bookkeeping.
! Also will do some dependent variable bookkeeping when a particular attribute has 
! been altered. Look at this routine's code for more details.
!
! set_flags_for_changed_attribute is an overloaded name for:
!   set_flags_for_changed_lat_attribute (lat, set_dependent)
!   set_flags_for_changed_real_attribute (ele, real_attrib, set_dependent)
!   set_flags_for_changed_inteter_attribute (ele, int_attrib, set_dependent)
!   set_flags_for_changed_logical_attribute (ele, logic_attrib, set_dependent)
!   set_flags_for_changed_all_attribute (ele, all_attrib, set_dependent)
!
! The set_flags_for_changed_lat_attribute (lat) routine is used when one
! does not know what has changed and wants a complete bookkeeping done.
!
! NOTE: The attribute argument MUST be the component that was changed. For example:
!     ele%value(x_offset$) = off_value
!     call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
! And NOT:
!     call set_flags_for_changed_attribute (ele, off_value)  ! WRONG
!
! Input:
!   lat           -- lat_struct: Lattice being modified.
!   ele           -- ele_struct, Element being modified.
!   real_attrib   -- real(rp), optional: Attribute that has been changed.
!                      For example: ele%value(hkick$).
!                      If not present then assume everything has potentially changed.
!   int_attrib    -- integer: Attribute that has been changed.
!                      For example: ele%mat6_calc_method.
!   logic_attrib  -- logical; Attribute that has been changed.
!                      For example: ele%is_on.
!   all_attrib    -- all_pointer_struct: Pointer to attribute.
!   set_dependent -- logical, optional: If False then dependent variable bookkeeping will not be done.
!                     Default is True. Do not set False unless you know what you are doing.
!
! Output:
!   lat  -- lat_struct: Lattice with appropriate changes.
!-

interface set_flags_for_changed_attribute
  module procedure set_flags_for_changed_real_attribute 
  module procedure set_flags_for_changed_integer_attribute 
  module procedure set_flags_for_changed_logical_attribute 
  module procedure set_flags_for_changed_all_attribute 
  module procedure set_flags_for_changed_lat_attribute 
end interface

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lattice_bookkeeper (lat, err_flag)
!
! Subroutine to do a "complete" bookkeeping job on a lattice:
!   lord/slave control
!   reference energy calc
!   s-position calc
!   geometry (floor position) calc
!
! Not done are Twiss, transfer matrices, and orbit calculations.  
!
! Note: This this routine does a complete job of bookking
! and could be unacceptably slow if lat%auto_bookkeeper = True.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat   -- lat_struct: Lattice needing bookkeeping.
!
! Output:
!   lat      -- lat_struct: Lattice with bookkeeping done.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine lattice_bookkeeper (lat, err_flag)

use precision_constants, only: e_muon

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (bookkeeping_state_struct), pointer :: stat

integer i, j

logical, optional :: err_flag
logical found, err, auto_saved

character(20), parameter :: r_name = 'lattice_bookkeeper'

! Set PTC E_MUON just to make sure it has the same value as bmad_com%electric_dipole_moment

E_MUON = bmad_com%electric_dipole_moment

! Turn on intelligent bookkeeping while this routine is running
! If bookkeeping has not been intelligent then mark everything as stale.

auto_saved = bmad_com%auto_bookkeeper
if (auto_saved) then
  bmad_com%auto_bookkeeper = .false.
  do i = 0, ubound(lat%branch, 1)
    branch => lat%branch(i)
    do j = 0, branch%n_ele_max
      call set_ele_status_stale (branch%ele(j), all_groups$, .false.) 
      call attributes_need_bookkeeping(branch%ele(j))
    enddo
  enddo
endif

if (present(err_flag)) err_flag = .true.

! Cases where bookkeeping routines have to be called multiple times:
!   * Control_bookkeeper to make sure that girders are correctly adjusted.
!       Especially when lattice has bends with field_master = True.
!   * lat_geometry in case the energy changes and there is a bend with field_master = T.
!   * Recompute ref energy for cases where a flexible patch has changed its geometry and this
!       affects the reference energy due to the presence of lcavity elements.

do i = 1, 3
  call control_bookkeeper (lat)
  call lat_geometry (lat)
  call s_calc (lat)
  call lat_compute_ref_energy_and_time (lat, err)
  if (err) return
enddo

call ptc_bookkeeper (lat)

! See if all status flags have been properly reset.
! Exception is mat6 flag since the bookkeeping routines do not touch this.

stat => lat%lord_state
if (stat%control == stale$ .or. stat%attributes == stale$ .or. stat%floor_position == stale$ .or. &
    stat%s_position == stale$ .or. stat%ref_energy == stale$) then
  call out_io (s_info$, r_name, 'Stale bookkeeping lord_status flags detected.', &
                                'Please contact DCS!', 'Status: \5i6\ ', &
          i_array = [stat%attributes, stat%control, stat%floor_position, stat%s_position, stat%ref_energy])
endif
call reset_status_flags_to_ok(stat)

do i = 0, ubound(lat%branch, 1)

  branch => lat%branch(i)
  stat => branch%param%bookkeeping_state
  if (stat%control == stale$ .or. stat%attributes == stale$ .or. stat%floor_position == stale$ .or. &
      stat%s_position == stale$ .or. stat%ref_energy == stale$) then
    call out_io (s_info$, r_name, 'Stale bookkeeping status flags detected at branch: \i0\.', &
                                  'Please contact DCS!', 'Status: \5i6\ ', &
            i_array = [i, stat%attributes, stat%control, stat%floor_position, stat%s_position, stat%ref_energy])
  endif
  call reset_status_flags_to_ok(stat)

  do j = 0, branch%n_ele_max
    ele => branch%ele(j)
    if (ele%key == null_ele$) cycle 
    stat => ele%bookkeeping_state
    if (stat%control == stale$ .or. stat%attributes == stale$ .or. stat%floor_position == stale$ .or. &
        stat%s_position == stale$ .or. stat%ref_energy == stale$) then
      call out_io (s_info$, r_name, &
            'Stale bookkeeping status flags detected at element: ' // trim(ele%name) // ' (\i0\>>\i0\).', &
            'Please contact DCS!', 'Status: \5i6\ ', &
            i_array = [i, j, stat%attributes, stat%control, stat%floor_position, stat%s_position, stat%ref_energy])
    endif
    call reset_status_flags_to_ok(stat)
  enddo

enddo

!

if (present(err_flag)) err_flag = .false.
bmad_com%auto_bookkeeper = auto_saved

!----------------------------------------------------------
contains

subroutine reset_status_flags_to_ok (stat)
  type (bookkeeping_state_struct) stat

  if (stat%control /= ok$        .and. stat%control /= super_ok$)        stat%control = ok$
  if (stat%attributes /= ok$     .and. stat%attributes /= super_ok$)     stat%attributes = ok$
  if (stat%floor_position /= ok$ .and. stat%floor_position /= super_ok$) stat%floor_position = ok$
  if (stat%s_position /= ok$     .and. stat%s_position /= super_ok$)     stat%s_position = ok$
  if (stat%ref_energy /= ok$     .and. stat%ref_energy /= super_ok$)     stat%ref_energy = ok$

end subroutine reset_status_flags_to_ok

end subroutine lattice_bookkeeper

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine control_bookkeeper (lat, ele, mark_eles_as_stale)
!
! Subroutine to transfer attibute information from lord to slave elements.
!
! If ele argument is present, bookkeeping will include all the slaves of ele
! but none of the lords.
!
! Note: This subroutine will call attribute_bookkeeper.
! Note: To do a complete bookkeeping job on a lattice use:
!   lattice_bookkeeper
!
! Modules needed:
!   use bmad
!
! Input:
!   lat    -- lat_struct: lattice to be used
!   ele    -- Ele_struct, optional: Element whose attribute values 
!               have been changed. If not present bookkeeping will be done 
!               for all elements.
!   mark_eles_as_stale
!          -- Logical, optional: If true then mark all elements for control bookkeeping
!               if using auto_bookkeeper. Default is True. Do not set this argument
!               unless you know what you are doing.
!-

recursive subroutine control_bookkeeper (lat, ele, mark_eles_as_stale)

type (lat_struct), target :: lat
type (ele_struct), optional :: ele
type (ele_struct), pointer :: slave, lord, branch_ele, ele2
type (branch_struct), pointer :: branch

integer i, j, ie, ib, n1, n2

logical, optional :: mark_eles_as_stale

character(*), parameter :: r_name = 'control_bookkeeper'

!----------------------------------------------------------------
! If ele is present we only do bookkeeping for this one element and its slaves

if (present(ele)) then
  call control_bookkeeper1 (lat, ele, .true.)
  return
endif

!----------------------------------------------------------------
! Else we need to make up all the lords...
! First mark all the elements needing bookkeeping

if (bmad_com%auto_bookkeeper .and. logic_option(.true., mark_eles_as_stale)) then
  lat%ele(:)%bookkeeping_state%control = stale$  ! Bookkeeping done on this element yet?
endif

! Bookkkeeping is done from the top level down.
! The top level elements are those lord elements that have no lords on top of them.

ie_loop: do ie = lat%n_ele_track+1, lat%n_ele_max
  ele2 => lat%ele(ie)
  if (ele2%key == null_ele$) then
    ele2%bookkeeping_state%control = ok$
    ele2%bookkeeping_state%attributes = ok$
    cycle
  endif
  if (ele2%n_lord > 0) cycle
  call control_bookkeeper1 (lat, ele2, .false.)
enddo ie_loop

! And now bookkeeping for the elements in the tracking lattice

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (.not. bmad_com%auto_bookkeeper .and. branch%param%bookkeeping_state%control /= stale$ .and. &
                                           branch%param%bookkeeping_state%attributes /= stale$) cycle

  do ie = 0, branch%n_ele_track
    ele2 => branch%ele(ie)
    if (ele2%bookkeeping_state%control /= stale$ .and. ele2%bookkeeping_state%attributes /= stale$) cycle
    call attribute_bookkeeper (ele2, branch%param)
    ele2%bookkeeping_state%control = ok$
  enddo

  branch%param%bookkeeping_state%attributes = ok$
  branch%param%bookkeeping_state%control = ok$

enddo

lat%lord_state%control = ok$
lat%lord_state%attributes = ok$

end subroutine control_bookkeeper

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine control_bookkeeper1 (lat, ele, force_bookkeeping)
!
! This routine is for control bookkeeping for a single element.
! This subroutine is only to be called from control_bookkeeper and is
! not meant for general use.
!-

recursive subroutine control_bookkeeper1 (lat, ele, force_bookkeeping)

type (lat_struct), target :: lat
type (ele_struct) ele
type (ele_struct), pointer :: slave

integer i
logical call_a_bookkeeper, force_bookkeeping

! Only do bookkeeping on this element if it is stale or bookkeeping is forced by the calling routine.

if (ele%bookkeeping_state%control == stale$ .or. ele%bookkeeping_state%attributes == stale$ .or. force_bookkeeping) then

  ! First make sure the attribute bookkeeping for this element is correct since
  ! the makeup_*_slave routines may need it.

  call attribute_bookkeeper (ele, lat%branch(ele%ix_branch)%param, force_bookkeeping)

  ! Slave bookkeeping

  call_a_bookkeeper = .false.

  if (ele%slave_status == super_slave$) then
    ! Attrubute bookkeeping is done in the makeup_super_slave
    call makeup_super_slave (lat, ele)

  elseif (ele%slave_status == multipass_slave$) then
    call makeup_multipass_slave (lat, ele)
    if (ele%n_lord > 1) call makeup_control_slave (lat, ele)
    call_a_bookkeeper = .true.

  elseif (ele%n_lord > 0 .and. ele%slave_status /= slice_slave$) then
    call makeup_control_slave (lat, ele)
    call_a_bookkeeper = .true.

  endif

  ! Lord bookkeeping

  if (ele%key == group$) then
    call makeup_group_lord (lat, ele)
    call_a_bookkeeper = .true.
  endif

  ! If bookkeeping has been done by a makeup_*_slave routine then
  ! attribute_bookkeeper must be called again.
  ! This is true even if the lattice is static since a slave element
  ! can have its lord's dependent attribute values.
  ! Example: super_slave will, at this point, have its lord's num_steps value but 
  ! num_steps in the slave is different from the lord due to differences in length.

  if (call_a_bookkeeper) call attribute_bookkeeper (ele, lat%branch(ele%ix_branch)%param, force_bookkeeping)

  ele%bookkeeping_state%control = ok$

endif

! Recursively call this routine on the slaves

do i = 1, ele%n_slave
  slave => pointer_to_slave (ele, i)
  call control_bookkeeper1 (lat, slave, force_bookkeeping)
enddo

end subroutine control_bookkeeper1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_group_lord (lat, lord)
!
! Subroutine to calculate the attributes of group slave elements.
! This routine is private to bookkeeper_mod.
!-

Subroutine makeup_group_lord (lat, lord)   

type (lat_struct), target :: lat
type (ele_struct) :: lord
type (ele_struct), pointer :: slave, slave2
type (control_struct), pointer :: control

real(rp), pointer :: r_ptr

integer ix_attrib, i, j

logical moved, err_flag

character(20) :: r_name = 'makeup_group_lord'

!

moved = .false.   ! have we longitudinally moved an element?

do i = 1, lord%n_slave
  slave => pointer_to_slave(lord, i, control)
  ix_attrib = control%ix_attrib
  if (ix_attrib == l$) moved = .true.

  select case (ix_attrib)

  !---------
  ! Edge: Varying lengths takes special code.

  case (start_edge$, end_edge$, s_position$, accordion_edge$, l$)

    if (slave%lord_status == multipass_lord$) then
      do j = 1, slave%n_slave
        slave2 => pointer_to_slave (slave, j)
        call change_this_edge (slave2, control)
      enddo
    else
      call change_this_edge (slave, control)
    endif

  !---------
  ! x_limit, y_limit, aperture

  case (x_limit$)
    call group_change_this (slave, x1_limit$, control, 1)
    call group_change_this (slave, x2_limit$, control, 1)

  case (y_limit$)
    call group_change_this (slave, y1_limit$, control, 1)
    call group_change_this (slave, y2_limit$, control, 1)

  case (aperture$) 
    call group_change_this (slave, x1_limit$, control, 1)
    call group_change_this (slave, x2_limit$, control, 1)
    call group_change_this (slave, y1_limit$, control, 1)
    call group_change_this (slave, y2_limit$, control, 1)

  !---------
  ! All else

  case default

    call group_change_this (slave, ix_attrib, control, 1)

  end select

enddo

!---------
! End stuff

if (moved) then
  call s_calc (lat)       ! recompute s distances
  call lat_geometry (lat)
endif

do i = 1, size(lord%control_var)
  lord%control_var(i)%old_value = lord%control_var(i)%value ! update old
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
elseif (this_slave%ix_ele < lat%n_ele_track) then
  ix_min = this_slave%ix_ele
  ix_max = this_slave%ix_ele
  branch => lat%branch(this_slave%ix_branch)
else
  call out_io (s_error$, r_name, &
                'A GROUP IS NOT ALLOWED TO CONTROL', &
                'A ' // control_name(slave%slave_status), &
                'YOU TRIED TO CONTROL: ' // slave%name)
  return
endif

! now that we have the ends we find the elements to either side whose length
! the group can adjust

if (ix_attrib /= end_edge$ .and. ix_attrib /= l$) then
  ix1 = ix_min - 1
  do
    if (attribute_name(branch%ele(ix1), l$) == 'L') exit  ! If has length attribute
    ix1 = ix1 - 1
    if (ix1 < 0) then
      call out_io (s_error$, r_name, &
                    'START_EDGE OF CONTROLED', &
                    'ELEMENT IS AT BEGINNING OF LAT AND CANNOT BE', &
                    'VARIED FOR GROUP: ' // lord%name)
      return
    endif
  enddo
endif

if (ix_attrib /= start_edge$ .and. ix_attrib /= l$) then
  ix2 = ix_max + 1 
  do
    if (attribute_name(branch%ele(ix2), l$) == 'L') exit  ! If has length attribute
    ix2 = ix2 + 1
    if (ix2 > branch%n_ele_track) then
      call out_io (s_error$, r_name, &
                    'END_EDGE OF CONTROLED', &
                    'ELEMENT IS AT END OF LAT AND CANNOT BE', &
                    'VARIED FOR GROUP: ' // lord%name)
      return
    endif
  enddo
endif

! put in changes

select case (ix_attrib)

case (l$)
  call group_change_this (branch%ele(ix_max), l$, ctl, 1)

case (start_edge$)
  call group_change_this (branch%ele(ix_min), l$, ctl, -1)
  call group_change_this (branch%ele(ix1), l$, ctl, 1)

case (end_edge$)
  call group_change_this (branch%ele(ix_max), l$, ctl, 1)
  call group_change_this (branch%ele(ix2), l$, ctl, -1)

case (accordion_edge$)
  call group_change_this (branch%ele(ix_min), l$, ctl, 1)
  call group_change_this (branch%ele(ix1), l$, ctl, -1)

  call group_change_this (branch%ele(ix_max), l$, ctl, 1)
  call group_change_this (branch%ele(ix2), l$, ctl, -1)

case (s_position$)
  call group_change_this (branch%ele(ix1), l$, ctl, 1)
  call group_change_this (branch%ele(ix2), l$, ctl, -1)

case (lord_pad1$)
  call group_change_this (branch%ele(ix1), l$, ctl, 1, this_slave, lord_pad1$)

case (lord_pad2$)
  call group_change_this (branch%ele(ix2), l$, ctl, 1, this_slave, lord_pad1$)

end select

end subroutine change_this_edge

!---------------------------------------------------------------------------------
! contains 
!+
! Note: It is assumed that edge, super_lord, and varied_length optional args are all present if
! any one of them is present.
!-

recursive subroutine group_change_this (ele, ix_attrib, ctl, dir, this_lord, this_pad)

type (ele_struct) ele
type (ele_struct), optional :: this_lord
type (ele_struct), pointer :: my_lord
type (control_struct) ctl
type (all_pointer_struct) :: a_ptr

integer dir
integer ix_attrib, il, ix_slave
integer, optional :: this_pad

real(rp) coef, new_val, val, val_old, delta

character(100) err_str

!

call pointer_to_indexed_attribute (ele, ix_attrib, .false., a_ptr, err_flag)
if (err_flag .and. global_com%exit_on_error) call err_exit

call evaluate_expression_stack (ctl%stack, val, err_flag, err_str, lord%control_var, .false.)
call evaluate_expression_stack (ctl%stack, val_old, err_flag, err_str, lord%control_var, .true.)
delta = val - val_old
if (err_flag) then
  call out_io (s_error$, r_name, err_str, 'FOR SLAVE: ' // slave%name, 'OF LORD: ' // lord%name)
  return
endif
a_ptr%r = a_ptr%r + delta * dir

call set_flags_for_changed_attribute (ele, a_ptr%r)
! super_slave length can be varied by a group so don't check this.
if (ele%slave_status /= super_slave$ .or. ix_attrib /= l$) then
  err_flag = attribute_free (ele, attribute_name(ele, ix_attrib), .true.)
endif

! Pad check

if (ele%lord_status == super_lord$ .and. a_ptr%r < 0) then
  if (ix_attrib == lord_pad1$) then
    call out_io (s_error$, r_name, 'GROUP ELEMENT: ' // lord%name, &
                                   'CONTROLS SUPER_LORD: ' // ele%name, &
                                   'AND LORD_PAD1 IS NOW NEGATIVE: \f8.3\ ', r_array = [a_ptr%r])
  elseif (ix_attrib == lord_pad2$) then
    call out_io (s_error$, r_name, 'GROUP ELEMENT: ' // lord%name, &
                                   'CONTROLS SUPER_LORD: ' // ele%name, &
                                   'AND LORD_PAD2 IS NOW NEGATIVE: \f8.3\ ', r_array = [a_ptr%r])
  endif
endif

! ele is a super_slave...

if (ele%slave_status == super_slave$) then
  if (ix_attrib /= l$) then
    call out_io (s_error$, r_name, &
                  'CONFUSED GROUP IS TRYING TO VARY SUPER_SLAVE ATTRIBUTE: ' // attribute_name(ele, ix_attrib))
    if (global_com%exit_on_error) call err_exit
    return
  endif

  do il = 1, ele%n_lord
    my_lord => pointer_to_lord(ele, il)
    if (my_lord%lord_status /= super_lord$) cycle
    call group_change_this (my_lord, ix_attrib, ctl, dir)
  enddo

  if (present(this_lord)) then
    call group_change_this (my_lord, ix_attrib, ctl, -dir)  ! Take out length change.
    call group_change_this (my_lord, this_pad, ctl, dir)    ! And change pad length instead.
  endif

endif

! ele is a multipass_slave...
! In the loop over all multipass_slaves, only modify the multipass_lord once

if (ele%slave_status == multipass_slave$) then
  my_lord => pointer_to_lord(ele, 1, ix_slave = ix_slave)
  if (ix_slave == 1) call group_change_this (my_lord, ix_attrib, ctl, 1)
endif

end subroutine group_change_this

end subroutine makeup_group_lord

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_multipass_slave (lat, slave)
!
! Subroutine to calcualte the attributes of multipass slave elements.
! This routine is not meant for general use.
!-

subroutine makeup_multipass_slave (lat, slave)

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

!

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
endif

! A slave's field_master = T irregardless of the lord's setting.
! This is to make attribute_bookkeeper compute the correct normalized field strength.

slave%value(E_tot_start$)    = slave_val(E_tot_start$)
slave%value(p0c_start$)      = slave_val(p0c_start$)
slave%value(e_tot$)          = slave_val(e_tot$)
slave%value(p0c$)            = slave_val(p0c$)
slave%value(delta_ref_time$) = slave_val(delta_ref_time$)
slave%value(ref_time_start$) = slave_val(ref_time_start$)
slave%value(n_ref_pass$)     = 0
if (attribute_index(slave, 'FIELD_MASTER') /= 0) slave%field_master = .true.

! A match element with match_end$: Restore initial Twiss parameters (which
! are calculated in twiss_propagate1).

if (lord%key == match$) then
  if (is_true(lord%value(match_end$))) then
    slave%value(beta_a0$)    = slave_val(beta_a0$)
    slave%value(beta_b0$)    = slave_val(beta_b0$)
    slave%value(alpha_a0$)   = slave_val(alpha_a0$)
    slave%value(alpha_b0$)   = slave_val(alpha_b0$)
    slave%value(eta_x0$)     = slave_val(eta_x0$)
    slave%value(eta_y0$)     = slave_val(eta_y0$)
    slave%value(etap_x0$)    = slave_val(etap_x0$)
    slave%value(etap_y0$)    = slave_val(etap_y0$)
    slave%value(c_11$:c_22$) = slave_val(c_11$:c_22$)
    slave%value(gamma_c$)    = slave_val(gamma_c$)
  endif

  if (is_true(lord%value(match_end_orbit$))) then
    slave%value(x0$)  = slave_val(x0$)
    slave%value(px0$) = slave_val(px0$)
    slave%value(y0$)  = slave_val(y0$)
    slave%value(py0$) = slave_val(py0$)
    slave%value(z0$)  = slave_val(z0$)
    slave%value(pz0$) = slave_val(pz0$)
  endif
endif

! Sbend field: The design bending strength is same for slave as lord. 
! So the error field must be adjusted so that total_field = design_field + err_field is the same. 
! Note: The lord's energy may not yet be set if bmad_parser is active. So only do calc if p0c is set.

if (slave%key == sbend$ .and. lord%value(p0c$) /= 0) then
  slave%value(b_field$) = lord%value(b_field$) * slave%value(p0c$) / lord%value(p0c$)
  slave%value(b_field_err$) = (lord%value(b_field$) + lord%value(b_field_err$)) - slave%value(b_field$)
  slave%value(g_err$) = (lord%value(g$) + lord%value(g_err$)) - slave%value(g$)
endif

! Multipoles. Note: p0c = 0 Can happen if not finished parsing lattice file.

if (associated (slave%a_pole) .and. slave%value(p0c$) /= 0 .and. .not. lord%field_master) then  
  slave%a_pole           = lord%a_pole * lord%value(p0c$) / slave%value(p0c$)
  slave%b_pole           = lord%b_pole * lord%value(p0c$) / slave%value(p0c$)
  slave%multipoles_on    = lord%multipoles_on
  slave%scale_multipoles = lord%scale_multipoles
endif

! Electric Multipoles

if (associated (slave%a_pole_elec)) then
  slave%a_pole_elec      = lord%a_pole_elec
  slave%b_pole_elec      = lord%b_pole_elec
  slave%multipoles_on    = lord%multipoles_on
endif

! RF wakes

call transfer_wake (lord%wake, slave%wake)

if (associated (slave%wake)) then
  do i = 1, size(lord%wake%lr)
    slave%wake%lr(i)%t_ref = lord%wake%lr(i)%t_ref - slave%ref_time
  enddo
endif

! methods

slave%taylor_map_includes_offsets = lord%taylor_map_includes_offsets
slave%is_on            = lord%is_on

! Handled by set_flags_for_changed_attribute

!! slave%aperture_at      = lord%aperture_at
!! slave%aperture_type    = lord%aperture_type
!! slave%mat6_calc_method = lord%mat6_calc_method
!! slave%tracking_method  = lord%tracking_method

end subroutine makeup_multipass_slave

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_super_slave (lat, slave)
!
! Subroutine to calcualte the attributes of superposition slave elements.
! This routine is not meant for general use.
!-
       
subroutine makeup_super_slave (lat, slave)

type (lat_struct), target :: lat
type (ele_struct) slave
type (ele_struct), pointer :: lord, slave0, lord1
type (ele_struct) :: sol_quad
type (branch_struct), pointer :: branch

integer i, j, ix, ix_lord, ix_order, ix_slave, n_major_lords, at

real(rp) tilt, k_x, k_y, x_kick, y_kick, ks, k1, coef
real(rp) x_o, y_o, x_p, y_p, s_slave, s_del, k2, k3, c, s
real(rp) sin_n, cos_n, a(0:n_pole_maxx), b(0:n_pole_maxx)
real(rp) knl(0:n_pole_maxx), t(0:n_pole_maxx), value(num_ele_attrib$)
real(rp) sum_1, sum_2, sum_3, sum_4, ks_sum, ks_xp_sum, ks_xo_sum
real(rp) ks_yp_sum, ks_yo_sum, l_slave, r_off(4), leng, offset
real(rp) t_1(4), t_2(4), T_end(4,4), mat4(4,4), mat4_inv(4,4), beta(4)
real(rp) T_tot(4,4), x_o_sol, x_p_sol, y_o_sol, y_p_sol

logical is_first, is_last, err_flag

character(20) :: r_name = 'makeup_super_slave'

! Super_slave:

branch => lat%branch(slave%ix_branch)
ix_slave = slave%ix_ele

call set_ele_status_stale (slave, attribute_group$)

if (slave%slave_status /= super_slave$) then
  call out_io(s_abort$, r_name, "ELEMENT IS NOT A SUPER SLAVE: " // slave%name)
  if (global_com%exit_on_error) call err_exit
  return
endif

slave%field_calc = refer_to_lords$

!-----------------------------------------------------------------------
! 1 super_lord for this super_slave: just transfer attributes except length

if (slave%n_lord == 1) then

  lord => pointer_to_lord (slave, 1, ix_slave = ix_order)

  is_first = (ix_order == 1)
  is_last  = (ix_order == lord%n_slave)

  ! If this is not the first slave: Transfer reference orbit from previous slave

  if (.not. is_first) then
    if (.not. all(slave%map_ref_orb_in%vec == branch%ele(ix_order-1)%map_ref_orb_out%vec)) then
      slave0 => pointer_to_slave(lord, ix_order-1)
      slave%map_ref_orb_in = slave0%map_ref_orb_out
      if (associated(slave%rad_int_cache)) slave%rad_int_cache%stale = .true. ! Forces recalc
    endif
  endif

  ! Find the offset from the longitudinal start of the lord to the start of the slave

  offset = 0 ! length of all slaves before this one
  do i = 1, ix_order-1
    slave0 => pointer_to_slave(lord, i)
    offset = offset + slave0%value(l$)
  enddo

  call makeup_super_slave1 (slave, lord, offset, branch%param, is_first, is_last, err_flag)

  return

endif

!-----------------------------------------------------------------------
! Multiple super_lords for this super_slave: 
! combine the lord elements.
                                         
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

value = 0
value(l$) = slave%value(l$)
value(E_tot_start$)    = slave%value(E_tot_start$)
value(p0c_start$)      = slave%value(p0c_start$)
value(E_tot$)          = slave%value(E_tot$)
value(p0c$)            = slave%value(p0c$)
value(delta_ref_time$) = slave%value(delta_ref_time$)
value(ref_time_start$) = slave%value(ref_time_start$)
value(fringe_at$)      = no_end$
value(fringe_type$)    = none$

slave%value(x1_limit$:y2_limit$) = 0

slave%aperture_at = no_end$
slave%is_on = .false.
s_slave = slave%s - value(l$)/2  ! center of slave

! A "major" element is something other than a pipe, monitor, etc.
! n_major_lords counts how many major lords there are.

n_major_lords = 0

! sum over all lords...

do j = 1, slave%n_lord

  lord => pointer_to_lord(slave, j, ix_slave = ix_order)
  if (lord%lord_status /= super_lord$) cycle

  is_first = (ix_order == 1)
  is_last  = (ix_order == lord%n_slave)

  ! Do some error checking.

  if (associated(lord%wake)) then
    call out_io (s_abort$, r_name, &
            'SUPERPOSITION OF ELEMENTS WITH WAKES NOT YET IMPLEMENTED!', &
            'SUPER_LORD: ' // lord%name)
    if (global_com%exit_on_error) call err_exit
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
      if (associated(slave%rad_int_cache)) slave%rad_int_cache%stale = .true. ! Forces recalc
    endif
  endif

  ! Choose the smallest ds_step of all the lords.

  if (value(ds_step$) == 0 .or. lord%value(ds_step$) < value(ds_step$)) &
                                        value(ds_step$) = lord%value(ds_step$)

  ! Methods...
  ! A "major" element is something other than a pipe.
  ! n_major_lords counts how many major lords there are.

  if (n_major_lords == 0) then
    slave%mat6_calc_method = lord%mat6_calc_method
    slave%tracking_method  = lord%tracking_method
    slave%taylor_map_includes_offsets = lord%taylor_map_includes_offsets
  endif

  if (ele_has(lord, 'FRINGE_TYPE')) then
    if (is_first .and. at_this_ele_end(entrance_end$, nint(lord%value(fringe_at$)))) then
      call set_fringe_on_off(value(fringe_at$), entrance_end$, on$)
      value(fringe_type$) = lord%value(fringe_type$)
    endif

    if (is_last .and. at_this_ele_end(exit_end$, nint(lord%value(fringe_at$)))) then
      call set_fringe_on_off(value(fringe_at$), exit_end$, on$)
      value(fringe_type$) = lord%value(fringe_type$)
    endif
  endif

  if (lord%key /= pipe$) then
    if (n_major_lords > 0) then
      if (slave%mat6_calc_method /= lord%mat6_calc_method) then
        call out_io(s_abort$, r_name, &
              'MAT6_CALC_METHOD DOES NOT AGREE FOR DIFFERENT SUPERPOSITION LORDS FOR SLAVE: ' // slave%name, &
              'Conflicting methods are: ' // trim(mat6_calc_method_name(lord%mat6_calc_method)) // ',  ' // & 
              mat6_calc_method_name(slave%mat6_calc_method))
        if (global_com%exit_on_error) call err_exit
      endif
      if (slave%tracking_method /= lord%tracking_method) then
        call out_io(s_abort$, r_name, &
             'TRACKING_METHOD DOES NOT AGREE FOR DIFFERENT SUPERPOSITION LORDS FOR SLAVE: ' // slave%name, &
             'Conflicting methods are: ' // trim(tracking_method_name(lord%tracking_method)) // ',  ' // & 
             tracking_method_name(slave%tracking_method))
        if (global_com%exit_on_error) call err_exit
      endif
      if (slave%taylor_map_includes_offsets .neqv. lord%taylor_map_includes_offsets) then
        call out_io(s_abort$, r_name, &
            'TAYLOR_MAP_INCLUDES_OFFSETS DOES NOT AGREE FOR DIFFERENT SUPERPOSITION LORDS FOR SLAVE: ' // slave%name)
        if (global_com%exit_on_error) call err_exit
      endif
      if ((is_first .or. is_last) .and. ele_has(lord, 'FRINGE_TYPE')) then
       if (value(fringe_type$) /= lord%value(fringe_type$)) then
         call out_io(s_abort$, r_name, &
            'FRINGE_TYPE DOES NOT AGREE FOR DIFFERENT SUPERPOSITION LORDS FOR SLAVE: ' // slave%name)
         if (global_com%exit_on_error) call err_exit
       endif
     endif
    endif

    n_major_lords = n_major_lords + 1
  endif

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
  elseif (lord%key == kicker$) then
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

    s_del = s_slave - (lord%s + lord%value(z_offset_tot$) - lord%value(l$)/2)
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

  ! bend_sol_quad

  case (bend_sol_quad$)
    call out_io (s_abort$, r_name, 'CODING NOT YET IMPLEMENTED FOR: ' // key_name(slave%key))
    if (global_com%exit_on_error) call err_exit

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
elseif (slave%key == kicker$) then
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

    call init_ele (sol_quad)
    sol_quad%key = sol_quad$
    sol_quad%value(ks$) = ks
    sol_quad%value(k1$) = k1
    sol_quad%value(l$)  = l_slave
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

! bend_sol_quad

case (bend_sol_quad$)
  call out_io (s_abort$, r_name, 'CODING NOT YET IMPLEMENTED FOR A: ' // key_name(slave%key))
  if (global_com%exit_on_error) call err_exit

end select

! If the slave has %field_master = T then we need to convert k1, etc values to field quantities.

8000 continue

! Coupler and aperture calc.

if (slave%key == lcavity$ .or. slave%key == rfcavity$) call compute_slave_coupler (slave)

if (slave%field_master) then
  slave%field_master = .false.   ! So attribute_bookkeeper will do the right thing.
  call attribute_bookkeeper (slave, branch%param)
  slave%field_master = .true.
else
  call attribute_bookkeeper (slave, branch%param)
endif

end subroutine makeup_super_slave

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine create_element_slice (sliced_ele, ele_in, l_slice, offset,
!                          param, include_upstream_end, include_downstream_end, err_flag, old_slice)
!
! Routine to create an element that represents a longitudinal slice of the original element.
!
! Note: This routine assumes that the following call has been made before hand:
!    call transfer_ele (ele_in, sliced_ele, .true.)
! transfer_ele only has to be done once as long as create_element_slice is repeatedly called 
! with the same ele_in.
!
! Note: To save tracking computation time, if ele_in has taylor, symp_lie_ptc, or symp_map 
! for tracking_method or mat6_calc_method, then this will be changed to symp_lie_bmad 
! for wigglers and bmad_standard for everything else.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele_in            -- Ele_struct: Original element to slice
!   l_slice           -- Real(rp): Length of the slice
!   offset            -- Real(rp): Offset of entrance end of sliced_ele from entrance end of ele_in.
!   param             -- Lat_param_struct: lattice paramters.
!   include_upstream_end   -- Logical: Sliced_ele contains the ele's entrance end?
!   include_downstream_end -- Logical: Sliced_ele contains the ele's exit end?
!   old_slice         -- Logical, optional: Previous slice. If present this saves computation
!                          time of the refernece energy and time at the start of the present slice.
!
! Output:
!   sliced_ele -- Ele_struct: Sliced_ele element with appropriate values set.
!   err_flag   -- Logical: Set True if there is an error. False otherwise.
!-

recursive subroutine create_element_slice (sliced_ele, ele_in, l_slice, offset, &
                             param, include_upstream_end, include_downstream_end, err_flag, old_slice)

type (ele_struct), target :: sliced_ele, ele_in, ele0
type (ele_struct), optional :: old_slice
type (ele_struct) :: ele2
type (lat_param_struct) param
type (coord_struct) time_ref_orb_out

real(rp) l_slice, offset, in_len, r
integer i
logical include_upstream_end, include_downstream_end, err_flag, err2_flag

character(24) :: r_name = 'create_element_slice'

! Init

err_flag = .true.
in_len = ele_in%value(l$)

sliced_ele%lord_status = not_a_lord$
sliced_ele%slave_status = slice_slave$
sliced_ele%ix_ele = -2  ! Indicate sliced ele is not an element in the lattice.
sliced_ele%value(l$) = l_slice
sliced_ele%s = ele_in%s - in_len + offset + sliced_ele%value(l$)

! A sad_mult with zero length is treated differently (edge fields ignored and
! multipoles are integrated multipoles instead of per unit length).
! Therefore, make sure a sad_mult has a finite length.

if (ele_in%key == sad_mult$ .and. l_slice == 0) sliced_ele%value(l$) = &
                                                            sign(bmad_com%significant_length/100, in_len)

! The sliced element is treated as a super_slave to the original element except
! if the original element is itself a super_slave in which case the sliced element 
! has the same lords as the original element.

if (ele_in%slave_status /= super_slave$) then
  sliced_ele%n_lord = 1
  sliced_ele%lord => ele_in
else
  nullify(sliced_ele%lord)
endif

! Err check. Remember: the element length may be negative

if (l_slice*in_len < 0 .or. abs(l_slice) > abs(in_len) + bmad_com%significant_length) then
  call out_io (s_fatal$, r_name, &
        'SLICE LENGTH IS OUT OF RANGE FOR ELEMENT: ' // ele_in%name, &
        'LENGTH: \2es12.3\ ', r_array = [l_slice, in_len])
  if (global_com%exit_on_error) call err_exit
  return
endif

if (ele_in%key == hybrid$) then
  call out_io (s_fatal$, r_name, &
        'CANNOT SLICE ELEMENT OF TYPE: ' // key_name(ele_in%key), &
        'CANNOT SLICE: ' // ele_in%name)
  if (global_com%exit_on_error) call err_exit
  return
endif

! Simple case where ele length is zero

if (in_len == 0) then
  err_flag = .false.
  return
endif

! Save values from old_slice if present

if (present(old_slice)) then
  ele0%value(p0c$)   = old_slice%value(p0c$)
  ele0%value(e_tot$) = old_slice%value(e_tot$)
  ele0%ref_time      = old_slice%ref_time
  time_ref_orb_out   = old_slice%time_ref_orb_out
endif

!

call makeup_super_slave1 (sliced_ele, ele_in, offset, param, include_upstream_end, include_downstream_end, err2_flag)
if (err2_flag) return

! For a sliced taylor element the %taylor%term components point to the lord components. 
! The routine deallocate_ele_pointers will only nullify and never deallocate these components of a slice_slave.

if (ele_in%key == taylor$) then
  do i = 1, 6
    sliced_ele%taylor(i)%term => ele_in%taylor(i)%term
  enddo
endif

! Use a speedier tracking method.

select case (sliced_ele%tracking_method)
case (taylor$, symp_map$, symp_lie_ptc$)
  select case (sliced_ele%key)
  case (wiggler$, undulator$); sliced_ele%tracking_method = symp_lie_bmad$
  case default;                sliced_ele%tracking_method = bmad_standard$
  end select
end select

select case (sliced_ele%mat6_calc_method)
case (taylor$, symp_map$, symp_lie_ptc$)
  select case (sliced_ele%key)
  case (wiggler$, undulator$); sliced_ele%mat6_calc_method = symp_lie_bmad$
  case default;                sliced_ele%mat6_calc_method = bmad_standard$
  end select
end select

sliced_ele%field_calc = refer_to_lords$

! Makeup_super_slave1 does not compute reference energy or time so need to do it here.

if (offset == 0) then
  ele0%value(p0c$)      = ele_in%value(p0c_start$)
  ele0%value(e_tot$)    = ele_in%value(e_tot_start$)
  ele0%ref_time         = ele_in%value(ref_time_start$)
  sliced_ele%time_ref_orb_in = ele_in%time_ref_orb_in

elseif (present(old_slice)) then
  sliced_ele%time_ref_orb_in = time_ref_orb_out

elseif (ele_has_constant_ds_dt_ref(ele_in)) then
  ele0%value(p0c$)      = ele_in%value(p0c$)
  ele0%value(e_tot$)    = ele_in%value(e_tot$)
  ele0%ref_time = ele_in%ref_time - ele_in%value(delta_ref_time$) * (ele_in%value(l$) - offset) / ele_in%value(l$)
  sliced_ele%time_ref_orb_in%vec = 0

else
  call transfer_ele (sliced_ele, ele2)
  call create_element_slice (ele2, ele_in, offset, 0.0_rp, param, .true., .false., err2_flag)
  if (err2_flag) return
  ele0%value(p0c$)      = ele2%value(p0c$)
  ele0%value(e_tot$)    = ele2%value(e_tot$)
  ele0%ref_time         = ele2%ref_time
  sliced_ele%time_ref_orb_in = ele2%time_ref_orb_out
endif

call ele_compute_ref_energy_and_time (ele0, sliced_ele, param, err2_flag)
if (err2_flag) return

if (.not. include_upstream_end) sliced_ele%time_ref_orb_in%location = inside$
if (.not. include_downstream_end) sliced_ele%time_ref_orb_out%location = inside$

! Round off can throw off the ending ref energy.
! This can be problematic when include_downstream_end = T since a particle traveling through and
! into the next element will have an energy mismatch and track1 flag this as an error.
! Solution is to just set the end ref energy to what it should be.

if (include_downstream_end) then
  sliced_ele%value(p0c$)      = ele_in%value(p0c$)
  sliced_ele%value(e_tot$)    = ele_in%value(e_tot$)
endif

err_flag = .false.

end subroutine create_element_slice

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine makeup_super_slave1 (slave, lord, offset, param, include_upstream_end, include_downstream_end)
!
! Routine to construct a super_slave from a super_lord when the slave has only one lord.
! Note: Reference energy and times are not computed in this routine.
!
! Modules needed:
!   use bmad
!
! Input:
!   slave  -- Ele_struct: Slave element.
!   lord   -- Ele_struct: Lord element.
!   offset -- Real(rp): offset of entrance end of slave from entrance end of the lord.
!   param  -- Lat_param_struct: lattice paramters.
!   include_upstream_end -- Logical: Slave contains the lord's entrance end?
!   include_downstream_end     -- Logical: Slave contains the lord's exit end?
!
! Output:
!   slave    -- Ele_struct: Slave element with appropriate values set.
!   err_flag -- Logical: Set true if there is an error. False otherwise.
!-

subroutine makeup_super_slave1 (slave, lord, offset, param, include_upstream_end, include_downstream_end, err_flag)

type (ele_struct), target :: slave, lord
type (ele_struct), pointer :: slave2
type (lat_param_struct) param

real(rp) offset, s_del, coef, lord_ang, slave_ang, angle
real(rp) value(num_ele_attrib$), cos_a, sin_a, dr(3)
real(rp) off(3), rot(3), cos_t, sin_t, m_trans(3,3)
real(rp) xp, yp, roll, r_roll, tilt, dx, dy
real(rp) w_inv(3,3), dl

integer i, ifr
logical include_upstream_end, include_downstream_end, err_flag, include_entrance, include_exit, err2_flag
logical has_fringe
character(24) :: r_name = 'makeup_super_slave1'

! Physically, the lord length cannot be less than the slave length.
! In case we are dealing with a non-physical situation, arbitrarily set coef = 1.

err_flag = .true.

if (abs(slave%value(l$)) >= abs(lord%value(l$))) then
  coef = 1
else
  coef = slave%value(l$) / lord%value(l$) 
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

value(l$)              = slave%value(l$)                ! do not change slave length, etc.
value(delta_ref_time$) = slave%value(delta_ref_time$)
value(ref_time_start$) = slave%value(ref_time_start$)
value(E_tot_start$)    = slave%value(E_tot_start$)
value(p0c_start$)      = slave%value(p0c_start$)
value(E_tot$)          = slave%value(E_tot$)
value(p0c$)            = slave%value(p0c$)
value(num_steps$)      = slave%value(num_steps$)

! Taylor element has a zero length map. Rule: The map gets applied at the entrance end.
! There is no reason why the entrance end was chosen over the exit end.
! Note: multipole, and ab_multipole have zero length maps but in this case multipole_ele_to_kt
! will scale the multipole strength proportional to the slave length.

slave%is_on = lord%is_on
has_fringe = .true.

if (lord%key == taylor$) then
  if (.not. include_entrance) slave%is_on = .false.
  has_fringe = .false.
elseif (lord%key == multipole$ .or. lord%key == ab_multipole$) then
  has_fringe = .false.
endif

! Ref energy shift for e_gun only happens at start of element.

if (lord%key == e_gun$ .and. .not. include_upstream_end) then
  value(e_tot_ref_init$) = lord%value(e_tot_start$)
  value(p0c_ref_init$)   = lord%value(p0c_start$)
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
  elseif (include_downstream_end) then
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

if (slave%key == rfcavity$) value(voltage$) = lord%value(voltage$) * coef

! s_del is the distance between lord and slave centers

if (has_orientation_attributes(slave)) then
  s_del = offset + slave%value(l$)/2 - lord%value(l$)/2

  if (slave%key == sbend$ .and. value(g$) /= 0) then

    roll = value(roll_tot$);      tilt = value(ref_tilt_tot$)
    off = [value(x_offset_tot$), value(y_offset_tot$), value(z_offset_tot$)]
    xp    = value(x_pitch_tot$);   yp    = value(y_pitch_tot$)

    value(ref_tilt$) = tilt

    if (any(off /= 0) .or. xp /= 0 .or. yp /= 0 .or. roll /= 0) then
      angle =  s_del * value(g$)
      cos_a = cos(angle); sin_a = sin(angle)
      dr = [2 * sin(angle/2)**2 / value(g$), 0.0_rp, sin_a / value(g$)]

      if (roll /= 0) then
        r_roll = value(rho$) * (1 - cos(lord%value(angle$)/2))
        off(1) = off(1) + r_roll*(cos(roll) - 1)
        off(2) = off(2) + r_roll*sin(roll)
      endif

      if (tilt == 0) then
        cos_t = 1; sin_t = 0
        off = [cos_a * off(1) + sin_a * off(3), off(2), -sin_a * off(1) + cos_a * off(3)]
        rot = [-cos_a * yp, xp, sin_a * yp + roll]
      else
        cos_t = cos(tilt);    sin_t = sin(tilt)
        m_trans(1,:) = [cos_a * cos_t**2 + sin_t**2, (cos_a - 1) * cos_t * sin_t, cos_t * sin_a]
        m_trans(2,:) = [(cos_a - 1) * cos_t * sin_t, cos_a * sin_t**2 + cos_t**2, sin_a * sin_t]
        m_trans(3,:) = [-cos_t * sin_a, -sin_a * sin_t, cos_a]
        rot = matmul(m_trans, [-yp, xp, roll])
        off = matmul(m_trans, off)
        dr = [cos_t * dr(1) + sin_t * dr(2), sin_t * dr(1) + cos_t * dr(2), dr(3)]
      endif

      if (any(rot /= 0)) then
        call axis_angle_to_w_mat (rot, norm2(rot), m_trans)
        off = off + matmul(m_trans, dr) - dr
      endif

      if (rot(3) /= 0) then
        r_roll = value(rho$) * (1 - cos(value(l$)*value(g$)/2))
        dx = r_roll * (1 - cos(rot(3)));  dy = -r_roll*sin(rot(3))
        off(1) = off(1) + dx * cos_t - dy * sin_t
        off(2) = off(2) + dx * sin_t + dy * cos_t
      endif

      value(x_offset$) = off(1)
      value(y_offset$) = off(2)
      value(z_offset$) = off(3)
      value(x_pitch$) =  rot(2)
      value(y_pitch$) = -rot(1)
      value(roll$)    =  rot(3)

    endif

  ! Not an sbend

  else
    value(tilt$)     = value(tilt_tot$)
    value(x_pitch$)  = value(x_pitch_tot$)
    value(y_pitch$)  = value(y_pitch_tot$)
    value(x_offset$) = value(x_offset_tot$) + s_del * value(x_pitch_tot$)
    value(y_offset$) = value(y_offset_tot$) + s_del * value(y_pitch_tot$)
    value(z_offset$) = value(z_offset_tot$)
  endif
endif

! Patch
! The rotation part of the patch is applied at the entrance end of the patch.
! Excluding the rotation, a patch is just a drift.

if (lord%key == patch$) then
  if ((include_upstream_end .and. lord%orientation == 1) .or. (include_downstream_end .and. lord%orientation == -1)) then
    call floor_angles_to_w_mat (lord%value(x_pitch$), lord%value(y_pitch$), lord%value(tilt$), w_mat_inv = w_inv)
    dl = lord%value(l$) - slave%value(l$)
    value(x_offset$)     = lord%value(x_offset$) - dl * w_inv(3,1)
    value(y_offset$)     = lord%value(y_offset$) - dl * w_inv(3,2)
    value(z_offset$)     = lord%value(z_offset$) - dl * w_inv(3,3)
    value(t_offset$)     = lord%value(t_offset$)
    value(e_tot_offset$) = lord%value(e_tot_offset$)
  else
    value(x_pitch$)      = 0
    value(y_pitch$)      = 0
    value(tilt$)         = 0
    value(x_offset$)     = 0
    value(y_offset$)     = 0
    value(z_offset$)     = slave%value(l$) ! L is set by create_element_slice
    value(t_offset$)     = 0
    value(e_tot_offset$) = 0
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
    call ele_compute_ref_energy_and_time (lord, slave, param, err2_flag)
  endif

endif

!

slave%value = value
slave%mat6_calc_method = lord%mat6_calc_method
slave%tracking_method  = lord%tracking_method
slave%taylor_map_includes_offsets = lord%taylor_map_includes_offsets

if (slave%tracking_method == bmad_standard$ .and. slave%key == em_field$) slave%tracking_method = runge_kutta$
if (slave%mat6_calc_method == bmad_standard$ .and. slave%key == em_field$) slave%mat6_calc_method = tracking$

! wiggler fields and electro-magnetic fields

if (slave%key == wiggler$ .or. slave%key == undulator$) slave%value(n_pole$) = lord%value(n_pole$) * coef

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
  slave%wake%lr%freq_in   = lord%wake%lr%freq_in
  slave%wake%lr%freq      = lord%wake%lr%freq
  slave%wake%lr%Q         = lord%wake%lr%Q
  slave%wake%lr%angle     = lord%wake%lr%angle
  slave%wake%lr%m         = lord%wake%lr%m
  slave%wake%lr%polarized = lord%wake%lr%polarized
  slave%wake%lr%r_over_q  = lord%wake%lr%r_over_q * coef
endif

!

if (slave%key == lcavity$ .or. slave%key == rfcavity$) call compute_slave_coupler (slave)

if (slave%key == lcavity$) then
  slave%value(coupler_at$) = no_end$
  slave%value(e_loss$) = lord%value(e_loss$) * coef
endif

slave%bookkeeping_state%attributes = stale$
call attribute_bookkeeper (slave, param)

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
! Subroutine makeup_control_slave (lat, slave)
!
! This routine is not meant for general use.
!-

subroutine makeup_control_slave (lat, slave)

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

character(40) a_name
character(*), parameter :: r_name = 'makeup_control_slave'
logical is_free

!

branch => slave%branch
n_attrib = 0
val_attrib = 0

call set_ele_status_stale (slave, attribute_group$)

l_stat = slave%lord_status
ix_slave = slave%ix_ele

on_an_offset_girder = .false.

do i = 1, slave%n_lord
  lord => pointer_to_lord(slave, i, control)

  if (lord%lord_status == multipass_lord$) cycle
  if (lord%key == group$) cycle

  if (lord%key == girder$ .and. has_orientation_attributes(slave)) then
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
    if (slave%key == sbend$ .and. vs(g$) /= 0) then
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

    on_an_offset_girder = .true.

    cycle
  endif

  if (lord%key /= overlay$) then
    call out_io (s_abort$, r_name, 'THE LORD IS NOT AN OVERLAY \i\ ', ix_slave)
    call type_ele (slave, .true., 0, .false., 0, .true.)
    if (global_com%exit_on_error) call err_exit
  endif     


  ! overlay lord

  iv = control%ix_attrib  
  select case (iv)

  case (x_limit$)
    call overlay_change_this(slave%value(x1_limit$), control)
    call overlay_change_this(slave%value(x2_limit$), control)
  case (y_limit$)
    call overlay_change_this(slave%value(y1_limit$), control)
    call overlay_change_this(slave%value(y2_limit$), control)
  case (aperture$)
    call overlay_change_this(slave%value(x1_limit$), control)
    call overlay_change_this(slave%value(x2_limit$), control)
    call overlay_change_this(slave%value(y1_limit$), control)
    call overlay_change_this(slave%value(y2_limit$), control)
  case default
    a_name = attribute_name(slave, iv)
    is_free = attribute_free (slave, a_name, .true., .true.)
    if (.not. is_free) then
      call out_io (s_abort$, r_name, 'OVERLAY LORD: ' // lord%name, &
           'IS TRYING TO VARY NON-FREE ATTRIBUTE: ' // trim(slave%name) // '[' // trim(a_name) // ']')
      err_flag = .true.
      return
    endif

    call pointer_to_indexed_attribute (slave, iv, .true., a_ptr, err_flag)
    if (err_flag) call err_exit
    call overlay_change_this(a_ptr%r, control)
  end select

enddo

! Transfer values from val_attrib to slave elements

do iv = 1, n_attrib

  a_ptr%r => ptr_attrib(iv)%r
  if (a_ptr%r == val_attrib(iv)) cycle
  a_ptr%r = val_attrib(iv)
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
  case (sbend$)
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

subroutine overlay_change_this (r_attrib, c)

type (ele_struct), pointer :: my_lord, my_slave
type (control_struct) c

real(rp), target :: r_attrib
real(rp) delta
integer iv
logical err_flag

character(100) err_str

!

call evaluate_expression_stack(c%stack, delta, err_flag, err_str, lord%control_var, .false.)
if (err_flag) then
  call out_io (s_error$, r_name, err_str, 'FOR SLAVE: ' // slave%name, 'OF LORD: ' // lord%name)
  return
endif

do iv = 1, n_attrib
  if (.not. associated(ptr_attrib(iv)%r, r_attrib)) cycle
  val_attrib(iv) = val_attrib(iv) + delta
  return
enddo

n_attrib = n_attrib + 1
ptr_attrib(n_attrib)%r => r_attrib
val_attrib(n_attrib) = val_attrib(n_attrib) + delta

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
type (surface_grid_struct), pointer :: grid
type (wall3d_section_struct), pointer :: sec

real(rp) angle, r_wall, dr_dtheta, x, y

integer i, j

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
  grid => ele%photon%surface%grid
  if (.not. allocated(grid%pt)) then
    call out_io (s_error$, r_name, 'ELEMENT APERTURE TYPE SET TO "SURFACE" BUT', &
                                   'NO SURFACE GRID IS DEFINED: ' // ele%name)
    return
  endif
  ele%value(x1_limit$) = -(grid%r0(1) + (lbound(grid%pt, 1) - 0.5) * grid%dr(1))
  ele%value(y1_limit$) = -(grid%r0(2) + (lbound(grid%pt, 2) - 0.5) * grid%dr(2))
  ele%value(x2_limit$) =  (grid%r0(1) + (ubound(grid%pt, 1) + 0.5) * grid%dr(1))
  ele%value(y2_limit$) =  (grid%r0(2) + (ubound(grid%pt, 2) + 0.5) * grid%dr(2))

end select

end subroutine aperture_bookkeeper

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

type (ele_struct), target :: ele
type (ele_struct), pointer :: slave

integer, target :: attrib
integer, pointer :: a_ptr
integer i

real(rp) dummy

logical, optional :: set_dependent

! This will set some generic flags

call set_flags_for_changed_real_attribute (ele, dummy, set_dependent)

!

a_ptr => attrib

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
    elseif (associated(a_ptr, ele%aperture_type)) then
      ele%aperture_type = a_ptr
    elseif (associated(a_ptr, ele%mat6_calc_method)) then
      ele%mat6_calc_method = a_ptr
    elseif (associated(a_ptr, ele%tracking_method)) then
      ele%tracking_method = a_ptr
    elseif (associated(a_ptr, ele%spin_tracking_method)) then
      ele%spin_tracking_method = a_ptr
    elseif (associated(a_ptr, ele%field_calc)) then
      ele%field_calc = a_ptr
    else
      exit
    endif
  enddo

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

type (ele_struct), target :: ele
type (ele_struct), pointer :: slave

integer i

real(rp) dummy

logical, target :: attrib
logical, pointer :: a_ptr
logical, optional :: set_dependent

! Call to set_flags_for_changed_real_attribute will set some generic flags

call set_flags_for_changed_real_attribute (ele, dummy, set_dependent)

a_ptr => attrib

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

type (ele_struct), target :: ele
type (branch_struct), pointer :: branch
type (grid_field_struct), pointer :: g_field
type (cylindrical_map_struct), pointer :: cl_map

real(rp), optional, target :: attrib
real(rp), pointer :: a_ptr
real(rp) v_mat(4,4), v_inv_mat(4,4), eta_vec(4), eta_xy_vec(4)
real(rp), target :: unknown_attrib

integer i

logical coupling_change, found, dep_set
logical, optional :: set_dependent

!-------------------
! For a particular elemement...

branch => ele%branch
dep_set = logic_option(.true., set_dependent)

! If a lord then set the control flag stale

if (ele%lord_status /= not_a_lord$) call set_ele_status_stale (ele, control_group$)

! Groups and overlays do not have any dependent attributes. 
! For all others set the attributes flag stale.

if (ele%key /= group$ .and. ele%key /= overlay$ .and. dep_set) then
  call set_ele_status_stale (ele, attribute_group$)
endif

! Transfer matrix calc needs to be flagged

if (ele%key /= overlay$ .and. ele%key /= group$ .and. &
    ele%lord_status /= multipass_lord$) then
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

! A length change involves changes in the floor position.

if (associated(a_ptr, ele%value(l$))) then
  if (ele%key /= overlay$ .and. ele%key /= group$) then
    call set_ele_status_stale (ele, s_and_floor_position_group$)
    call set_ele_status_stale (ele, floor_position_group$)
    call set_ele_status_stale (ele, ref_energy_group$)
  endif
  if (ele%value(p0c$) /= ele%value(p0c_start$)) call set_ele_status_stale (ele, ref_energy_group$)
endif

! E_tot and p0c can be varied in an init_ele or a multipass lord with n_ref_pass = 0.
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

  if (dep_set) then
    call convert_total_energy_to (ele%value(e_tot$), branch%param%particle, pc = ele%value(p0c$))
    if (ele%key == beginning_ele$) then
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

  if (dep_set) then
    call convert_pc_to (ele%value(p0c$), branch%param%particle, e_tot = ele%value(e_tot$))
    if (ele%key == beginning_ele$) then
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
  if (dep_set) then
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
  if (dep_set) then
    call convert_pc_to (ele%value(p0c_start$), branch%param%particle, e_tot = ele%value(e_tot_start$))
  endif
  return
endif

!

select case (ele%key)
case (beginning_ele$) 
  coupling_change = .false.

  if (associated(a_ptr, ele%a%beta) .or. associated(a_ptr, ele%a%alpha)) then
    if (dep_set) then
      if (ele%a%beta /= 0) ele%a%gamma = (1 + ele%a%alpha**2) / ele%a%beta
    endif
    return
  endif

  if (dep_set) then
    if (associated(a_ptr, ele%b%beta) .or. associated(a_ptr, ele%b%alpha)) then
      if (ele%b%beta /= 0) ele%b%gamma = (1 + ele%b%alpha**2) / ele%b%beta
    endif
    return
  endif

  if (dep_set) then
    if (associated(a_ptr, ele%c_mat(1,1)) .or. associated(a_ptr, ele%c_mat(1,2)) .or. & 
            associated(a_ptr, ele%c_mat(2,1)) .or. associated(a_ptr, ele%c_mat(2,2))) then
      ele%gamma_c = sqrt(1 - ele%c_mat(1,1)*ele%c_mat(2,2) + ele%c_mat(1,2)*ele%c_mat(2,1))
      coupling_change = .true.
    endif
  endif

  if (dep_set) then
    if (associated(a_ptr, ele%x%eta) .or. associated(a_ptr, ele%x%etap) .or. &
        associated(a_ptr, ele%y%eta) .or. associated(a_ptr, ele%y%etap) .or. &
        coupling_change) then 
      call make_v_mats (ele, v_mat, v_inv_mat)
      eta_xy_vec = [ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap]
      eta_vec = matmul (v_inv_mat, eta_xy_vec)
      ele%a%eta  = eta_vec(1)
      ele%a%etap = eta_vec(2)
      ele%b%eta  = eta_vec(3)
      ele%b%etap = eta_vec(4)
      return
    endif
  endif

  if (dep_set) then
    if (associated(a_ptr, ele%a%eta) .or. associated(a_ptr, ele%a%etap) .or. &
        associated(a_ptr, ele%b%eta) .or. associated(a_ptr, ele%b%etap)) then 
      call make_v_mats (ele, v_mat, v_inv_mat)
      eta_vec = [ele%a%eta, ele%a%etap, ele%b%eta, ele%b%etap]
      eta_xy_vec = matmul (v_mat, eta_vec)
      ele%x%eta  = eta_xy_vec(1)
      ele%x%etap = eta_xy_vec(2)
      ele%y%eta  = eta_xy_vec(3)
      ele%y%etap = eta_xy_vec(4)
      return
    endif
  endif

  if (associated(a_ptr, ele%floor%r(1)) .or. associated(a_ptr, ele%floor%r(2)) .or. &
      associated(a_ptr, ele%floor%r(3)) .or. associated(a_ptr, ele%floor%theta) .or. &
      associated(a_ptr, ele%floor%phi) .or. associated(a_ptr, ele%floor%psi)) then
    call set_ele_status_stale (ele, floor_position_group$)
    return
  endif

case (fork$, photon_fork$)

case (rfcavity$)
  if (associated(a_ptr, ele%value(voltage$)) .and. ele%value(l$) /= 0) ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)

case (lcavity$, e_gun$)

  if (associated(a_ptr, ele%value(gradient$)) .or. associated(a_ptr, ele%value(phi0$)) .or. &
      associated(a_ptr, ele%value(voltage$)) .or. associated(a_ptr, ele%value(rf_frequency$)) .or. &
      associated(a_ptr, ele%value(phi0_autoscale$))) then
    call set_ele_status_stale (ele, ref_energy_group$)
  endif

  if (associated(a_ptr, ele%value(voltage$)) .and. ele%value(l$) /= 0) ele%value(gradient$) = ele%value(voltage$) / ele%value(l$)
  if (associated(a_ptr, ele%value(voltage_err$)) .and. ele%value(l$) /= 0) ele%value(gradient_err$) = ele%value(voltage_err$) / ele%value(l$)

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

case (patch$)
  ! Any attribute change will shift the reference time.
  call set_ele_status_stale (ele, ref_energy_group$)
  call set_ele_status_stale (ele, floor_position_group$)

case (sbend$)
  if (associated(a_ptr, ele%value(angle$)) .or. associated(a_ptr, ele%value(g$)) .or. &
      associated(a_ptr, ele%value(rho$))) then
    call set_ele_status_stale (ele, floor_position_group$)
  endif

end select

end subroutine set_flags_for_changed_real_attribute

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine set_on_off (key, lat, switch, orb, use_ref_orb, ix_branch, saved_values, ix_attrib)
!
! Routine to turn on or off a set of elements (quadrupoles, rfcavities, etc.) in a lattice. 
! An element that is turned off acts like a drift.
! lat_make_mat6 will be called to remake lat%ele()%mat6.
!
! This routine can also be used to turn on/off a specific attribute in a set of elements using the
! ix_attrib and saved_values arguments. For example, to turn on/off the K2 component of all bends. 
! For real and integer attributes, switch = off$ means to set to zero.
! For real and integer attributes, switch = on$ and switch = restore_state$ have the same meaning 
! and the value that is set is optained from the saved_values argument (which must be present.)
!
! Modules needed:
!   use bmad
!
! Input:
!   key             -- integer: Class name of elements to be turned on or off. [quadrupole$, etc.]
!   lat             -- lat_struct: lattice structure holding the elements.
!   switch          -- integer: 
!                       on$            => Turn elements on. If saved_values argument is present, use this.
!                                          If not present (only for logical attributes), set to True.
!                       off$           => Turn elements off (but will not store the present state).
!                       off_and_save$  => Save on/off state and then turn elements off.
!                       save_state$    => Save present on/off state. No turning on or off is done.
!                       restore_state$ => Restore saved on/off state from saved_values argument.
!   orb(0:)         -- coord_struct, optional: Needed for lat_make_mat6
!   use_ref_orb     -- logical, optional: If present and true then use 
!                        ele%map_ref_orb for the reference orbit for
!                        calculating %mat6. Default is false.
!   ix_branch       -- integer, optional: If present then only set for 
!                        this lattice branch.
!   saved_values(:) -- real(rp), allocatable, optional: Saved values of the component.
!                       Must be present if needed (EG if switch = restore_state$, etc.).
!   ix_attrib       -- integer, optional: Attribute to turn on/off. Eg: k2$, multipoles_on$, etc.
!                       Default is is_on$.
!                   
!
! Output:
!   lat             -- lat_struct: Modified lattice.
!   saved_values(:) -- real(rp), allocatable, optional: Saved values of the component.
!                       Must be present if needed (EG if switch = off_and_save$, etc.).
!-

subroutine set_on_off (key, lat, switch, orb, use_ref_orb, ix_branch, saved_values, ix_attrib)

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct), optional :: orb(0:)
type (all_pointer_struct) a_ptr

real(rp), optional, allocatable :: saved_values(:)
real(rp) old_state

integer :: key, switch
integer, optional :: ix_branch, ix_attrib
integer i, ib, ixa, n_set, n_save_old

logical, optional :: use_ref_orb
logical err_flag

character(20) :: r_name = 'set_on_off'

!

ixa = integer_option(is_on$, ix_attrib)
n_save_old = 0
if (present(saved_values)) n_save_old = size(saved_values)

select case (switch)
case (restore_state$, off_and_save$, save_state$)
  if (present(saved_values)) then
    if (.not. allocated(saved_values)) allocate(saved_values(20))
  else
    call out_io (s_fatal$, r_name, 'SAVED_VALUES ARGUMENT NOT PRESENT.')
    if (global_com%exit_on_error) call err_exit
    return
  endif
end select

!  

do ib = 0, ubound(lat%branch, 1)

  if (present(ix_branch)) then
    if (ix_branch /= ib) cycle
  endif

  branch => lat%branch(ib)

  n_set = 0
  do i = 1, branch%n_ele_max

    ele => branch%ele(i)
    if (ele%key /= key) cycle

    n_set = n_set + 1

    call pointer_to_indexed_attribute (ele, ixa, .true., a_ptr, err_flag)
    if (err_flag) return

    if (present(saved_values)) then
      if (n_set > size(saved_values)) call re_allocate(saved_values, 2*n_set, .false.)
    endif

    old_state = value_of_all_ptr(a_ptr)

    select case (switch)
    case (on$)
      if (present(saved_values)) then
        call set_value_of(a_ptr, saved_values(n_set))
      else
        call set_value_of(a_ptr, true$)
      endif

    case (off$)
      call set_value_of(a_ptr, false$)  ! false$ = 0

    case (save_state$)
      saved_values(n_set) = value_of_all_ptr(a_ptr)
      cycle

    case (restore_state$)
      call set_value_of(a_ptr, saved_values(n_set))

    case (off_and_save$)
      saved_values(n_set) = value_of_all_ptr(a_ptr)
      call set_value_of(a_ptr, false$)

    case default
      call out_io (s_abort$, r_name, 'BAD SWITCH: \i\ ', switch)
      if (global_com%exit_on_error) call err_exit
    end select

    if (old_state == value_of_all_ptr(a_ptr)) cycle

    call set_flags_for_changed_attribute (ele, a_ptr)

    if (logic_option (.false., use_ref_orb)) then
      call make_mat6(ele, branch%param, ele%map_ref_orb_in)
    else
      call lat_make_mat6(lat, i, orb, ib)
    endif

  enddo

  if (key == rfcavity$) then   ! Reset 1-turn maps
    branch%param%t1_with_rf = 0  
    branch%param%t1_no_rf = 0    
  endif

enddo

! Error check

if (n_set == 0) return ! No setting done so no error checking needed

if ((switch == on$ .and. (associated(a_ptr%r) .or. associated(a_ptr%i))) .or. &
    (switch == on$ .and. associated(a_ptr%l) .and. present(saved_values)) .or. &
     switch == restore_state$) then
  if (n_save_old < n_set) then
    call out_io (s_fatal$, r_name, 'SAVED_VALUES IS BEING USED BEFORE BEING SET!')
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

!-------------------------------------------------------------------------
contains

subroutine set_value_of (a_ptr, val)

type (all_pointer_struct) a_ptr
real(rp) val

!

if (associated(a_ptr%r)) then
  a_ptr%r = val
elseif (associated(a_ptr%i)) then
  a_ptr%i = nint(val)
else
  a_ptr%l = is_true(val)
endif

end subroutine set_value_of

end subroutine set_on_off

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
  return
end select

if (ele%key == sad_mult$) then
  ele0 => ele
  if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) ele0 => pointer_to_lord(ele0, 1)
  ele%value(check_sum$) = 0
  do i = 0, ubound(ele%a_pole, 1)
    ele%value(check_sum$) = ele%value(check_sum$) + fraction(ele0%a_pole(i)) + fraction(ele0%b_pole(i))
    ele%value(check_sum$) = ele%value(check_sum$) + (exponent(ele0%a_pole(i)) + exponent(fraction(ele0%b_pole(i)))) / 10
  enddo
else
  ele%value(check_sum$) = ele%tracking_method
  if (ele%is_on) ele%value(check_sum$) = ele%value(check_sum$) + 1
endif

dv = abs(ele%value - ele%old_value)
dv(x1_limit$:y2_limit$) = 0  ! Limit changes do not need bookkeeping
dv(custom_attribute1$:custom_attribute_max$) = 0
if (present(dval)) dval = dv

if (all(dv == 0) .and. ele%key /= capillary$) then
  ele%bookkeeping_state%attributes = ok$
else
  ele%bookkeeping_state%attributes = stale$
endif

end subroutine attributes_need_bookkeeping

end module
