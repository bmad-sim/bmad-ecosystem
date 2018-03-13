!+
! Subroutine set_ele_status_stale (ele, status_group, set_slaves)
!
! Routine to set a status flags to stale in an element and the corresponding 
! ones for any slaves the element has.
!
! Also the branch%param structure of the branch the element is in is set.
!
! For example: status_group = ref_energy_group$ sets stale:
!   ele%bookkeeping_state%ref_energy 
!   ele%bookkeeping_state%floor_position
!   ele%bookkeeping_state%mat6
! See the code for more details.
! 
! Output:
!   ele           -- ele_struct: Element.
!     %bookkeeping_state   -- Status block to set.
!   status_group  -- Integer: Which flag groups to set. Possibilities are:
!                      attribute_group$, control_group$, floor_position_group$, s_position_group$, 
!                      s_and_floor_position_group$, ref_energy_group$, or mat6_group$, all_groups$
!   set_slaves    -- Logical, optional: If present and False then do not set
!                      the status for any slaves. Default is True.
!-

recursive subroutine set_ele_status_stale (ele, status_group, set_slaves)

use bmad_routine_interface, dummy => set_ele_status_stale

implicit none

type (bookkeeping_state_struct), pointer :: state
type (ele_struct), target :: ele
type (ele_struct), pointer :: slave, ele2
integer status_group, i
logical, optional :: set_slaves

! Only set overall lattice status flags if the element is part of a lattice.

if (ele%ix_ele > -1 .and. associated(ele%branch)) then
  ! If a lord
  if (ele%ix_branch == 0 .and. ele%ix_ele > ele%branch%n_ele_track) then
    state => ele%branch%lat%lord_state
  else
    state => ele%branch%param%bookkeeping_state
  endif
else
  nullify(state)
endif

!

select case (status_group)

case (attribute_group$)
  call set_attributes_status
  call set_mat6_status
  call set_ptc_status

case (control_group$)
  call set_control_status

case (floor_position_group$)
  call set_floor_position_status
  call set_mat6_status
  call set_ptc_status

case (s_and_floor_position_group$)
  call set_s_position_status
  call set_floor_position_status
  call set_mat6_status
  call set_ptc_status

case (s_position_group$)
  call set_s_position_status
  call set_mat6_status
  call set_ptc_status

case (ref_energy_group$)
  call set_ref_energy_status
  call set_mat6_status
  call set_attributes_status ! EG: k1 <--> b1_gradient calc needed 
  call set_ptc_status

case (mat6_group$)
  call set_mat6_status

case (rad_int_group$)
  call set_rad_int_status

case (all_groups$)
  call set_attributes_status
  call set_control_status
  call set_ref_energy_status
  call set_floor_position_status
  call set_s_position_status
  call set_mat6_status
  call set_rad_int_status
  call set_ptc_status

case default
   if (global_com%exit_on_error) call err_exit   ! Should not be here

end select

! Set slave

if (logic_option(.true., set_slaves)) then
  do i = 1, ele%n_slave
    slave => pointer_to_slave (ele, i)
    call set_ele_status_stale (slave, status_group)
  enddo
endif

!----------------------------------------------------------------------------
contains

subroutine set_attributes_status
  if (ele%key == overlay$) return
  if (ele%key == group$) return
  ele%bookkeeping_state%attributes = stale$
  if (associated(state)) state%attributes = stale$
end subroutine set_attributes_status

!----------------------------------------------------------------------------
! contains

subroutine set_control_status
  if (ele%lord_status == not_a_lord$ .and. ele%n_lord == 0) return
  ele%bookkeeping_state%control = stale$
  if (associated(state)) state%control = stale$
end subroutine set_control_status

!----------------------------------------------------------------------------
! contains

subroutine set_floor_position_status
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%floor_position = stale$
  if (associated(state)) state%floor_position = stale$
  ! If there is a flexible patch that depends upon the position of this element then
  ! set the patch's floor_position status stale as well as all elements in between.
  if (ele%ix_ele > 1 .and. associated(ele%branch)) then
    if (ele%ix_ele > ele%branch%n_ele_track) return  ! Do not need to check lord elements.
    ele2 => pointer_to_next_ele(ele, -1)
    do
      if (ele2%ix_ele == 0) return
      if (ele2%key /= patch$) then
        if (ele2%value(l$) /= 0) return
        ele2 => pointer_to_next_ele(ele2, -1)
        cycle
      endif
      if (is_false(ele2%value(flexible$))) return
      ele%branch%ele(ele2%ix_ele:ele%ix_ele-1)%bookkeeping_state%floor_position = stale$
      return
    enddo
  endif
end subroutine set_floor_position_status

!----------------------------------------------------------------------------
! contains

subroutine set_s_position_status
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%s_position = stale$
  if (associated(state)) state%s_position = stale$
end subroutine set_s_position_status

!----------------------------------------------------------------------------
! contains

subroutine set_ref_energy_status
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%ref_energy = stale$
  if (associated(state)) state%ref_energy = stale$
end subroutine set_ref_energy_status

!----------------------------------------------------------------------------
! contains

subroutine set_rad_int_status
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%rad_int = stale$
  if (associated(state)) state%rad_int = stale$
end subroutine set_rad_int_status

!----------------------------------------------------------------------------
! contains

subroutine set_ptc_status
  if (ele%key == overlay$ .or. ele%key == group$) return
  ele%bookkeeping_state%ptc = stale$
  if (associated(state)) state%ptc = stale$
end subroutine set_ptc_status

!----------------------------------------------------------------------------
! contains

! Ignore if the element does not have an associated linear transfer map
! Also the branch status does not get set since the transfer map calc
! must always check a branch due to possible reference orbit shifts.

subroutine set_mat6_status
  if (ele%key == overlay$) return
  if (ele%key == group$) return
  if (ele%lord_status == multipass_lord$) return
  ele%bookkeeping_state%mat6 = stale$
end subroutine set_mat6_status

end subroutine set_ele_status_stale 

