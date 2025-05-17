!+
! Function pointer_to_super_lord (slave, control, ix_slave_back, ix_control, ix_ic) result (lord_ptr)
!
! Function to point to the super_lord of a super_slave or the lord of a slick_slave.
!
! If slave is a slice_slave with a control chain:
!   super_lord -> super_slave -> slice_slave
! Then pointer_to_super_lord will point to the super_lord and not the super_slave. 
! Access to the super_slave is via slave%lord.
!
! Also see:
!   pointer_to_lord
!   pointer_to_slave
!   pointer_to_ele
!   num_lords
!
! Input:
!   slave            -- ele_struct: Slave element.
!
! Output:
!   lord_ptr        -- ele_struct, pointer: Pointer to the lord.
!                        Nullified if there is an error.
!   control         -- control_struct, pointer, optional: Pointer to control info for this lord/slave relationship.
!                        Nullified if there is an error.
!   ix_slave_back   -- integer, optional: Index back to the slave. That is, pointer_to_slave(lord_ptr, ix_slave_back) 
!                        will point back to slave. Set to -1 if there is an error or the slave is a slice_slave.
!   ix_control      -- integer, optional: Index in lat%control(:) array the control argument is at.
!                        For ramper lord elements, ix_control is index for the lord%control%ramper(:) array.
!   ix_ic           -- integer, optional: Index of the lat%ic(:) element associated with the control argument.
!-

function pointer_to_super_lord (slave, control, ix_slave_back, ix_control, ix_ic) result (lord_ptr)

use bmad_routine_interface, except_dummy => pointer_to_super_lord

implicit none

type (ele_struct), target :: slave
type (control_struct), pointer, optional :: control
type (control_struct), pointer :: ctl
type (ele_struct), pointer :: lord_ptr, ele_ptr
type (lat_struct), pointer :: lat

integer, optional :: ix_slave_back, ix_control, ix_ic
integer icon, ix

character(*), parameter :: r_name = 'pointer_to_super_lord'

! Case where there is no lord

if (present(control)) nullify(control)
if (present(ix_slave_back)) ix_slave_back = -1
if (present(ix_control)) ix_control = -1
if (present(ix_ic)) ix_ic = -1
nullify(lord_ptr)

lat => slave%branch%lat

! If a slice_ele is a slave of a super_slave, return a lord of the super_slave.

if (slave%slave_status == slice_slave$) then
  ele_ptr => slave%lord
  if (ele_ptr%slave_status /= super_slave$) then
    lord_ptr => ele_ptr
    return
  endif
else
  ele_ptr => slave
endif

!

if (ele_ptr%slave_status /= super_slave$) then
  call out_io (s_fatal$, r_name, 'Routine is being called with lattice element ' // ele_full_name(slave), &
                                 'that is not a super_slave nor a slice_slave!')
  if (global_com%exit_on_error) call err_exit
  return
endif

do ix = 1, ele_ptr%n_lord
  lord_ptr => pointer_to_lord(ele_ptr, ix)
  if (lord_ptr%key /= pipe$) exit
  if (ix == ele_ptr%n_lord) exit
enddo

if (lord_ptr%lord_status /= super_lord$) then  ! Can happen with group or overlay lords present.
  ix = ix - 1
  lord_ptr => pointer_to_lord(ele_ptr, ix)
endif  


icon = lat%ic(slave%ic1_lord + ix - 1)
ctl => lat%control(icon)
lord_ptr => lat%branch(ctl%lord%ix_branch)%ele(ctl%lord%ix_ele)

if (present(control)) control => ctl
if (present(ix_control)) ix_control = icon
if (present(ix_ic)) ix_ic = slave%ic1_lord + ix - 1
if (present(ix_slave_back)) ix_slave_back = icon - lord_ptr%ix1_slave + 1

end function pointer_to_super_lord

