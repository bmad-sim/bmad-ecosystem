!+
! Function pointer_to_lord (slave, ix_lord, control, ix_slave_back, lord_type, ix_control, ix_ic) result (lord_ptr)
!
! Function to point to a lord of a slave.
!
! If lord_type = all$ (the default), the range for ix_lord is:
!   1 to slave%n_lord                                 for "regular" lords.
!   slave%n_lord+1 to slave%n_lord+slave%n_lord_field for field overlap lords.
!
! If lord_type = field_lord$, only the field overlap lords may be accessed and the range for ix_lord is:
!   1 to slave%n_lord_field  
!
! If slave arg is a slice_slave with control chain:
!   super_lord -> super_slave -> slice_slave
! Then pointer_to_lord will point to the super_lord and not the super_slave. 
! Access to the super_slave is via slave%lord.
!
! Also see:
!   pointer_to_slave
!   pointer_to_ele
!   num_lords
!
! Input:
!   slave            -- ele_struct: Slave element.
!   ix_lord          -- integer: Index of the lord.
!   lord_type        -- integer, optional: See above.
!
!
! Output:
!   lord_ptr        -- ele_struct, pointer: Pointer to the lord.
!                        Nullified if there is an error.
!   control         -- control_struct, pointer, optional: Pointer to control info for this lord/slave relationship.
!                        Nullified if there is an error.
!   ix_slave_back   -- integer, optional: Index back to the slave. That is, pointer_to_slave(lord_ptr, ix_slave_back) 
!                        will point back to slave. Set to -1 if there is an error or the slave is a slice_slave.
!   ix_control      -- integer, optional: Index in lat%control(:) array the control argument is at.
!   ix_ic           -- integer, optional: Index of the lat%ic(:) element associated with the control argument.
!-

function pointer_to_lord (slave, ix_lord, control, ix_slave_back, lord_type, ix_control, ix_ic) result (lord_ptr)

use equal_mod, except_dummy => pointer_to_lord

implicit none

type (ele_struct), target :: slave
type (control_struct), pointer, optional :: control
type (control_struct), pointer :: ctl
type (ele_struct), pointer :: lord_ptr
type (lat_struct), pointer :: lat

integer, optional :: ix_slave_back, lord_type, ix_control, ix_ic
integer i, ix_lord, icon, ixl

character(*), parameter :: r_name = 'pointer_to_lord'

! Case where there is no lord

if (present(control)) nullify(control)
if (present(ix_slave_back)) ix_slave_back = -1
if (present(ix_control)) ix_control = -1
if (present(ix_ic)) ix_ic = -1

ixl = ix_lord
if (integer_option(all$, lord_type) == field_lord$) ixl = ixl + slave%n_lord

if (ixl > slave%n_lord+slave%n_lord_field .or. ix_lord < 1) then
  nullify(lord_ptr)
  return
endif

! If a slice_ele is a slave of a super_slave, return a lord of the super_slave.

lat => slave%branch%lat

if (slave%slave_status == slice_slave$) then
  lord_ptr => slave%lord
  if (lord_ptr%slave_status == super_slave$) then
    icon = lat%ic(slave%ic1_lord + ixl - 1)
    ctl => lat%control(icon)
    lord_ptr => lat%branch(ctl%lord%ix_branch)%ele(ctl%lord%ix_ele)
  endif

  return
endif

! Point to the lord

icon = lat%ic(slave%ic1_lord + ixl - 1)
ctl => lat%control(icon)
lord_ptr => lat%branch(ctl%lord%ix_branch)%ele(ctl%lord%ix_ele)

if (present(control)) control => ctl
if (present(ix_control)) ix_control = icon
if (present(ix_ic)) ix_ic = slave%ic1_lord + ixl - 1
if (present(ix_slave_back)) ix_slave_back = icon - lord_ptr%ix1_slave + 1

end function pointer_to_lord

