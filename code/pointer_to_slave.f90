!+
! Function pointer_to_slave (lord, ix_slave, control, field_overlap_ptr) result (slave_ptr)
!
! Function to point to a slave of a lord.
!
! Note: Ramper lords do not have any associated slaves (slaves are assigned dynamically at run time).
!
! If field_overlap_ptr = False (default), the range for ix_slave is:
!   1 to lord%n_slave                                 for "regular" slaves.
!   lord%n_slave+1 to lord%n_slave+lord%n_slave_field for field overlap slaves.
!
! If field_overlap_ptr = True, only the field overlap slaves may be accessed and the range for ix_slave is:
!   1 to lord%n_slave_field  
!
! Also see:
!   pointer_to_lord
!   pointer_to_ele
!   num_lords
!
! Input:
!   lord               -- Ele_struct: Lord element
!   ix_slave           -- Integer: Index of the slave. 
!   field_overlap_ptr  -- logical, optional: Slave pointed to restricted to be a field overlap slave?
!                           Default is False.
!
! Output:
!   slave_ptr  -- Ele_struct, pointer: Pointer to the slave.
!                   Nullified if there is an error.
!   control    -- control_struct, pointer, optional: Pointer to control info for this lord/slave relationship.
!                   Nullified if there is an error.
!-

function pointer_to_slave (lord, ix_slave, control, field_overlap_ptr) result (slave_ptr)

use equal_mod, except_dummy => pointer_to_slave

implicit none

type (ele_struct), target :: lord
type (control_struct), pointer, optional :: control
type (ele_struct), pointer :: slave_ptr
type (control_struct), pointer :: con
type (lat_struct), pointer :: lat

integer ix_slave, icon, ixs
logical, optional :: field_overlap_ptr

!

ixs = ix_slave
if (logic_option(.false., field_overlap_ptr)) ixs = ixs + lord%n_slave

if (ixs > lord%n_slave+lord%n_slave_field .or. ix_slave < 1) then
  nullify(slave_ptr)
  if (present(control)) nullify(control)
  return
endif

lat => lord%branch%lat
icon = lord%ix1_slave + ixs - 1
con => lat%control(icon)
slave_ptr => lat%branch(con%slave%ix_branch)%ele(con%slave%ix_ele)
if (present(control)) control => con

end function pointer_to_slave

