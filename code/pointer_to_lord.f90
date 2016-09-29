!+
! Function pointer_to_lord (slave, ix_lord, control, ix_slave, field_overlap_ptr) result (lord_ptr)
!
! Function to point to a lord of a slave.
!
! If field_overlap_ptr = False (default), the range for ix_lord is:
!   1 to slave%n_lord                                 for "regular" lords.
!   slave%n_lord+1 to slave%n_lord+slave%n_lord_field for field overlap lords.
!
! If field_overlap_ptr = True, only the field overlap lords may be accessed and the range for ix_lord is:
!   1 to slave%n_lord_field  
!
! Also see:
!   pointer_to_slave
!   pointer_to_ele
!
! Modules Needed:
!   use bmad_utils_mod
!
! Input:
!   slave              -- Ele_struct: Slave element.
!   ix_lord            -- Integer: Index of the lord.
!   field_overlap_ptr  -- logical, optional: Slave pointed to restricted to be a field overlap slave?
!                           Default is False.
!
! Output:
!   lord_ptr   -- Ele_struct, pointer: Pointer to the lord.
!                   Nullified if there is an error.
!   control    -- control_struct, pointer, optional: Pointer to control info for this lord/slave relationship.
!                   Nullified if there is an error.
!   ix_slave   -- Integer, optional: Index back to the slave. That is, 
!                   pointer_to_slave(lord_ptr, ix_slave) will point back to slave. 
!                   Set to -1 is there is an error or the slave is a slice_slave.
!-

function pointer_to_lord (slave, ix_lord, control, ix_slave, field_overlap_ptr) result (lord_ptr)

use basic_bmad_interface, except_dummy => pointer_to_lord

implicit none

type (ele_struct), target :: slave
type (control_struct), pointer, optional :: control
type (control_struct), pointer :: ctl
type (ele_struct), pointer :: lord_ptr
type (lat_struct), pointer :: lat

integer, optional :: ix_slave
integer i, ix_lord, icon, ixl

logical, optional :: field_overlap_ptr
character(*), parameter :: r_name = 'pointer_to_lord'

! Case where there is no lord

ixl = ix_lord
if (logic_option(.false., field_overlap_ptr)) ixl = ixl + slave%n_lord

if (ixl > slave%n_lord+slave%n_lord_field .or. ix_lord < 1) then
  nullify(lord_ptr)
  if (present(control)) nullify(control)
  if (present(ix_slave)) ix_slave = -1
  return
endif

! If a slice_ele is a slave of a super_slave, return a lord of the super_slave.

lat => slave%branch%lat

if (slave%slave_status == slice_slave$) then
  if (present(control)) nullify(control)
  if (present(ix_slave)) ix_slave = -1

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

! There must be a corresponding ix_slave value such that
!   pointer_to_slave(lord_ptr, ix_slave) => slave 

if (present(ix_slave)) then

  do i = 1, lord_ptr%n_slave
    if (associated (pointer_to_slave(lord_ptr, i), slave)) then
      ix_slave = i
      return
    endif
  enddo

  ! If ix_slave not found then this is an error

  call out_io (s_fatal$, r_name, 'CANNOT FIND SLAVE INDEX FOR LORD!')
  if (global_com%exit_on_error) call err_exit   

endif

end function pointer_to_lord

