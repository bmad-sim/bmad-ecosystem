!+
! Subroutine apply_rampers_to_slave (slave, err_flag)
!
! Routine to apply the ramper elements of a lattice to a particular element.
! Also see: apply_all_rampers.
!
! Input:
!   slave       -- ele_struct: Element to apply ramper elements to.
!
! Output:
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!-

subroutine apply_rampers_to_slave (slave, err_flag)

use bmad, except_dummy => apply_rampers_to_slave

implicit none

type (ele_struct), target :: slave
type (ele_struct), pointer :: ramper
type (lat_struct), pointer :: lat
type (control_ramp1_struct), pointer :: r1

integer iv, key, ix, ir, n, is, ix_con
logical err_flag, ok, found

character(40) name
character(*), parameter :: r_name = 'apply_rampers_to_slave'

! Init

err_flag = .false.
lat => slave%branch%lat

! Calculate slave values

do ir = 1, slave%n_lord_ramper
  ramper => pointer_to_lord(slave, ir, lord_type = ramper_lord$, ix_control = ix_con)
  if (.not. ramper%is_on) cycle
  call this_slave_bookkeeper(ramper, slave, ramper%control%ramp(ix_con))
enddo

call attribute_bookkeeper(slave, .true.)

!-------------------------------------------------------------------
contains


subroutine this_slave_bookkeeper (this_ramp, slave, r1)

type (ele_struct), target :: this_ramp, slave
type (ele_struct), pointer :: slave2
type (control_ramp1_struct) r1
type (all_pointer_struct) a_ptr

real(rp) value
logical err_flag

!

value = ramper_value(this_ramp, r1, err_flag)

if (r1%is_controller) then
  slave2 => pointer_to_ele(lat, r1%slave)
else
  slave2 => slave  
endif

call pointer_to_attribute (slave2, r1%attribute, .true., a_ptr, err_flag, .false.)
if (err_flag .or. .not. associated(a_ptr%r)) then
  return
endif

a_ptr%r = value
call set_flags_for_changed_attribute (slave2, a_ptr%r, .true.)
if (r1%is_controller) call control_bookkeeper(lat, slave2)

end subroutine this_slave_bookkeeper

end subroutine
