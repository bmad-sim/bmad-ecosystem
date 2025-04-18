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

recursive subroutine apply_rampers_to_slave (slave, err_flag)

use bmad_routine_interface, except_dummy => apply_rampers_to_slave

implicit none

type (ele_struct), target :: slave
type (ele_struct), pointer :: ramper, lord
type (lat_struct), pointer :: lat
type (control_ramp1_struct), pointer :: r1

integer iv, ie, key, ix, ir, n, is, ix_con
logical err_flag, ok, found

character(40) name
character(*), parameter :: r_name = 'apply_rampers_to_slave'

! Init

err_flag = .false.
lat => slave%branch%lat

! Apply to lords first

do ie = 1, slave%n_lord
  lord => pointer_to_lord(slave, ie)
  call apply_rampers_to_slave(lord, err_flag)
enddo

! Calculate slave values

do ir = 1, slave%n_lord_ramper
  ramper => pointer_to_lord(slave, ir, lord_type = ramper_lord$, ix_control = ix_con)
  if (.not. ramper%is_on) cycle
  call this_slave_bookkeeper(ramper, slave, ramper%control%ramp(ix_con), slave%control%ramper_lord(ir)%attrib_ptr)
enddo

call attribute_bookkeeper(slave, .true.)

!-------------------------------------------------------------------
contains

subroutine this_slave_bookkeeper (this_ramp, slave, r1, attrib_ptr)

type (ele_struct), target :: this_ramp, slave
type (control_ramp1_struct) r1

real(rp) value, attrib_ptr
logical err_flag, has_wild

!

value = ramper_value(this_ramp, r1, err_flag)
attrib_ptr = value
call set_flags_for_changed_attribute (slave, attrib_ptr, .true.)
if (r1%is_controller) call control_bookkeeper(lat, slave)

end subroutine this_slave_bookkeeper

end subroutine
