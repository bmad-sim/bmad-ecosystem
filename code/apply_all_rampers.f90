!+
! Subroutine apply_all_rampers (lat, err_flag)
!
! Routine to apply all the rampers in a lattice.
! Also see: apply_rampers_to_slave.
!
! Input:
!   lat         -- lat_struct: Lattice.
!
! Output:
!   lat         -- lat_struct: Lattice with rampers applied.
!   err_flag    -- logical: Set True if there is an error. False otherwise.
!-

subroutine apply_all_rampers (lat, err_flag)

use bmad, except_dummy => apply_all_rampers

implicit none

type (lat_struct), target :: lat
type (ele_pointer_struct), allocatable :: rampers(:)
type (branch_struct), pointer :: branch

integer ie, ib, n_ramp
logical err_flag

!

call lat_ele_locator('RAMPER::*', lat, rampers, n_ramp)
if (n_ramp == 0) return

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 0, branch%n_ele_max
    call apply_rampers_to_slave(branch%ele(ie), rampers(1:n_ramp), err_flag)
    if (err_flag) return
  enddo
enddo

call set_flags_for_changed_attribute (lat)

end subroutine
