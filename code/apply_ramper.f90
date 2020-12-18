!+
! Subroutine apply_ramper (ele, ramper, err_flag)
!
! Input:
!   ele         -- ele_struct: Element to apply ramper to.
!   ramper      -- ele_struct: Ramper element.
!
! Output:
!   start_orb   -- coord_struct: Modified starting position for track1 to use.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!-

subroutine apply_ramper (ele, ramper, err_flag)

use bmad, except_dummy => apply_ramper

implicit none

type (ele_struct) :: ele, ramper

integer iv
logical err_flag

character(*), parameter :: r_name = 'apply_ramper'

!

err_flag = .false.

do iv = 1, size(ramper%control%ramp)

enddo


end subroutine
