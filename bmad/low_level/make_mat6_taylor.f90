!+
! Subroutine make_mat6_taylor (ele, param, start_orb, end_orb, err_flag)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Input:
!   ele      -- Ele_struct: Element to track through.
!   param     -- lat_param_struct: Parameters are needed for some elements.
!   start_orb -- coord_struct: Starting coords.
!
! Output:
!   ele       -- Ele_struct: Element with transfer matrix.
!     %vec0     -- 0th order map component
!     %mat6     -- 6x6 transfer matrix.
!   end_orb   -- Coord_struct: Coordinates at the end of element.
!   err       -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine make_mat6_taylor (ele, param, start_orb, end_orb, err_flag)

use bmad_interface, dummy => make_mat6_taylor
use ptc_interface_mod, only: ele_to_taylor

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: start_orb, end_orb
type (lat_param_struct)  param

logical, optional :: err_flag

!

if (present(err_flag)) err_flag = .false.

if (.not. associated(ele%taylor(1)%term)) call ele_to_taylor(ele, param, start_orb)

call mat_make_unit (ele%mat6)

end_orb = start_orb
call track1_taylor (end_orb, ele, param, mat6 = ele%mat6, make_matrix = .true.)

ele%vec0 = end_orb%vec - matmul(ele%mat6, start_orb%vec)

end subroutine

