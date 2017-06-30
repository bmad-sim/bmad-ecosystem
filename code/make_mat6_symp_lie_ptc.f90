!+
! Subroutine make_mat6_symp_lie_ptc (ele, param, c0, c1)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix
!   param  -- lat_param_struct: Parameters are needed for some elements.
!   c0     -- Coord_struct: Coordinates at the beginning of element. 
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %vec0  -- 0th order map component
!     %mat6  -- 6x6 transfer matrix.
!   c1     -- Coord_struct: Coordinates at end of element.
!-

subroutine make_mat6_symp_lie_ptc (ele, param, c0, c1)

use ptc_interface_mod, except_dummy => make_mat6_symp_lie_ptc

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: c0, c1
type (lat_param_struct)  param
type (taylor_struct) bmad_taylor(6)

logical st

!

st = bmad_com%spin_tracking_on
bmad_com%spin_tracking_on = .false.

call ele_to_taylor(ele, param, bmad_taylor, c0, .true.)
call taylor_to_mat6 (bmad_taylor, c0%vec, ele%vec0, ele%mat6)
call track1_taylor (c0, ele, param, c1, bmad_taylor)
call kill_taylor (bmad_taylor)

bmad_com%spin_tracking_on = st

end subroutine

