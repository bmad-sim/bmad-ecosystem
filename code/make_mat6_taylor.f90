!+
! Subroutine make_mat6_taylor (ele, param, orb_in)
!
! Subroutine to make the 6x6 transfer matrix for an element. 
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix
!   param  -- lat_param_struct: Parameters are needed for some elements.
!   orb_in -- Coord_struct: Coordinates at the beginning of element. 
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %vec0  -- 0th order map component
!     %mat6  -- 1st order map component (6x6 transfer matrix).
!-

#include "CESR_platform.inc"

subroutine make_mat6_taylor (ele, param, orb_in)

use ptc_interface_mod, except => make_mat6_taylor
use make_mat6_mod, only: tilt_mat6

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: orb0, orb_in, orb_out
type (lat_param_struct)  param

!

if (.not. associated(ele%taylor(1)%term)) call ele_to_taylor(ele, param, orb_in)

if (ele%map_with_offsets) then
  call taylor_to_mat6 (ele%taylor, orb_in%vec, ele%vec0, ele%mat6)

else
  orb0 = orb_in
  call offset_particle (ele, param, orb0, set$, set_canonical = .false.)
  call taylor_to_mat6 (ele%taylor, orb0%vec, ele%vec0, ele%mat6, orb_out%vec)
  call offset_particle (ele, param, orb_out, unset$, set_canonical = .false.)

  if (ele%value(tilt_tot$) /= 0) call tilt_mat6 (ele%mat6, ele%value(tilt_tot$))

  call mat6_add_pitch (ele, ele%mat6)

  ele%vec0 = orb_out%vec - matmul(ele%mat6, orb_in%vec)

endif

end subroutine

