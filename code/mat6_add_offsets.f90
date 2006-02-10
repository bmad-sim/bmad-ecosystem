#include "CESR_platform.inc"

!+
! Subroutine mat6_add_offsets (ele)
!
! Subroutine to add in the affect of an element's orientation in space to
! to the computed Jacobian matrix.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element with given orientation.
!     %vec0(6)   -- Real(rp): 0th order part of the transfer map.
!     %mat6(6,6) -- Real(rp): 1st order part of the transfer map (Jacobian).
!
! Output:
!   ele       -- Ele_struct: Element with given orientation.
!     %vec0(6)   -- Real(rp): 0th order part of the transfer map.
!     %mat6(6,6) -- Real(rp): 1st order xfer map.
!-

subroutine mat6_add_offsets (ele)

use bmad_struct
use bmad_interface
use bmad_utils_mod

implicit none

type (ele_struct) ele
type (coord_struct) orb
type (param_struct) param

! the new vec0 is obtained by just tracking through the element

orb%vec = 0
call offset_particle (orb, param, set$, set_canonical = .false., set_hvkicks = .false.)
orb%vec = ele%vec0 + matmul (ele%mat6, orb%vec)
call offset_particle (orb, param, unset$, set_canonical = .false., set_hvkicks = .false.)
ele%vec0 = orb%vec

! calculate the new Jacobian.

if (ele%value(tilt_tot$) /= 0) call tilt_mat6 (ele%mat6, ele%value(tilt_tot$))
call mat6_add_pitch (ele, ele%mat6)

end subroutine

