!+
! Subroutine make_mat6_tracking (ele, param, c0, c1)
!
! Subroutine to make the 6x6 transfer matrix for an element using the
! Present tracking method.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele    -- Ele_struct: Element with transfer matrix
!   param  -- Param_struct: Parameters are needed for some elements.
!   c0     -- Coord_struct: Coordinates at the beginning of element. 
!   bmad_com%d_orb(6) -- Vector of offsets to use. 
!               Default if d_orb = 0 is to set the offset to 1e-5
!
! Output:
!   ele    -- Ele_struct: Element with transfer matrix.
!     %mat6  -- 6x6 transfer matrix.
!   c1     -- Coord_struct: Coordinates at the end of element.
!-

#include "CESR_platform.inc"

subroutine make_mat6_tracking (ele, param, c0, c1)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct), target :: ele
  type (coord_struct) :: c0, c1, start1, end1
  type (param_struct)  param

  real (rp) error, del_orb(6)
  integer i

!

  call track1 (c0, ele, param, c1)

  if (all(bmad_com%d_orb == 0)) then
    del_orb = 1e-5
  else
    del_orb = bmad_com%d_orb
  endif

  do i = 1, 6
    start1 = c0
    start1%vec(i) = start1%vec(i) + del_orb(i)
    call track1 (start1, ele, param, end1)
    ele%mat6(1:6, i) = (end1%vec - c1%vec) / del_orb(i)
  enddo

end subroutine

