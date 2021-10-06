!+
! Subroutine spin_map1_normalize (spin1)
!
! Routine to normalize a spin linear map so that the 0th order part has unit mangitude and
! the 1st order part is normal to the 0th order part.
!
! Input:
!   spin1(0:3,0:6)    -- real(rp): Unnormalized spin map.
!
! Output:
!   spin1(0:3,0:6)    -- real(rp): Normalized spin map.
!-

subroutine spin_map1_normalize (spin1)

use bmad_struct

implicit none

real(rp) spin1(0:3,0:6), f
integer i

!

spin1(:,0) = spin1(:,0) / norm2(spin1(:,0))

do i = 1, 6
  f = dot_product(spin1(:,0), spin1(:,i))
  spin1(:,i) = spin1(:,i) - f * spin1(:,0)
enddo

end subroutine
