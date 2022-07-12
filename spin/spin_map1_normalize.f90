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
character(*), parameter :: r_name = 'spin_map1_normalize'

!

f = norm2(spin1(:,0))
if (f < 0.5 .or. f > 2) then
  call out_io (s_warn$, r_name, 'Spin map norm is: ' // real_str(f, 3) // ' which is very different from 1.', &
                                '  This will make spin tracking unreliable')
endif

spin1(:,0) = spin1(:,0) / f

do i = 1, 6
  f = dot_product(spin1(:,0), spin1(:,i))
  spin1(:,i) = spin1(:,i) - f * spin1(:,0)
enddo

end subroutine
