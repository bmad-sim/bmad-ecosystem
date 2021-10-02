!+
! Function spin_taylor_to_linear (spin_taylor) result (spin_map1)
!
! Routine to truncate a Taylor spin map to order 1.
!
! Input:
!   spin_taylor(0:3)    -- taylor_struct: Taylor spin map.
!
! Output:
!   spin_map1(0:3,0:6)  -- real(rp): First otfrt spin_map.
!-

function spin_taylor_to_linear (spin_taylor) result (spin_map1)

use bmad_routine_interface, dummy => spin_taylor_to_linear

implicit none

type (taylor_struct), target :: spin_taylor(0:3)
type (taylor_struct), pointer :: st

real(rp) spin_map1(0:3,0:6)
integer i, k, n, p

!

spin_map1 = 0

do i = 0, 3
  st => spin_taylor(i)
  do k = 1, size(st%term)
    n = sum(st%term(k)%expn)
    select case (n)
    case (0)
      spin_map1(i,0) = st%term(k)%coef
    case (1)
      do p = 1, 6
        if (st%term(k)%expn(p) == 0) cycle
        spin_map1(i,p) = st%term(k)%coef
        exit
      enddo
    end select
  enddo
enddo

end function
