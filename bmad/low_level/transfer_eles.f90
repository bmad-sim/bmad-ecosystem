!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!+
! Subroutine transfer_eles (ele1, ele2)
!
! Subroutine to set ele2 = ele1. 
! This is a plain transfer of information not using the overloaded equal sign.
! Thus, at the end, ele2's pointers point to the same memory as ele1's.
!
! NOTE: Do not use this routine unless you know what you are doing!
!
! Input:
!   ele1(:) -- Ele_struct:
!
! Output:
!   ele2(:) -- Ele_struct:
!-

subroutine transfer_eles (ele1, ele2)

use bmad_struct

implicit none

type (ele_struct), intent(inout) :: ele1(:)
type (ele_struct), intent(inout) :: ele2(:)

ele2 = ele1

end subroutine transfer_eles

