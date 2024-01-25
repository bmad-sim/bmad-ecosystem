!+
! Subroutine kill_taylor (bmad_taylor)
!
! Subroutine to deallocate a taylor map.
! It is OK if the taylor has already been deallocated.
!
! Input:
!   bmad_taylor(:)   -- Taylor_struct, optional: Taylor to be deallocated. 
!
! Output:
!   bmad_taylor(:)   -- Taylor_struct, optional: deallocated Taylor structure.
!-

subroutine kill_taylor (bmad_taylor)

use bmad_struct

implicit none

type (taylor_struct) :: bmad_taylor(:)

integer i

!

do i = 1, size(bmad_taylor)
  if (associated(bmad_taylor(i)%term)) deallocate (bmad_taylor(i)%term)
enddo

end subroutine kill_taylor

