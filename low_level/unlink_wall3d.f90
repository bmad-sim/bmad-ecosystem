!+
! Subroutine unlink_wall3d (wall3d)
!
! Routine to deallocate a wall3d pointer.
!
! Input:
!   wall3d(:) -- wall3d_struct, pointer: Pointer to wall3d structure.
!
! Output:
!   wall3d(:) -- wall3d_struct, pointer: deallocated
!-

subroutine unlink_wall3d (wall3d)

use bmad_struct

implicit none

type (wall3d_struct), pointer :: wall3d(:)
integer i

!

if (associated (wall3d)) then
  wall3d%n_link = wall3d%n_link - 1
  if (wall3d(1)%n_link == 0) then
    do i = 1, size(wall3d)
      deallocate (wall3d(i)%section)
    enddo
    deallocate (wall3d)
  else
    nullify(wall3d)
  endif
endif

end subroutine unlink_wall3d

