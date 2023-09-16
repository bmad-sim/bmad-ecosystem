!+
! Subroutine transfer_wall3d (wall3d_in, wall3d_out)
!
! Subroutine to point wall3d_out => wall3d_in
!
! Input:
!   wall3d_in(:)  -- Wall3d_struct, pointer: Input wall3dgler field.
!
! Output:
!   wall3d_out(:) -- Wall3d_struct, pointer: Output wall3dgler field.
!-

subroutine transfer_wall3d (wall3d_in, wall3d_out)

use bmad_routine_interface, dummy => transfer_wall3d

implicit none

type (wall3d_struct), pointer :: wall3d_in(:), wall3d_out(:)

!

if (.not. associated(wall3d_in) .and. .not. associated(wall3d_out)) return
if (associated(wall3d_in, wall3d_out)) return

! If both associated must be pointing to different memory locations

if (associated(wall3d_out)) call unlink_wall3d(wall3d_out)

if (associated(wall3d_in)) then 
  wall3d_out => wall3d_in
  wall3d_out%n_link = wall3d_out%n_link + 1
endif

end subroutine transfer_wall3d

