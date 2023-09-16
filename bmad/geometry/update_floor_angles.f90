!+
! Subroutine update_floor_angles (floor, floor0)
!
! Routine to calculate floor angles from its W matrix.
!
! Input:
!   floor -- floor_position_struct: Position with input w matrix.
!     %w      -- w matrix.
!   floor0 -- floor_position_struct, optional: Reference position. There are two solutions related by:
!                   [theta, phi, psi] & [pi+theta, pi-phi, pi+psi]
!                 If floor0 is present, choose the solution "nearest" the angles in floor0.
!
! Output:
!   floor -- floor_position_struct: Position with output angles.
!     %theta, %phi, %psi -- Orientation angles.
!-

subroutine update_floor_angles (floor, floor0)

use bmad_interface, dummy => update_floor_angles

implicit none

type(floor_position_struct) :: floor
type(floor_position_struct), optional :: floor0

!

if (present(floor0)) then
  call floor_w_mat_to_angles (floor%W, floor%theta, floor%phi, floor%psi, floor0)
else
  call floor_w_mat_to_angles (floor%W, floor%theta, floor%phi, floor%psi)
endif

end subroutine update_floor_angles
