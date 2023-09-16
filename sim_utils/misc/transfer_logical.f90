!-
! Subroutine transfer_logical (from_logic, to_logic)
!
! Routine to transfer the logical value from one logical variable to another.
! This routine is used in conjunction with the routines:
!   set_logical_to_garbage
!   logical_is_garbage
! See these routines for more details.
!
! The idea of this routine is that the simple set:
!   to_logic = from_logic
! might possibly result in a run-time error if from_logic is in a garbage state.
! Whether such a run-time error would be generated has not been tested though.
!
! Input:
!   from_logic -- Logical: Input logical.
!
! Output:
!   to_logci -- Logical: Logical to set.
!-

subroutine transfer_logical (from_logic, to_logic)

use physical_constants
implicit none
integer from_logic, to_logic

! This routine works by pretending that from_logic and to_logic arguments are integers.

to_logic = from_logic

end subroutine
