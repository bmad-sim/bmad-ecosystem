!+
! Subroutine set_logical_to_garbage (logic)
!
! Routine to set the bit representation of a logical to a value which is different from the bit
! representation of .true. and .false.
!
! This routine is useful when using Fortran namelist input and one wants to determine if a logical is
! being set. The idea is to use this routine to set the bit pattern of the logical before the namelist
! read and then use the routine logical_is_garbage to test if the logical has been set.
!
! Also see the transfer_logical routine
!
! Output:
!   logic -- Logical: Variable with bit representation set to a "garbage" value.
!-

subroutine set_logical_to_garbage (logic)

use sim_utils_struct
implicit none
integer logic

! This routine works by pretending that the logic argument is an integer

logic = int_garbage$

end subroutine
