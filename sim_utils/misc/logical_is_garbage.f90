!+
! Function logical_is_garbage (logic) result (is_garbage)
!
! Routine to test the bit representation of a logical to see if it has the same bit representation
! as set by the set_logical_to_garbage routine.
!
! This routine is useful when using Fortran namelist input and one wants to determine if a logical is
! being set. The idea is to use set_logical_to_garbage to set the bit pattern of the logical before the 
! namelist read and then use this routine to test if the logical has been set.
!
! Also see the transfer_logical routine
!
! Input:
!   logic -- Logical: Variable with bit representation set to a "garbage" value.
!
! Outpu:
!   is_garbage -- Logical: Set True if the bit pattern of the logic argument is garbage. False otherwise
!-

function logical_is_garbage (logic) result (is_garbage)

use sim_utils_struct
implicit none
integer logic
logical is_garbage

! This routine works by pretending that the logic argument is an integer

is_garbage = (logic == int_garbage$)

end function
