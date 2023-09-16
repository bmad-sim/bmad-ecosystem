!+
! Function ele_loc_name (ele, show_branch0, parens) result (str)
!
! #####################################################################
! THIS ROUTINE IS DEPRECATED. PLEASE USE ELE_FULL_NAME ROUTINE INSTEAD!
! #####################################################################
!
! Routine to encode an element's location into a string.
! Example output:
!   "34"     ! Input: lat%ele(34) which is equivalent to lat%branch(0)%ele(34)
!   "0>>34"  ! Same as above if show_branch0 is set to True.
!   "1>>56"  ! Input: lat%branch(1)%ele(56).
!   "(34)"   ! If parens = "()".
!
! Input:
!   ele           -- Ele_struct: Element in a lattice
!   show_branch0  -- Logical, optional: Explicitly show branch for main 
!                     lattice elements? Default is False.
!   parens        -- character(2), optional: If present, enclose location string using the two characters
!                      supplied. Typically parens will be set to "()" or "[]".
!
! Output:
!   str           -- character(10): Output string. Left justified.
!-

function ele_loc_name (ele, show_branch0, parens) result (str)

use bmad_struct

implicit none

type (ele_struct) ele
logical, optional :: show_branch0

character(10) str
character(2), optional :: parens

!

if (ele%ix_branch == 0 .and. .not. logic_option(.false., show_branch0)) then
  write (str, '(i0)') ele%ix_ele
else
  write (str, '(i0, a, i0)') ele%ix_branch, '>>', ele%ix_ele
endif

if (present(parens)) str = parens(1:1) // trim(str) // parens(2:2)

end function ele_loc_name 

