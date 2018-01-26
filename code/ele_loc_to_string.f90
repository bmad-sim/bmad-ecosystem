!+
! Function ele_loc_to_string (ele, show_branch0) result (str)
!
! Routine to encode an element's location into a string.
! Example output:
!   "34"     ! Input: lat%ele(34) which is equivalent to lat%branch(0)%ele(34)
!   "0>>34"  ! Same as above if show_branch0 is set to True.
!   "1>>56"  ! Input: lat%branch(1)%ele(56).
!
! Modules needed:
!   use bmad
!
! Input:
!   ele          -- Ele_struct: Element in a lattice
!   show_branch0 -- Logical, optional: Explicitly show branch for main 
!                     lattice elements? Default is False.
!
! Output:
!   str(10)     -- Character: Output string. Left justified.
!-

function ele_loc_to_string (ele, show_branch0) result (str)

use bmad_struct

implicit none

type (ele_struct) ele
logical, optional :: show_branch0

character(10) str

!

if (ele%ix_branch == 0 .and. .not. logic_option(.false., show_branch0)) then
  write (str, '(i0)') ele%ix_ele
else
  write (str, '(i0, a, i0)') ele%ix_branch, '>>', ele%ix_ele
endif

end function ele_loc_to_string 

