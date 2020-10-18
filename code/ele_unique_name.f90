!+
! Function ele_unique_name (ele, order) result (unique_name)
!
! Routine to construct a name that is unique for the given element.
!
! For example, if there are multiple elements named "Q" in multiple branches a
! unique name may be something like "2>>QQ##5" which would designate the 5th element
! named "QQ" in branch #2.
!
! Note: The super_lord elements, which live in branch #0, are associated with the branch that its super_slaves live in.
! That is, if a marker is superimpsed upon "2>>QQ##5", the "QQ" lord element will be placed in the lord
! section of branch 0 but this lord's order index, associated branch, and unique name are not changed.
!
! Also see the routine: create_unique_ele_names
!
! Input:
!   ele     -- ele_struct: Element to construct a unique name for.
!   order   -- lat_ele_order_struct: Information on element ordering. 
!               Before calling this routine, use the routine ele_order_calc to compute this argument.
!
! Output:
!   unique_name   -- character(40): Unique name that can can be used to identify ele.
!                     The simplist name will be constructed. For example, if the element name is
!                     unique, unique_name will be set to the element name.
!-

function ele_unique_name (ele, order) result (unique_name)

use bmad_struct

implicit none

type (ele_struct), target :: ele
type (lat_struct), pointer :: lat
type (lat_ele_order_struct), target :: order
type (lat_ele_order1_struct), pointer :: ord1

character(40) unique_name

!

lat => ele%branch%lat
ord1 => order%branch(ele%ix_branch)%ele(ele%ix_ele)
select case (ord1%ix_order)
case (-1)
  unique_name = ele%name
case (0)
  unique_name = int_str(ord1%ix_branch) // '>>' // trim(ele%name)
case default
  if (ubound(lat%branch, 1) == 0) then
    unique_name = trim(ele%name) // '##' // int_str(ord1%ix_order)
  else
    unique_name = int_str(ord1%ix_branch) // '>>' // trim(ele%name) // '##' // int_str(ord1%ix_order)
  endif
end select

end function ele_unique_name
