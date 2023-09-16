!+
! Function ele_full_name (ele, template) result (str)
!
! Routine to encode an element's name and/or location into a string.
!
! The template uses the tokens:
!   "@N"   Translates to element name.
!   "%#"   Translates to "ix_ele" for elements in branch 0 and "ix_branch>>ix_ele" for all others.
!   "&#"   Translates to "ix_ele" if the lattice has only one branch. Otherwise translates to "ix_ele>>ix_branch".
!   "!#"   Always translates to "ix_ele>>ix_branch".
!
! Examples:
!   template     str
!   --------     ------------
!   "@N (%#)"    "Q7 (34)"        ! If Q7 is in branch 0
!   "@N (%#)"    "Q7 (2>>34)"     ! With Q7 in branch 2
!   "[&#]"       "[34]"           ! If lattice has only one branch.
!   "[&#]"       "[0>>34]"        ! If lattice has multiple branches.
!   "{!#}"       "{0>>34}"        ! "!#" will always translate to "ix_ele>>ix_branch"
! Input:
!   ele           -- ele_struct: Element in a lattice
!   template      -- character(*), optional: Encoding template. Default is "@N (&#)".
!
! Output:
!   str           -- character(:), allocatable: Name/location string. 
!-

function ele_full_name (ele, template) result (str)

use bmad_struct

implicit none

type (ele_struct) ele
character(*), optional :: template
character(:), allocatable :: str
integer ix

!

if (present(template)) then
  str = template
else
  str = '@N (&#)'
endif

!

ix = index(str, '@N')
if (ix == 0) then
  str = str
else
  str = str(1:ix-1) // trim(ele%name) // str(ix+2:)
endif

!

ix = index(str, '!#')
if (ix /= 0) str = str(1:ix-1) // int_str(ele%ix_branch) // '>>' // int_str(ele%ix_ele) // str(ix+2:)

!

ix = index(str, '&#')
if (ix /= 0) then
  if (size(ele%branch%lat%branch) == 1) then
    str = str(1:ix-1) // int_str(ele%ix_ele) // str(ix+2:)
  else
    str = str(1:ix-1) // int_str(ele%ix_branch) // '>>' // int_str(ele%ix_ele) // str(ix+2:)
  endif
endif

!

ix = index(str, '%#')
if (ix /= 0) then
  if (ele%ix_branch == 0) then
    str = str(1:ix-1) // int_str(ele%ix_ele) // str(ix+2:)
  else
    str = str(1:ix-1) // int_str(ele%ix_branch) // '>>' // int_str(ele%ix_ele) // str(ix+2:)
  endif
endif

end function ele_full_name 

