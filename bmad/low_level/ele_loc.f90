!+
! Function ele_loc (ele) result (loc)
!
! Function to return an lat_ele_loc_struct identifying where an element is in the lattice.
!
! Input:
!   ele     -- Ele_struct: Element to be identified
!
! Output:
!   loc     -- Lat_ele_loc_struct: Element identifier.
!-

function ele_loc (ele) result (loc)

use bmad_struct

implicit none

type (ele_struct) ele
type (lat_ele_loc_struct) loc

!

loc%ix_ele = ele%ix_ele
loc%ix_branch = ele%ix_branch

end function ele_loc
