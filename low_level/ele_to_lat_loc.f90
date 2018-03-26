!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function ele_to_lat_loc (ele) result (ele_loc)
!
! Function to return an lat_ele_loc_struct identifying where an element is in the lattice.
!
! Input:
!   ele -- Ele_struct: Element to be identified
!
! Output:
!   ele_loc -- Lat_ele_loc_struct: Element identifier.
!-

function ele_to_lat_loc (ele) result (ele_loc)

use bmad_struct

type (ele_struct) ele
type (lat_ele_loc_struct) ele_loc

!

ele_loc%ix_ele = ele%ix_ele
ele_loc%ix_branch = ele%ix_branch

end function ele_to_lat_loc
