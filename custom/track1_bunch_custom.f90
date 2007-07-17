!+
! Subroutine track1_custom (start, ele, param, end)
!
! Dummy routine for custom tracking. 
! If called, this routine will generate an error message and quit.
! This routine needs to be replaced for a custom calculation.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position.
!   ele    -- Ele_struct: Element.
!   param  -- lat_param_struct: Lattice parameters.
!
! Output:
!   end   -- Coord_struct: End position.
!   param -- lat_param_struct: Lattice parameters.
!     %lost -- Logical. Set to true if a particle is lost.
!-

#include "CESR_platform.inc"

subroutine track1_bunch_custom (bunch_start, lat, ix_ele, bunch_end)

use bmad_interface, except => track1_bunch_custom
use bmad_struct
use beam_def_struct

implicit none

type (bunch_struct) bunch_start, bunch_end
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele

integer ix_ele

!

print *, 'ERROR: DUMMY TRACK1_BUNCH_CUSTOM CALLED FOR: ', lat%ele(ix_ele)%name
print *, '       EITHER CUSTOM TRACKING_METHOD WAS CALLED BY MISTAKE,'
print *, '       OR THE CORRECT ROUTINE WAS NOT LINKED IN!'
call err_exit

end subroutine
