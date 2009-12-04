!+
! Subroutine track1_bunch_custom (bunch_start, lat, ele, bunch_end)
!
! Dummy routine for custom bunch tracking. 
! If called, this routine will generate an error message and quit.
! This routine needs to be replaced for a custom calculation.
!
! Modules Needed:
!   use beam_def_struct
!
! Input:
!   bunch_start -- bunch_struct: Starting bunch position.
!   lat         -- lat_struct: Lattice containing element to be tracked through.
!   ele         -- Ele_struct: Element to track through.
!
! Output:
!   bunch_end -- Bunch_struct: Ending bunch position.
!-

#include "CESR_platform.inc"

subroutine track1_bunch_custom (bunch_start, lat, ele, bunch_end)

use bmad_interface, except_dummy => track1_bunch_custom
use bmad_struct
use beam_def_struct

implicit none

type (bunch_struct) bunch_start, bunch_end
type (lat_struct), target :: lat
type (ele_struct) :: ele

!

print *, 'ERROR: DUMMY TRACK1_BUNCH_CUSTOM CALLED FOR: ', ele%name
print *, '       EITHER CUSTOM TRACKING_METHOD WAS CALLED BY MISTAKE,'
print *, '       OR THE CORRECT ROUTINE WAS NOT LINKED IN!'
call err_exit

end subroutine
