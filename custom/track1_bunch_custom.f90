!+
! Subroutine track1_bunch_custom (bunch_start, lat, ele, bunch_end, err_flag)
!
! Dummy routine for custom bunch tracking. 
! If called, this routine will generate an error message and quit.
! This routine needs to be replaced for a custom calculation.
!
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below."
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
!   err_flag  -- Logical: Set true if there is an error. False otherwise.
!-

#include "CESR_platform.inc"

subroutine track1_bunch_custom (bunch_start, lat, ele, bunch_end, err_flag)

use bmad_interface, except_dummy => track1_bunch_custom
use bmad_struct
use beam_def_struct

implicit none

type (bunch_struct) bunch_start, bunch_end
type (lat_struct), target :: lat
type (ele_struct) :: ele
logical err_flag

character(32) :: r_name = 'track1_bunch_custom'

!

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
err_flag = .true.

end subroutine
