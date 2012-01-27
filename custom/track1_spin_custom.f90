!+
! Subroutine track1_spin_custom (start, ele, param, end, track)
!
! Dummy routine for custom spin tracking. 
! If called, this routine will generate an error message and quit.
! This routine needs to be replaced for a custom calculation.
!
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below."
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
!   track -- track_struct, optional: Structure holding the track information if the 
!             tracking method does tracking step-by-step.
!-

#include "CESR_platform.inc"

subroutine track1_spin_custom (start, ele, param, end, track)

use bmad_interface, except_dummy => track1_spin_custom

implicit none

type (coord_struct) :: start
type (coord_struct) :: end
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

!

print *, 'ERROR: DUMMY TRACK1_SPIN_CUSTOM CALLED FOR: ', ele%name
print *, '       EITHER CUSTOM TRACKING_METHOD WAS CALLED BY MISTAKE,'
print *, '       OR THE CORRECT ROUTINE WAS NOT LINKED IN!'
call err_exit

end%vec = 0  ! so compiler will not complain

end subroutine
