!+
! Subroutine radiation_integrals_custom (lat, ir, orb, err_flag)
!
! Dummy routine for custom elements. Will generate an error if called.
! A valid radiation_integrals_custom is needed only if the 
! radiation_integrals routine is being used.
!
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below."
!
! Modules needed:
!   use rad_int_common
!
! Input:
!   lat    -- lat_struct: Lattice with the custom element.
!   ir     -- Integer: lat%ele(ir) is the custom element.
!   orb(:) -- Coord_struct: Orbit around which integrals are to be evaluated.
!
! Output:
!   ric      -- Rad_int_all_ele_struct: Common block for storing the results.
!   err_flag -- Logical: Set true if there is an error. False otherwise.
!-

subroutine radiation_integrals_custom (lat, ir, orb, err_flag)

use bmad_struct
use bmad_interface

implicit none

type (lat_struct) lat
type (coord_struct) orb(0:)
integer ir
logical err_flag
character(32) :: r_name = 'radiation_integrals_custom'

!

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
err_flag = .true.

end subroutine
