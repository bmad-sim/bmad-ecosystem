!+
! Subroutine radiation_integrals_custom (lat, ir, orb, rad_int1, err_flag)
!
! Dummy routine for custom elements. Will generate an error if called.
! A valid radiation_integrals_custom is needed only if the 
! radiation_integrals routine is being used.
!
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below.
!
! Input:
!   lat         -- lat_struct: Lattice with the custom element.
!   ir          -- integer: lat%ele(ir) is the custom element.
!   orb(:)      -- coord_struct: Orbit around which integrals are to be evaluated.
!
! Output:
!   rad_int1    -- rad_int1_struct: Structure for storing the results.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!-

subroutine radiation_integrals_custom (lat, ir, orb, rad_int1, err_flag)

use bmad_interface, dummy => radiation_integrals_custom

implicit none

type (lat_struct) lat
type (coord_struct) orb(0:)
type (rad_int1_struct) rad_int1
integer ir
logical err_flag
character(32) :: r_name = 'radiation_integrals_custom'

!

call out_io (s_fatal$, r_name, 'THIS DUMMY ROUTINE SHOULD NOT HAVE BEEN CALLED IN THE FIRST PLACE.')
err_flag = .true.

end subroutine
