!+
! Subroutine custom_radiation_integrals (lat, ir, orb)
!
! Dummy routine for custom elements. Will generate an error if called.
! A valid custom_radiation_integrals is needed only if the 
! radiation_integrals routine is being used.
!
! Modules and arguments needed to construct a custom_emitt_calc routine:
!
! Modules needed:
!   use rad_int_common
!
! Input:
!   lat    -- lat_struct: Lattice with the custom element.
!   ir      -- Integer: lat%ele(ir) is the custom element.
!   orb(:) -- Coord_struct: Orbit around which integrals are to be evaluated.
!
! Output:
!   ric  -- Rad_int_common_struct: Common block for stolat the results.
!     %i1(ir)  -- I1 integral.
!     %i2(ir)  -- I2 integral.
!     %i3(ir)  -- I3 integral.
!     %i4a(ir) -- I4a integral.
!     %i4b(ir) -- I4b integral.
!     %i5a(ir) -- I5a integral.
!     %i5b(ir) -- I5b integral.
!-

#include "CESR_platform.inc"

subroutine custom_radiation_integrals (lat, ir, orb)

  use bmad_struct, only: lat_struct, coord_struct

  implicit none

  type (lat_struct) lat
  type (coord_struct) orb(0:)
  integer ir

  print *, 'ERROR IN CUSTOM_RADIATION_INTEGRALS: THIS DUMMY ROUTINE SHOULD NOT' 
  print *, '      HAVE BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
