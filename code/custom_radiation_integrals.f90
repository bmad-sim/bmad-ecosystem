!+
! Subroutine custom_radiation_integrals (ring, ir, orb_)
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
!   ring    -- Ring_struct: Lattice with the custom element.
!   ir      -- Integer: ring%ele_(ir) is the custom element.
!   orb_(:) -- Coord_struct: Orbit around which integrals are to be evaluated.
!
! Output:
!   ric  -- Rad_int_common_struct: Common block for storing the results.
!     %i1_(ir)  -- I1 integral.
!     %i2_(ir)  -- I2 integral.
!     %i3_(ir)  -- I3 integral.
!     %i4a_(ir) -- I4a integral.
!     %i4b_(ir) -- I4b integral.
!     %i5a_(ir) -- I5a integral.
!     %i5b_(ir) -- I5b integral.
!-

#include "CESR_platform.inc"

subroutine custom_radiation_integrals

  print *, 'ERROR IN CUSTOM_RADIATION_INTEGRALS: THIS DUMMY ROUTINE SHOULD NOT' 
  print *, '      HAVE BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
