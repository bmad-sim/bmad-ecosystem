!+
! Subroutine custom_emitt_calc (ring, ir, i2, i3, i5a, i5b) 
!
! Dummy routine for custom elements. Will generate an error if called.
! A valid custom_emitt_calc is needed only if the emitt_calc routine 
! is being used.
!
! Arguments needed to construct a custom_emitt_calc routine:
!
! Input:
!   ring -- Ring_struct: Lattice with the custom element.
!   ir   -- Integer: ring%ele_(ir) is the custom element.
!
! Output:
!   i2   -- Real(rp): Integral of g^2 (g = 1/bending_radius)
!   i3   -- Real(rp): Integral of g^3
!   i5a  -- Real(rp): a mode i5 integral
!   i5b  -- Real(rp): b mode i5 integral
!-

#include "CESR_platform.inc"

subroutine custom_emitt_calc

  print *, 'ERROR IN CUSTOM_EMITT_CALC: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
