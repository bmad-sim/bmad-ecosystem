!+
! Subroutine em_field_custom (ele, param, s_pos, here, field, calc_dfield)
!
! Dummy routine. Will generate an error if called.
! A valid em_field_custom is needed only if the em_field routine is being used.
!
! Arguments needed to construct a custom_emitt_calc routine:
!
! Input:
!   ele         -- Ele_struct: Custom element.
!   param       -- Param_struct: Lattice parameters.
!   s_pos       -- Real(rp): Longitudinal position of the reference particle.
!   here        -- Coord_struct: Coords with respect to the reference particle.
!   calc_dfield -- Logical, optional: If present and True then the field 
!                     derivative matrix is wanted by the calling program.
!
! Output:
!   field -- Em_field_struct: Structure hoding the field values.
!-

#include "CESR_platform.inc"

subroutine em_field_custom

  print *, 'ERROR IN EM_FIELD_CUSTOM: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
