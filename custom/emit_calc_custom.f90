!+
! Subroutine emit_calc_custom (lat, ir, i2, i3, i5a, i5b) 
!
! Dummy routine for custom elements. Will generate an error if called.
! A valid emit_calc_custom is needed only if the emit_calc routine 
! is being used with custom elements.
!
! Arguments needed to construct a emit_calc_custom routine:
!
! Input:
!   lat -- lat_struct: Lattice with the custom element.
!   ir  -- Integer: lat%ele(ir) is the custom element.
!
! Output:
!   i2   -- Real(rp): Integral of g^2 (g = 1/bending_radius)
!   i3   -- Real(rp): Integral of g^3
!   i5a  -- Real(rp): a mode i5 integral
!   i5b  -- Real(rp): b mode i5 integral
!-

subroutine emit_calc_custom

  print *, 'ERROR IN EMIT_CALC_CUSTOM: THIS DUMMY ROUTINE SHOULD NOT HAVE'
  print *, '      BEEN CALLED IN THE FIRST PLACE.'
  call err_exit

end subroutine
