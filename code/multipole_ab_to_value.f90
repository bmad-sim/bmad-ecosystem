!+
! Subroutine multipole_ab_to_value (a, b, value)
!
! Subroutine to transform a, b multipole value arrays to a ele%value array.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   a(0:n_pole_maxx) -- Real: Array of An values.
!   b(0:n_pole_maxx) -- Real: Array of Bn values.
!
! Output:
!   value(n_attrib_maxx) -- Real: Array of An/Bn values.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:54  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine multipole_ab_to_value (a, b, value)

  use bmad_struct

  real a(0:n_pole_maxx), b(0:n_pole_maxx), value(n_attrib_maxx)

!

  value(ix1_m$:ix2_m$-1:2) = a
  value(ix1_m$+1:ix2_m$:2) = b

end subroutine
