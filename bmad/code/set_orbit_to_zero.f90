!+
! Subroutine set_orbit_to_zero (orbit, n1, n2, ix_noset)
!
! Routine to set the orbit of a subset orbit(n1:n2) of a coord_struct array to zero.
! This routine actually initalizes each orbit(n) so that, in particular, orbit(n)%state = not_set$
!
! Input:
!   n1        -- integer: Lower bound of orbit(:) array subset.
!   n2        -- integer: Upper bound of orbit(:) array subset.
!   ix_noset  -- integer, optional: If present then orbit(ix_noset) will not be zeroed.
!
! Output:
!   orbit(:)  -- coord_struct: Array with particle positions in the range orbit(n1:n2)
!                   set to zero except for orbit(ix_noset).
!-

subroutine set_orbit_to_zero (orbit, n1, n2, ix_noset)

use bmad_struct

implicit none

type (coord_struct) orbit(0:)
integer n, n1, n2
integer, optional :: ix_noset

! Note: n will never equal -1 so the if statement does what it should do.

do n = n1, n2
  if (n == integer_option(-1, ix_noset)) cycle 
  orbit(n) = coord_struct()
enddo

end subroutine set_orbit_to_zero
