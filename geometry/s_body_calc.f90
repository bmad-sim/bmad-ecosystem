!+
! Function s_body_calc (orbit, ele) result (s_body)
!
! Routine to calculate the longitudinal distance from an elements entrance face to a particle.
!
! Input:
!   orbit     -- coord_struct: Particle coordinates.
!   ele       -- ele_struct: Lattice element
!
! Output:
!   s_body    -- real(rp): Body postion.
!-

function s_body_calc (orbit, ele) result (s_body)

use bmad_struct
implicit none

type (coord_struct) orbit
type (ele_struct) ele
real(rp) s_body

!

if (ele%orientation == 1) then
  s_body = orbit%s - ele%s_start
else
  s_body = ele%s - orbit%s
endif

end function
