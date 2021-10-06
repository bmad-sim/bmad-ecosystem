!+
! Subroutine map1_make_unit(map1)
!
! Routine to make the unit spin/orbit map
!
! Output:
!   map1    -- spin_orbit_map1_struct: Unit map.
!-

subroutine map1_make_unit(map1)

use bmad_struct

implicit none

type (spin_orbit_map1_struct) map1

!

map1 = spin_orbit_map1_struct()
call mat_make_unit (map1%orb_mat)
map1%spin_q(0,0) = 1

end subroutine
