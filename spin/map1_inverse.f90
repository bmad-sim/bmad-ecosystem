!+
! Function map1_inverse (map1) result (inv_map1)
!
! Routine to form the inverse of a spin/orbit linear map
!
! Input:
!   map1      -- spin_orbit_map1_struct: Input map.
!
! Output:
!   inv_map1  -- spin_orbit_map1_struct: Inverse map.
!-

function map1_inverse (map1) result (inv_map1)

use bmad_struct
implicit none

type (spin_orbit_map1_struct) map1, inv_map1

!

call mat_inverse(map1%orb_mat, inv_map1%orb_mat)
inv_map1%vec0 = -matmul(inv_map1%orb_mat, map1%vec0)

inv_map1%spin_q(0,:)  =  map1%spin_q(0,:)
inv_map1%spin_q(1:,:) = -map1%spin_q(1:,:)

end function

