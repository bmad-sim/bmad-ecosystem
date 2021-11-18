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

real(rp) mq(0:3,6)
integer i

!

call mat_inverse(map1%orb_mat, inv_map1%orb_mat)
inv_map1%vec0 = -matmul(inv_map1%orb_mat, map1%vec0)

inv_map1%spin_q(0,0)   =  map1%spin_q(0,0)
inv_map1%spin_q(1:3,0) = -map1%spin_q(1:3,0)

mq = matmul(map1%spin_q(:,1:6), inv_map1%orb_mat)
do i = 1, 6
  inv_map1%spin_q(:,i) = -quat_mul(quat_mul(inv_map1%spin_q(:,0), mq(:,i)), inv_map1%spin_q(:,0))
enddo

end function

