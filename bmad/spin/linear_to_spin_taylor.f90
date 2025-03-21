!+
! Subroutine linear_to_spin_taylor(q_map, spin_taylor)
!
! Routine to create a 1st order spin Taylor map from a linear quaternion map
!
! Input:
!   q_map(0:3,0:6)    -- real(rp): Linear quaternion map.
!
! Output:
!   spin_taylor(0:3)  -- taylor_struct: Taylor map
!-

subroutine linear_to_spin_taylor(q_map, spin_taylor)

use bmad_interface, dummy => linear_to_spin_taylor

type (taylor_struct) spin_taylor(0:3)
real(rp) q_map(0:3, 0:6)
integer i, j

!

do i = 0, 3
  call init_taylor_series(spin_taylor(i), 7)
  spin_taylor(i)%term(1) = taylor_term_struct(q_map(i,0), [0,0,0,0,0,0])
  do j = 1, 6
    spin_taylor(i)%term(j+1) = taylor_term_struct(q_map(i,j), taylor_expn([j]))
  enddo
enddo

end subroutine linear_to_spin_taylor

