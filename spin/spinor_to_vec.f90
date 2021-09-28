!+
! Function spinor_to_vec (spinor) result (vec)
!
! Converts a spinor to a spin cartesian vector
!
! Input:
!   spinor(2)  -- complex(rp): Spinor
!
! Output
!   vec(3)     -- Real(rp): spin vector in cartesian coordinates
!-

function spinor_to_vec (spinor) result (vec)

use equal_mod, dummy_except => spinor_to_vec

implicit none

complex(rp) spinor(2)
real(rp) vec(3)

!

! vec = conjg(spinor) * pauli(i)%sigma * spinor done explicitly
vec(1) = 2 * (real(spinor(1))*real(spinor(2)) + aimag(spinor(1))*aimag(spinor(2)))
vec(2) = 2 * (real(spinor(1))*aimag(spinor(2)) - aimag(spinor(1))*real(spinor(2)))
vec(3) = real(spinor(1))**2 + aimag(spinor(1))**2 - real(spinor(2))**2 - aimag(spinor(2))**2

end function spinor_to_vec


