!+
! Function polar_to_spinor (polar) result (spinor)
!
! Converts a spin vector in polar coordinates to a spinor
!
! Input:
!   polar     -- spin_polar_struct: includes polar phase
!
! Output:
!   spinor(2)   -- complex(rp): Spinor
!-

function polar_to_spinor (polar) result (spinor)

use equal_mod, dummy_except => polar_to_spinor

implicit none

type (spin_polar_struct) polar
complex(rp) :: spinor(2)

!

spinor(1) = sqrt(polar%polarization) * Exp(i_imag * polar%xi) * cos(polar%theta / 2.0d0)
spinor(2) = sqrt(polar%polarization) * Exp(i_imag * (polar%xi+polar%phi)) * sin(polar%theta / 2.0d0)

end function polar_to_spinor


