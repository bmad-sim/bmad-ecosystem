!+
! Function tao_curve_ix_uni (curve) result (ix_uni)
!
! Routine to return the universe index associated with a plot curve.
!
! Input:
!   curve     -- tao_curve_struct: Curve.
!
! Output:
!   ix_uni    -- integer: Universe index.
!-

function tao_curve_ix_uni (curve) result (ix_uni)

use tao_interface, dummy => tao_curve_ix_uni

implicit none

type (tao_curve_struct) curve
integer ix_uni

!

if (curve%ix_universe == -1) then
  ix_uni = curve%g%ix_universe
else
  ix_uni = curve%ix_universe
endif

end function tao_curve_ix_uni
