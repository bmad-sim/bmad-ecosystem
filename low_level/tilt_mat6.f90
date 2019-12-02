!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+      
! Subroutine tilt_mat6 (mat6, tilt)
!
! Subroutine to transform a 6x6 transfer matrix to a new reference frame
! that is tilted in (x, Px, y, Py) with respect to the old reference frame.
!     mat6 -> tilt_mat * mat6 * tilt_mat_inverse
!
! Input:
!   mat6(6,6) -- Real(rp): Untilted matrix.
!   tilt      -- Real(rp): Tilt angle.
!
! Output:
!   mat6(6,6) -- Real(rp): Tilted matrix.
!-

subroutine tilt_mat6 (mat6, tilt)

use bmad_struct

implicit none

real(rp) tilt, mat6(6,6), mm(6,6)
real(rp) c, s

!

if (tilt == 0) return

c = cos(tilt)
s = sin(tilt)

mm(1,:) = c * mat6(1,:) - s * mat6(3,:)
mm(2,:) = c * mat6(2,:) - s * mat6(4,:)
mm(3,:) = c * mat6(3,:) + s * mat6(1,:)
mm(4,:) = c * mat6(4,:) + s * mat6(2,:)
mm(5,:) =     mat6(5,:)
mm(6,:) =     mat6(6,:)

mat6(:,1) = mm(:,1) * c - mm(:,3) * s
mat6(:,2) = mm(:,2) * c - mm(:,4) * s
mat6(:,3) = mm(:,3) * c + mm(:,1) * s
mat6(:,4) = mm(:,4) * c + mm(:,2) * s
mat6(:,5) = mm(:,5)
mat6(:,6) = mm(:,6)

end subroutine tilt_mat6

