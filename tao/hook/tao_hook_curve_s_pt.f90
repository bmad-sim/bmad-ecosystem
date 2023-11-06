!+
! Function tao_hook_curve_s_pt (s_default, ix_now, x1, x2, n_pts, tao_lat, curve) result (s_pt)
!
! Routine to calculate custom s-values for evaluating points on a plot curve.
! Tao by default will pick evenally spaced points but this is sometimes not wanted if 
! the lattice length is varied.
!
! Input:
!   s_default -- real(rp): The default evaluation point value.
!   ix_now    -- integer: The index of the evaluation point. Runs from 1 to n_pts.
!   x1        -- real(rp): Lower bound evaluation point value when ix_now = 1.
!   x2        -- real(rp): Upper bound evaluation point value when ix_now = n_pts.
!   tao_lat   -- tao_latticc_struct: The lattice used.
!   curve     -- tao_curve_struct: Curve under consideration.
!
! Output:   
!   s_pt      -- real(rp): The evaluation point value. The behavior of the unaltered tao_hook_curve_s_pt
!                 routine will be to set s_pt = s_default.
!-

function tao_hook_curve_s_pt (s_default, ix_now, x1, x2, n_pts, tao_lat, curve) result (s_pt)

use tao_interface

implicit none

type (tao_lattice_struct) tao_lat
type (tao_curve_struct) curve

real(rp) s_default, x1, x2, s_pt
integer ix_now, n_pts

! 

s_pt = s_default

end function
