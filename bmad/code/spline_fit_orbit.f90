!+
! Subroutine spline_fit_orbit (start_orb, end_orb, spline_x, spline_y)
!
! Routine to fit the orbit through an element to a cubic spline.
! s = 0                        -> particle at start_orb
! s = end_orb%s - start_orb%s  -> particle at end_orb
!
! When interpolating:
!   px = (1+orb%vec(6)) * orb%direction * dspline_x/ds
!
! Input:
!   start_orb       -- coord_struct: Starting coords.
!   end_orb         -- coord_struct: Ending coords.
!
!   spline_x(0:3)   -- real(rp): Spline coefs for the horizontal trajectory.
!   spline_y(0:3)   -- real(rp): Spline coefs for vertical trajectory.
!-

subroutine spline_fit_orbit (start_orb, end_orb, spline_x, spline_y)

use bmad_routine_interface, dummy_except => spline_fit_orbit

implicit none

type (coord_struct) start_orb, end_orb
real(rp) spline_x(0:3), spline_y(0:3)
real(rp) ds, alpha, beta

!

ds = end_orb%s - start_orb%s

! X

spline_x(0) = start_orb%vec(1)
spline_x(1) = start_orb%direction * start_orb%vec(2) / (1 + start_orb%vec(6))

alpha = end_orb%vec(1) - spline_x(0) - spline_x(1) * ds
beta = end_orb%direction * end_orb%vec(2) / (1 + end_orb%vec(6)) - spline_x(1)

spline_x(2) = 3 * alpha / ds**2 - beta / ds
spline_x(3) = beta / ds**2 - 2 * alpha / ds**3

! Y

spline_y(0) = start_orb%vec(3)
spline_y(1) = start_orb%direction * start_orb%vec(4) / (1 + start_orb%vec(6))

alpha = end_orb%vec(3) - spline_y(0) - spline_y(1) * ds
beta = end_orb%direction * end_orb%vec(4) / (1 + end_orb%vec(6)) - spline_y(1)

spline_y(2) = 3 * alpha / ds**2 - beta / ds
spline_y(3) = beta / ds**2 - 2 * alpha / ds**3


end subroutine spline_fit_orbit


