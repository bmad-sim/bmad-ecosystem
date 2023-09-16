!+
! Subroutine quad_mat2_calc (k1, length, rel_p, mat2, z_coef, dz_dpz_coef)
!
! Subroutine to calculate the 2x2 transfer matrix for a quad for one plane. 
!
! Input:
!   k1       -- Real(rp): Quad strength: k1 > 0 ==> defocus.
!   length   -- Real(rp): Quad length
!   rel_p    -- Real(rp), Relative momentum P/P0.      
!
! Output:
!   mat2(2,2)      -- Real(rp): Transfer matrix.
!   z_coef(3)      -- Real(rp), optional: Coefficients for calculating the
!                       the change in z position:
!                          z = Integral [-(px/(1+pz))^2/2 ds]
!                            = c(1) * x_0^2 + c(2) * x_0 * px_0 + c(3) * px_0^2 
!   dz_dpz_coef(3) -- Real(rp), optional: Coefficients for calculating the
!                       the mat6(5,6) Jacobian matrix element:
!                         dz_dpz = c(1) * x_0^2 + c(2) * x_0 * px_0 + c(3) * px_0^2 
!-

subroutine quad_mat2_calc (k1, length, rel_p, mat2, z_coef, dz_dpz_coef)

use bmad_struct

implicit none

real(rp) length, mat2(:,:), cx, sx
real(rp) k1, sqrt_k, sk_l, k_l2, zc(3), dsx, dcx, rel_p
real(rp), optional :: z_coef(3), dz_dpz_coef(3)

!

sqrt_k = sqrt(abs(k1))
sk_l = sqrt_k * length

if (abs(sk_l) < 1d-10) then
  k_l2 = k1 * length**2
  cx = 1 + k_l2 / 2
  sx = (1 + k_l2 / 6) * length
elseif (k1 < 0) then       ! focus
  cx = cos(sk_l)
  sx = sin(sk_l) / sqrt_k
else                       ! defocus
  cx = cosh(sk_l)
  sx = sinh(sk_l) / sqrt_k
endif

mat2(1,1) = cx
mat2(1,2) = sx / rel_p
mat2(2,1) = k1 * sx * rel_p
mat2(2,2) = cx

!

if (present(z_coef) .or. present(dz_dpz_coef)) then
  zc(1) = k1 * (-cx * sx + length) / 4
  zc(2) = -k1 * sx**2 / (2 * rel_p)
  zc(3) = -(cx * sx + length) / (4 * rel_p**2)
  if (present(z_coef)) z_coef = zc
endif

! dz_dpz_coef

if (present(dz_dpz_coef)) then

  if (abs(sk_l) < 1d-10) then
    dcx = -k_l2 / (2 * rel_p)
    dsx = -k_l2 * length / (6 * rel_p)
  else
    dcx = -k1 * sx * length / (2 * rel_p)
    dsx = (sx - length * cx) / (2 * rel_p)
  endif

  dz_dpz_coef(1) =   -zc(1)/rel_p - k1 * (cx * dsx + dcx * sx) / 4
  dz_dpz_coef(2) = -2*zc(2)/rel_p - k1 * sx * dsx/rel_p
  dz_dpz_coef(3) = -2*zc(3)/rel_p - (cx * dsx + dcx * sx) / (4 * rel_p**2)

endif

end subroutine quad_mat2_calc

