!+
! Subroutine gen_grad1_to_gg_taylor (ele, gen_grad, iz, gg_taylor)
!
! Routine to take gen_gradients coefs at a given z-plane and return the equivalent Taylor field map.
! The monomial field-expansion coefficients are obtained from the generalized-gradient coefficient
! table (gen_gradients_mod / gg_coef_table_mod), which includes the reference-frame curvature (g_ref).
! Also see: gen_grad_at_s_to_gg_taylor.
!
! Input:
!   ele           -- ele_struct: Element containing the map.
!   gen_grad      -- gen_gradients_struct: Gen_gradients map.
!   iz            -- integer: z-plane index to evaluate.
!
! Output:
!   gg_taylor(3)  -- gg_taylor_struct: Map for (Bx, By, Bs) or (Ex, Ey, Es) fields.
!-

subroutine gen_grad1_to_gg_taylor (ele, gen_grad, iz, gg_taylor)

use bmad_interface, dummy => gen_grad1_to_gg_taylor
use gen_gradients_mod

implicit none

type (ele_struct) ele
type (gen_gradients_struct), target :: gen_grad
type (gg_taylor_struct), target :: gg_taylor(3)
type (gen_grad_curve_struct), pointer :: crv

real(rp) scale
real(rp) :: gval(3, 0:gg_coef_max_n$, 0:gg_coef_max_m$+1)
real(rp) :: kcoef(6, 0:gg_coef_max_pq$, 0:gg_coef_max_pq$)
real(rp) :: kbump(4:6, 0:gg_coef_max_pq$, 0:gg_coef_max_pq$)

integer iz, i, ip, iq, m, n

! Build the GG derivative-value array at the plane (z_rel = 0 => the stored derivatives at the plane).

gval = 0
do i = 1, size(gen_grad%curve)
  crv => gen_grad%curve(i)
  do m = 0, crv%m_max
    gval(crv%kind, crv%n, m) = crv%deriv(iz, m)
  enddo
enddo

call gg_accumulate_coefs (gen_grad%g_ref, gval, .true., .false., kcoef, kbump)

! Emit the field monomials (Bx, By, Bs) as Taylor terms coef * x^ip * y^iq.

scale = gen_grad%field_scale * master_parameter_value(gen_grad%master_parameter, ele)

do i = 1, 3
  n = count(kcoef(i,:,:) /= 0)
  call init_gg_taylor_series(gg_taylor(i), n)

  n = 0
  do ip = 0, gg_coef_max_pq$      ! x power
    do iq = 0, gg_coef_max_pq$    ! y power
      if (kcoef(i,ip,iq) == 0) cycle
      n = n + 1
      gg_taylor(i)%term(n)%coef = kcoef(i,ip,iq) * scale
      gg_taylor(i)%term(n)%expn = [ip, iq]
    enddo
  enddo
enddo

end subroutine
