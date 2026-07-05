!+
! Subroutine gen_grad_at_s_to_gg_taylor (ele, gen_grad, s_pos, gg_taylor)
!
! Routine to return the equivalent Taylor field map at a point s_pos. The monomial field-expansion
! coefficients are obtained from the generalized-gradient coefficient table (gen_gradients_mod /
! gg_coef_table_mod), which includes the reference-frame curvature (g_ref).
! Also see gen_grad1_to_gg_taylor.
!
! Input:
!   ele           -- ele_struct: Element containing the map.
!   gen_grad      -- gen_gradients_struct: Gen_gradients map.
!   s_pos         -- real(rp): Position to evaluate gg_taylor at.
!
! Output:
!   gg_taylor(3)  -- gg_taylor_struct: Map for (Bx, By, Bs) or (Ex, Ey, Es) fields.
!-

subroutine gen_grad_at_s_to_gg_taylor (ele, gen_grad, s_pos, gg_taylor)

use bmad_interface, dummy => gen_grad_at_s_to_gg_taylor
use gen_gradients_mod

implicit none

type (ele_struct) ele
type (gen_gradients_struct), target :: gen_grad
type (gg_taylor_struct), target :: gg_taylor(3)
type (gen_grad_curve_struct), pointer :: crv

real(rp) s_pos, s0, scale, z_rel, s_here
real(rp) :: gval(3, 0:gg_coef_max_n$, 0:gg_coef_max_m$+1)
real(rp) :: kcoef(6, 0:gg_coef_max_pq$, 0:gg_coef_max_pq$)
real(rp) :: kbump(4:6, 0:gg_coef_max_pq$, 0:gg_coef_max_pq$)

integer iz0, i, ip, iq, m, n
character(*), parameter :: r_name = 'gen_grad_at_s_to_gg_taylor'

! Find where to interpolate.

select case (gen_grad%ele_anchor_pt)
case (anchor_beginning$); s0 = 0
case (anchor_center$);    s0 = ele%value(l$) / 2
case (anchor_end$);       s0 = ele%value(l$)
end select

s_here = s_pos - s0 - gen_grad%r0(3)

iz0 = floor(s_here / gen_grad%dz)
if (iz0 < gen_grad%iz0) iz0 = iz0 + 1 ! Allow one dz width out-of-bounds.
if (iz0 >= gen_grad%iz1) iz0 = iz0 - 1 ! Allow one dz width out-of-bounds.

if (iz0 < gen_grad%iz0 .or. iz0 >= gen_grad%iz1) then
  call out_io (s_error$, r_name, 'PARTICLE Z  \F10.3\ POSITION OUT OF BOUNDS.', &
                                 'FOR GEN_GRADIENTS IN ELEMENT: ' // ele%name, r_array = [s_pos])
  return
endif

z_rel = s_here - iz0 * gen_grad%dz

! Interpolate each GG derivative tower to (iz0, z_rel) to build the value array.

gval = 0
do i = 1, size(gen_grad%curve)
  crv => gen_grad%curve(i)
  do m = 0, crv%m_max
    gval(crv%kind, crv%n, m) = poly_eval(crv%deriv(iz0, m:), z_rel, diff_coef=.true.)
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
