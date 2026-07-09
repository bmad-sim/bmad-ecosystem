!+
! Subroutine gen_grad_at_s_to_gg_a_taylor (ele, gen_grad, s_pos, gg_a, gg_da_ds)
!
! Routine to return the equivalent Taylor map for the magnetic (or electric) vector potential
! A = (Ax, Ay, As), along with its s-derivative dA/ds, at a point s_pos. The monomial
! expansion coefficients are obtained from the generalized-gradient coefficient table
! (gen_gradients_mod / gg_coef_table_mod), which includes the reference-frame curvature (g_ref).
!
! This is the vector-potential analog of gen_grad_at_s_to_gg_taylor (which returns the field
! B = (Bx, By, Bs)). Since gg_a holds A as a polynomial in (x, y), the transverse derivatives
! dA/dx and dA/dy are obtained by differentiating gg_a directly. The s-derivative dA/ds cannot
! be recovered from gg_a and so is returned separately in gg_da_ds.
!
! Input:
!   ele           -- ele_struct: Element containing the map.
!   gen_grad      -- gen_gradients_struct: Gen_gradients map.
!   s_pos         -- real(rp): Position to evaluate the map at.
!
! Output:
!   gg_a(3)          -- gg_taylor_struct: Map for the (Ax, Ay, As) vector potential.
!   gg_da_ds(3)      -- gg_taylor_struct, optional: Map for the s-derivative (dAx/ds, dAy/ds, dAs/ds).
!-

subroutine gen_grad_at_s_to_gg_a_taylor (ele, gen_grad, s_pos, gg_a, gg_da_ds)

use bmad_interface, dummy => gen_grad_at_s_to_gg_a_taylor
use gen_gradients_mod

implicit none

type (ele_struct) ele
type (gen_gradients_struct), target :: gen_grad
type (gg_taylor_struct), target :: gg_a(3)
type (gg_taylor_struct), optional, target :: gg_da_ds(3)
type (gen_grad_curve_struct), pointer :: crv

logical do_ds

real(rp) s_pos, s0, scale, z_rel, s_here
real(rp) :: gval(3, 0:gg_coef_max_n$, 0:gg_coef_max_m$+1)
real(rp) :: kcoef(6, 0:gg_coef_max_pq$, 0:gg_coef_max_pq$)
real(rp) :: kbump(4:6, 0:gg_coef_max_pq$, 0:gg_coef_max_pq$)

integer iz0, i, ip, iq, m, n, comp
character(*), parameter :: r_name = 'gen_grad_at_s_to_gg_a_taylor'

do_ds = present(gg_da_ds)

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

call gg_accumulate_coefs (gen_grad%g_ref, gval, .false., .true., kcoef, kbump)

! Emit the vector-potential monomials (Ax, Ay, As) as Taylor terms coef * x^ip * y^iq, and
! likewise the s-derivative (dA/ds) from the bumped-order coefficient array.

scale = gen_grad%field_scale * master_parameter_value(gen_grad%master_parameter, ele)

do i = 1, 3
  comp = i + 3   ! A components are 4:6 in kcoef/kbump.

  ! Vector potential A.

  n = count(kcoef(comp,:,:) /= 0)
  call init_gg_taylor_series(gg_a(i), n)

  n = 0
  do ip = 0, gg_coef_max_pq$      ! x power
    do iq = 0, gg_coef_max_pq$    ! y power
      if (kcoef(comp,ip,iq) == 0) cycle
      n = n + 1
      gg_a(i)%term(n)%coef = kcoef(comp,ip,iq) * scale
      gg_a(i)%term(n)%expn = [ip, iq]
    enddo
  enddo

  ! s-derivative dA/ds.

  if (.not. do_ds) cycle

  n = count(kbump(comp,:,:) /= 0)
  call init_gg_taylor_series(gg_da_ds(i), n)

  n = 0
  do ip = 0, gg_coef_max_pq$      ! x power
    do iq = 0, gg_coef_max_pq$    ! y power
      if (kbump(comp,ip,iq) == 0) cycle
      n = n + 1
      gg_da_ds(i)%term(n)%coef = kbump(comp,ip,iq) * scale
      gg_da_ds(i)%term(n)%expn = [ip, iq]
    enddo
  enddo
enddo

end subroutine
