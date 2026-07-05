!+
! Module gen_gradients_mod
!
! Routines to evaluate the magnetic field B = (Bx,By,Bs), the vector potential
! A = (Ax,Ay,As), and the Jacobian dA(i,j) = dA_i/du_j (u = x,y,s) from a set of
! generalized-gradient (GG) derivative values, in a straight OR constant-curvature
! (g_ref = 1/rho) reference frame. The monomial coefficient table lives in
! gg_coef_table_mod and is ported from GeneralizedGradients.jl (issue #2086).
!
! The low-level kernel gg_field_potential_calc is decoupled from any lattice
! structure: it takes the GG values already interpolated to the evaluation
! s-plane. Higher-level routines (in em_field_calc, etc.) build that value array
! from a gen_gradients_struct map and call the kernel.
!-

module gen_gradients_mod

use gg_coef_table_mod

implicit none

contains

!+
! Subroutine gg_field_potential_calc (g_ref, gval, x, y, calc_field, calc_potential, B, A, dA)
!
! Evaluate B, A, and dA at transverse position (x, y) from GG derivative values.
!
! Each table entry contributes  comp += coef * g_ref**k * x**p * y**q * G(kind,n,m)
! to output component `comp` (1=Bx 2=By 3=Bs 4=Ax 5=Ay 6=As). The (x,y) partials
! of A come from the monomial partials; the d/ds partial comes from bumping the GG
! derivative order (G(kind,n,m) -> G(kind,n,m+1)), so gval must be supplied up to
! order gg_coef_max_m$+1.
!
! Input:
!   g_ref           -- real(rp): Reference-frame curvature 1/rho (0 => straight).
!   gval(:,:,:)     -- real(rp): GG derivative values, indexed
!                         (kind = gg_a$/gg_b$/gg_bs$, n = 0:gg_coef_max_n$, m = 0:gg_coef_max_m$+1).
!                         Absent orders/harmonics are zero. For gg_bs$ use n = 0.
!   x, y            -- real(rp): Transverse position relative to the GG expansion axis.
!   calc_field      -- logical: Compute B?
!   calc_potential  -- logical: Compute A and dA?
!
! Output:
!   B(3)            -- real(rp): (Bx, By, Bs) if calc_field, else 0.
!   A(3)            -- real(rp): (Ax, Ay, As) if calc_potential, else 0.
!   dA(3,3)         -- real(rp): dA(i,j) = dA_i/du_j, (u1,u2,u3)=(x,y,s), if calc_potential, else 0.
!-

subroutine gg_field_potential_calc (g_ref, gval, x, y, calc_field, calc_potential, B, A, dA)

real(rp), intent(in) :: g_ref, x, y
real(rp), intent(in) :: gval(3, 0:gg_coef_max_n$, 0:gg_coef_max_m$+1)
logical, intent(in) :: calc_field, calc_potential
real(rp), intent(out) :: B(3), A(3), dA(3,3)

real(rp) :: kcoef(6, 0:gg_coef_max_pq$, 0:gg_coef_max_pq$)   ! (x,y) coef array per comp
real(rp) :: kbump(4:6, 0:gg_coef_max_pq$, 0:gg_coef_max_pq$) ! bumped-order coefs for dA/ds
real(rp) :: xp(0:gg_coef_max_pq$), yq(0:gg_coef_max_pq$)
real(rp) :: c, val, dvx, dvy
integer comp, i, j

!

B = 0; A = 0; dA = 0
if (.not. calc_field .and. .not. calc_potential) return

call gg_accumulate_coefs (g_ref, gval, calc_field, calc_potential, kcoef, kbump)

! Powers of x, y.

xp(0) = 1; yq(0) = 1
do i = 1, gg_coef_max_pq$
  xp(i) = xp(i-1) * x
  yq(i) = yq(i-1) * y
enddo

! Field: value only.

if (calc_field) then
  do comp = 1, 3
    val = 0
    do i = 0, gg_coef_max_pq$
      do j = 0, gg_coef_max_pq$
        if (kcoef(comp,i,j) /= 0) val = val + kcoef(comp,i,j) * xp(i) * yq(j)
      enddo
    enddo
    B(comp) = val
  enddo
endif

! Vector potential: value plus (x,y) partials from the same array, d/ds from bumped array.

if (calc_potential) then
  do comp = 4, 6
    val = 0; dvx = 0; dvy = 0
    do i = 0, gg_coef_max_pq$
      do j = 0, gg_coef_max_pq$
        c = kcoef(comp,i,j)
        if (c /= 0) then
          val = val + c * xp(i) * yq(j)
          if (i > 0) dvx = dvx + c * i * xp(i-1) * yq(j)
          if (j > 0) dvy = dvy + c * j * xp(i) * yq(j-1)
        endif
      enddo
    enddo
    A(comp-3) = val
    dA(comp-3, 1) = dvx
    dA(comp-3, 2) = dvy

    val = 0
    do i = 0, gg_coef_max_pq$
      do j = 0, gg_coef_max_pq$
        if (kbump(comp,i,j) /= 0) val = val + kbump(comp,i,j) * xp(i) * yq(j)
      enddo
    enddo
    dA(comp-3, 3) = val
  enddo
endif

end subroutine gg_field_potential_calc

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine gg_accumulate_coefs (g_ref, gval, do_field, do_potential, kcoef, kbump)
!
! Accumulate the monomial (x,y) coefficient arrays for the field and/or vector potential from a
! set of GG derivative values. For output component comp (1=Bx 2=By 3=Bs 4=Ax 5=Ay 6=As):
!   comp-polynomial = Sum_{p,q} kcoef(comp,p,q) * x^p * y^q
! and, for the A components (4:6), kbump holds the same sum with GG orders bumped by one (used to
! form dA/ds). This is the shared core of gg_field_potential_calc and the gg_taylor builders.
!
! Input:
!   g_ref           -- real(rp): Reference-frame curvature 1/rho.
!   gval(:,:,:)     -- real(rp): GG derivative values (kind, n, m); m up to gg_coef_max_m$+1.
!   do_field        -- logical: Accumulate the B components (1:3)?
!   do_potential    -- logical: Accumulate the A components (4:6) and kbump?
!
! Output:
!   kcoef(6,0:,0:)  -- real(rp): Monomial coefficient array per component.
!   kbump(4:6,0:,0:)-- real(rp): Bumped-order (dA/ds) coefficient array for the A components.
!-

subroutine gg_accumulate_coefs (g_ref, gval, do_field, do_potential, kcoef, kbump)

real(rp), intent(in) :: g_ref
real(rp), intent(in) :: gval(3, 0:gg_coef_max_n$, 0:gg_coef_max_m$+1)
logical, intent(in) :: do_field, do_potential
real(rp), intent(out) :: kcoef(6, 0:gg_coef_max_pq$, 0:gg_coef_max_pq$)
real(rp), intent(out) :: kbump(4:6, 0:gg_coef_max_pq$, 0:gg_coef_max_pq$)

real(rp) :: gpow(0:gg_coef_max_k$), c, g
integer it, comp, knd, n, m, p, q, kk, i

!

call gg_coef_table_init()

kcoef = 0; kbump = 0

gpow(0) = 1
do i = 1, gg_coef_max_k$
  gpow(i) = gpow(i-1) * g_ref
enddo

do it = 1, gg_n_coef_term$
  comp = gg_coef_table(it)%comp
  if (comp <= 3 .and. .not. do_field) cycle
  if (comp >= 4 .and. .not. do_potential) cycle

  knd = gg_coef_table(it)%kind
  n   = gg_coef_table(it)%n
  m   = gg_coef_table(it)%m
  p   = gg_coef_table(it)%p
  q   = gg_coef_table(it)%q
  kk  = gg_coef_table(it)%k
  c   = gg_coef_table(it)%coef * gpow(kk)

  g = gval(knd, n, m)
  if (g /= 0) kcoef(comp, p, q) = kcoef(comp, p, q) + c * g

  if (comp >= 4) then
    g = gval(knd, n, m+1)
    if (g /= 0) kbump(comp, p, q) = kbump(comp, p, q) + c * g
  endif
enddo

end subroutine gg_accumulate_coefs

end module gen_gradients_mod
