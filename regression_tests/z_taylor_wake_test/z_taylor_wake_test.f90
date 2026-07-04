program z_taylor_wake_test

! Regression test for the Taylor-expanded 3D short-range wake (sr_wake z_taylor).
!
! Every term h_ab of the test lattice is a constant wake A_k = k*1e15 (k = term
! number) over z in [-100um, 0]. For a constant kernel and particles placed on
! the wake table grid points, the binned FFT convolution is exact, so the
! tracked kicks must match the analytic second-order Taylor sums:
!
!   dpz_i = -(q/p0c) * Sum_k c_k A_k T_k(i) * (S_k(i)/2 + Sum_{j ahead} S_k(j))
!   dpx_i = +(q/p0c) * Sum_k' c_k A_k U_k(i) * (Sum_{j ahead} S_k(j) s_ij + S_k(i) dz/4)
!
! where S_k is the source monomial (1, x, y, xy, x^2-y^2), T_k the witness
! monomial, U_k the witness-derivative monomial, c_k = 2 for cross terms with
! both indices transverse, s_ij = z_j - z_i, and the dz/4 self term is the
! discrete Panofsky-Wenzel integral of the halved z = 0 kernel bin.
! dpy is analogous. Particle z values are grid-aligned with adjacent occupied
! bins so bin-accumulation errors cannot hide (see z_wake_compare_test).

use beam_mod
use bmad
use wake_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (bunch_struct), target :: bunch
type (coord_struct), pointer :: p

integer, parameter :: np = 5
real(rp), parameter :: a0 = 1e15_rp

real(rp) :: zz(np), xx(np), yy(np), q, p0c, dz, amp(13)
real(rp) :: dpx(np), dpy(np), dpz(np), ex_pz(np), ex_px(np), ex_py(np)
real(rp) :: s_mono(np,5), sij, wsum, tsum
integer :: i, j, k, nargs

! Term k -> (source monomial index, witness code, multiplicity c_k).
! Source monomials: 1=const, 2=x, 3=y, 4=xy, 5=x^2-y^2.
! Witness code: 0=const, 1=x, 2=y, 3=xy, 4=x^2-y^2 (for the longitudinal kick).
integer, parameter :: k_src(13) = [1, 2, 3, 1, 1, 5, 4, 2, 2, 3, 3, 1, 1]
integer, parameter :: k_wit(13) = [0, 0, 0, 1, 2, 0, 0, 1, 2, 1, 2, 4, 3]
real(rp), parameter :: k_c(13) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 2]

character(200) :: lat_file

! Setup ----------------------------------------------------------------

lat_file = 'z_taylor_wake_test.bmad'
nargs = command_argument_count()
if (nargs == 1) call get_command_argument(1, lat_file)

call bmad_parser (lat_file, lat)
ele => lat%ele(1)
p0c = ele%value(p0c$)
dz = ele%wake%sr%z_taylor%dz

q = 1.0e-12_rp
do k = 1, 13
  amp(k) = k * a0
enddo

! Bunch: head first (largest z), grid-aligned (10um grid, grid-aligned mean),
! leading three particles in adjacent bins. Mixed-sign transverse offsets.

zz = [30e-6_rp, 20e-6_rp, 10e-6_rp, -20e-6_rp, -40e-6_rp]
xx = [0.7e-3_rp, -1.1e-3_rp, 0.4e-3_rp,  1.6e-3_rp, -0.5e-3_rp]
yy = [-0.9e-3_rp, 0.6e-3_rp, 1.3e-3_rp, -0.2e-3_rp,  0.8e-3_rp]

call reallocate_bunch(bunch, np)
do i = 1, np
  call init_coord(bunch%particle(i), [xx(i), 0.0_rp, yy(i), 0.0_rp, zz(i), 0.0_rp], &
                  lat%ele(0), element_end = downstream_end$)
  bunch%particle(i)%charge = q
enddo
bunch%n_live = np
bunch%charge_live = np*q
bunch%ix_z = 0
call order_particles_in_z(bunch)

! Track ----------------------------------------------------------------

bmad_com%sr_wakes_on = .true.
call track1_sr_wake(bunch, ele)

do i = 1, np
  p => bunch%particle(i)
  dpx(i) = p%vec(2)
  dpy(i) = p%vec(4)
  dpz(i) = p%vec(6)
enddo

! Analytic expectation -------------------------------------------------

do i = 1, np
  s_mono(i,:) = [1.0_rp, xx(i), yy(i), xx(i)*yy(i), xx(i)**2 - yy(i)**2]
enddo

ex_pz = 0;  ex_px = 0;  ex_py = 0

do i = 1, np
  do k = 1, 13
    ! Longitudinal: leading sum with half self term.
    wsum = 0.5_rp * s_mono(i, k_src(k))
    do j = 1, np
      if (zz(j) > zz(i)) wsum = wsum + s_mono(j, k_src(k))
    enddo
    select case (k_wit(k))
    case (0);  ex_pz(i) = ex_pz(i) - k_c(k) * amp(k) * wsum * q / p0c
    case (1);  ex_pz(i) = ex_pz(i) - k_c(k) * amp(k) * wsum * xx(i) * q / p0c
    case (2);  ex_pz(i) = ex_pz(i) - k_c(k) * amp(k) * wsum * yy(i) * q / p0c
    case (3);  ex_pz(i) = ex_pz(i) - k_c(k) * amp(k) * wsum * xx(i) * yy(i) * q / p0c
    case (4);  ex_pz(i) = ex_pz(i) - k_c(k) * amp(k) * wsum * (xx(i)**2 - yy(i)**2) * q / p0c
    end select

    ! Transverse: integrated (Panofsky-Wenzel) leading sum.
    if (k_wit(k) == 0) cycle
    tsum = s_mono(i, k_src(k)) * dz / 4
    do j = 1, np
      sij = zz(j) - zz(i)
      if (sij > 0) tsum = tsum + s_mono(j, k_src(k)) * sij
    enddo
    tsum = k_c(k) * amp(k) * tsum * q / p0c

    select case (k)
    case (sr_z_taylor_w03$, sr_z_taylor_w13$, sr_z_taylor_w23$)  ! witness x
      ex_px(i) = ex_px(i) + tsum
    case (sr_z_taylor_w04$, sr_z_taylor_w14$, sr_z_taylor_w24$)  ! witness y
      ex_py(i) = ex_py(i) + tsum
    case (sr_z_taylor_w34$)   ! witness x*y: kicks in both planes
      ex_px(i) = ex_px(i) + tsum * yy(i)
      ex_py(i) = ex_py(i) + tsum * xx(i)
    case (sr_z_taylor_w33$)   ! witness x^2 - y^2: quadrupole-like
      ex_px(i) = ex_px(i) + 2 * tsum * xx(i)
      ex_py(i) = ex_py(i) - 2 * tsum * yy(i)
    end select
  enddo
enddo

! Report ---------------------------------------------------------------

print '(a)', ' '
print '(a)', '   i      tracked dpz         exact dpz        tracked dpx         exact dpx        tracked dpy         exact dpy'
do i = 1, np
  print '(i4, 6es19.10)', i, dpz(i), ex_pz(i), dpx(i), ex_px(i), dpy(i), ex_py(i)
enddo
print '(a, es14.6)', 'max |dpz diff| = ', maxval(abs(dpz - ex_pz))
print '(a, es14.6)', 'max |dpx diff| = ', maxval(abs(dpx - ex_px))
print '(a, es14.6)', 'max |dpy diff| = ', maxval(abs(dpy - ex_py))

open (1, file = 'output.now')
do i = 1, np
  write (1, '(a, i0, a, es20.10)') '"dpz-', i, '" REL 1E-10  ', dpz(i)
  write (1, '(a, i0, a, es20.10)') '"dpx-', i, '" REL 1E-10  ', dpx(i)
  write (1, '(a, i0, a, es20.10)') '"dpy-', i, '" REL 1E-10  ', dpy(i)
enddo
do i = 1, np
  write (1, '(a, i0, a, es20.10)') '"diff-pz-', i, '" ABS 1E-14  ', dpz(i) - ex_pz(i)
  write (1, '(a, i0, a, es20.10)') '"diff-px-', i, '" ABS 1E-14  ', dpx(i) - ex_px(i)
  write (1, '(a, i0, a, es20.10)') '"diff-py-', i, '" ABS 1E-14  ', dpy(i) - ex_py(i)
enddo
close (1)

end program z_taylor_wake_test
