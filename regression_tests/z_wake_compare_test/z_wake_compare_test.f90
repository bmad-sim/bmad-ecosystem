program z_wake_compare_test

! Regression test: the z_long table wake and the pseudomode sr longitudinal
! wake must give identical particle kicks when the table is constructed to
! match the analytic mode wake function.
!
! Mode wake:  W(z) = A*sin(2*pi*phi + k*z)*exp(damp*z) for z <= 0,  0 for z > 0.
!
! This test uses the simplest case: damp = 0, k = 0, phi = 0.25 so that
! W(z) = A (constant) for z <= 0. With the fundamental theorem of beam
! loading correction (W(0) -> W(0)/2 in the convolution kernel) the two
! methods must agree to machine precision.
!
! The x_leading / y_leading cases compare the same constant wake with the
! deposited charge weighted by the leading-particle transverse offset.
! These guard the moment-weighted binning in sr_z_long_wake (a y_leading
! bin-accumulation bug was fixed there in 2026). The offset bunches use
! z values on the table grid (the mode/table comparison is only exact for
! grid-aligned particles) but in ADJACENT bins: the historic bug read the
! neighboring bin while depositing, which is only detectable when that
! bin is already occupied. With a constant wake the kicks depend only on
! the z ordering, so the expected values match the simple leading-sum
! formula regardless of spacing.

use beam_mod
use bmad
use wake_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_mode, ele_tab, ele_tab_ss
type (ele_struct), pointer :: ele_mode_xl, ele_tab_xl, ele_mode_yl, ele_tab_yl
type (bunch_struct), target :: bunch_mode, bunch_tab, bunch_tab_ss
type (bunch_struct), target :: bunch_mode_xl, bunch_tab_xl, bunch_mode_yl, bunch_tab_yl
type (coord_struct), pointer :: p

real(rp) :: zlist(5), zlist_off(5), xoff(5), yoff(5), q
real(rp) :: dpz_mode(5), dpz_tab(5), dpz_tab_ss(5)
real(rp) :: dpz_mode_xl(5), dpz_tab_xl(5), dpz_mode_yl(5), dpz_tab_yl(5)
real(rp) :: diff(5), diff_ss(5), diff_xl(5), diff_yl(5)
integer :: i, n, nargs

character(200) :: lat_file

! Lattice ----------------------------------------------------------------

lat_file = 'z_wake_compare_test.bmad'
nargs = command_argument_count()
if (nargs == 1) call get_command_argument(1, lat_file)

call bmad_parser (lat_file, lat)

ele_mode    => lat%ele(1)
ele_tab     => lat%ele(2)
ele_tab_ss  => lat%ele(3)
ele_mode_xl => lat%ele(4)
ele_tab_xl  => lat%ele(5)
ele_mode_yl => lat%ele(6)
ele_tab_yl  => lat%ele(7)

q = 1.0e-12_rp

! Bunch ------------------------------------------------------------------
! Convention: large z = head of bunch (leader), small z = tail (trailer).

zlist = [40e-6_rp, 20e-6_rp, 0.0_rp, -20e-6_rp, -40e-6_rp]
n = size(zlist)

! Transverse offsets for the x_leading / y_leading checks. Mixed signs so
! that a binning error cannot cancel. The z values keep the same ordering
! as zlist and are grid-aligned (multiples of the 10 um table spacing with
! a grid-aligned mean), but the leading three particles occupy adjacent
! bins so the deposit loop encounters already-occupied neighbor bins.

zlist_off = [30e-6_rp, 20e-6_rp, 10e-6_rp, -20e-6_rp, -40e-6_rp]
xoff = [0.7e-3_rp, -1.1e-3_rp, 0.4e-3_rp,  1.6e-3_rp, -0.5e-3_rp]
yoff = [-0.9e-3_rp, 0.6e-3_rp, 1.3e-3_rp, -0.2e-3_rp,  0.8e-3_rp]

call reallocate_bunch(bunch_mode,   n)
call reallocate_bunch(bunch_tab,    n)
call reallocate_bunch(bunch_tab_ss, n)
call reallocate_bunch(bunch_mode_xl, n)
call reallocate_bunch(bunch_tab_xl,  n)
call reallocate_bunch(bunch_mode_yl, n)
call reallocate_bunch(bunch_tab_yl,  n)

do i = 1, n
  call init_coord(bunch_mode%particle(i), [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, zlist(i), 0.0_rp], &
                  lat%ele(0), element_end = downstream_end$)
  bunch_mode%particle(i)%charge = q
  bunch_tab%particle(i)    = bunch_mode%particle(i)
  bunch_tab_ss%particle(i) = bunch_mode%particle(i)

  call init_coord(bunch_mode_xl%particle(i), [xoff(i), 0.0_rp, yoff(i), 0.0_rp, zlist_off(i), 0.0_rp], &
                  lat%ele(0), element_end = downstream_end$)
  bunch_mode_xl%particle(i)%charge = q
  bunch_tab_xl%particle(i)  = bunch_mode_xl%particle(i)
  bunch_mode_yl%particle(i) = bunch_mode_xl%particle(i)
  bunch_tab_yl%particle(i)  = bunch_mode_xl%particle(i)
enddo
bunch_mode%n_live   = n
bunch_tab%n_live    = n
bunch_tab_ss%n_live = n
bunch_mode%charge_live   = n*q
bunch_tab%charge_live    = n*q
bunch_tab_ss%charge_live = n*q
bunch_mode%ix_z   = 0
bunch_tab%ix_z    = 0
bunch_tab_ss%ix_z = 0
bunch_mode_xl%n_live = n;  bunch_mode_xl%charge_live = n*q;  bunch_mode_xl%ix_z = 0
bunch_tab_xl%n_live  = n;  bunch_tab_xl%charge_live  = n*q;  bunch_tab_xl%ix_z  = 0
bunch_mode_yl%n_live = n;  bunch_mode_yl%charge_live = n*q;  bunch_mode_yl%ix_z = 0
bunch_tab_yl%n_live  = n;  bunch_tab_yl%charge_live  = n*q;  bunch_tab_yl%ix_z  = 0

call order_particles_in_z(bunch_mode)
call order_particles_in_z(bunch_tab)
call order_particles_in_z(bunch_tab_ss)
call order_particles_in_z(bunch_mode_xl)
call order_particles_in_z(bunch_tab_xl)
call order_particles_in_z(bunch_mode_yl)
call order_particles_in_z(bunch_tab_yl)

! Apply wakes ------------------------------------------------------------

bmad_com%sr_wakes_on = .true.

call track1_sr_wake(bunch_mode,   ele_mode)
call track1_sr_wake(bunch_tab,    ele_tab)
call track1_sr_wake(bunch_tab_ss, ele_tab_ss)
call track1_sr_wake(bunch_mode_xl, ele_mode_xl)
call track1_sr_wake(bunch_tab_xl,  ele_tab_xl)
call track1_sr_wake(bunch_mode_yl, ele_mode_yl)
call track1_sr_wake(bunch_tab_yl,  ele_tab_yl)

do i = 1, n
  p => bunch_mode%particle(i);   dpz_mode(i)   = p%vec(6)
  p => bunch_tab%particle(i);    dpz_tab(i)    = p%vec(6)
  p => bunch_tab_ss%particle(i); dpz_tab_ss(i) = p%vec(6)
  diff(i)    = dpz_tab(i)    - dpz_mode(i)
  diff_ss(i) = dpz_tab_ss(i) - dpz_mode(i)

  p => bunch_mode_xl%particle(i); dpz_mode_xl(i) = p%vec(6)
  p => bunch_tab_xl%particle(i);  dpz_tab_xl(i)  = p%vec(6)
  p => bunch_mode_yl%particle(i); dpz_mode_yl(i) = p%vec(6)
  p => bunch_tab_yl%particle(i);  dpz_tab_yl(i)  = p%vec(6)
  diff_xl(i) = dpz_tab_xl(i) - dpz_mode_xl(i)
  diff_yl(i) = dpz_tab_yl(i) - dpz_mode_yl(i)
enddo

! Echo to stdout for human reading.

print '(a)', ' '
print '(a)', '   i     z [m]        pseudomode dpz       z_long dpz         diff      z_long(small_sig)    diff_ss'
do i = 1, n
  print '(i4, es14.6, 2es20.10, es12.2, es20.10, es12.2)', &
        i, zlist(i), dpz_mode(i), dpz_tab(i), diff(i), dpz_tab_ss(i), diff_ss(i)
enddo
print '(a, es14.6)', 'max |diff|    = ', maxval(abs(diff))
print '(a, es14.6)', 'max |diff_ss| = ', maxval(abs(diff_ss))

print '(a)', ' '
print '(a)', '   i     z [m]      mode dpz (x_lead)   z_long dpz (x_lead)   diff_xl   mode dpz (y_lead)   z_long dpz (y_lead)   diff_yl'
do i = 1, n
  print '(i4, es14.6, 2es20.10, es12.2, 2es20.10, es12.2)', &
        i, zlist_off(i), dpz_mode_xl(i), dpz_tab_xl(i), diff_xl(i), dpz_mode_yl(i), dpz_tab_yl(i), diff_yl(i)
enddo
print '(a, es14.6)', 'max |diff_xl| = ', maxval(abs(diff_xl))
print '(a, es14.6)', 'max |diff_yl| = ', maxval(abs(diff_yl))

! Regression test output.

open (1, file = 'output.now')

do i = 1, n
  write (1, '(a, i0, a, es20.10)') '"dpz-mode-',   i, '" REL 1E-10  ', dpz_mode(i)
  write (1, '(a, i0, a, es20.10)') '"dpz-tab-',    i, '" REL 1E-10  ', dpz_tab(i)
  write (1, '(a, i0, a, es20.10)') '"dpz-tab-ss-', i, '" REL 1E-10  ', dpz_tab_ss(i)
enddo

do i = 1, n
  write (1, '(a, i0, a, es20.10)') '"dpz-mode-xl-', i, '" REL 1E-10  ', dpz_mode_xl(i)
  write (1, '(a, i0, a, es20.10)') '"dpz-tab-xl-',  i, '" REL 1E-10  ', dpz_tab_xl(i)
  write (1, '(a, i0, a, es20.10)') '"dpz-mode-yl-', i, '" REL 1E-10  ', dpz_mode_yl(i)
  write (1, '(a, i0, a, es20.10)') '"dpz-tab-yl-',  i, '" REL 1E-10  ', dpz_tab_yl(i)
enddo

! The two methods must agree well below the wake kick magnitude.
do i = 1, n
  write (1, '(a, i0, a, es20.10)') '"diff-',    i, '" ABS 1E-15  ', diff(i)
  write (1, '(a, i0, a, es20.10)') '"diff-ss-', i, '" ABS 1E-15  ', diff_ss(i)
  write (1, '(a, i0, a, es20.10)') '"diff-xl-', i, '" ABS 1E-15  ', diff_xl(i)
  write (1, '(a, i0, a, es20.10)') '"diff-yl-', i, '" ABS 1E-15  ', diff_yl(i)
enddo

close (1)

end program z_wake_compare_test
