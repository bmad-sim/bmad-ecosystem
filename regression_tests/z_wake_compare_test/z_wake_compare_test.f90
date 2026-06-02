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

use beam_mod
use bmad
use wake_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_mode, ele_tab, ele_tab_ss
type (bunch_struct), target :: bunch_mode, bunch_tab, bunch_tab_ss
type (coord_struct), pointer :: p

real(rp) :: zlist(5), p0c, q
real(rp) :: dpz_mode(5), dpz_tab(5), dpz_tab_ss(5)
real(rp) :: diff(5), diff_ss(5)
integer :: i, n, nargs

character(200) :: lat_file

! Lattice ----------------------------------------------------------------

lat_file = 'z_wake_compare_test.bmad'
nargs = command_argument_count()
if (nargs == 1) call get_command_argument(1, lat_file)

call bmad_parser (lat_file, lat)

ele_mode   => lat%ele(1)
ele_tab    => lat%ele(2)
ele_tab_ss => lat%ele(3)

p0c = lat%ele(0)%value(p0c$)
q = 1.0e-12_rp

! Bunch ------------------------------------------------------------------
! Convention: large z = head of bunch (leader), small z = tail (trailer).

zlist = [40e-6_rp, 20e-6_rp, 0.0_rp, -20e-6_rp, -40e-6_rp]
n = size(zlist)

call reallocate_bunch(bunch_mode,   n)
call reallocate_bunch(bunch_tab,    n)
call reallocate_bunch(bunch_tab_ss, n)

do i = 1, n
  call init_coord(bunch_mode%particle(i), [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, zlist(i), 0.0_rp], &
                  lat%ele(0), element_end = downstream_end$)
  bunch_mode%particle(i)%charge = q
  bunch_tab%particle(i)    = bunch_mode%particle(i)
  bunch_tab_ss%particle(i) = bunch_mode%particle(i)
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

call order_particles_in_z(bunch_mode)
call order_particles_in_z(bunch_tab)
call order_particles_in_z(bunch_tab_ss)

! Apply wakes ------------------------------------------------------------

bmad_com%sr_wakes_on = .true.

call track1_sr_wake(bunch_mode,   ele_mode)
call track1_sr_wake(bunch_tab,    ele_tab)
call track1_sr_wake(bunch_tab_ss, ele_tab_ss)

do i = 1, n
  p => bunch_mode%particle(i);   dpz_mode(i)   = p%vec(6)
  p => bunch_tab%particle(i);    dpz_tab(i)    = p%vec(6)
  p => bunch_tab_ss%particle(i); dpz_tab_ss(i) = p%vec(6)
  diff(i)    = dpz_tab(i)    - dpz_mode(i)
  diff_ss(i) = dpz_tab_ss(i) - dpz_mode(i)
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

! Regression test output.

open (1, file = 'output.now')

do i = 1, n
  write (1, '(a, i0, a, es20.10)') '"dpz-mode-',   i, '" REL 1E-10  ', dpz_mode(i)
  write (1, '(a, i0, a, es20.10)') '"dpz-tab-',    i, '" REL 1E-10  ', dpz_tab(i)
  write (1, '(a, i0, a, es20.10)') '"dpz-tab-ss-', i, '" REL 1E-10  ', dpz_tab_ss(i)
enddo

! The two methods must agree well below the wake kick magnitude.
do i = 1, n
  write (1, '(a, i0, a, es20.10)') '"diff-',    i, '" ABS 1E-15  ', diff(i)
  write (1, '(a, i0, a, es20.10)') '"diff-ss-', i, '" ABS 1E-15  ', diff_ss(i)
enddo

close (1)

end program z_wake_compare_test
