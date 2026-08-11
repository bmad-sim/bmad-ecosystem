!+
! Program exact_bend_edge_test
!
! Regression test for the Bmad-coordinate exact_bend_edge_kick.
!
! For a set of starting phase-space vectors, both edges, and both time
! directions, the program compares the Bmad-coordinate exact_bend_edge_kick
! against the reference PTC-coordinate exact_bend_edge_kick_ptc (orbit, arrival
! time, and transfer matrix), and checks that the new transfer matrix is
! symplectic. It also benchmarks the run time of the two implementations
! (printed to stdout, not part of the regression output).
!-

program exact_bend_edge_test

use bmad
use fringe_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct) orb0, orbn, orbo

real(rp) matn(6,6), mato(6,6), mat0(6,6)
real(rp) start(6,6), orb_d, mat_d, t_d, symp
real(rp) t1, t2

integer i, pa, td, j, n_iter
integer edge(2)

!

open (1, file = 'output.now')

call bmad_parser ('exact_bend_edge_test.bmad', lat)
ele => lat%ele(1)
call mat_make_unit (mat0)
edge = [first_track_edge$, second_track_edge$]

! Representative starting phase-space vectors.

start(:,1) = [ 0.010_rp,  0.020_rp,  0.030_rp,  0.040_rp,  0.000_rp,  0.500_rp]
start(:,2) = [-0.008_rp,  0.015_rp,  0.006_rp, -0.012_rp,  0.010_rp, -0.050_rp]
start(:,3) = [ 0.002_rp, -0.030_rp, -0.020_rp,  0.005_rp, -0.004_rp,  0.100_rp]
start(:,4) = [ 0.000_rp,  0.000_rp,  0.000_rp,  0.000_rp,  0.000_rp,  0.000_rp]
start(:,5) = [ 0.025_rp, -0.010_rp,  0.018_rp,  0.022_rp,  0.003_rp, -0.200_rp]
start(:,6) = [-0.015_rp,  0.028_rp, -0.026_rp, -0.019_rp,  0.007_rp,  0.300_rp]

orb_d = 0;  mat_d = 0;  t_d = 0;  symp = 0

do i = 1, 6
  do pa = 1, 2
    do td = -1, 1, 2
      call init_coord (orb0, lat%particle_start, ele, upstream_end$)
      orb0%vec = start(:,i)
      orb0%time_dir = td

      orbn = orb0;  matn = mat0
      call exact_bend_edge_kick (ele, lat%param, edge(pa), orbn, matn, .true.)

      orbo = orb0;  mato = mat0
      call exact_bend_edge_kick_ptc (ele, lat%param, edge(pa), orbo, mato, .true.)

      orb_d = max(orb_d, maxval(abs(orbn%vec - orbo%vec)))
      t_d   = max(t_d,   abs(orbn%t - orbo%t))
      mat_d = max(mat_d, maxval(abs(matn - mato)))
      symp  = max(symp,  mat_symp_error(matn))
    enddo
  enddo
enddo

write (1, '(a, es14.6)') '"orbit-new-minus-old"  ABS 1E-11', orb_d
write (1, '(a, es14.6)') '"matrix-new-minus-old" ABS 1E-10', mat_d
write (1, '(a, es14.6)') '"time-new-minus-old"   ABS 1E-16', t_d
write (1, '(a, es14.6)') '"new-symplectic-error" ABS 1E-9 ', symp

close (1)

! Benchmark (printed to stdout only). Skipped under the regression harness (which
! runs the program with no arguments); run e.g. "exact_bend_edge_test bench" to time.

if (command_argument_count() == 0) stop

n_iter = 2000000

call init_coord (orb0, lat%particle_start, ele, upstream_end$)
orb0%vec = start(:,1)

call cpu_time (t1)
do j = 1, n_iter
  orbn = orb0
  call exact_bend_edge_kick (ele, lat%param, first_track_edge$, orbn)
enddo
call cpu_time (t2)
print '(a, f9.2, a)', 'Bmad exact_bend_edge_kick      (orbit)  : ', 1e9_rp*(t2-t1)/n_iter, ' ns/call'

call cpu_time (t1)
do j = 1, n_iter
  orbo = orb0
  call exact_bend_edge_kick_ptc (ele, lat%param, first_track_edge$, orbo)
enddo
call cpu_time (t2)
print '(a, f9.2, a)', 'PTC  exact_bend_edge_kick_ptc  (orbit)  : ', 1e9_rp*(t2-t1)/n_iter, ' ns/call'

call cpu_time (t1)
do j = 1, n_iter
  orbn = orb0;  matn = mat0
  call exact_bend_edge_kick (ele, lat%param, first_track_edge$, orbn, matn, .true.)
enddo
call cpu_time (t2)
print '(a, f9.2, a)', 'Bmad exact_bend_edge_kick      (matrix) : ', 1e9_rp*(t2-t1)/n_iter, ' ns/call'

call cpu_time (t1)
do j = 1, n_iter
  orbo = orb0;  mato = mat0
  call exact_bend_edge_kick_ptc (ele, lat%param, first_track_edge$, orbo, mato, .true.)
enddo
call cpu_time (t2)
print '(a, f9.2, a)', 'PTC  exact_bend_edge_kick_ptc  (matrix) : ', 1e9_rp*(t2-t1)/n_iter, ' ns/call'

end program
