!+
! Benchmark: LSC (Longitudinal Space Charge) kick computation
!
! Exercises the vectorized lsc_kick_params_calc through the full
! track1_bunch_csr tracking path. Uses a simple bend with varying
! n_bin to demonstrate scaling.
!
! The benchmark sweeps n_bin = 40, 80, 160, 320 and reports timing
! for each. Larger n_bin values emphasize the O(n_bin^2) LSC cost
! where vectorization has the biggest impact.
!-

program perf_lsc

use beam_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (bunch_struct) :: bunch_init, bunch
type (coord_struct), allocatable :: centroid(:)
type (beam_init_struct) :: beam_init
type (branch_struct), pointer :: branch

real(rp) :: t_start, t_end, t_total
integer :: i, iter, n_iter
integer :: n_bin_vals(6), ib
logical :: err

!

n_bin_vals = [40, 80, 160, 320, 640, 1280]
n_iter = 5

call bmad_parser ('perf_lsc.bmad', lat)

branch => lat%branch(0)
bmad_com%csr_and_space_charge_on = .true.

call ran_seed_put(1)
beam_init%n_particle = 2000
beam_init%a_norm_emit = 1.0e-6
beam_init%b_norm_emit = 1.0e-6
beam_init%dPz_dz = 0.0
beam_init%bunch_charge = 77.0e-12
beam_init%sig_pz = 0e-9
beam_init%sig_z = 0.000899377  ! 3 ps * cLight

call init_bunch_distribution (branch%ele(0), branch%param, beam_init, 0, bunch_init)

call reallocate_coord(centroid, lat, 0)
call init_coord (centroid(0), lat%particle_start, branch%ele(0), downstream_end$)
call track_all (lat, centroid, 0)

ele => branch%ele(1)
ele%csr_method = off$
ele%space_charge_method = slice$

open (1, file = 'output.now')
write(*, '(a)') '=== LSC Benchmark ==='
write(*, '(a, i6)') 'Particles: ', beam_init%n_particle
write(*, '(a, i6)') 'Iterations:', n_iter
write(*, '(a)') ''

do ib = 1, size(n_bin_vals)
  space_charge_com%n_bin = n_bin_vals(ib)
  space_charge_com%ds_track_step = 0.01
  space_charge_com%particle_bin_span = 2
  space_charge_com%lsc_sigma_cutoff = 0.1
  space_charge_com%beam_chamber_height = 0.0254
  space_charge_com%n_shield_images = 0

  ! Warm-up
  bunch = bunch_init
  call track1_bunch(bunch, ele, err, centroid)

  ! Timed iterations
  call cpu_time(t_start)
  do iter = 1, n_iter
    bunch = bunch_init
    ! Reset kick_lsc for clean timing
    call track1_bunch(bunch, ele, err, centroid)
  enddo
  call cpu_time(t_end)

  t_total = (t_end - t_start) / n_iter

  write(*, '(a, i4, a, f10.4, a)') 'n_bin=', n_bin_vals(ib), '  time=', t_total * 1000, ' ms/track'
enddo

write(*, '(a)') ''
write(*, '(a)') 'Benchmark complete.'

write (1, '(a)') '"benchmark_complete" STR "T"'
close (1)

end program
