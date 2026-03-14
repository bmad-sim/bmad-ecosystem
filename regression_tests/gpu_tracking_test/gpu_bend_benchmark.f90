!+
! gpu_bend_benchmark
!
! Quick 1M particle benchmark: CPU vs GPU through 6 bend elements
! exercising all GPU-accelerated bend parameters.
!-

program gpu_bend_benchmark

use bmad
use beam_mod
use gpu_tracking_mod

implicit none

type (lat_struct), target :: lat
type (beam_init_struct) :: beam_init
type (beam_struct) :: beam_cpu, beam_gpu
type (branch_struct), pointer :: branch

integer :: j, k, np
integer :: n_state_mismatch, n_lost_cpu, n_lost_gpu, n_alive
real(rp) :: max_diff, tol
real(rp) :: t_cpu_start, t_cpu_end, t_gpu_start, t_gpu_end
real(rp) :: dt_cpu, dt_gpu
logical :: err, pass

tol = 1d-12

! ---- Beam parameters ----
beam_init%n_particle = 1000000
beam_init%random_engine = 'quasi'
beam_init%a_emit = 1e-9
beam_init%b_emit = 1e-9
beam_init%dPz_dz = 0
beam_init%n_bunch = 1
beam_init%bunch_charge = 1e-9
beam_init%sig_pz = 1e-3
beam_init%sig_z = 1e-4
beam_init%random_sigma_cutoff = 4

! ---- Initialize GPU ----
call gpu_tracking_init()
if (.not. bmad_com%gpu_tracking_on) then
  print *, 'FATAL: GPU tracking not available'
  stop
endif

! ---- Warmup GPU ----
beam_init%n_particle = 10
call bmad_parser('lat_bend_only.bmad', lat)
branch => lat%branch(0)
call init_beam_distribution(branch%ele(0), branch%param, beam_init, beam_gpu, err)
call track_beam(lat, beam_gpu, err=err)

! ---- Parse bend benchmark lattice ----
beam_init%n_particle = 1000000
call bmad_parser('lat_bend_benchmark.bmad', lat)
branch => lat%branch(0)

print *
print *, '=================================================================='
print *, '  Bend Benchmark — 1M particles, 6 bend elements'
print *, '=================================================================='
print *

do j = 1, branch%n_ele_track
  print '(A,I3,A,A,A,A)', '   ', j, '. ', trim(key_name(branch%ele(j)%key)), &
    '  ', trim(branch%ele(j)%name)
enddo
print *

! ---- Initialize identical beams ----
print *, 'Initializing 1M particle beams...'
call init_beam_distribution(branch%ele(0), branch%param, beam_init, beam_cpu, err)
call init_beam_distribution(branch%ele(0), branch%param, beam_init, beam_gpu, err)
beam_gpu%bunch(1)%particle = beam_cpu%bunch(1)%particle
np = size(beam_cpu%bunch(1)%particle)
print '(A,I10,A)', '  ', np, ' particles initialized'
print *

! ---- CPU tracking ----
print *, 'Tracking CPU...'
bmad_com%gpu_tracking_on = .false.
call cpu_time(t_cpu_start)
call track_beam(lat, beam_cpu, err=err)
call cpu_time(t_cpu_end)
dt_cpu = t_cpu_end - t_cpu_start
print '(A,F10.3,A)', '  CPU time: ', dt_cpu, ' s'

! ---- GPU tracking ----
print *, 'Tracking GPU...'
bmad_com%gpu_tracking_on = .true.
call cpu_time(t_gpu_start)
call track_beam(lat, beam_gpu, err=err)
call cpu_time(t_gpu_end)
dt_gpu = t_gpu_end - t_gpu_start
print '(A,F10.3,A)', '  GPU time: ', dt_gpu, ' s'
print *

! ---- Compare results ----
n_state_mismatch = 0
n_lost_cpu = 0
n_lost_gpu = 0
n_alive = 0
max_diff = 0

do j = 1, np
  if (beam_cpu%bunch(1)%particle(j)%state /= alive$) n_lost_cpu = n_lost_cpu + 1
  if (beam_gpu%bunch(1)%particle(j)%state /= alive$) n_lost_gpu = n_lost_gpu + 1
  if (beam_cpu%bunch(1)%particle(j)%state /= beam_gpu%bunch(1)%particle(j)%state) then
    n_state_mismatch = n_state_mismatch + 1
  endif
  if (beam_cpu%bunch(1)%particle(j)%state == alive$ .and. &
      beam_gpu%bunch(1)%particle(j)%state == alive$) then
    n_alive = n_alive + 1
    do k = 1, 6
      max_diff = max(max_diff, abs(beam_cpu%bunch(1)%particle(j)%vec(k) &
                                  - beam_gpu%bunch(1)%particle(j)%vec(k)))
    enddo
    max_diff = max(max_diff, abs(beam_cpu%bunch(1)%particle(j)%t &
                                - beam_gpu%bunch(1)%particle(j)%t))
  endif
enddo

! ---- Report ----
print *, '=================================================================='
print *, '  Results'
print *, '=================================================================='
print '(A,I10)',    '  Total particles:      ', np
print '(A,I10)',    '  Lost (CPU):           ', n_lost_cpu
print '(A,I10)',    '  Lost (GPU):           ', n_lost_gpu
print '(A,I10)',    '  State mismatches:     ', n_state_mismatch
print '(A,I10)',    '  Alive in both:        ', n_alive
print '(A,ES12.4)', '  Max coord difference: ', max_diff
print '(A,ES12.4)', '  Tolerance:            ', tol
print *
print '(A,F10.3,A)', '  CPU time:  ', dt_cpu, ' s'
print '(A,F10.3,A)', '  GPU time:  ', dt_gpu, ' s'
if (dt_gpu > 0) then
  print '(A,F10.2,A)', '  Speedup:   ', dt_cpu / dt_gpu, 'x'
endif
print *

pass = (n_state_mismatch == 0) .and. (max_diff < tol)
if (pass) then
  print *, '  PASS'
else
  print *, '  FAIL'
endif
print *, '=================================================================='

if (.not. pass) stop 1

end program gpu_bend_benchmark
