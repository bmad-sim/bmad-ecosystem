!+
! gpu_tracking_stress_test
!
! Comprehensive GPU tracking test: 1M particles through a 21-element lattice
! exercising every parameter that affects drift, quadrupole, bend, lcavity,
! and pipe tracking.
!
! Lattice (lat_kitchen_sink.bmad):
!   1. Plain drift                    10. Mid drift
!   2. Plain quad                     11. Plain bend
!   3. Drift with aperture            12. Bend with k1
!   4. Quad with misalignment         13. Bend fringe+misalign
!   5. Quad with fringe               14. Bend with multipoles
!   6. Quad with mag multipoles       15. Plain lcavity
!   7. Quad with elec multipoles      16. Lcavity misalign+phase
!   8. Quad with aperture             17. Plain pipe
!   9. Quad combo                     18. Pipe with misalignment
!                                     19. Pipe with aperture
!                                     20. Pipe combo (misalign+aperture)
!                                     21. Final drift
!
! Compares CPU vs GPU tracking of all 1M particles, reports:
!   - max coordinate difference (alive particles)
!   - particle state agreement (lost vs alive)
!   - wall-clock time for each
!-

program gpu_tracking_stress_test

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
real(rp) :: dt_cpu, dt_gpu
integer(8) :: clock_start, clock_end, clock_rate
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
  print *, 'FATAL: GPU tracking not available (set ACC_ENABLE_GPU_TRACKING=Y)'
  stop
endif

! ---- Warmup GPU (small beam to trigger CUDA context init) ----
beam_init%n_particle = 10
call bmad_parser('lat_drift_only.bmad', lat)
branch => lat%branch(0)
call init_beam_distribution(branch%ele(0), branch%param, beam_init, beam_gpu, err)
call track_beam(lat, beam_gpu, err=err)

! ---- Parse the kitchen sink lattice ----
beam_init%n_particle = 1000000
call bmad_parser('lat_kitchen_sink.bmad', lat)
branch => lat%branch(0)

print *
print *, '=================================================================='
print *, '  GPU Tracking Stress Test — 1M particles, 21 elements'
print *, '=================================================================='
print *

! Print lattice summary
print *, 'Lattice elements:'
do j = 1, branch%n_ele_track
  print '(A,I3,A,A,A,A)', '   ', j, '. ', trim(key_name(branch%ele(j)%key)), &
    '  ', trim(branch%ele(j)%name)
enddo
print *

! ---- Initialize identical beams ----
print *, 'Initializing 1M particle beams...'
call init_beam_distribution(branch%ele(0), branch%param, beam_init, beam_cpu, err)
if (err) then
  print *, 'FATAL: beam init failed'
  stop
endif
call init_beam_distribution(branch%ele(0), branch%param, beam_init, beam_gpu, err)
beam_gpu%bunch(1)%particle = beam_cpu%bunch(1)%particle

np = size(beam_cpu%bunch(1)%particle)
print '(A,I10,A)', '  ', np, ' particles initialized'
print *

! ---- CPU tracking ----
print *, 'Tracking CPU...'
bmad_com%gpu_tracking_on = .false.
call system_clock(clock_start, clock_rate)
call track_beam(lat, beam_cpu, err=err)
call system_clock(clock_end)
dt_cpu = real(clock_end - clock_start, rp) / real(clock_rate, rp)
print '(A,F10.3,A)', '  CPU time: ', dt_cpu, ' s'

! ---- GPU tracking ----
print *, 'Tracking GPU...'
bmad_com%gpu_tracking_on = .true.
call system_clock(clock_start, clock_rate)
call track_beam(lat, beam_gpu, err=err)
call system_clock(clock_end)
dt_gpu = real(clock_end - clock_start, rp) / real(clock_rate, rp)
print '(A,F10.3,A)', '  GPU time: ', dt_gpu, ' s'
print *

! ---- Compare results ----
print *, 'Comparing results...'

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

  ! Compare coordinates for particles alive in both
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
print *
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
  print *, '  PASS — GPU matches CPU'
else
  print *, '  FAIL — GPU does not match CPU'
  if (n_state_mismatch > 0) print '(A,I10,A)', '    ', n_state_mismatch, ' particles have different alive/lost state'
  if (max_diff >= tol) print '(A,ES12.4,A,ES12.4)', '    max_diff=', max_diff, ' exceeds tol=', tol
endif
print *, '=================================================================='

if (.not. pass) stop 1

end program gpu_tracking_stress_test
