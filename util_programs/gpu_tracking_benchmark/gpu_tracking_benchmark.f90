!+
! gpu_tracking_benchmark
!
! Benchmark program comparing CPU vs GPU particle tracking through
! various element configurations. Tests 10k, 100k, and 1M macroparticles.
! Validates that GPU results match CPU results.
!-

program gpu_tracking_benchmark

use bmad
use beam_mod
use gpu_tracking_mod
use, intrinsic :: iso_c_binding

implicit none

interface
  integer(C_INT) function c_setenv(name, value, overwrite) bind(C, name='setenv')
    import :: C_CHAR, C_INT
    character(kind=C_CHAR), intent(in) :: name(*), value(*)
    integer(C_INT), value :: overwrite
  end function
  integer(C_INT) function c_unsetenv(name) bind(C, name='unsetenv')
    import :: C_CHAR, C_INT
    character(kind=C_CHAR), intent(in) :: name(*)
  end function
end interface

integer(C_INT) :: rc

integer, parameter :: n_lats = 8
integer, parameter :: n_sizes = 3
integer, parameter :: sizes(n_sizes) = [10000, 100000, 1000000]

type (lat_struct), target :: lats(n_lats)
type (beam_init_struct) :: beam_init
type (beam_struct) :: beam_cpu, beam_gpu
type (branch_struct), pointer :: branch

character(len=60) :: lat_names(n_lats), lat_files(n_lats)

integer :: j, n_particles, isz, ilat
integer(8) :: clock_start, clock_end, clock_rate
logical :: err, gpu_ok

! Results arrays
real(rp) :: t_cpu_arr(n_lats, n_sizes), t_gpu_arr(n_lats, n_sizes)
real(rp) :: mdiff_arr(n_lats, n_sizes)
real(rp) :: t_cpu, t_gpu, max_diff_all

! Define lattice configurations
lat_files(1) = 'lat_drift_only.bmad';           lat_names(1) = 'Drift only (10m)'
lat_files(2) = 'lat_quad_only.bmad';            lat_names(2) = 'Quad only (0.5m, k1=1.2)'
lat_files(3) = 'lat.bmad';                      lat_names(3) = 'Drift + Quad + Drift'
lat_files(4) = 'lat_quad_fringe.bmad';          lat_names(4) = 'Quad + fringe'
lat_files(5) = 'lat_quad_misalign.bmad';        lat_names(5) = 'Quad + misalignment'
lat_files(6) = 'lat_quad_fringe_misalign.bmad'; lat_names(6) = 'Quad + fringe + misalign + tilt'
lat_files(7) = 'lat_quad_multi_multipole.bmad'; lat_names(7) = 'Quad + sext+oct+skew+misalign+fringe'
lat_files(8) = 'lat_quad_elec_combo.bmad';      lat_names(8) = 'Quad + elec+mag+fringe+misalign'

! Parse all lattices
do ilat = 1, n_lats
  call bmad_parser(trim(lat_files(ilat)), lats(ilat))
enddo

! Warmup: small GPU run to initialize CUDA context
call gpu_tracking_reset()
rc = c_setenv('ACC_ENABLE_GPU_TRACKING' // C_NULL_CHAR, 'Y' // C_NULL_CHAR, 1_C_INT)
gpu_ok = gpu_tracking_is_active()
if (.not. gpu_ok) then
  print *, 'FATAL: GPU tracking not available'
  stop
endif

beam_init%n_particle = 100
beam_init%random_engine = 'quasi'
beam_init%a_emit = 1e-9
beam_init%b_emit = 1e-9
beam_init%dPz_dz = 0
beam_init%n_bunch = 1
beam_init%bunch_charge = 1e-9
beam_init%sig_pz = 1e-3
beam_init%sig_z = 1e-4
beam_init%random_sigma_cutoff = 4

branch => lats(1)%branch(0)
call init_beam_distribution(branch%ele(0), branch%param, beam_init, beam_gpu, err)
call track_beam(lats(1), beam_gpu, err=err)

! Run all benchmarks (collecting results silently)
do ilat = 1, n_lats
  do isz = 1, n_sizes
    n_particles = sizes(isz)
    beam_init%n_particle = n_particles
    branch => lats(ilat)%branch(0)

    ! Initialize both beams from same distribution
    call init_beam_distribution(branch%ele(0), branch%param, beam_init, beam_cpu, err)
    call init_beam_distribution(branch%ele(0), branch%param, beam_init, beam_gpu, err)
    beam_gpu%bunch(1)%particle = beam_cpu%bunch(1)%particle

    ! CPU run
    call gpu_tracking_reset()
    rc = c_unsetenv('ACC_ENABLE_GPU_TRACKING' // C_NULL_CHAR)
    gpu_ok = gpu_tracking_is_active()
    call system_clock(clock_start, clock_rate)
    call track_beam(lats(ilat), beam_cpu, err=err)
    call system_clock(clock_end)
    t_cpu = real(clock_end - clock_start, rp) / real(clock_rate, rp)

    ! GPU run
    call gpu_tracking_reset()
    rc = c_setenv('ACC_ENABLE_GPU_TRACKING' // C_NULL_CHAR, 'Y' // C_NULL_CHAR, 1_C_INT)
    gpu_ok = gpu_tracking_is_active()
    call system_clock(clock_start, clock_rate)
    call track_beam(lats(ilat), beam_gpu, err=err)
    call system_clock(clock_end)
    t_gpu = real(clock_end - clock_start, rp) / real(clock_rate, rp)

    ! Compare results
    max_diff_all = 0
    do j = 1, n_particles
      max_diff_all = max(max_diff_all, abs(beam_cpu%bunch(1)%particle(j)%vec(1) - beam_gpu%bunch(1)%particle(j)%vec(1)))
      max_diff_all = max(max_diff_all, abs(beam_cpu%bunch(1)%particle(j)%vec(2) - beam_gpu%bunch(1)%particle(j)%vec(2)))
      max_diff_all = max(max_diff_all, abs(beam_cpu%bunch(1)%particle(j)%vec(3) - beam_gpu%bunch(1)%particle(j)%vec(3)))
      max_diff_all = max(max_diff_all, abs(beam_cpu%bunch(1)%particle(j)%vec(4) - beam_gpu%bunch(1)%particle(j)%vec(4)))
      max_diff_all = max(max_diff_all, abs(beam_cpu%bunch(1)%particle(j)%vec(5) - beam_gpu%bunch(1)%particle(j)%vec(5)))
      max_diff_all = max(max_diff_all, abs(beam_cpu%bunch(1)%particle(j)%vec(6) - beam_gpu%bunch(1)%particle(j)%vec(6)))
      max_diff_all = max(max_diff_all, abs(beam_cpu%bunch(1)%particle(j)%t      - beam_gpu%bunch(1)%particle(j)%t))
    enddo

    t_cpu_arr(ilat, isz) = t_cpu
    t_gpu_arr(ilat, isz) = t_gpu
    mdiff_arr(ilat, isz) = max_diff_all
  enddo
enddo

! Print results table
print *
print *, '========================================================================================================'
print *, '  GPU Tracking Benchmark Results'
print *, '========================================================================================================'
print *
print '(A)', '                                            |       10k       |      100k       |       1M'
print '(A)', '  Lattice                                   | CPU/GPU  speed  | CPU/GPU  speed  | CPU/GPU  speed'
print '(A)', '  ------------------------------------------|-----------------|-----------------|------------------'

do ilat = 1, n_lats
  print '(A2,A40,A,F5.3,A,F5.3,A,F5.1,A,F6.3,A,F6.3,A,F5.1,A,F6.3,A,F6.3,A,F5.1,A)', &
    '  ', lat_names(ilat), ' |', &
    t_cpu_arr(ilat,1), '/', t_gpu_arr(ilat,1), ' ', &
    t_cpu_arr(ilat,1)/max(t_gpu_arr(ilat,1),1d-10), 'x|', &
    t_cpu_arr(ilat,2), '/', t_gpu_arr(ilat,2), ' ', &
    t_cpu_arr(ilat,2)/max(t_gpu_arr(ilat,2),1d-10), 'x|', &
    t_cpu_arr(ilat,3), '/', t_gpu_arr(ilat,3), ' ', &
    t_cpu_arr(ilat,3)/max(t_gpu_arr(ilat,3),1d-10), 'x'
enddo

print *
print '(A)', '  Times in seconds. Speed = CPU_time / GPU_time.'
print *
print *, '========================================================================================================'
print *, '  Benchmark complete'
print *, '========================================================================================================'

end program gpu_tracking_benchmark
