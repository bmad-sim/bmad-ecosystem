module ts_mod

use bsim_interface
use mode3_mod
 
implicit none

integer, parameter :: master_rank$   = 0
integer, parameter :: job_tag$       = 1000
integer, parameter :: have_data_tag$ = 1001
integer, parameter :: results_tag$   = 1002

type ts_params_struct
  character(100) :: lat_file = '', dat_out_file = '', quad_mask = ''
  real(rp) :: Q_a0 = 0, Q_a1 = 0, dQ_a = 0
  real(rp) :: Q_b0 = 0, Q_b1 = 0, dQ_b = 0
  real(rp) :: Q_z0 = 0, Q_z1 = 0, dQ_z = 0
  real(rp) :: pz0 = 0, pz1 = 0, dpz = 0
  real(rp) :: a_emit = 0, b_emit = 0, sig_pz
  real(rp) :: a0_amp = 0, b0_amp = 0, pz0_amp = 0
  real(rp) :: timer_print_dtime = 120
  integer :: n_turn = 0
  integer :: ix_branch = 0
  logical :: use_phase_trombone = .false.
  logical :: debug = .false.
  logical :: rf_on = .false.
end type

type ts_com_struct
  character(100) master_input_file
  type (lat_struct) ring
  type (coord_struct), allocatable :: closed_orb(:)
  real(rp) sig_a, sig_b, sig_pz, a_emit, b_emit
  integer n_a, n_b, n_z
  integer int_Qa, int_Qb
  integer :: mpi_rank = master_rank$
  real(rp) :: time_start = 0
  logical :: using_mpi = .false.
end type

type ts_data_struct
  integer :: ix_q(3) = 0
  real(rp) :: tune(3) = 0   ! Either (Qa, Qb, Qz) with RF on or (Qa, Qb, pz) with RF off.
  real(rp) :: amp_max(3) = -1
  real(rp) :: amp_ave(3) = -1
  real(rp) :: amp_rms(3) = -1
  integer :: n_data_track = 0
end type

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ts_init_params (ts, ts_com)

type (ts_params_struct), target :: ts
type (ts_com_struct), target :: ts_com
type (normal_modes_struct) mode
type (ele_struct), pointer :: ele

integer n_arg
logical err

namelist / params / bmad_com, ts

!---------------------------------------------
! Read in the parameters

n_arg = command_argument_count()
if (n_arg > 1) then
  print *, 'Usage: tune_scan <input_file>'
  print *, 'Default: <input_file> = tune_scan.init'
  stop
endif

ts_com%master_input_file = 'tune_scan.init'
if (n_arg == 1) call get_command_argument(1, ts_com%master_input_file)

open (unit= 1, file = ts_com%master_input_file, status = 'old', action = 'read')
read(1, nml = params)
close (1)

ts%Q_z0 = abs(ts%Q_z0)
ts%Q_z1 = abs(ts%Q_z1)
ts%dQ_z = abs(ts%dQ_z)

if (ts%dat_out_file == '') call file_suffixer(ts_com%master_input_file, ts%dat_out_file, 'dat', .true.)

!---------------------------------------------
! Calculate number of steps to take from range and step size

ts_com%n_a = 0
ts_com%n_b = 0
ts_com%n_z = 0

if (ts%dQ_a > 0) ts_com%n_a = nint(abs((ts%Q_a1 - ts%Q_a0) / ts%dQ_a))
if (ts%dQ_b > 0) ts_com%n_b = nint(abs((ts%Q_b1 - ts%Q_b0) / ts%dQ_b))

if (ts%rf_on) then
  if (ts%dQ_z > 0) ts_com%n_z = nint(abs((ts%Q_z1 - ts%Q_z0) / ts%dQ_z))
else
  if (ts%dpz > 0)  ts_com%n_z = nint(abs((ts%pz1 - ts%pz0) / ts%dpz))
endif

if (ts_com%mpi_rank == master_rank$) then
  print '(3(a, f10.6) a, i0, a)', 'ts%dQ_a, ts%Q_a0, ts%Q_a1 = [', ts%dQ_a, ', ', ts%Q_a0, ', ', ts%Q_a1, '],   n_a = [0, ', ts_com%n_a, ']'
  print '(3(a, f10.6) a, i0, a)', 'ts%dQ_b, ts%Q_b0, ts%Q_b1 = [', ts%dQ_b, ', ', ts%Q_b0, ', ', ts%Q_b1, '],   n_b = [0, ', ts_com%n_b, ']'
  if (ts%rf_on) then
    print '(3(a, f10.6) a, i0, a)', 'ts%dQ_z, ts%Q_z0, ts%Q_z1 = [', ts%dQ_z, ', ', ts%Q_z0, ', ', ts%Q_z1, '],   n_z = [0, ', ts_com%n_z, ']'
  else
    print '(3(a, f10.6) a, i0, a)', 'ts%dpz, ts%pz0, ts%pz1    = [', ts%dpz,  ', ', ts%pz0,  ', ', ts%pz1,  '],   n_z = [0, ', ts_com%n_z, ']'
  endif
endif

!---------------------------------------------
! Initialize lattice

bmad_com%auto_bookkeeper = .false.
global_com%exit_on_error = .false.

call bmad_parser(ts%lat_file, ts_com%ring, err_flag = err)
if (err) stop

if (.not. ts%rf_on) call set_on_off(rfcavity$, ts_com%ring, off$)
if (ts%use_phase_trombone) call insert_phase_trombone(ts_com%ring%branch(0))

allocate(ts_com%closed_orb(0:ts_com%ring%n_ele_max))
bmad_com%aperture_limit_on = .true.

if (ts%rf_on) then
  call closed_orbit_calc(ts_com%ring, ts_com%closed_orb, 6)
else
  call closed_orbit_calc(ts_com%ring, ts_com%closed_orb, 4)
endif

call lat_make_mat6(ts_com%ring, -1, ts_com%closed_orb)
call twiss_at_start(ts_com%ring)
call twiss_propagate_all(ts_com%ring)
call radiation_integrals (ts_com%ring, ts_com%closed_orb, mode)
call calc_z_tune (ts_com%ring%branch(0))

!---------------------------------------------

ts_com%int_Qa = int(ts_com%ring%ele(ts_com%ring%n_ele_track)%a%phi / twopi)
ts_com%int_Qb = int(ts_com%ring%ele(ts_com%ring%n_ele_track)%b%phi / twopi)
ele => ts_com%ring%ele(0)

ts_com%a_emit = ts%a_emit
if (ts_com%a_emit == 0) ts_com%a_emit = mode%a%emittance
ts_com%sig_a = sqrt(ts_com%a_emit * ele%a%beta)

ts_com%b_emit = ts%b_emit
if (ts_com%b_emit == 0) ts_com%b_emit = mode%b%emittance
ts_com%sig_b = sqrt(ts_com%b_emit * ele%b%beta)

ts_com%sig_pz = ts%sig_pz
if (ts_com%sig_pz == 0) ts_com%sig_pz = mode%sigE_E

end subroutine ts_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

recursive subroutine ts_track_particle(ts, ts_com, ja, jb, jz, ts_dat)

type (ts_params_struct) ts
type (ts_com_struct) ts_com
type (ts_data_struct) ts_dat
type (coord_struct) coords
type (lat_struct), target :: ring
type (coord_struct), allocatable :: closed_orb(:), orbit(:)
type (ele_struct), pointer :: ele

real(rp) Jvec0(1:6), Jvec(1:6), init_vec(6), amp(3), r(3), v_mat(4,4)
integer ja, jb, jz, nt, track_state, status
logical ok, error

!

ts_dat = ts_data_struct()
ts_dat%ix_q = [ja, jb, jz]

ts_dat%tune(1) = ts_com%int_Qa + ts%Q_a0 + ja*ts%dQ_a
ts_dat%tune(2) = ts_com%int_Qb + ts%Q_b0 + jb*ts%dQ_b
if (ts%rf_on) then
  ts_dat%tune(3) = -(ts%Q_z0 + jz*ts%dQ_z)
else
  ts_dat%tune(3) = ts%pz0 + jz*ts%dpz
endif

ring = ts_com%ring              ! Use copy in case tune setting fails, which may garble the lattice
closed_orb = ts_com%closed_orb  ! Use copy in case tune setting fails, which may garble the closed orbit
ele => ring%ele(0)

ok = set_tune_3d (ring%branch(0), ts_dat%tune, ts%quad_mask, ts%use_phase_trombone, ts%rf_on)  ! Tunes in radians.
if (.not. ok) return    ! Tunes could not be set, probably on a resonance.

if (ts%rf_on) then
  call closed_orbit_calc(ring, closed_orb, 6)
else
  closed_orb(0)%vec(6) = ts_dat%tune(3)
  call closed_orbit_calc(ring, closed_orb, 4)
endif

call lat_make_mat6(ring, -1, closed_orb)
call twiss_at_start(ring, status, ts%ix_branch)
if (status /= ok$) then
  print '(a, 3f10.4)', 'Twiss calc fail at tunes:', ts_dat%tune/twopi
  return
endif

call twiss_propagate_all(ring)
call calc_z_tune (ring%branch(0))

!---------------------------------------------
! Calculate initial vector from particle actions.
! Note: Use ts_com%ring for converting between action/xyz which makes the conversion consistent between tune points.

if (ts%rf_on) then
  Jvec0 = [ts%a0_amp * sqrt(ts_com%a_emit), 0.0_rp, &
           ts%b0_amp * sqrt(ts_com%b_emit), 0.0_rp, &
           0.0_rp, ts%pz0_amp * ts_com%sig_pz]
  call action_to_xyz(ring, 0, Jvec0, init_vec, error)
  if (error) return  ! decomposition of the 1-turn matrix failed.
else
  call make_v_mats(ele, v_mat)
  init_vec = 0
  init_vec(1) = ts%a0_amp * sqrt(ts_com%a_emit / ele%a%gamma)
  init_vec(3) = ts%b0_amp * sqrt(ts_com%b_emit / ele%b%gamma)
  init_vec(1:4) = matmul(v_mat, init_vec(1:4))
endif

call reallocate_coord(orbit, ring)
orbit(0)%vec(:) = closed_orb(0)%vec(:) + init_vec(:)

do nt = 1, ts%n_turn
  call track_all(ring, orbit, 0, track_state)
  orbit(0) = orbit(ring%n_ele_track)
  if (track_state /= moving_forward$) exit

  coords%vec(:) = orbit(0)%vec(:) - closed_orb(0)%vec(:)
  if (ts%rf_on) then
    call xyz_to_action(ring, 0, coords%vec, Jvec, error)
    amp = [sqrt(Jvec(1)**2 + Jvec(2)**2)/sqrt(ts_com%a_emit), sqrt(Jvec(3)**2 + Jvec(4)**2)/sqrt(ts_com%b_emit), &
                                                                         sqrt(Jvec(5)**2 + Jvec(6)**2)/ts_com%sig_pz]
    if (error) then
       print *, "BAD: xyz_to_action returned error."
    endif
  else
    coords%vec(6) = 0  ! The dispersion has already been subtracted off
    call orbit_amplitude_calc(ele, coords, amp(1), amp(2))
    amp = [2*amp(1) / ts_com%a_emit, 2*amp(2) / ts_com%b_emit, 0.0_rp]
  endif

  ts_dat%amp_max = max(ts_dat%amp_max, amp)
  ts_dat%amp_rms = ts_dat%amp_rms + amp**2
  ts_dat%amp_ave = ts_dat%amp_ave + amp
  ts_dat%n_data_track = ts_dat%n_data_track + 1
enddo

ts_dat%amp_ave = ts_dat%amp_ave/ts_dat%n_data_track
r = ts_dat%amp_rms/ts_dat%n_data_track - ts_dat%amp_ave**2
ts_dat%amp_rms = sqrt(max(0.0_rp, r))

end subroutine ts_track_particle

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ts_write_results (ts, ts_com, ts_dat)

type (ts_params_struct) ts
type (ts_com_struct) ts_com
type (ts_data_struct), target :: ts_dat(0:,0:,0:)
type (ts_data_struct), pointer :: t

integer ja, jb, jz

!

open(unit = 23, file = ts%dat_out_file)

write (23, '(a, a)')         '# lat_file                   = ', quote(ts%lat_file)
write (23, '(a, a)')         '# quad_mask                  = ', quote(ts%quad_mask)
write (23, '(a, es12.4)')    '# Q_a0                       = ', ts%Q_a0
write (23, '(a, es12.4)')    '# Q_a1                       = ', ts%Q_a1
write (23, '(a, es12.4)')    '# dQ_a                       = ', ts%dQ_a
write (23, '(a, es12.4)')    '# a0_amp                     = ', ts%a0_amp
write (23, '(a, es12.4)')    '# Q_b0                       = ', ts%Q_b0
write (23, '(a, es12.4)')    '# Q_b1                       = ', ts%Q_b1
write (23, '(a, es12.4)')    '# dQ_b                       = ', ts%dQ_b
write (23, '(a, es12.4)')    '# b0_amp                     = ', ts%b0_amp
if (ts%rf_on) then
  write (23, '(a, es12.4)')  '# Q_z0                       = ', ts%Q_z0
  write (23, '(a, es12.4)')  '# Q_z1                       = ', ts%Q_z1
  write (23, '(a, es12.4)')  '# dQ_z                       = ', ts%dQ_z
  write (23, '(a, es12.4)')  '# pz0_amp                    = ', ts%pz0_amp
else
  write (23, '(a, es12.4)')  '# pz0                        = ', ts%pz0
  write (23, '(a, es12.4)')  '# pz1                        = ', ts%pz1
  write (23, '(a, es12.4)')  '# dpz                        = ', ts%dpz
endif
write (23, '(a, es12.4)')    '# a_emit                     = ', ts%a_emit
write (23, '(a, es12.4)')    '# b_emit                     = ', ts%b_emit
write (23, '(a, es12.4)')    '# sig_pz                     = ', ts%sig_pz
write (23, '(a, i8)')        '# na_max                     = ', ts_com%n_a
write (23, '(a, i8)')        '# nb_max                     = ', ts_com%n_b
write (23, '(a, i8)')        '# nz_max                     = ', ts_com%n_z
write (23, '(a, i8)')        '# n_turn                     = ', ts%n_turn
write (23, '(a, es12.4, a)') '# sigma_a                    = ', ts_com%sig_a,  '  # Used in calculation'
write (23, '(a, es12.4, a)') '# sigma_b                    = ', ts_com%sig_b,  '  # Used in calculation'
write (23, '(a, es12.4, a)') '# sigma_pz                   = ', ts_com%sig_pz, '  # Used in calculation'
write (23, '(a, es12.4, a)') '# emittance_a                = ', ts_com%a_emit, '  # Used in calculation'
write (23, '(a, es12.4, a)') '# emittance_b                = ', ts_com%b_emit, '  # Used in calculation'
write (23, '(a, l4)')        '# radiation_damping_on       = ', bmad_com%radiation_damping_on
write (23, '(a, l4)')        '# radiation_fluctuations_on  = ', bmad_com%radiation_fluctuations_on
write (23, '(a, l4)')        '# rf_on                      = ', ts%rf_on
write (23, '(a, l4)')        '# use_phase_trombone         = ', ts%use_phase_trombone

if (ts%rf_on) then
  write (23, '(a, a4, 2a6, 3a10, a12, 3a13, 3a13)') '#-', 'ja', 'jb', 'jz', 'Q_a', 'Q_b', 'Q_z', 'data_turns', &
                              'Amp_a_rms', 'Amp_b_rms', 'Amp_z_rms', 'Amp_a_max', 'Amp_b_max', 'Amp_z_max'
else
  write (23, '(a, a4, 2a6, 3a10, a12, 3a13, 3a13)') '#-', 'ja', 'jb', 'jz', 'Q_a', 'Q_b', 'pz', 'data_turns', &
                              'Amp_a_rms', 'Amp_b_rms', 'Amp_z_rms', 'Amp_a_max', 'Amp_b_max', 'Amp_z_max'
endif

do jz = 0, ts_com%n_z
do jb = 0, ts_com%n_b
do ja = 0, ts_com%n_a
  t => ts_dat(ja, jb, jz)
  if (ts%rf_on) then
    write(23, '(3i6, 3f10.5, i12, 3es13.4, 3es13.4)') ja, jb, jz, &
                  t%tune(1)-ts_com%int_Qa, t%tune(2)-ts_com%int_Qb, -t%tune(3), t%n_data_track, t%amp_rms, t%amp_max
  else
    write(23, '(3i6, 3f10.5, i12, 3es13.4, 3es13.4)') ja, jb, jz, &
                  t%tune(1)-ts_com%int_Qa, t%tune(2)-ts_com%int_Qb, t%tune(3), t%n_data_track, t%amp_rms, t%amp_max
  endif
enddo
enddo
enddo

close(23)

end subroutine ts_write_results

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ts_print_mpi_info (ts, ts_com, line, do_print)

type (ts_params_struct) ts
type (ts_com_struct) ts_com

real(rp) time_now
character(*) line
character(20) time_str
logical, optional :: do_print

!

if (.not. logic_option(ts%debug, do_print)) return

call run_timer ('ABS', time_now)
call date_and_time_stamp (time_str)
print '(a, f8.2, 2a, 2x, i0, 2a)', 'dTime:', (time_now-ts_com%time_start)/60, &
                                        ' Now: ', time_str, ts_com%mpi_rank, ': ', trim(line)

end subroutine ts_print_mpi_info

end module
