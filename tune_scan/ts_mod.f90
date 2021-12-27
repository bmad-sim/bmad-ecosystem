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
  real(rp) :: a_emit = 0, b_emit = 0
  real(rp) :: a0_amp = 0, b0_amp = 0, pz0_amp = 0
  integer :: n_turn_data = 0, n_turn_init = 0
  logical :: use_phase_trombone = .false.
  logical :: debug = .false.
end type

type ts_com_struct
  character(100) master_input_file
  type (lat_struct) ring
  type (coord_struct), allocatable :: closed_orb(:)
  real(rp) sig_a, sig_b, sig_pz
  integer n_a, n_b, n_z
  integer int_Qa, int_Qb
  integer :: mpi_rank = master_rank$
  real(rp) :: time_start = 0
  logical :: using_mpi = .false.
end type

type ts_data_struct
  real(rp) tune(3)
  real(rp) A_max(3)
  real(rp) A_rms(3)
  real(rp) n_turn_track
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
print *, 'Opening: ', trim(ts_com%master_input_file)

open (unit= 1, file = ts_com%master_input_file, status = 'old')
read(1, nml = params)
close (1)

ts%Q_z0 = abs(ts%Q_z0)
ts%Q_z1 = abs(ts%Q_z1)
ts%dQ_z = abs(ts%dQ_z)

if (ts%dat_out_file == '') call file_suffixer(ts_com%master_input_file, ts%dat_out_file, 'dat', .true.)

!---------------------------------------------
! Calculate number of steps to take from range and step size

if (ts%dQ_a > 0) then
   ts_com%n_a = nint(abs((ts%Q_a1 - ts%Q_a0) / ts%dQ_a))
else
   ts_com%n_a = 0
endif
if (ts%dQ_b > 0) then
   ts_com%n_b = nint(abs((ts%Q_b1 - ts%Q_b0) / ts%dQ_b))
else
   ts_com%n_b = 0
endif
if (ts%dQ_z > 0) then
   ts_com%n_z = nint(abs((ts%Q_z1 - ts%Q_z0) / ts%dQ_z))
else
   ts_com%n_z = 0
endif

if (ts_com%mpi_rank == master_rank$) then
  print '(3(a, f10.6) a, i0, a)', 'ts%dQ_a, ts%Q_a0, ts%Q_a1 = [', ts%dQ_a, ', ', ts%Q_a0, ', ', ts%Q_a1, '],   n_a = [0, ', ts_com%n_a, ']'
  print '(3(a, f10.6) a, i0, a)', 'ts%dQ_b, ts%Q_b0, ts%Q_b1 = [', ts%dQ_b, ', ', ts%Q_b0, ', ', ts%Q_b1, '],   n_b = [0, ', ts_com%n_b, ']'
  print '(3(a, f10.6) a, i0, a)', 'ts%dQ_z, ts%Q_z0, ts%Q_z1 = [', ts%dQ_z, ', ', ts%Q_z0, ', ', ts%Q_z1, '],   n_z = [0, ', ts_com%n_z, ']'
endif

!---------------------------------------------
! Initialize lattice

bmad_com%auto_bookkeeper = .false.
global_com%exit_on_error = .false.
call bmad_parser(ts%lat_file, ts_com%ring)
if (ts%use_phase_trombone) call insert_phase_trombone(ts_com%ring%branch(0))

allocate(ts_com%closed_orb(0:ts_com%ring%n_ele_max))
bmad_com%aperture_limit_on = .true.

print *, "Calculating closed orbit, transfer matrices, and Twiss parameters..."
call set_on_off(rfcavity$, ts_com%ring, on$)

call closed_orbit_calc(ts_com%ring, ts_com%closed_orb, 6)
call lat_make_mat6(ts_com%ring, -1, ts_com%closed_orb)
call twiss_at_start(ts_com%ring)
call twiss_propagate_all(ts_com%ring)
call radiation_integrals (ts_com%ring, ts_com%closed_orb, mode)
call calc_z_tune (ts_com%ring)

!---------------------------------------------

ts_com%int_Qa = int(ts_com%ring%ele(ts_com%ring%n_ele_track)%a%phi / twopi)
ts_com%int_Qb = int(ts_com%ring%ele(ts_com%ring%n_ele_track)%b%phi / twopi)
ele => ts_com%ring%ele(0)

if (ts%a_emit == 0) then
  ts_com%sig_a = sqrt(mode%a%emittance * ele%a%beta)
else
  ts_com%sig_a = sqrt(ts%a_emit * ele%a%beta)
endif

if (ts%b_emit == 0) then
  ts_com%sig_b = sqrt(mode%b%emittance * ele%b%beta)
else
  ts_com%sig_b = sqrt(ts%b_emit * ele%b%beta)
endif

ts_com%sig_pz = mode%sigE_E

end subroutine ts_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ts_track_particle(ts, ts_com, ja, jb, jz, ts_dat)

type (ts_params_struct) ts
type (ts_com_struct) ts_com
type (ts_data_struct) ts_dat
type (coord_struct) coords
type (lat_struct) ring
type (coord_struct), allocatable :: closed_orb(:), orbit(:)

real(rp) Jvec0(1:6), Jvec(1:6), init_vec(6), amp(3)
integer ja, jb, jz, nt, track_state
logical ok, error

!

ts_dat%tune(1) = ts_com%int_Qa + ts%Q_a0 + ja*ts%dQ_a
ts_dat%tune(2) = ts_com%int_Qb + ts%Q_b0 + jb*ts%dQ_b
ts_dat%tune(3) = -1.0*(ts%Q_z0 + jz*ts%dQ_z)

ts_dat%a_max = 1000.0_rp
ts_dat%a_rms = 1000.0_rp
ts_dat%n_turn_track = 0

ring = ts_com%ring              ! Use copy in case tune setting fails, which may garble the lattice
closed_orb = ts_com%closed_orb  ! Use copy in case tune setting fails, which may garble the closed orbit

call set_tune3 (ring, ts_dat%tune, ts%use_phase_trombone, ts%quad_mask, ok)  ! Takes tunes that have not not been multiplied by 2pi.

if (.not. ok) return    ! Tunes could not be set, probably on a resonance.

call closed_orbit_calc(ring, closed_orb,6)
call lat_make_mat6(ring, -1, closed_orb)
call twiss_at_start(ring)
call twiss_propagate_all(ring)
call calc_z_tune (ring)

!---------------------------------------------
! Calculate initial vector from particle actions.

Jvec0 = [ts%a0_amp * ts_com%sig_a, 0.0_rp, &
         ts%b0_amp * ts_com%sig_b, 0.0_rp, &
         0.0_rp, ts%pz0_amp * ts_com%sig_pz]
call action_to_xyz(ts_com%ring, 0, Jvec0, init_vec, error)

if (error) return  ! decomposition of the 1-turn matrix failed.

call reallocate_coord(orbit, ring)
orbit(0)%vec(:) = closed_orb(0)%vec(:) + init_vec(:)
ts_dat%a_max = 0
ts_dat%a_rms = 0

do nt = 1, ts%n_turn_init + ts%n_turn_data
  call track_all(ring, orbit, 0, track_state)
  orbit(0) = orbit(ring%n_ele_track)

  if (track_state /= moving_forward$) exit

  if (nt > ts%n_turn_init) then ! save this data
    coords%vec(:) = orbit(0)%vec(:) - closed_orb(0)%vec(:)
    call xyz_to_action(ring, 0, coords%vec, Jvec, error)
    if (error) then
       print *, "BAD: xyz_to_action returned error."
    endif
    amp = [(Jvec(1)**2 + Jvec(2)**2)/2.0d0, (Jvec(3)**2 + Jvec(4)**2)/2.0d0, (Jvec(5)**2 + Jvec(6)**2)/2.0d0]
    ts_dat%a_max = max(ts_dat%a_max, amp)
    ts_dat%a_rms = ts_dat%a_rms + amp**2
  endif
enddo

ts_dat%n_turn_track = nt
ts_dat%a_rms = sqrt(ts_dat%a_rms/ts%n_turn_data)

end subroutine ts_track_particle

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine ts_write_results (ts, ts_com, ts_dat)

type (ts_params_struct) ts
type (ts_com_struct) ts_com
type (ts_data_struct), target :: ts_dat(:,:,:)
type (ts_data_struct), pointer :: t

integer ja, jb, jz

!

open(unit = 23, file = ts%dat_out_file)
write(23, '(3a10,3a13,a15, 3a13)') '#      Q_a', 'Q_b', 'Q_z', 'amp_a_max', 'amp_b_max', 'amp_z_max', 'turns', &
                                                                              'amp_a_rms', 'amp_b_rms', 'amp_z_rms'

do jz = 0, ts_com%n_z
do jb = 0, ts_com%n_b
do ja = 0, ts_com%n_a
  t => ts_dat(ja, jb, jz)
  write(23, '(3f10.5, 3es13.4, i15, 3es13.4)') t%tune(1)-ts_com%int_Qa, t%tune(2)-ts_com%int_Qb, t%tune(3), &
                                                                                t%a_max, t%n_turn_track, t%a_rms
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
