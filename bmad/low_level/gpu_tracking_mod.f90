module gpu_tracking_mod

use bmad_struct
use bmad_routine_interface, dummy_gtm => track_a_drift

implicit none
private

public :: gpu_tracking_is_active
public :: gpu_tracking_reset
public :: track_bunch_thru_drift_gpu
public :: track_bunch_thru_quad_gpu

! Cached runtime dispatch decision (checked once per process)
logical, save :: gpu_trk_checked  = .false.
logical, save :: gpu_trk_active   = .false.

#ifdef USE_GPU_TRACKING
! ----- C interfaces (gpu_tracking_kernels.cu) ---------------------------------
interface
  subroutine gpu_track_drift(vx, vpx, vy, vpy, vz, vpz, &
                             state, beta, p0c, s_pos, t_time, &
                             mc2, length, n) bind(C, name='gpu_track_drift_')
    use, intrinsic :: iso_c_binding
    real(C_DOUBLE), intent(inout) :: vx(*), vpx(*), vy(*), vpy(*), vz(*), vpz(*)
    integer(C_INT), intent(inout) :: state(*)
    real(C_DOUBLE), intent(inout) :: beta(*), p0c(*), s_pos(*), t_time(*)
    real(C_DOUBLE), value, intent(in) :: mc2, length
    integer(C_INT), value, intent(in) :: n
  end subroutine

  subroutine gpu_track_quad(vx, vpx, vy, vpy, vz, vpz, &
                            state, beta, p0c, t_time, &
                            mc2, b1, ele_length, delta_ref_time, &
                            e_tot_ele, charge_dir, n) bind(C, name='gpu_track_quad_')
    use, intrinsic :: iso_c_binding
    real(C_DOUBLE), intent(inout) :: vx(*), vpx(*), vy(*), vpy(*), vz(*), vpz(*)
    integer(C_INT), intent(inout) :: state(*)
    real(C_DOUBLE), intent(inout) :: beta(*), p0c(*), t_time(*)
    real(C_DOUBLE), value, intent(in) :: mc2, b1, ele_length
    real(C_DOUBLE), value, intent(in) :: delta_ref_time, e_tot_ele, charge_dir
    integer(C_INT), value, intent(in) :: n
  end subroutine

  function gpu_tracking_available() result(avail) bind(C, name='gpu_tracking_available_')
    use, intrinsic :: iso_c_binding
    integer(C_INT) :: avail
  end function
end interface
#endif

contains

!------------------------------------------------------------------------
! gpu_tracking_is_active — check env var + GPU presence (once per process)
!------------------------------------------------------------------------
function gpu_tracking_is_active() result (is_active)
logical :: is_active
character(len=32) :: env_val
integer :: env_len, env_stat

if (.not. gpu_trk_checked) then
  gpu_trk_checked = .true.
#ifdef USE_GPU_TRACKING
  call get_environment_variable('ACC_ENABLE_GPU_TRACKING', env_val, env_len, env_stat)
  if (env_stat == 0 .and. trim(env_val) == 'Y') then
    if (gpu_tracking_available() == 1) then
      gpu_trk_active = .true.
      print *, 'gpu_tracking: ACC_ENABLE_GPU_TRACKING=Y and CUDA GPU detected — using GPU tracking'
    else
      print *, 'gpu_tracking: ACC_ENABLE_GPU_TRACKING=Y but no CUDA GPU found — using CPU tracking'
    endif
  endif
#endif
endif

is_active = gpu_trk_active
end function

!------------------------------------------------------------------------
! gpu_tracking_reset — force re-check of env var on next call
!------------------------------------------------------------------------
subroutine gpu_tracking_reset()
gpu_trk_checked = .false.
gpu_trk_active  = .false.
end subroutine

!------------------------------------------------------------------------
! track_bunch_thru_drift_gpu
!
! GPU batch tracking of all particles in a bunch through a drift.
! Extracts particle data into SoA arrays, calls CUDA kernel, writes back.
!------------------------------------------------------------------------
subroutine track_bunch_thru_drift_gpu (bunch, ele)

use, intrinsic :: iso_c_binding

type (bunch_struct), intent(inout) :: bunch
type (ele_struct),   intent(in)    :: ele

#ifdef USE_GPU_TRACKING
integer(C_INT) :: n
integer :: j
real(rp) :: length, mc2

real(C_DOUBLE), allocatable :: vx(:), vpx(:), vy(:), vpy(:), vz(:), vpz(:)
real(C_DOUBLE), allocatable :: beta_a(:), p0c_a(:), s_a(:), t_a(:)
integer(C_INT), allocatable :: state_a(:)

n = size(bunch%particle)
length = ele%value(l$)
if (length == 0) return

mc2 = mass_of(bunch%particle(1)%species)

! Allocate SoA arrays
allocate(vx(n), vpx(n), vy(n), vpy(n), vz(n), vpz(n))
allocate(state_a(n), beta_a(n), p0c_a(n), s_a(n), t_a(n))

! AoS -> SoA extraction
do j = 1, n
  vx(j)      = bunch%particle(j)%vec(1)
  vpx(j)     = bunch%particle(j)%vec(2)
  vy(j)      = bunch%particle(j)%vec(3)
  vpy(j)     = bunch%particle(j)%vec(4)
  vz(j)      = bunch%particle(j)%vec(5)
  vpz(j)     = bunch%particle(j)%vec(6)
  state_a(j) = bunch%particle(j)%state
  beta_a(j)  = bunch%particle(j)%beta
  p0c_a(j)   = bunch%particle(j)%p0c
  s_a(j)     = bunch%particle(j)%s
  t_a(j)     = bunch%particle(j)%t
enddo

! Call CUDA kernel
call gpu_track_drift(vx, vpx, vy, vpy, vz, vpz, &
                     state_a, beta_a, p0c_a, s_a, t_a, &
                     mc2, length, n)

! SoA -> AoS write-back
do j = 1, n
  bunch%particle(j)%vec(1) = vx(j)
  bunch%particle(j)%vec(2) = vpx(j)
  bunch%particle(j)%vec(3) = vy(j)
  bunch%particle(j)%vec(4) = vpy(j)
  bunch%particle(j)%vec(5) = vz(j)
  bunch%particle(j)%vec(6) = vpz(j)
  bunch%particle(j)%state  = state_a(j)
  bunch%particle(j)%s      = s_a(j)
  bunch%particle(j)%t      = t_a(j)
  ! Set location to downstream end
  bunch%particle(j)%location = downstream_end$
  bunch%particle(j)%ix_ele   = ele%ix_ele
  bunch%particle(j)%ix_branch = ele%ix_branch
enddo

deallocate(vx, vpx, vy, vpy, vz, vpz)
deallocate(state_a, beta_a, p0c_a, s_a, t_a)
#endif

end subroutine track_bunch_thru_drift_gpu

!------------------------------------------------------------------------
! track_bunch_thru_quad_gpu
!
! GPU batch tracking through a quadrupole.  Only handles the simple case:
! no fringe fields, no multipoles beyond b1, no electric multipoles.
! Sets did_track = .false. and returns immediately if the element has
! features the GPU kernel cannot handle, so the caller can fall back to
! CPU tracking.
!------------------------------------------------------------------------
subroutine track_bunch_thru_quad_gpu (bunch, ele, param, did_track)

use, intrinsic :: iso_c_binding

type (bunch_struct),     intent(inout) :: bunch
type (ele_struct),       intent(in)    :: ele
type (lat_param_struct), intent(in)    :: param
logical,                 intent(out)   :: did_track

#ifdef USE_GPU_TRACKING
integer(C_INT) :: n
integer :: j, ix_mag_max, ix_elec_max
real(rp) :: ele_length, mc2, b1, delta_ref_time, e_tot_ele
real(rp) :: charge_dir, rel_tracking_charge
real(rp) :: an(0:n_pole_maxx), bn(0:n_pole_maxx)
real(rp) :: an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)
type (fringe_field_info_struct) :: fringe_info
logical :: has_misalign

real(C_DOUBLE), allocatable :: vx(:), vpx(:), vy(:), vpy(:), vz(:), vpz(:)
real(C_DOUBLE), allocatable :: beta_a(:), p0c_a(:), t_a(:)
integer(C_INT), allocatable :: state_a(:)

did_track = .false.

n = size(bunch%particle)
ele_length = ele%value(l$)
if (ele_length == 0) then
  did_track = .true.
  return
endif

mc2 = mass_of(bunch%particle(1)%species)
delta_ref_time = ele%value(delta_ref_time$)
e_tot_ele = ele%value(e_tot$)

! --- Safety checks: bail out to CPU if element has unsupported features ---

has_misalign = ele%bookkeeping_state%has_misalign

! Bail out if fringe fields are active
call init_fringe_info(fringe_info, ele)
if (fringe_info%has_fringe) return

! Get the quad gradient b1 and check for extra multipoles
call multipole_ele_to_ab(ele, .false., ix_mag_max, an, bn, magnetic$, include_kicks$, b1)

! Bail out if higher magnetic multipoles are present (GPU kernel only handles b1)
if (ix_mag_max > -1) return

! Bail out if electric multipoles are present
call multipole_ele_to_ab(ele, .false., ix_elec_max, an_elec, bn_elec, electric$)
if (ix_elec_max > -1) return

! Compute charge_dir: rel_charge * orientation (forward tracking: direction=1, time_dir=1)
rel_tracking_charge = rel_tracking_charge_to_mass(bunch%particle(1), param%particle)
charge_dir = rel_tracking_charge * ele%orientation

! Allocate SoA arrays
allocate(vx(n), vpx(n), vy(n), vpy(n), vz(n), vpz(n))
allocate(state_a(n), beta_a(n), p0c_a(n), t_a(n))

! Apply entrance misalignment transform (lab → body frame) on CPU
if (has_misalign) then
  do j = 1, n
    if (bunch%particle(j)%state == alive$) then
      call offset_particle(ele, set$, bunch%particle(j), set_hvkicks = .false.)
    endif
  enddo
endif

! AoS -> SoA extraction (now in body frame if misaligned)
do j = 1, n
  vx(j)      = bunch%particle(j)%vec(1)
  vpx(j)     = bunch%particle(j)%vec(2)
  vy(j)      = bunch%particle(j)%vec(3)
  vpy(j)     = bunch%particle(j)%vec(4)
  vz(j)      = bunch%particle(j)%vec(5)
  vpz(j)     = bunch%particle(j)%vec(6)
  state_a(j) = bunch%particle(j)%state
  beta_a(j)  = bunch%particle(j)%beta
  p0c_a(j)   = bunch%particle(j)%p0c
  t_a(j)     = bunch%particle(j)%t
enddo

! All checks passed — proceed with GPU tracking
did_track = .true.

! Call CUDA kernel
call gpu_track_quad(vx, vpx, vy, vpy, vz, vpz, &
                    state_a, beta_a, p0c_a, t_a, &
                    mc2, b1, ele_length, delta_ref_time, &
                    e_tot_ele, charge_dir, n)

! SoA -> AoS write-back (still in body frame if misaligned)
do j = 1, n
  bunch%particle(j)%vec(1) = vx(j)
  bunch%particle(j)%vec(2) = vpx(j)
  bunch%particle(j)%vec(3) = vy(j)
  bunch%particle(j)%vec(4) = vpy(j)
  bunch%particle(j)%vec(5) = vz(j)
  bunch%particle(j)%vec(6) = vpz(j)
  bunch%particle(j)%state  = state_a(j)
  bunch%particle(j)%t      = t_a(j)
  bunch%particle(j)%location = downstream_end$
  bunch%particle(j)%ix_ele   = ele%ix_ele
  bunch%particle(j)%ix_branch = ele%ix_branch
enddo

deallocate(vx, vpx, vy, vpy, vz, vpz)
deallocate(state_a, beta_a, p0c_a, t_a)

! Apply exit misalignment transform (body → lab frame) on CPU
if (has_misalign) then
  do j = 1, n
    if (bunch%particle(j)%state == alive$) then
      call offset_particle(ele, unset$, bunch%particle(j), set_hvkicks = .false.)
    endif
  enddo
endif
#endif

end subroutine track_bunch_thru_quad_gpu

end module gpu_tracking_mod
