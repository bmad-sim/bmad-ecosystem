module gpu_tracking_mod

use bmad_struct
use bmad_routine_interface, dummy_gtm => track_a_drift

implicit none
private

public :: gpu_tracking_init
public :: ele_gpu_eligible
public :: track_bunch_thru_drift_gpu
public :: track_bunch_thru_quad_gpu

! Whether gpu_tracking_init has been called
logical, save :: gpu_trk_initialized = .false.
! Whether CUDA hardware is present (checked once, does not change)
logical, save :: gpu_hw_available = .false.

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
                            e_tot_ele, charge_dir, n_particles, &
                            a2_arr, b2_arr, cm_arr, &
                            ix_mag_max, n_step, &
                            ea2_arr, eb2_arr, ix_elec_max) bind(C, name='gpu_track_quad_')
    use, intrinsic :: iso_c_binding
    real(C_DOUBLE), intent(inout) :: vx(*), vpx(*), vy(*), vpy(*), vz(*), vpz(*)
    integer(C_INT), intent(inout) :: state(*)
    real(C_DOUBLE), intent(inout) :: beta(*), p0c(*), t_time(*)
    real(C_DOUBLE), value, intent(in) :: mc2, b1, ele_length
    real(C_DOUBLE), value, intent(in) :: delta_ref_time, e_tot_ele, charge_dir
    integer(C_INT), value, intent(in) :: n_particles
    real(C_DOUBLE), intent(in) :: a2_arr(*), b2_arr(*), cm_arr(*)
    integer(C_INT), value, intent(in) :: ix_mag_max, n_step
    real(C_DOUBLE), intent(in) :: ea2_arr(*), eb2_arr(*)
    integer(C_INT), value, intent(in) :: ix_elec_max
  end subroutine

  function gpu_tracking_available() result(avail) bind(C, name='gpu_tracking_available_')
    use, intrinsic :: iso_c_binding
    integer(C_INT) :: avail
  end function
end interface
#endif

contains

!------------------------------------------------------------------------
! gpu_tracking_init — initialize GPU tracking from env var (call once)
!
! Reads ACC_ENABLE_GPU_TRACKING env var and checks for CUDA hardware.
! Sets bmad_com%gpu_tracking_on = .true. if env var is 'Y' and GPU is
! present. After this, callers can toggle bmad_com%gpu_tracking_on
! directly. Does nothing if called more than once.
!------------------------------------------------------------------------
subroutine gpu_tracking_init()
character(len=32) :: env_val
integer :: env_len, env_stat

if (gpu_trk_initialized) return
gpu_trk_initialized = .true.

#ifdef USE_GPU_TRACKING
if (gpu_tracking_available() == 1) then
  gpu_hw_available = .true.
endif

call get_environment_variable('ACC_ENABLE_GPU_TRACKING', env_val, env_len, env_stat)
if (env_stat == 0 .and. trim(env_val) == 'Y') then
  if (gpu_hw_available) then
    bmad_com%gpu_tracking_on = .true.
    print *, 'gpu_tracking: GPU tracking enabled (ACC_ENABLE_GPU_TRACKING=Y, CUDA GPU detected)'
  else
    print *, 'gpu_tracking: ACC_ENABLE_GPU_TRACKING=Y but no CUDA GPU found — using CPU tracking'
  endif
endif
#endif

end subroutine

!------------------------------------------------------------------------
! ele_gpu_eligible — check if an element can be GPU-tracked
!
! Returns .true. if the element's intrinsic properties allow GPU tracking.
! This checks element type, tracking method, and on/off state.
! Runtime conditions (bmad_com flags, particle direction, wakefields)
! are NOT checked here — those are evaluated at dispatch time.
!
! Currently supported element types: drift, quadrupole.
!------------------------------------------------------------------------
function ele_gpu_eligible(ele) result (eligible)
type (ele_struct), intent(in) :: ele
logical :: eligible

eligible = .false.

! Must use bmad_standard tracking
if (ele%tracking_method /= bmad_standard$) return

! Must be turned on
if (.not. ele%is_on) return

! Check supported element types
select case (ele%key)
case (drift$, quadrupole$)
  eligible = .true.
end select

end function ele_gpu_eligible

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

use multipole_mod, only: ab_multipole_kicks
use, intrinsic :: iso_c_binding

type (bunch_struct),     intent(inout) :: bunch
type (ele_struct),       intent(in)    :: ele
type (lat_param_struct), intent(in)    :: param
logical,                 intent(out)   :: did_track

#ifdef USE_GPU_TRACKING
integer, parameter :: n_multi = n_pole_maxx + 1  ! = 22
integer(C_INT) :: n
integer :: j, nn, mm, ix_mag_max, ix_elec_max, n_step
real(rp) :: ele_length, mc2, b1, delta_ref_time, e_tot_ele
real(rp) :: charge_dir, rel_tracking_charge, r_step, length
real(rp) :: an(0:n_pole_maxx), bn(0:n_pole_maxx)
real(rp) :: an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)
real(rp) :: f_charge, r_ratio
type (fringe_field_info_struct) :: fringe_info
logical :: has_misalign, has_mag_multipoles

! Precomputed scaled multipole arrays and c_multi coefficients for CUDA
real(C_DOUBLE) :: a2_arr(0:n_pole_maxx), b2_arr(0:n_pole_maxx)
real(C_DOUBLE) :: ea2_arr(0:n_pole_maxx), eb2_arr(0:n_pole_maxx)
real(C_DOUBLE) :: cm_arr(0:n_pole_maxx, 0:n_pole_maxx)
real(rp) :: f_elec, step_len_val
logical :: has_elec_multipoles

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

! Fringe fields are handled on CPU (before/after GPU body tracking)
call init_fringe_info(fringe_info, ele)

! Get the quad gradient b1 and check for extra multipoles
call multipole_ele_to_ab(ele, .false., ix_mag_max, an, bn, magnetic$, include_kicks$, b1)
call multipole_ele_to_ab(ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

! Determine n_step for split-step integration
has_mag_multipoles = (ix_mag_max > -1)
length = bunch%particle(1)%time_dir * ele_length
n_step = 1
if (has_mag_multipoles .or. ix_elec_max > -1) &
  n_step = max(nint(abs(length) / ele%value(ds_step$)), 1)

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

! Apply entrance fringe kick on CPU (after misalignment, before body tracking)
if (fringe_info%has_fringe) then
  fringe_info%particle_at = first_track_edge$
  do j = 1, n
    if (bunch%particle(j)%state == alive$) then
      call apply_element_edge_kick(bunch%particle(j), fringe_info, ele, param, .false.)
      if (bunch%particle(j)%state /= alive$) cycle
    endif
  enddo
endif

! Precompute scaled multipole arrays for CUDA kernel
! The kernel uses a2[n], b2[n] which include charge/scale factors.
! Scale factor = r_step for full kicks (half is applied inside kernel).
! The c_multi coefficients (with no_n_fact=.true.) are precomputed.
a2_arr = 0
b2_arr = 0
ea2_arr = 0
eb2_arr = 0
cm_arr = 0
has_elec_multipoles = (ix_elec_max > -1)
r_ratio = ele%value(p0c$) / bunch%particle(1)%p0c
r_step = real(bunch%particle(1)%time_dir, rp) / n_step
step_len_val = ele_length * r_step

! Precompute scaled magnetic multipole arrays for CUDA kernel
if (has_mag_multipoles) then
  f_charge = bunch%particle(1)%direction * ele%orientation * &
             charge_to_mass_of(bunch%particle(1)%species) / charge_to_mass_of(ele%ref_species)
  do nn = 0, ix_mag_max
    a2_arr(nn) = r_ratio * an(nn) * f_charge * r_step
    b2_arr(nn) = r_ratio * bn(nn) * f_charge * r_step
  enddo
endif

! Precompute scaled electric multipole arrays for CUDA kernel
! Electric scale: charge_of(species) / (beta * p0c) * step_len
! Beta varies per particle, so we factor out 1/beta here and apply it in the kernel.
! p0c is the same for all particles at a given element.
if (has_elec_multipoles) then
  f_elec = charge_of(bunch%particle(1)%species) / bunch%particle(1)%p0c
  do nn = 0, ix_elec_max
    ea2_arr(nn) =  r_ratio * an_elec(nn) * f_elec * step_len_val
    eb2_arr(nn) = -r_ratio * bn_elec(nn) * f_elec * step_len_val
  enddo
endif

! Precompute c_multi table (shared by magnetic and electric)
nn = max(ix_mag_max, ix_elec_max)
do j = 0, nn
  do mm = 0, j
    cm_arr(j, mm) = c_multi(j, mm, .true.)
  enddo
enddo

! AoS -> SoA extraction (now in body frame, after fringe)
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

! Call CUDA kernel (with multipole and step parameters)
call gpu_track_quad(vx, vpx, vy, vpy, vz, vpz, &
                    state_a, beta_a, p0c_a, t_a, &
                    mc2, b1, ele_length, delta_ref_time, &
                    e_tot_ele, charge_dir, n, &
                    a2_arr, b2_arr, cm_arr, &
                    int(ix_mag_max, C_INT), int(n_step, C_INT), &
                    ea2_arr, eb2_arr, int(ix_elec_max, C_INT))

! SoA -> AoS write-back (still in body frame if misaligned)
do j = 1, n
  bunch%particle(j)%vec(1) = vx(j)
  bunch%particle(j)%vec(2) = vpx(j)
  bunch%particle(j)%vec(3) = vy(j)
  bunch%particle(j)%vec(4) = vpy(j)
  bunch%particle(j)%vec(5) = vz(j)
  bunch%particle(j)%vec(6) = vpz(j)
  bunch%particle(j)%state  = state_a(j)
  bunch%particle(j)%beta   = beta_a(j)
  bunch%particle(j)%t      = t_a(j)
  bunch%particle(j)%location = downstream_end$
  bunch%particle(j)%ix_ele   = ele%ix_ele
  bunch%particle(j)%ix_branch = ele%ix_branch
enddo

deallocate(vx, vpx, vy, vpy, vz, vpz)
deallocate(state_a, beta_a, p0c_a, t_a)

! Apply exit fringe kick on CPU (after body tracking, before misalignment undo)
if (fringe_info%has_fringe) then
  fringe_info%particle_at = second_track_edge$
  do j = 1, n
    if (bunch%particle(j)%state == alive$) then
      call apply_element_edge_kick(bunch%particle(j), fringe_info, ele, param, .false.)
      if (bunch%particle(j)%state /= alive$) cycle
    endif
  enddo
endif

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
