module gpu_tracking_mod

use bmad_struct
use bmad_routine_interface, dummy_gtm => track_a_drift

implicit none
private

public :: gpu_tracking_init
public :: ele_gpu_eligible
public :: track_bunch_thru_drift_gpu
public :: track_bunch_thru_quad_gpu
public :: track_bunch_thru_bend_gpu
public :: track_bunch_thru_lcavity_gpu

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

  subroutine gpu_track_bend(vx, vpx, vy, vpy, vz, vpz, &
                            state, beta, p0c, t_time, &
                            mc2, g, g_tot, dg, b1, &
                            ele_length, delta_ref_time, e_tot_ele, &
                            rel_charge_dir, charge_dir_for_multipole, &
                            p0c_ele, n_particles, &
                            a2_arr, b2_arr, cm_arr, &
                            ix_mag_max, n_step, &
                            ea2_arr, eb2_arr, ix_elec_max) bind(C, name='gpu_track_bend_')
    use, intrinsic :: iso_c_binding
    real(C_DOUBLE), intent(inout) :: vx(*), vpx(*), vy(*), vpy(*), vz(*), vpz(*)
    integer(C_INT), intent(inout) :: state(*)
    real(C_DOUBLE), intent(inout) :: beta(*), p0c(*), t_time(*)
    real(C_DOUBLE), value, intent(in) :: mc2, g, g_tot, dg, b1
    real(C_DOUBLE), value, intent(in) :: ele_length, delta_ref_time, e_tot_ele
    real(C_DOUBLE), value, intent(in) :: rel_charge_dir, charge_dir_for_multipole, p0c_ele
    integer(C_INT), value, intent(in) :: n_particles
    real(C_DOUBLE), intent(in) :: a2_arr(*), b2_arr(*), cm_arr(*)
    integer(C_INT), value, intent(in) :: ix_mag_max, n_step
    real(C_DOUBLE), intent(in) :: ea2_arr(*), eb2_arr(*)
    integer(C_INT), value, intent(in) :: ix_elec_max
  end subroutine

  subroutine gpu_track_lcavity(vx, vpx, vy, vpy, vz, vpz, &
                               state, beta, p0c, t_time, &
                               mc2, &
                               step_s0, step_s, step_p0c, step_p1c, &
                               step_scale, step_time, &
                               n_rf_steps, &
                               voltage, voltage_err, field_autoscale, &
                               rf_frequency, phi0_total, &
                               voltage_tot, l_active, &
                               cavity_type, &
                               fringe_at, charge_ratio, &
                               n_particles) bind(C, name='gpu_track_lcavity_')
    use, intrinsic :: iso_c_binding
    real(C_DOUBLE), intent(inout) :: vx(*), vpx(*), vy(*), vpy(*), vz(*), vpz(*)
    integer(C_INT), intent(inout) :: state(*)
    real(C_DOUBLE), intent(inout) :: beta(*), p0c(*), t_time(*)
    real(C_DOUBLE), value, intent(in) :: mc2
    real(C_DOUBLE), intent(in) :: step_s0(*), step_s(*), step_p0c(*), step_p1c(*)
    real(C_DOUBLE), intent(in) :: step_scale(*), step_time(*)
    integer(C_INT), value, intent(in) :: n_rf_steps
    real(C_DOUBLE), value, intent(in) :: voltage, voltage_err, field_autoscale
    real(C_DOUBLE), value, intent(in) :: rf_frequency, phi0_total
    real(C_DOUBLE), value, intent(in) :: voltage_tot, l_active
    integer(C_INT), value, intent(in) :: cavity_type
    integer(C_INT), value, intent(in) :: fringe_at
    real(C_DOUBLE), value, intent(in) :: charge_ratio
    integer(C_INT), value, intent(in) :: n_particles
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
! Currently supported element types: drift, quadrupole, sbend, lcavity.
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
case (drift$, quadrupole$, sbend$, lcavity$)
  eligible = .true.
end select

end function ele_gpu_eligible

!------------------------------------------------------------------------
! bunch_to_soa — extract particle data from AoS bunch to SoA arrays
!------------------------------------------------------------------------
subroutine bunch_to_soa(bunch, n, vx, vpx, vy, vpy, vz, vpz, state_a, beta_a, p0c_a, t_a)

use, intrinsic :: iso_c_binding

type (bunch_struct), intent(in) :: bunch
integer,             intent(in) :: n
real(C_DOUBLE),      intent(out) :: vx(n), vpx(n), vy(n), vpy(n), vz(n), vpz(n)
real(C_DOUBLE),      intent(out) :: beta_a(n), p0c_a(n), t_a(n)
integer(C_INT),      intent(out) :: state_a(n)

integer :: j

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

end subroutine bunch_to_soa

!------------------------------------------------------------------------
! soa_to_bunch — write SoA arrays back to AoS bunch structure
!
! copy_beta/copy_p0c control optional write-back of beta and p0c.
! lcavity: both true. quad/bend: copy_beta true. drift: both false.
!------------------------------------------------------------------------
subroutine soa_to_bunch(bunch, ele, n, vx, vpx, vy, vpy, vz, vpz, state_a, beta_a, p0c_a, t_a, &
                         copy_beta, copy_p0c)

use, intrinsic :: iso_c_binding

type (bunch_struct), intent(inout) :: bunch
type (ele_struct),   intent(in)    :: ele
integer,             intent(in)    :: n
real(C_DOUBLE),      intent(in)    :: vx(n), vpx(n), vy(n), vpy(n), vz(n), vpz(n)
real(C_DOUBLE),      intent(in)    :: beta_a(n), p0c_a(n), t_a(n)
integer(C_INT),      intent(in)    :: state_a(n)
logical,             intent(in)    :: copy_beta, copy_p0c

integer :: j

do j = 1, n
  bunch%particle(j)%vec(1)    = vx(j)
  bunch%particle(j)%vec(2)    = vpx(j)
  bunch%particle(j)%vec(3)    = vy(j)
  bunch%particle(j)%vec(4)    = vpy(j)
  bunch%particle(j)%vec(5)    = vz(j)
  bunch%particle(j)%vec(6)    = vpz(j)
  bunch%particle(j)%state     = state_a(j)
  bunch%particle(j)%t         = t_a(j)
  if (copy_beta) bunch%particle(j)%beta = beta_a(j)
  if (copy_p0c)  bunch%particle(j)%p0c  = p0c_a(j)
  bunch%particle(j)%location  = downstream_end$
  bunch%particle(j)%ix_ele    = ele%ix_ele
  bunch%particle(j)%ix_branch = ele%ix_branch
enddo

end subroutine soa_to_bunch

!------------------------------------------------------------------------
! apply_misalign_to_bunch — apply misalignment transform to all alive particles
!------------------------------------------------------------------------
subroutine apply_misalign_to_bunch(bunch, ele, n, set_or_unset)

type (bunch_struct), intent(inout) :: bunch
type (ele_struct),   intent(in)    :: ele
integer,             intent(in)    :: n
integer,             intent(in)    :: set_or_unset  ! set$ or unset$

integer :: j

do j = 1, n
  if (bunch%particle(j)%state == alive$) then
    call offset_particle(ele, set_or_unset, bunch%particle(j), set_hvkicks = .false.)
  endif
enddo

end subroutine apply_misalign_to_bunch

!------------------------------------------------------------------------
! apply_fringe_to_bunch — apply fringe kicks to all alive particles
!------------------------------------------------------------------------
subroutine apply_fringe_to_bunch(bunch, ele, param, n, fringe_info, particle_at)

type (bunch_struct),              intent(inout) :: bunch
type (ele_struct),                intent(in)    :: ele
type (lat_param_struct),          intent(in)    :: param
integer,                          intent(in)    :: n
type (fringe_field_info_struct),  intent(inout) :: fringe_info
integer,                          intent(in)    :: particle_at

integer :: j

fringe_info%particle_at = particle_at
do j = 1, n
  if (bunch%particle(j)%state == alive$) then
    call apply_element_edge_kick(bunch%particle(j), fringe_info, ele, param, .false.)
    if (bunch%particle(j)%state /= alive$) cycle
  endif
enddo

end subroutine apply_fringe_to_bunch

!------------------------------------------------------------------------
! precompute_multipole_arrays — compute scaled multipole coefficients for CUDA
!
! Computes a2/b2 (magnetic), ea2/eb2 (electric), and cm (c_multi) arrays
! from the raw multipole data.  These include all element-level scaling
! factors; per-particle factors (1/beta for electric, (1+g*x) for bends)
! are applied in the CUDA kernels.
!------------------------------------------------------------------------
subroutine precompute_multipole_arrays(particle1, ele, &
    ix_mag_max, an, bn, ix_elec_max, an_elec, bn_elec, &
    ele_length, n_step, &
    a2_arr, b2_arr, ea2_arr, eb2_arr, cm_arr)

use, intrinsic :: iso_c_binding

type (coord_struct),  intent(in)  :: particle1
type (ele_struct),    intent(in)  :: ele
integer,              intent(in)  :: ix_mag_max, ix_elec_max, n_step
real(rp),             intent(in)  :: an(0:), bn(0:)
real(rp),             intent(in)  :: an_elec(0:), bn_elec(0:)
real(rp),             intent(in)  :: ele_length
real(C_DOUBLE),       intent(out) :: a2_arr(0:n_pole_maxx), b2_arr(0:n_pole_maxx)
real(C_DOUBLE),       intent(out) :: ea2_arr(0:n_pole_maxx), eb2_arr(0:n_pole_maxx)
real(C_DOUBLE),       intent(out) :: cm_arr(0:n_pole_maxx, 0:n_pole_maxx)

integer :: nn, mm
real(rp) :: r_ratio, r_step, step_len_val, f_charge, f_elec

a2_arr = 0; b2_arr = 0; ea2_arr = 0; eb2_arr = 0; cm_arr = 0

r_ratio = ele%value(p0c$) / particle1%p0c
r_step = 1.0_rp / n_step
step_len_val = ele_length / n_step

! Magnetic multipoles
if (ix_mag_max > -1) then
  f_charge = particle1%direction * ele%orientation * &
             charge_to_mass_of(particle1%species) / charge_to_mass_of(ele%ref_species)
  do nn = 0, ix_mag_max
    a2_arr(nn) = r_ratio * an(nn) * f_charge * r_step
    b2_arr(nn) = r_ratio * bn(nn) * f_charge * r_step
  enddo
endif

! Electric multipoles (1/beta applied per-particle in kernel)
if (ix_elec_max > -1) then
  f_elec = charge_of(particle1%species) / particle1%p0c
  do nn = 0, ix_elec_max
    ea2_arr(nn) =  r_ratio * an_elec(nn) * f_elec * step_len_val
    eb2_arr(nn) = -r_ratio * bn_elec(nn) * f_elec * step_len_val
  enddo
endif

! c_multi coefficient table
do nn = 0, max(ix_mag_max, ix_elec_max)
  do mm = 0, nn
    cm_arr(nn, mm) = c_multi(nn, mm, .true.)
  enddo
enddo

end subroutine precompute_multipole_arrays

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
call bunch_to_soa(bunch, n, vx, vpx, vy, vpy, vz, vpz, state_a, beta_a, p0c_a, t_a)
do j = 1, n
  s_a(j) = bunch%particle(j)%s
enddo

! Call CUDA kernel
call gpu_track_drift(vx, vpx, vy, vpy, vz, vpz, &
                     state_a, beta_a, p0c_a, s_a, t_a, &
                     mc2, length, n)

! SoA -> AoS write-back
call soa_to_bunch(bunch, ele, n, vx, vpx, vy, vpy, vz, vpz, state_a, beta_a, p0c_a, t_a, &
                   .false., .false.)
do j = 1, n
  bunch%particle(j)%s = s_a(j)
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

! Apply entrance misalignment and fringe on CPU
if (has_misalign) call apply_misalign_to_bunch(bunch, ele, n, set$)
if (fringe_info%has_fringe) call apply_fringe_to_bunch(bunch, ele, param, n, fringe_info, first_track_edge$)

! Precompute scaled multipole arrays for CUDA kernel
call precompute_multipole_arrays(bunch%particle(1), ele, &
    ix_mag_max, an, bn, ix_elec_max, an_elec, bn_elec, &
    ele_length, n_step, a2_arr, b2_arr, ea2_arr, eb2_arr, cm_arr)

! AoS -> SoA extraction (now in body frame, after fringe)
call bunch_to_soa(bunch, n, vx, vpx, vy, vpy, vz, vpz, state_a, beta_a, p0c_a, t_a)

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
call soa_to_bunch(bunch, ele, n, vx, vpx, vy, vpy, vz, vpz, state_a, beta_a, p0c_a, t_a, &
                   .true., .false.)

deallocate(vx, vpx, vy, vpy, vz, vpz)
deallocate(state_a, beta_a, p0c_a, t_a)

! Apply exit fringe and misalignment on CPU
if (fringe_info%has_fringe) call apply_fringe_to_bunch(bunch, ele, param, n, fringe_info, second_track_edge$)
if (has_misalign) call apply_misalign_to_bunch(bunch, ele, n, unset$)
#endif

end subroutine track_bunch_thru_quad_gpu

!------------------------------------------------------------------------
! track_bunch_thru_bend_gpu
!
! GPU batch tracking through a bend (sbend).
! CPU sandwich pattern: misalignment + fringe on CPU, body on GPU.
! Handles all three body paths: general bend, k1 map, and drift fallback.
!------------------------------------------------------------------------
subroutine track_bunch_thru_bend_gpu (bunch, ele, param, did_track)

use multipole_mod, only: ab_multipole_kicks
use, intrinsic :: iso_c_binding

type (bunch_struct),     intent(inout) :: bunch
type (ele_struct), target, intent(inout) :: ele
type (lat_param_struct), intent(inout) :: param
logical,                 intent(out)   :: did_track

#ifdef USE_GPU_TRACKING
integer, parameter :: n_multi = n_pole_maxx + 1
integer(C_INT) :: n
integer :: j, nn, mm, ix_mag_max, ix_elec_max, n_step
real(rp) :: ele_length, mc2, b1, delta_ref_time, e_tot_ele, p0c_ele
real(rp) :: g, g_tot, dg, rel_charge_dir, c_dir
real(rp) :: r_step, length, step_len_val
real(rp) :: an(0:n_pole_maxx), bn(0:n_pole_maxx)
real(rp) :: an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)
real(rp) :: f_charge, r_ratio, f_elec, f_p0c
type (fringe_field_info_struct) :: fringe_info
logical :: has_misalign, has_mag_multipoles, has_elec_multipoles
type (coord_struct) :: start_orb

real(C_DOUBLE) :: a2_arr(0:n_pole_maxx), b2_arr(0:n_pole_maxx)
real(C_DOUBLE) :: ea2_arr(0:n_pole_maxx), eb2_arr(0:n_pole_maxx)
real(C_DOUBLE) :: cm_arr(0:n_pole_maxx, 0:n_pole_maxx)

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
p0c_ele = ele%value(p0c$)

! Bail if exact_multipoles is on — too complex for GPU
if (nint(ele%value(exact_multipoles$)) /= off$) return

has_misalign = ele%bookkeeping_state%has_misalign

! Compute charge/direction factors
rel_charge_dir = ele%orientation * bunch%particle(1)%direction * &
                 rel_tracking_charge_to_mass(bunch%particle(1), param%particle)
c_dir = ele%orientation * bunch%particle(1)%direction * charge_of(bunch%particle(1)%species)

! Get multipoles
call multipole_ele_to_ab(ele, .false., ix_mag_max, an, bn, magnetic$, include_kicks$, b1)
b1 = b1 * rel_charge_dir
if (abs(b1) < 1d-10) then
  bn(1) = b1
  b1 = 0
endif

call multipole_ele_to_ab(ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

! Compute g and g_tot
g = ele%value(g$)
length = bunch%particle(1)%time_dir * ele_length
if (length == 0) then
  dg = 0
else
  dg = bn(0) / ele_length
  bn(0) = 0
endif
g_tot = (g + dg) * rel_charge_dir

! Determine n_step
has_mag_multipoles = (ix_mag_max > -1)
has_elec_multipoles = (ix_elec_max > -1)
n_step = 1
if (has_mag_multipoles .or. has_elec_multipoles) &
  n_step = max(nint(abs(length) / ele%value(ds_step$)), 1)
r_step = real(bunch%particle(1)%time_dir, rp) / n_step
step_len_val = ele_length / n_step

! Fringe info
call init_fringe_info(fringe_info, ele)

! Allocate SoA arrays
allocate(vx(n), vpx(n), vy(n), vpy(n), vz(n), vpz(n))
allocate(state_a(n), beta_a(n), p0c_a(n), t_a(n))

! Apply entrance misalignment and fringe on CPU
if (has_misalign) call apply_misalign_to_bunch(bunch, ele, n, set$)
if (fringe_info%has_fringe) call apply_fringe_to_bunch(bunch, ele, param, n, fringe_info, first_track_edge$)

! Precompute scaled multipole arrays for CUDA kernel
! Note: the (1+g*x) curvature factor for bends is applied per-particle in the kernel.
call precompute_multipole_arrays(bunch%particle(1), ele, &
    ix_mag_max, an, bn, ix_elec_max, an_elec, bn_elec, &
    ele_length, n_step, a2_arr, b2_arr, ea2_arr, eb2_arr, cm_arr)

! AoS -> SoA
call bunch_to_soa(bunch, n, vx, vpx, vy, vpy, vz, vpz, state_a, beta_a, p0c_a, t_a)

did_track = .true.

! Call CUDA kernel
call gpu_track_bend(vx, vpx, vy, vpy, vz, vpz, &
                    state_a, beta_a, p0c_a, t_a, &
                    mc2, g, g_tot, dg, b1, &
                    ele_length, delta_ref_time, e_tot_ele, &
                    rel_charge_dir, c_dir, p0c_ele, n, &
                    a2_arr, b2_arr, cm_arr, &
                    int(ix_mag_max, C_INT), int(n_step, C_INT), &
                    ea2_arr, eb2_arr, int(ix_elec_max, C_INT))

! SoA -> AoS
call soa_to_bunch(bunch, ele, n, vx, vpx, vy, vpy, vz, vpz, state_a, beta_a, p0c_a, t_a, &
                   .true., .false.)

deallocate(vx, vpx, vy, vpy, vz, vpz)
deallocate(state_a, beta_a, p0c_a, t_a)

! Apply exit fringe and misalignment on CPU
if (fringe_info%has_fringe) call apply_fringe_to_bunch(bunch, ele, param, n, fringe_info, second_track_edge$)
if (has_misalign) call apply_misalign_to_bunch(bunch, ele, n, unset$)

! Update s position
do j = 1, n
  if (bunch%particle(j)%state == alive$) then
    bunch%particle(j)%s = ele%s
  endif
enddo
#endif

end subroutine track_bunch_thru_bend_gpu

!------------------------------------------------------------------------
! track_bunch_thru_lcavity_gpu
!
! GPU batch tracking through an lcavity (linac cavity).
! Handles the stair-step RF approximation with energy kicks,
! ponderomotive transverse kicks (standing wave), and coordinate
! transformations.
!
! CPU sandwich: misalignment. Fringe kicks handled on GPU.
! Falls back to CPU if: multipoles present, solenoid (ks/=0),
!   zero length, zero rf_frequency, absolute time tracking,
!   coupler kicks (coupler_strength /= 0).
!------------------------------------------------------------------------
subroutine track_bunch_thru_lcavity_gpu (bunch, ele, param, did_track)

use, intrinsic :: iso_c_binding

type (bunch_struct),     intent(inout) :: bunch
type (ele_struct), target, intent(inout) :: ele
type (lat_param_struct), intent(inout) :: param
logical,                 intent(out)   :: did_track

#ifdef USE_GPU_TRACKING
integer(C_INT) :: n
integer :: j, nn, ix_mag_max, ix_elec_max, n_steps
real(rp) :: mc2, phi0_total
real(rp) :: an(0:n_pole_maxx), bn(0:n_pole_maxx)
real(rp) :: an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)
logical :: has_misalign
integer :: i_fringe_at
real(rp) :: charge_ratio_val
type (ele_struct), pointer :: lord
type (rf_stair_step_struct), pointer :: step

real(C_DOUBLE), allocatable :: vx(:), vpx(:), vy(:), vpy(:), vz(:), vpz(:)
real(C_DOUBLE), allocatable :: beta_a(:), p0c_a(:), t_a(:)
integer(C_INT), allocatable :: state_a(:)
real(C_DOUBLE), allocatable :: h_step_s0(:), h_step_s(:)
real(C_DOUBLE), allocatable :: h_step_p0c(:), h_step_p1c(:)
real(C_DOUBLE), allocatable :: h_step_scale(:), h_step_time(:)

did_track = .false.

n = size(bunch%particle)

! --- Bail out conditions ---

! Zero length → track on CPU
if (ele%value(l$) == 0) return

! Zero rf_frequency with non-zero voltage → CPU
if (ele%value(rf_frequency$) == 0) return

! Absolute time tracking → CPU (phase computation differs)
if (bmad_com%absolute_time_tracking) return

! Get the super lord for RF step data
lord => pointer_to_super_lord(ele)

! Solenoid → CPU (step_drift uses solenoid_track_and_mat)
if (lord%value(ks$) /= 0) return

! Multipoles → CPU (scale varies per step, complex interaction)
call multipole_ele_to_ab(ele, .false., ix_mag_max, an, bn, magnetic$, include_kicks$)
call multipole_ele_to_ab(ele, .false., ix_elec_max, an_elec, bn_elec, electric$)
if (ix_mag_max > -1) return
if (ix_elec_max > -1) return

! Check for coupler kicks → CPU
if (ele%value(coupler_strength$) /= 0) return

! Compute fringe parameters (fringe is handled on GPU)
if (nint(lord%value(fringe_type$)) == none$) then
  i_fringe_at = 0
else
  i_fringe_at = nint(lord%value(fringe_at$))
  if (i_fringe_at < 1 .or. i_fringe_at > 3) i_fringe_at = 0
endif
charge_ratio_val = charge_of(bunch%particle(1)%species) / (2.0_rp * charge_of(lord%ref_species))

mc2 = mass_of(bunch%particle(1)%species)
n_steps = nint(lord%value(n_rf_steps$))
has_misalign = ele%bookkeeping_state%has_misalign

! Compute total phase offset (relative time tracking, non-multipass)
phi0_total = lord%value(phi0$) + lord%value(phi0_err$) + lord%value(phi0_multipass$)

! Extract step data from the lord's RF step array (indices 0..n_steps+1)
allocate(h_step_s0(n_steps+2), h_step_s(n_steps+2))
allocate(h_step_p0c(n_steps+2), h_step_p1c(n_steps+2))
allocate(h_step_scale(n_steps+2), h_step_time(n_steps+2))

do j = 0, n_steps + 1
  step => lord%rf%steps(j)
  h_step_s0(j+1)    = step%s0
  h_step_s(j+1)     = step%s
  h_step_p0c(j+1)   = step%p0c
  h_step_p1c(j+1)   = step%p1c
  h_step_scale(j+1) = step%scale
  h_step_time(j+1)  = step%time
enddo

! Allocate SoA arrays
allocate(vx(n), vpx(n), vy(n), vpy(n), vz(n), vpz(n))
allocate(state_a(n), beta_a(n), p0c_a(n), t_a(n))

! Apply entrance misalignment (lab → body frame) on CPU
if (has_misalign) call apply_misalign_to_bunch(bunch, ele, n, set$)

! AoS -> SoA extraction
call bunch_to_soa(bunch, n, vx, vpx, vy, vpy, vz, vpz, state_a, beta_a, p0c_a, t_a)

did_track = .true.

! Call CUDA kernel
call gpu_track_lcavity(vx, vpx, vy, vpy, vz, vpz, &
                       state_a, beta_a, p0c_a, t_a, &
                       mc2, &
                       h_step_s0, h_step_s, h_step_p0c, h_step_p1c, &
                       h_step_scale, h_step_time, &
                       int(n_steps, C_INT), &
                       lord%value(voltage$), lord%value(voltage_err$), &
                       lord%value(field_autoscale$), &
                       ele%value(rf_frequency$), phi0_total, &
                       lord%value(voltage_tot$), lord%value(l_active$), &
                       int(nint(lord%value(cavity_type$)), C_INT), &
                       int(i_fringe_at, C_INT), charge_ratio_val, &
                       int(n, C_INT))

! SoA -> AoS write-back
call soa_to_bunch(bunch, ele, n, vx, vpx, vy, vpy, vz, vpz, state_a, beta_a, p0c_a, t_a, &
                   .true., .true.)

deallocate(vx, vpx, vy, vpy, vz, vpz)
deallocate(state_a, beta_a, p0c_a, t_a)
deallocate(h_step_s0, h_step_s, h_step_p0c, h_step_p1c, h_step_scale, h_step_time)

! Apply exit misalignment (body → lab frame) on CPU
if (has_misalign) call apply_misalign_to_bunch(bunch, ele, n, unset$)

! Update s position
do j = 1, n
  if (bunch%particle(j)%state == alive$) then
    bunch%particle(j)%s = ele%s
  endif
enddo
#endif

end subroutine track_bunch_thru_lcavity_gpu

end module gpu_tracking_mod
