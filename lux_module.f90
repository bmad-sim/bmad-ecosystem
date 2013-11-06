module lux_module

use random_mod
use photon_init_mod
use photon_target_mod
use track1_photon_mod
use em_field_mod
use bmad

!

type lux_bend_slice_struct
  type (twiss_struct) a, b
  type (em_field_struct) field
  type (coord_struct) orbit
  real(rp) integrated_emit_prob    ! Probability of photon emission from element start to slice
  real(rp) g_x, g_y
  real(rp) c_mat(2,2), gamma_c
end type

! %bend_slice(i) are the emission parameters at the end of the i^th slice

type lux_source_struct
  type (ele_struct), pointer :: source_ele, branch_ele, det_ele
  integer n_bend_slice                                       ! Number of slices
  type (lux_bend_slice_struct), allocatable :: bend_slice(:) ! Size: (0:n_bend_slice)
  real(rp) E_min, E_max                                      ! Photon energy range 
end type

type lux_params_struct
  real(rp) sig_E, dE_center
  real(rp) :: y_max = 0.5, phi_max = 0.7
  real(rp) :: Intensity_min_det_pixel_cutoff = 1e-6
  real(rp) :: Intensity_min_photon1_cutoff = 1e-6
  real(rp) :: stop_total_intensity = 10
  real(rp) :: emit_x, emit_y, frac_sig_E          ! Charged particle beam emittances.
  integer :: n_energy_pts = 1
  integer :: stop_num_photons = 0
  real(rp) :: window_width = 800.0_rp, window_height = 400.0_rp
  real(rp) :: e_field_x, e_field_y
  character(16) :: source_type = 'SPHERICAL'  ! or 'PLANAR' or 'BEND'
  character(16) energy_spectrum               ! 'UNIFORM' or 'GAUSSIAN'
  character(16) transverse_distribution       ! 'UNIFORM' or 'GAUSSIAN'
  logical :: use_tile = .false.
  logical :: dE_relative_to_ref = .true.
  integer :: ix_tracking_mode
  logical :: debug = .false.     ! For debugging
end type

type lux_photon_struct
  type (coord_struct), allocatable :: orb(:)
end type

!type lux_detector_pixel_struct
!  integer :: n_photon = 0
!  real(rp) :: intensity = 0
!  real(rp) :: E_sum = 0, E2_sum = 0  ! Scratch space
!  real(rp) :: E_ave = 0, E_rms = 0   ! Average and rms photon energy
!end type
!
!type lux_detector_struct
!  type (lux_detector_pixel_struct), allocatable :: pixel(:,:)
!end type

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine generate_photon (photon, lat, ix_energy, lux_param, source)
!
! Routine to generate the starting photon coordinates
!
! Input:
!   lat           -- lat_struct: Lattice.
!   ix_energy     -- Integer: Energy slice index.
!   source        -- lux_source_struct: Source parameters
!
! Ouput:
!   photon     -- lux_photon_struct: Initialized starting coords.
!-

subroutine generate_photon (photon, lat, ix_energy, lux_param, source)

use nr

implicit none

type (lat_struct) lat
type (lux_photon_struct), target :: photon
type (coord_struct) charged_orb
type (coord_struct), pointer :: orb
type (lux_params_struct) lux_param
type (lux_source_struct), target :: source
type (lux_bend_slice_struct), pointer :: sl(:)
type (ele_struct) ele
type (ele_struct), pointer :: source_ele

real(rp) x, y, phi, r(3), dir(2), ds, rr, r_emit(5), prob, f
real(rp) v_mat(4,4), v_inv_mat(4,4), vec(4), dE, e(3), b(3)
real(rp) g_bend(3), gamma_electron

integer ix, n_slice, ix_energy

!

orb => photon%orb(0)
source_ele => source%source_ele

!-----------------------------------------------------
! Bend

select case (lux_param%source_type)
case ('BEND')
  sl => source%bend_slice
  n_slice = ubound(sl, 1)

  ! Find where photon emitted

  call ran_uniform(rr)  ! longitudinal position
  call bracket_index (sl%integrated_emit_prob, 0, n_slice, rr, ix)
  if (ix == n_slice) ix = n_slice - 1
  f = (rr - sl(ix)%integrated_emit_prob) / (sl(ix+1)%integrated_emit_prob - sl(ix)%integrated_emit_prob)

  ! Calculate electron average position

  charged_orb = sl(ix)%orbit
  charged_orb%vec = (1-f) * sl(ix)%orbit%vec + f * sl(ix+1)%orbit%vec

  ele%a = average_twiss(1-f, sl(ix)%a, sl(ix+1)%a)
  ele%b = average_twiss(1-f, sl(ix)%b, sl(ix+1)%b)
  ele%c_mat   = (1-f) * sl(ix)%c_mat   + f * sl(ix+1)%c_mat
  ele%gamma_c = (1-f) * sl(ix)%gamma_c + f * sl(ix+1)%gamma_c
  call make_v_mats (ele, v_mat, v_inv_mat)

  ! Add offsets due to finite bunch size to the electron position.
  ! To do this must transform to the normal mode coords

  call ran_gauss(r_emit)  ! electron momentum offset.
  dE = r_emit(5) * lux_param%frac_sig_E
  charged_orb%vec(6) = charged_orb%vec(6) + dE

  vec = matmul (v_inv_mat, charged_orb%vec(1:4))
  vec(1:2) = vec(1:2) + charged_orb%vec(1:2) + [ele%a%eta, ele%a%etap] * dE
  vec(3:4) = vec(3:4) + charged_orb%vec(1:2) + [ele%b%eta, ele%b%etap] * dE

  vec(1) = vec(1) + sqrt(lux_param%emit_x * ele%a%beta) * r_emit(1)
  vec(2) = vec(2) + sqrt(lux_param%emit_x / ele%a%beta) * (r_emit(2) - ele%a%alpha * r_emit(1))

  vec(3) = vec(3) + sqrt(lux_param%emit_y * ele%b%beta) * r_emit(3)
  vec(4) = vec(4) + sqrt(lux_param%emit_y / ele%b%beta) * (r_emit(4) - ele%b%alpha * r_emit(3))

  charged_orb%vec(1:4) = matmul(v_mat, vec)

  ! Calculate bending strength

  B = (1-f) * sl(ix)%field%b + f * sl(ix+1)%field%b
  E = 0
  g_bend = g_bend_from_em_field (B, E, charged_orb)
  
  ! Init photon

  gamma_electron = source_ele%value(p0c$) * &
                      (1 + sl(ix)%orbit%vec(6)) / sl(ix)%orbit%beta / mass_of(sl(ix)%orbit%species)
  if (lux_param%ix_tracking_mode == coherent$ .and. ix_energy > 0) then
    rr = (ix_energy - 0.5_rp) / lux_param%n_energy_pts
    call photon_init (g_bend(1), g_bend(2), gamma_electron, orb, source%E_min, source%E_max, rr)
  else
    call photon_init (g_bend(1), g_bend(2), gamma_electron, orb, source%E_min, source%E_max, -1.0_rp)
  endif
  call absolute_photon_position (charged_orb, orb)
  orb%s = source_ele%s - (1 - rr) * source_ele%value(l$) 

  ! Track to branch element.

  ds = source%branch_ele%s - orb%s  

  if (source_ele%key == sbend$) then
    call track_a_bend_photon (orb, ele, ds)
  else
    call track_a_drift_photon (orb, ds)
  endif

  return

end select

!-----------------------------------------------------
! Spherical & Planar

! Set position

if (lux_param%transverse_distribution == 'UNIFORM') then
  call ran_uniform(r)
  rr = (2 * rr - 1) / 2.0
elseif (lux_param%transverse_distribution == 'GAUSSIAN') then
  call ran_gauss(r)
else
  print *, 'BAD LUX_PARAM%TRANSVERSE_DISTRIBUTION: ', lux_param%transverse_distribution
  call err_exit
endif

orb%vec(1) = source%source_ele%value(x_half_length$) * r(1)
orb%vec(3) = source%source_ele%value(y_half_length$) * r(2)
orb%vec(5) = source%source_ele%value(l$) * (r(3) + 0.5)

select case (lux_param%source_type)
case ('SPHERICAL')

  ! Set direction

  if (lux_param%use_tile) then
    dir = 2*dir - 1    ! Renormalize to be in range [-1, 1]

    y   = dir(1) * lux_param%y_max
    phi = dir(2) * lux_param%phi_max

    orb%vec(2) = sqrt(1 - y**2) * sin(phi)
    orb%vec(4) = y

    orb%vec(6) = sqrt(1 - orb%vec(2)**2 - orb%vec(4)**2)

  elseif (source_ele%key == x_ray_init$) then
    call isotropic_photon_emission (source_ele, lat%param, orb, +1, twopi)

  else
    print *, 'NOT YET IMPLEMENTED...'
    stop
  endif

case ('PLANAR')
  orb%vec(2:4:2) = lat%beam_start%vec(2:4:2)
  orb%vec(6) = sqrt(1 - orb%vec(2)**2 - orb%vec(4)**2)

end select

! Set energy

if (lux_param%ix_tracking_mode == coherent$ .and. ix_energy > 0) then
  if (lux_param%energy_spectrum == 'UNIFORM') then
    rr = (ix_energy - 0.5_rp) / lux_param%n_energy_pts
  else if (lux_param%energy_spectrum == 'GAUSSIAN') then
    rr = sqrt_2 * erfc((ix_energy - 0.5_rp) / lux_param%n_energy_pts) 
  else
    print *, 'BAD LUX_PARAM%ENERGY_SPECTRUM SETTING: ', lux_param%energy_spectrum 
  endif

else
  if (lux_param%energy_spectrum == 'UNIFORM') then
    call ran_uniform(rr)
    rr = (2 * rr - 1) / 2.0
  else if (lux_param%energy_spectrum == 'GAUSSIAN') then
    call ran_gauss(rr)
  else
    print *, 'BAD LUX_PARAM%ENERGY_SPECTRUM SETTING: ', lux_param%energy_spectrum 
  endif

endif

orb%p0c = lux_param%sig_E * rr + lux_param%dE_center
if (lux_param%dE_relative_to_ref) orb%p0c = orb%p0c + source%source_ele%value(p0c$) 

call init_coord (orb, orb%vec, source%source_ele, .false., photon$, 1, orb%p0c) 
orb%s = orb%vec(5) + orb%s - source%source_ele%value(l$) + source%source_ele%value(z_offset_tot$)
orb%t = 0

x = lux_param%e_field_x; y = lux_param%e_field_y
if (x == 0 .and. y == 0) then
  call ran_uniform(rr)
  orb%field(1) = cos(twopi * rr)
  orb%field(2) = sin(twopi * rr)
else
 orb%field(1) = x / sqrt(x**2 + y**2)
 orb%field(2) = y / sqrt(x**2 + y**2)
endif

! Translate from element to lab coordinates
! and track to entrance end of source%source_ele

call offset_photon (source%source_ele, orb, unset$)

call track_a_drift_photon (orb, -orb%s)

end subroutine generate_photon

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine lux_setup (photon, lat, lux_param, source)
!
! Routine 
!
! Input:
!   photon      -- lux_photon_struct
!   lat         -- lat_struct
!   lux_param   -- lux_params_struct
!
! Output:
!   source      -- lux_source_struct
!-

subroutine lux_setup (photon, lat, lux_param, source)

implicit none

type (lat_struct), target :: lat
type (lux_params_struct) lux_param
type (lux_source_struct), target :: source
type (lux_photon_struct), target :: photon
type (ele_struct) twiss_ele
type (ele_struct), pointer :: ele
type (lux_bend_slice_struct), pointer :: sl(:)
type (coord_struct) orb
type (coord_struct), pointer :: orbit
type (branch_struct), pointer :: branch

real(rp) vz, rho
real(rp) phi, y, ds, s_now, prob1, prob2, g_bend(3), g_abs
real(rp) emit_prob, old_emit_prob, gamma

integer i, j, k, n, n_phi, n_y, ip, iy, iz, ie, iz2
integer track_state, iy0, iy1, ip0, ip1, iz0, iz1
integer n_slice, n_z, ix

logical err

!-------------------------------------------------------------
! Init

branch => lat%branch(source%det_ele%ix_branch)

!-------------------------------------------------------------
! Lattice has sample element.

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)
  select case (ele%key)
  case (sample$, diffraction_plate$) 
    call photon_target_setup (ele)
  case (x_ray_init$)
    if (.not. lux_param%use_tile) call photon_target_setup (ele)
  end select
enddo

!-------------------------------------------------------------
! Sbend source

if (lux_param%source_type == 'BEND') then
  n_slice = max(1, nint(source%source_ele%value(l$) / source%source_ele%value(ds_step$)))
  source%n_bend_slice = n_slice
  allocate (source%bend_slice(0:n_slice))

  source%E_min = lux_param%dE_center - lux_param%sig_E
  if (lux_param%dE_relative_to_ref) source%E_min = source%E_min + source%det_ele%value(p0c$) 
  if (lux_param%sig_E == 0) then
    source%E_max = source%E_min + 1d-10  ! Need some small offset for the calculation
  else
    source%E_max = source%E_min + 2 * lux_param%sig_E
  endif

  ! Track through source ele and gather data

  call init_coord (orb, lat%beam_start, source%source_ele, .false.)
  twiss_ele = pointer_to_next_ele (source%source_ele, -1)
  ds = source%source_ele%value(l$) / n_slice
  s_now = 0
  gamma = (orb%p0c / orb%beta) / mass_of(orb%species)
  sl => source%bend_slice
  sl(0)%integrated_emit_prob = 0

  do i = 0, n_slice
    sl(i)%a = twiss_ele%a
    sl(i)%b = twiss_ele%b
    sl(i)%c_mat   = twiss_ele%c_mat
    sl(i)%gamma_c = twiss_ele%gamma_c
    sl(i)%orbit   = orb
    call em_field_calc (source%source_ele, lat%param, s_now, 0.0_rp, orb, .false., sl(i)%field)
    g_bend = g_bend_from_em_field (sl(i)%field%b, sl(i)%field%e, orb)
    g_abs = norm2(g_bend)
    prob1 = photon_energy_integ_prob (source%E_min, g_abs, gamma)
    prob2 = photon_energy_integ_prob (source%E_max, g_abs, gamma)
    emit_prob = g_abs * (prob2 - prob1)
    if (i /= 0) sl(i)%integrated_emit_prob = sl(i)%integrated_emit_prob + (old_emit_prob + emit_prob) / 2
    old_emit_prob = emit_prob 

    if (i == n_slice) exit
    call twiss_and_track_intra_ele (source%source_ele, lat%param, s_now, s_now+ds, &
              .true., .true., orb, orb, twiss_ele, twiss_ele, err)
    if (err) call err_exit
    s_now = s_now + ds
  enddo

  sl%integrated_emit_prob = sl%integrated_emit_prob / sl(n_slice)%integrated_emit_prob

  return
endif

end subroutine lux_setup 

end module
