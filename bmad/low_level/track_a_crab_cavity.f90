!+
! Subroutine track_a_crab_cavity (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through an crab_cavity element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: crab_cavity element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_crab_cavity (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_crab_cavity

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (ele_struct), pointer :: ref_ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) voltage, phase0, phase, t0, length, charge_dir, dt_length, beta_ref
real(rp) k_rf, dl, beta_old, pz_old, h, pc, s_here, dt_ref
real(rp) mat_2(6,6), E_old, E_new
real(rp) an(0:n_pole_maxx), bn(0:n_pole_maxx), an_elec(0:n_pole_maxx), bn_elec(0:n_pole_maxx)

integer i, n_slice, orientation
integer ix_mag_max, ix_elec_max

logical, optional :: make_matrix
logical err

character(*), parameter :: r_name = 'track_a_crab_cavity'

!

call multipole_ele_to_ab (ele, .false., ix_mag_max, an,      bn,      magnetic$, include_kicks$)
call multipole_ele_to_ab (ele, .false., ix_elec_max, an_elec, bn_elec, electric$)

ref_ele => ele
if (ref_ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) ref_ele => pointer_to_super_lord(ele)
dt_ref = ele%s_start - ref_ele%s_start

!

call offset_particle (ele, set$, orbit, mat6 = mat6, make_matrix = make_matrix)

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, 1.0_rp/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, ele%value(l$)/2, mat6, make_matrix)

length = ele%value(l$) * orbit%time_dir
n_slice = max(1, nint(length / ele%value(ds_step$))) 
!n_slice = 1
dl = length / n_slice
charge_dir = rel_tracking_charge_to_mass(orbit, param%particle) * ele%orientation
voltage = orbit%time_dir * e_accel_field(ele, voltage$, .true.) * charge_dir / (ele%value(p0c$) * n_slice)
beta_ref = ele%value(p0c$) / ele%value(e_tot$)
dt_length = length / (c_light * beta_ref)
k_rf = twopi * ele%value(rf_frequency$) / c_light
h = mass_of(orbit%species)/ele%value(p0c$)

! Track through slices.

call track_this_drift(orbit, dl/2, ele, phase, mat6, make_matrix)

do i = 1, n_slice

  ! Note: particle_rf_time is referenced to the lord element if ele is a super or slice slave.
  ! Thus we need dt_ref.
  s_here = (i - 0.5_rp) * ele%value(l$) / n_slice
  phase0 = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) + ele%value(phi0_autoscale$) - &
          (particle_rf_time (orbit, ele, .false., s_here) - rf_ref_time_offset(ele) - &
          (s_here+dt_ref)/(orbit%beta * c_light)) * ele%value(rf_frequency$))
  if (ele%orientation == -1) phase0 = phase0 + twopi * ele%value(rf_frequency$) * dt_length
  phase = phase0

  orbit%vec(2) = orbit%vec(2) + voltage * sin(phase)
  pz_old = orbit%vec(6)
  beta_old = orbit%beta

  E_old = orbit%p0c * (1.0_rp + orbit%vec(6)) / beta_old
  E_new = E_old + voltage * cos(phase) * k_rf * orbit%vec(1) * orbit%p0c
  call convert_total_energy_to (E_new, orbit%species, beta = orbit%beta, pc = pc, err_flag = err, print_err = .false.)
  if (err) then
    orbit%state = lost_pz$
    return
  endif

  orbit%vec(6) = (pc - orbit%p0c) / orbit%p0c

  if (logic_option(.false., make_matrix)) then
    mat_2 = mat6

    mat_2(2,:) = mat6(2,:) + voltage*k_rf*cos(phase)*(mat6(5,:)/beta_old - &
                 mat6(6,:)*orbit%vec(5)/beta_old**2 * h**2/(h**2+(1+pz_old)**2)**(1.5))

    mat_2(6,:) = E_new/pc * ( mat6(6,:)*beta_old + voltage*k_rf*(cos(phase)*mat6(1,:) - & 
                 sin(phase)*k_rf*orbit%vec(1)*(mat6(5,:)/beta_old - &
                 mat6(6,:)*orbit%vec(5)/beta_old**2 * h**2/(h**2+(1+pz_old)**2)**(1.5))))
    
    mat_2(5,:) = mat6(5,:)*orbit%beta/beta_old + orbit%vec(5)*(mat_2(6,:)/beta_old * h**2/(h**2+(1+orbit%vec(6))**2)**(1.5) - &
                 mat6(6,:)*orbit%beta/beta_old**2 * h**2/(h**2+(1+pz_old)**2)**(1.5))
    
    mat6 = mat_2
  endif
  
  orbit%vec(5) = orbit%vec(5) * orbit%beta / beta_old   !! Reversed betas
  
  if (i == n_slice) exit
  call track_this_drift(orbit, dl, ele, phase, mat6, make_matrix)
enddo

call track_this_drift(orbit, dl/2, ele, phase, mat6, make_matrix)

! coupler kick, multipoles, back to lab coords.

if (ix_mag_max > -1)  call ab_multipole_kicks (an,      bn,      ix_mag_max,  ele, orbit, magnetic$, 1.0_rp/2,   mat6, make_matrix)
if (ix_elec_max > -1) call ab_multipole_kicks (an_elec, bn_elec, ix_elec_max, ele, orbit, electric$, ele%value(l$)/2, mat6, make_matrix)

call offset_particle (ele, unset$, orbit, mat6 = mat6, make_matrix = make_matrix)

!-------------------------
contains

subroutine track_this_drift (orbit, dl, ele, phase, mat6, make_matrix)

type (coord_struct) orbit
type (ele_struct) ele
real(rp) mat6(6,6)
real(rp) z, dl, phase
logical make_matrix

!

z = orbit%vec(5)
call track_a_drift (orbit, dl, mat6, make_matrix)
!! phase = phase + twopi * ele%value(rf_frequency$) * (orbit%vec(5) - z) / (c_light * orbit%beta)

end subroutine track_this_drift

end subroutine
