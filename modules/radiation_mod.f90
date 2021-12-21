module radiation_mod

use rad_int_common
use symp_lie_mod

contains

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!+
! Subroutine release_rad_int_cache (ix_cache)
!
! Subroutine to release the memory associated with caching wiggler values.
! See the radiation_integrals routine for further details.
!
! Input:
!   ix_cache -- Integer: Cache number.
!
! Output:
!   ix_cache -- Integer: Cache number set to 0,
!-

subroutine release_rad_int_cache (ix_cache)

implicit none

integer i, ix_cache

!

rad_int_cache_common(ix_cache)%in_use = .false.
ix_cache = 0

end subroutine release_rad_int_cache 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine calc_radiation_tracking_integrals (ele, orbit, edge, use_half_length, int_gx, int_gy, int_g2, int_g3)
!
! Routine to calculate the integrated g bending strength parameters for half the element. 
!
! g = 1/rho where rho is the radius of curvature. g points radially outward in the bending plane.
! If the particle is at the starting edge then the calculation is over the first half of the element.
! If the particle is at the exit edge then the calculation is over the 2nd half of the element.
! The calculation is in laboratory (not element body) coords.
!
! Input:
!   orbit             -- coord_struct: Particle position at the entrance edge if edge = start_edge$ and the 
!                         exit edge if edge = end_edge$.
!   ele               -- ele_struct: Element that causes radiation.
!   edge              -- integer: Where the particle is: start_edge$ or end_edge$.
!   use_half_length   -- logical: Use only half the element length in the calc?
!
! Output:
!   int_gx    -- real(rp): Integral of x-component of g.
!   int_gy    -- real(rp): Integral of y-component of g.
!   int_g2    -- real(rp): Integral of g^2.
!   int_g3    -- real(rp): Integral of g^3.
!-

subroutine calc_radiation_tracking_integrals (ele, orbit, edge, use_half_length, int_gx, int_gy, int_g2, int_g3)

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (ele_struct), pointer :: field_ele
type (coord_struct) :: orbit2

real(rp) len2, int_gx, int_gy, int_g2, int_g3, kx, ky, kx_tot, ky_tot, s_here
real(rp) eff_len, g2, g3, gx, gy, cos_t, sin_t, tilt
real(rp) a_pole_mag(0:n_pole_maxx), b_pole_mag(0:n_pole_maxx)
real(rp) a_pole_elec(0:n_pole_maxx), b_pole_elec(0:n_pole_maxx)

integer edge, direc
integer i, j, ix_mag_max, ix_elec_max

logical use_half_length, err_flag

character(*), parameter :: r_name = 'calc_radiation_tracking_integrals'

! Init

int_gx = 0; int_gy = 0
int_g2 = 0; int_g3 = 0

select case (ele%key)
case (quadrupole$, sextupole$, octupole$, sbend$, sol_quad$, wiggler$, undulator$, em_field$)
! All other types ignored.
case default
  return
end select

! The total radiation length is the element length + any change in path length.
! If entering the element then the length over which radiation is generated
! is taken to be 1/2 the element length.
! If leaving the element the radiation length is taken to be 1/2 the element length + delta_Z

if (edge == start_edge$) then
  direc = +1
  s_here = 0
elseif (edge == end_edge$) then
  direc = -1
  s_here = ele%value(l$)
else
  call out_io (s_fatal$, r_name, 'BAD EDGE ARGUMENT:', edge)
  if (global_com%exit_on_error) call err_exit
endif

! The problem with a negative element length is that it is not possible to undo the stochastic part of the radiation kick.
! In this case the best thing is to just set everything to zero

if (use_half_length) then
  eff_len = ele%value(l$) / 2
else
  eff_len = ele%value(l$)
endif

if (eff_len <= 0) return

!---------------------------------
! Calculate the radius of curvature for an on-energy particle
! Wiggler, undulator, em_field case

if (ele%key == wiggler$ .or. ele%key == undulator$ .or. ele%key == em_field$) then
  int_gx = 0
  int_gy = 0
  field_ele => pointer_to_field_ele (ele, 1)

  if (field_ele%field_calc == planar_model$) then
    g2 = abs(field_ele%value(k1y$))
    g3 = 4 * sqrt(2*g2)**3 / (3 * pi)  
    int_g2 = eff_len * g2
    int_g3 = eff_len * g3

  elseif (field_ele%field_calc == helical_model$) then
    g2 = abs(field_ele%value(k1y$))
    g3 = sqrt(g2)**3
    int_g2 = eff_len * g2
    int_g3 = eff_len * g3

  else
    if (associated(ele%rad_int_cache)) field_ele => ele
    int_g2 = eff_len * (field_ele%rad_int_cache%g2_0 + dot_product(orbit%vec(1:4)-field_ele%rad_int_cache%orb0(1:4), field_ele%rad_int_cache%dg2_dorb(1:4)))
    int_g3 = eff_len * (field_ele%rad_int_cache%g3_0 + dot_product(orbit%vec(1:4)-field_ele%rad_int_cache%orb0(1:4), field_ele%rad_int_cache%dg3_dorb(1:4)))
    if (int_g3 < 0) int_g3 = 0
  endif

  return
endif

!---------------------------------------------------------
! Everything else but wiggler, undulator, em_field

! First calculate things in the element body coords.

orbit2 = orbit
call offset_particle (ele, set$, orbit2, s_pos = s_here)
call canonical_to_angle_coords (orbit2)
orbit2%vec(1) = orbit2%vec(1) + direc * orbit2%vec(2) * eff_len / 2.0_rp ! Extrapolate to center of region 1/4 of way into element.
orbit2%vec(3) = orbit2%vec(3) + direc * orbit2%vec(4) * eff_len / 2.0_rp

call multipole_ele_to_ab (ele, .false., ix_mag_max, a_pole_mag, b_pole_mag, magnetic$, include_kicks$)
call multipole_ele_to_ab (ele, .false., ix_elec_max, a_pole_elec, b_pole_elec, electric$, include_kicks$)

kx_tot = 0
ky_tot = 0

do i = 0, ix_mag_max
  call ab_multipole_kick (a_pole_mag(i), b_pole_mag(i), i, orbit2%species, ele%orientation, orbit2, kx, ky, pole_type = magnetic$)
  kx_tot = kx_tot + kx
  ky_tot = ky_tot + ky
enddo

do i = 0, ix_elec_max
  call ab_multipole_kick (a_pole_elec(i), b_pole_elec(i), i, orbit2%species, ele%orientation, orbit2, kx, ky, pole_type = electric$)
  kx_tot = kx_tot + kx
  ky_tot = ky_tot + ky
enddo

! A positive kick means that g is negative.

select case (ele%key)
case (sbend$)
  gx = -kx_tot/ele%value(l$) + ele%value(g$)
  gy = -ky_tot/ele%value(l$)
  g2 = gx**2 + gy**2
  g3 = sqrt(g2)**3
  tilt = ele%value(ref_tilt$)

case default
  gx = -kx_tot/ele%value(l$)
  gy = -ky_tot/ele%value(l$)
  g2 = gx**2 + gy**2
  g3 = sqrt(g2)**3
  tilt = ele%value(tilt$)
end select

! Now rotate to laboratory coords

len2 = eff_len * (1.0_rp + ele%value(g$) * orbit2%vec(1))
if (len2 < 0) then   ! Can happen if g * x < -1. Effectively the particle is lost.
  orbit%state = lost_neg_x_aperture$
  return
endif

int_g2 = len2 * g2
int_g3 = len2 * g3

if (tilt == 0) then
  int_gx = len2 * gx
  int_gy = len2 * gy
else
  cos_t = cos(tilt)
  sin_t = sin(tilt)
  int_gx = len2 * (gx * cos_t - gy * sin_t)
  int_gy = len2 * (gx * sin_t + gy * cos_t)
endif

end subroutine calc_radiation_tracking_integrals

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track1_radiation (orbit, ele, edge)
!
! Subroutine to apply a kick to a particle to account for radiation dampling and/or fluctuations.
!
! For tracking through a given element, this routine should be called initially when
! the particle is at the entrance end and at the end when the particle is at the exit end, when
! the orbit is with respect to laboratory (not element body) coordinates.
! That is, each time this routine is called it applies half the radiation kick for the entire element.
!
! Also see: track1_radiation_center
! Note: This routine is called by track1.
!
! Input:
!   orbit     -- coord_struct: Particle position before radiation applied.
!   ele       -- ele_struct: Element generating radiation.
!   edge      -- integer: Where the particle is: start_edge$ or end_edge$.
!
! Output:
!   orbit     -- coord_struct: Particle position after radiation has been applied.
!-

subroutine track1_radiation (orbit, ele, edge)

use random_mod

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele

integer :: edge

real(rp) int_gx, int_gy, this_ran, mc2, int_g2, int_g3
real(rp) gamma_0, dE_p, fact_d, fact_f, q_charge2, p_spin, spin_norm(3), norm, rel_p
real(rp), parameter :: rad_fluct_const = 55.0_rp * classical_radius_factor * h_bar_planck * c_light / (24.0_rp * sqrt_3)
real(rp), parameter :: spin_const = 5.0_rp * sqrt_3 * classical_radius_factor * h_bar_planck * c_light / 16
real(rp), parameter :: damp_const = 2 * classical_radius_factor / 3
real(rp), parameter :: c1_spin = 2.0_rp / 9.0_rp, c2_spin = 8.0_rp / (5.0_rp * sqrt_3)

character(*), parameter :: r_name = 'track1_radiation'

!

if (.not. bmad_com%radiation_damping_on .and. .not. bmad_com%radiation_fluctuations_on) return

call calc_radiation_tracking_integrals (ele, orbit, edge, .true., int_gx, int_gy, int_g2, int_g3)
if (int_g2 == 0 .or. orbit%state /= alive$) return

! Apply the radiation kicks
! Basic equation is E_radiated = xi * (dE/dt) * sqrt(L) / c_light
! where xi is a random number with sigma = 1.

mc2 = mass_of(ele%ref_species)
q_charge2 = charge_of(orbit%species)**2
gamma_0 = ele%value(e_tot$) / mc2

fact_d = 0
if (bmad_com%radiation_damping_on) then
  fact_d = damp_const * q_charge2 * gamma_0**3 * int_g2 / mc2
  if (bmad_com%backwards_time_tracking_on) fact_d = -fact_d
endif

fact_f = 0
if (bmad_com%radiation_fluctuations_on) then
  call ran_gauss (this_ran)
  fact_f = sqrt(rad_fluct_const * q_charge2 * gamma_0**5 * int_g3) * this_ran / mc2
endif

rel_p = 1 + orbit%vec(6)
dE_p = rel_p * (fact_d + fact_f)
if (bmad_com%radiation_zero_average) then
  if (ele%key == sbend$ .or. ele%key == wiggler$ .or. ele%key == undulator$) dE_p = &
                                                  dE_p + ele%value(dpz_rad_damp_ave$) / (2 * rel_p)
endif
dE_p = dE_p * synch_rad_com%scale 

orbit%vec(2) = orbit%vec(2) * (1 - dE_p)
orbit%vec(4) = orbit%vec(4) * (1 - dE_p)
orbit%vec(6) = orbit%vec(6) - dE_p * rel_p

! Sokolov-Ternov Spin flip
! The equation is not correct

!if (bmad_com%spin_tracking_on .and. bmad_com%spin_sokolov_ternov_flipping_on) then
!  norm = norm2(orbit%spin)
!  if (norm /= 0) then
!    spin_norm = orbit%spin / norm
!    call ran_uniform (this_ran)
!    p_spin = (spin_const * q_charge2 * gamma_0**5 / (orbit%beta * mc2**2)) * &
!          (int_g3 - c1_spin * int_g3 * (spin_norm(3))**2 + c2_spin * int_g2 * dot_product([g_y, -g_x], spin_norm(1:2))) 
!    if (this_ran < p_spin) orbit%spin = -orbit%spin  ! spin flip
!  endif
!endif

end subroutine track1_radiation 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track1_radiation_center (orbit, ele1, ele2, rad_damp, rad_fluct)
!
! Used for elements that have been split in half: This routine applies a kick to a particle 
! to account for radiation dampling and/or fluctuations.
!
! Also see: track1_radiation.
!
! Input:
!   orbit     -- coord_struct: Particle at center of element before radiation applied.
!   ele1      -- ele_struct: First half of the split element.
!   ele2      -- ele_struct: Second half of the split element.
!   rad_damp  -- logical, optional: If present, override setting of bmad_com%radiation_damping_on.
!   rad_fluct -- logical, optional: If present, override setting of bmad_com%radiation_fluctuations_on.
!
! Output:
!   orbit     -- coord_struct: Particle position after radiation has been applied.
!-

subroutine track1_radiation_center (orbit, ele1, ele2, rad_damp, rad_fluct)

use random_mod

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele1, ele2
type (lat_param_struct) :: param

real(rp) int_gx, int_gy, this_ran, mc2, int_g2, int_g3, gxi, gyi, g2i, g3i
real(rp) gamma_0, dE_p, fact_d, fact_f, q_charge2, p_spin, spin_norm(3), norm, rel_p
real(rp), parameter :: rad_fluct_const = 55.0_rp * classical_radius_factor * h_bar_planck * c_light / (24.0_rp * sqrt_3)
real(rp), parameter :: spin_const = 5.0_rp * sqrt_3 * classical_radius_factor * h_bar_planck * c_light / 16
real(rp), parameter :: damp_const = 2 * classical_radius_factor / 3
real(rp), parameter :: c1_spin = 2.0_rp / 9.0_rp, c2_spin = 8.0_rp / (5.0_rp * sqrt_3)

logical, optional :: rad_damp, rad_fluct
logical r_damp, r_fluct
character(*), parameter :: r_name = 'track1_radiation_center'

!

r_damp  = logic_option(bmad_com%radiation_damping_on, rad_damp)
r_fluct = logic_option(bmad_com%radiation_fluctuations_on, rad_fluct)
if (.not. r_damp .and. .not. r_fluct) return

call calc_radiation_tracking_integrals (ele1, orbit, start_edge$, .false., int_gx, int_gy, int_g2, int_g3)
call calc_radiation_tracking_integrals (ele2, orbit, end_edge$, .false., gxi, gyi, g2i, g3i)
int_gx = int_gx + gxi
int_gy = int_gy + gyi
int_g2 = int_g2 + g2i
int_g3 = int_g3 + g3i

if (int_g2 == 0) return

! Apply the radiation kicks.
! Basic equation is E_radiated = xi * (dE/dt) * sqrt(L) / c_light.
! where xi is a Gaussian random number with sigma = 1.

mc2 = mass_of(ele1%ref_species)
q_charge2 = charge_of(orbit%species)**2
gamma_0 = ele1%value(e_tot$) / mc2

fact_d = 0
if (r_damp) then
  fact_d = damp_const * q_charge2 * gamma_0**3 * int_g2 / mc2
  if (bmad_com%backwards_time_tracking_on) fact_d = -fact_d
endif

fact_f = 0
if (r_fluct) then
  call ran_gauss (this_ran)
  fact_f = sqrt(rad_fluct_const * q_charge2 * gamma_0**5 * int_g3) * this_ran / mc2
endif

rel_p = 1 + orbit%vec(6) 
dE_p = rel_p * (fact_d + fact_f)
if (bmad_com%radiation_zero_average) then
  if (ele1%key == sbend$ .or. ele1%key == wiggler$ .or. ele1%key == undulator$) dE_p = &
                                dE_p + (ele1%value(dpz_rad_damp_ave$) + ele2%value(dpz_rad_damp_ave$)) / rel_p
endif
dE_p = dE_p * synch_rad_com%scale 

orbit%vec(2) = orbit%vec(2) * (1 - dE_p)
orbit%vec(4) = orbit%vec(4) * (1 - dE_p)
orbit%vec(6) = orbit%vec(6) - dE_p * rel_p

end subroutine track1_radiation_center 

end module
