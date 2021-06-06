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
! Subroutine calc_radiation_tracking_integrals (ele, orbit, param, edge, use_half_length, int_gx, int_gy, int_g2, int_g3)
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

subroutine calc_radiation_tracking_integrals (ele, orbit, param, edge, use_half_length, int_gx, int_gy, int_g2, int_g3)

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (coord_struct) :: orbit2
type (coord_struct) start0_orb, start_orb, end_orb
type (track_struct), save, target :: track_save
type (track_struct), pointer :: track

real(rp) len2, int_gx, int_gy, int_g2, int_g3, kx, ky, kx_tot, ky_tot, s_here
real(rp) eff_len, g2, g3, gx, gy, cos_t, sin_t, tilt
real(rp) a_pole_mag(0:n_pole_maxx), b_pole_mag(0:n_pole_maxx)
real(rp) a_pole_elec(0:n_pole_maxx), b_pole_elec(0:n_pole_maxx)
real(rp), parameter :: del_orb = 1d-4

integer edge, direc, track_method_saved
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

  if (ele%field_calc == planar_model$) then
    g2 = abs(ele%value(k1y$))
    g3 = 4 * sqrt(2*g2)**3 / (3 * pi)  
    int_g2 = eff_len * g2
    int_g3 = eff_len * g3

  elseif (ele%field_calc == helical_model$) then
    g2 = abs(ele%value(k1y$))
    g3 = sqrt(g2)**3
    int_g2 = eff_len * g2
    int_g3 = eff_len * g3

  else
    if (.not. associated(ele%rad_int_cache) .or. ele%rad_int_cache%stale) then
      if (.not. associated(ele%rad_int_cache)) allocate (ele%rad_int_cache)

      ! If %map_ref_orb_in has not been set (can happen before the closed orbit is 
      ! calculated or during optimization), use %time_ref_orb_in.
      if (ele%map_ref_orb_in%state == not_set$ .or. ele%map_ref_orb_in%p0c /= ele%value(p0c$)) then
        start0_orb = ele%time_ref_orb_in
      else
        start0_orb = ele%map_ref_orb_in
      endif
      ele%rad_int_cache%orb0 = start0_orb%vec

      if (global_com%be_thread_safe) then
        allocate(track)
      else
        track => track_save
      endif

      track%n_pt = -1
      track_method_saved = ele%tracking_method
      if (ele%tracking_method == taylor$) ele%tracking_method = runge_kutta$
      call track1 (start0_orb, ele, param, end_orb, track, err_flag, .true.)
      call calc_g (track, ele%rad_int_cache%g2_0, ele%rad_int_cache%g3_0)

      do j = 1, 4
        start_orb = start0_orb
        start_orb%vec(j) = start_orb%vec(j) + del_orb
        track%n_pt = -1
        call track1 (start_orb, ele, param, end_orb, track, err_flag, .true.)
        call calc_g (track, g2, g3)
        ele%rad_int_cache%dg2_dorb(j) = (g2 - ele%rad_int_cache%g2_0) / del_orb
        ele%rad_int_cache%dg3_dorb(j) = (g3 - ele%rad_int_cache%g3_0) / del_orb
      enddo

      ele%rad_int_cache%stale = .false.
      ele%tracking_method = track_method_saved

      if (global_com%be_thread_safe) then
        deallocate(track)
      endif
    endif

    int_g2 = eff_len * (ele%rad_int_cache%g2_0 + dot_product(orbit%vec(1:4)-ele%rad_int_cache%orb0(1:4), ele%rad_int_cache%dg2_dorb(1:4)))
    int_g3 = eff_len * (ele%rad_int_cache%g3_0 + dot_product(orbit%vec(1:4)-ele%rad_int_cache%orb0(1:4), ele%rad_int_cache%dg3_dorb(1:4)))
    if (int_g3 < 0) int_g3 = 0
  endif

  return
endif

!---------------------------------------------------------
! Everything else but wiggler, undulator, em_field

! First calculate things in the element body coords.

orbit2 = orbit
call offset_particle (ele, param, set$, orbit2, s_pos = s_here)
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
  gx = -kx_tot/ele%value(l$) + ele%value(g$) + ele%value(dg$) 
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

!-------------------------------------------------------
contains

subroutine calc_g (track, g2, g3)

type (track_struct) track
real(rp) g2, g3, g2_here, g3_here, g(3), f, s0
integer j, n1

! g2 is the average g^2 over the element for an on-energy particle.

track%pt(:)%orb%vec(6) = 0  ! on-energy

g2 = 0; g3 = 0

n1 = track%n_pt
s0 = ele%s_start

do j = 0, n1

  call g_bending_strength_from_em_field (ele, param, track%pt(j)%orb%s - s0, track%pt(j)%orb, .false., g)

  g2_here = g(1)**2 + g(2)**2 ! = g_x^2 + g_y^2
  g3_here = sqrt(g2_here)**3

  if (j == 0) then
    g2 = g2 + g2_here * track%pt(1)%s_body / 2.0_rp
    g3 = g3 + g3_here * track%pt(1)%s_body / 2.0_rp
  elseif (j == n1) then
    g2 = g2 + g2_here * (track%pt(n1)%s_body - track%pt(n1-1)%s_body) / 2.0_rp
    g3 = g3 + g3_here * (track%pt(n1)%s_body - track%pt(n1-1)%s_body) / 2.0_rp
  else
    g2 = g2 + g2_here * (track%pt(j-1)%s_body - track%pt(j+1)%s_body) / 2.0_rp
    g3 = g3 + g3_here * (track%pt(j-1)%s_body - track%pt(j+1)%s_body) / 2.0_rp
  endif


enddo

g2 = g2 / track%pt(n1)%s_body
g3 = g3 / track%pt(n1)%s_body

end subroutine calc_g

end subroutine calc_radiation_tracking_integrals

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track1_radiation (orbit, ele, param, edge)
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

subroutine track1_radiation (orbit, ele, param, edge)

use random_mod

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

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

call calc_radiation_tracking_integrals (ele, orbit, param, edge, .true., int_gx, int_gy, int_g2, int_g3)
if (int_g2 == 0 .or. orbit%state /= alive$) return

! Apply the radiation kicks
! Basic equation is E_radiated = xi * (dE/dt) * sqrt(L) / c_light
! where xi is a random number with sigma = 1.

mc2 = mass_of(param%particle)
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
! Subroutine track1_radiation_center (orbit, ele1, ele2, param, edge)
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
!
! Output:
!   orbit     -- coord_struct: Particle position after radiation has been applied.
!-

subroutine track1_radiation_center (orbit, ele1, ele2, param)

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

character(*), parameter :: r_name = 'track1_radiation_center'

!

if (.not. bmad_com%radiation_damping_on .and. .not. bmad_com%radiation_fluctuations_on) return

call calc_radiation_tracking_integrals (ele1, orbit, param, start_edge$, .false., int_gx, int_gy, int_g2, int_g3)
call calc_radiation_tracking_integrals (ele2, orbit, param, end_edge$, .false., gxi, gyi, g2i, g3i)
int_gx = int_gx + gxi
int_gy = int_gy + gyi
int_g2 = int_g2 + g2i
int_g3 = int_g3 + g3i

if (int_g2 == 0) return

! Apply the radiation kicks
! Basic equation is E_radiated = xi * (dE/dt) * sqrt(L) / c_light
! where xi is a random number with sigma = 1.

mc2 = mass_of(param%particle)
q_charge2 = charge_of(orbit%species)**2
gamma_0 = ele1%value(e_tot$) / mc2

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
  if (ele1%key == sbend$ .or. ele1%key == wiggler$ .or. ele1%key == undulator$) dE_p = &
                                dE_p + (ele1%value(dpz_rad_damp_ave$) + ele2%value(dpz_rad_damp_ave$)) / rel_p
endif
dE_p = dE_p * synch_rad_com%scale 

orbit%vec(2) = orbit%vec(2) * (1 - dE_p)
orbit%vec(4) = orbit%vec(4) * (1 - dE_p)
orbit%vec(6) = orbit%vec(6) - dE_p * rel_p

end subroutine track1_radiation_center 

end module
