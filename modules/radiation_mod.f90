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
! Modules needed:
!   use radiation_mod
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
! Subroutine calc_radiation_tracking_g_factors (ele, orbit, param, edge, len2, g_x, g_y, g, g2, g3)
!
! Routine to calculate the g bending strength parameters for half the element. 
! g = 1/rho where rho is the radius of curvature. g points raially outward in the plane of bending.
! If the particle is at the starting edge then the calculation is over the first half of the element.
! If the particle is at the exit edge then the calculation is over the 2nd half of the element.
!
! Input:
!   orbit     -- coord_struct: Particle position.
!   ele       -- ele_struct: Element that causes radiation.
!   edge      -- integer: Where the particle is: start_edge$ or end_edge$.
!     
!
! Output:
!   len2      -- real(rp): Half length traveled by particle. Set to zero if no bending present.
!   g_x       -- real(rp): Average of x-component of g.
!   g_y       -- real(rp): Average of y-component of g.
!   g2        -- real(rp): Average of g^2.
!   g3        -- real(rp): Average of g^3.
!-

subroutine calc_radiation_tracking_g_factors (ele, orbit, param, edge, len2, g_x, g_y, g2, g3)

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (coord_struct) :: orbit2
type (coord_struct) start0_orb, start_orb, end_orb
type (track_struct), save, target :: track_save
type (track_struct), pointer :: track

real(rp) len2, g_x, g_y, g2, g3, kx, ky, kx_tot, ky_tot, s_here
real(rp) a_pole_mag(0:n_pole_maxx), b_pole_mag(0:n_pole_maxx), a_pole_elec(0:n_pole_maxx), b_pole_elec(0:n_pole_maxx)
real(rp), parameter :: del_orb = 1d-4

integer edge, direc
integer i, j, ix_mag_max, ix_elec_max

character(*), parameter :: r_name = 'calc_radiation_tracking_g_factors'

!

select case (ele%key)
case (quadrupole$, sextupole$, octupole$, sbend$, sol_quad$, wiggler$, undulator$, em_field$)
case default
  len2 = 0
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

len2 = ele%value(l$) / 2
if (len2 < 0) len2 = 0

! Get the coords in the frame of reference of the element

if (ele%key /= wiggler$ .and. ele%key /= undulator$) then
  orbit2 = orbit
  call offset_particle (ele, param, set$, orbit2, s_pos = s_here)
  call canonical_to_angle_coords (orbit2)
  orbit2%vec(1) = orbit2%vec(1) + direc * orbit2%vec(2) * len2 / 2.0_rp ! Extrapolate to center of region 1/4 of way into element.
  orbit2%vec(3) = orbit2%vec(3) + direc * orbit2%vec(4) * len2 / 2.0_rp

  call multipole_ele_to_ab (ele, .false., ix_mag_max, a_pole_mag, b_pole_mag, magnetic$, .true.)
  call multipole_ele_to_ab (ele, .false., ix_elec_max, a_pole_elec, b_pole_elec, electric$, .true.)

  kx_tot = 0
  ky_tot = 0

  do i = 0, ix_mag_max
    call ab_multipole_kick (a_pole_mag(i), b_pole_mag(i), i, orbit%species, ele%orientation, orbit, kx, ky, pole_type = magnetic$)
    kx_tot = kx_tot + kx
    ky_tot = ky_tot + ky
  enddo

  do i = 0, ix_mag_max
    call ab_multipole_kick (a_pole_elec(i), b_pole_elec(i), i, orbit%species, ele%orientation, orbit, kx, ky, pole_type = electric$)
    kx_tot = kx_tot + kx
    ky_tot = ky_tot + ky
  enddo
endif

! Calculate the radius of curvature for an on-energy particle

select case (ele%key)

case (quadrupole$, sol_quad$)
  g_x = -ele%value(k1$) * orbit2%vec(1) - kx_tot
  g_y =  ele%value(k1$) * orbit2%vec(3) - ky_tot
  g2 = g_x**2 + g_y**2
  g3 = sqrt(g2)**3

case (sextupole$, octupole$)
  g_x = -kx_tot
  g_y = -ky_tot
  g2 = g_x**2 + g_y**2
  g3 = sqrt(g2)**3

case (sbend$)
  g_x = ele%value(g$) + ele%value(g_err$) + ele%value(k1$) * orbit2%vec(1) - kx_tot
  g_y = ele%value(k1$) * orbit2%vec(3) - ky_tot
  g2 = g_x**2 + g_y**2
  g3 = sqrt(g2)**3
  len2 = len2 * (1.0_rp + ele%value(g$) * orbit2%vec(1))

case (wiggler$, undulator$, em_field$)
  g_x = 0
  g_y = 0
  if (ele%sub_key == map_type$) then
    if (.not. associated(ele%rad_int_cache) .or. ele%rad_int_cache%stale) then


      if (.not. associated(ele%rad_int_cache)) allocate (ele%rad_int_cache)
      ele%rad_int_cache%orb0 = ele%map_ref_orb_in%vec

      if (global_com%be_thread_safe) then
        allocate(track)
      else
        track => track_save
      endif

      track%n_pt = -1
      call symp_lie_bmad (ele, param, ele%map_ref_orb_in, end_orb, .false., track)
      call calc_g (track, ele%rad_int_cache%g2_0, ele%rad_int_cache%g3_0)

      do j = 1, 4
        start_orb = ele%map_ref_orb_in
        start_orb%vec(j) = start_orb%vec(j) + del_orb
        track%n_pt = -1
        call symp_lie_bmad (ele, param, start_orb, end_orb, .false., track)
        call calc_g (track, g2, g3)
        ele%rad_int_cache%dg2_dorb(j) = (g2 - ele%rad_int_cache%g2_0) / del_orb
        ele%rad_int_cache%dg3_dorb(j) = (g3 - ele%rad_int_cache%g3_0) / del_orb
      enddo

      ele%rad_int_cache%stale = .false.

      if (global_com%be_thread_safe) then
        deallocate(track%orb, track%field, track%map)
        deallocate(track)
      endif
    endif



    g2 = ele%rad_int_cache%g2_0 + dot_product(orbit%vec(1:4)-ele%rad_int_cache%orb0(1:4), ele%rad_int_cache%dg2_dorb(1:4))
    g3 = ele%rad_int_cache%g3_0 + dot_product(orbit%vec(1:4)-ele%rad_int_cache%orb0(1:4), ele%rad_int_cache%dg3_dorb(1:4))
    if (g3 < 0) g3 = 0
  elseif (ele%sub_key == periodic_type$) then
    g2 = abs(ele%value(k1$))
    g3 = 4 * sqrt(2*g2)**3 / (3 * pi)  
  endif

end select

!-------------------------------------------------------
contains

subroutine calc_g (track, g2, g3)

type (track_struct) track
real(rp) g2, g3, g2_here, g3_here, g(3), f, s0
integer j, n0, n1

! g2 is the average g^2 over the element for an on-energy particle.
! Note: em_field_g_bend assumes orb is lab (not local) coords.

track%orb(:)%vec(6) = 0  ! on-energy

g2 = 0; g3 = 0

n0 = lbound(track%orb, 1)
n1 = track%n_pt
s0 = ele%s_start

do j = n0, n1

  call em_field_g_bend (ele, param, track%orb(j)%s - s0, track%orb(j), g)

  g2_here = g(1)**2 + g(2)**2 ! = g_x^2 + g_y^2
  g3_here = sqrt(g2_here)**3

  if (j == n0 .or. j == n1) then
    g2_here = g2_here / 2
    g3_here = g3_here / 2
  endif

  g2 = g2 + g2_here
  g3 = g3 + g3_here

enddo

g2 = g2 / (n1 - n0 + 1)
g3 = g3 / (n1 - n0 + 1)

end subroutine calc_g

end subroutine calc_radiation_tracking_g_factors

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track1_radiation (orbit, ele, param, edge)
!
! Subroutine to apply a kick to a particle to account for radiation dampling and/or fluctuations.
! "Baier-Katkov" spin flips are included.
!
! For tracking through a given element, this routine should be called initially when
! the particle is at the entrance end and at the end when the particle is at the exit end.
! That is, each time this routine is called it applies half the radiation kick for the entire element.
!
! Note: If both bmad_com%radiation_damping_on and not bmad_com%radiation_fluctuations_on are
! False then no spin flipping is done.  
!
! Note: This routine is called by track1.
!
! Input:
!   orbit     -- coord_struct: Particle position before radiation applied.
!   ele       -- ele_struct: Element that causes radiation.
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
type (ele_struct), pointer :: ele0
type (lat_param_struct) :: param

integer :: edge

real(rp) s_len, g_x, g_y, this_ran, mc2, g2, g3
real(rp) gamma_0, dE_p, fact_d, fact_f, q_charge2, p_spin, spin_norm(3), norm
real(rp), parameter :: rad_fluct_const = 55.0_rp * classical_radius_factor * h_bar_planck * c_light / (24.0_rp * sqrt_3)
real(rp), parameter :: spin_const = 5.0_rp * sqrt_3 * classical_radius_factor * h_bar_planck * c_light / 16
real(rp), parameter :: c1_spin = 2.0_rp / 9.0_rp, c2_spin = 8.0_rp / (5.0_rp * sqrt_3)

character(*), parameter :: r_name = 'track1_radiation'

!

if (.not. bmad_com%radiation_damping_on .and. .not. bmad_com%radiation_fluctuations_on) return

call calc_radiation_tracking_g_factors (ele, orbit, param, edge, s_len, g_x, g_y, g2, g3)
if (s_len == 0) return

! Apply the radiation kicks
! Basic equation is E_radiated = xi * (dE/dt) * sqrt(L) / c_light
! where xi is a random number with sigma = 1.

mc2 = mass_of(param%particle)
q_charge2 = charge_of(orbit%species)**2
gamma_0 = ele%value(e_tot$) / mc2

fact_d = 0
if (bmad_com%radiation_damping_on) then
  fact_d = 2 * classical_radius_factor * q_charge2 * gamma_0**3 * g2 * s_len / (3 * mc2)
  if (param%backwards_time_tracking) fact_d = -fact_d
endif

fact_f = 0
if (bmad_com%radiation_fluctuations_on) then
  call ran_gauss (this_ran)
  fact_f = sqrt(rad_fluct_const * q_charge2 * s_len * gamma_0**5 * g3) * this_ran / mc2
endif

dE_p = (1 + orbit%vec(6)) * (fact_d + fact_f) * synch_rad_com%scale 

orbit%vec(2) = orbit%vec(2) * (1 - dE_p)
orbit%vec(4) = orbit%vec(4) * (1 - dE_p)
orbit%vec(6) = orbit%vec(6)  - dE_p * (1 + orbit%vec(6))

! This is useful for debugging.

if (synch_rad_com%i_calc_on) then
  synch_rad_com%i2 = synch_rad_com%i2 + g2 * s_len
  synch_rad_com%i3 = synch_rad_com%i3 + g3 * s_len
  if (associated(ele%branch)) then
    if (edge == start_edge$) then
      ele0 => ele%branch%ele(ele%ix_ele-1)
    else
      ele0 => ele
    endif
    synch_rad_com%i5a = synch_rad_com%i5a + g3 * s_len * (ele0%a%gamma * ele0%a%eta**2 + &
                  2 * ele0%a%alpha * ele0%a%eta * ele0%a%etap + ele0%a%beta * ele0%a%etap**2)
    synch_rad_com%i5b = synch_rad_com%i5b + g3 * s_len * (ele0%b%gamma * ele0%b%eta**2 + &
                  2 * ele0%b%alpha * ele0%b%eta * ele0%b%etap + ele0%b%beta * ele0%b%etap**2)
  endif
endif

! Sokolov-Ternov Spin flip?

if (bmad_com%spin_tracking_on .and. bmad_com%spin_sokolov_ternov_flipping_on) then
  norm = norm2(orbit%spin)
  if (norm /= 0) then
    spin_norm = orbit%spin / norm
    call ran_uniform (this_ran)
    p_spin = (spin_const * q_charge2 * gamma_0**5 * s_len / (orbit%beta * mc2**2)) * &
          (g3 - c1_spin * g3 * (spin_norm(3))**2 + c2_spin * g2 * dot_product([g_y, -g_x], spin_norm(1:2))) 
    if (this_ran < p_spin) orbit%spin = -orbit%spin  ! spin flip
  endif
endif

end subroutine track1_radiation 

end module
