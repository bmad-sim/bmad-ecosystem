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
! Subroutine track1_radiation (orbit, ele, edge)
!
! Subroutine to apply a kick to a particle to account for radiation dampling and/or fluctuations.
!
! For tracking through a given element, this routine should be called initially when
! the particle is at the entrance end and at the end when the particle is at the exit end, when
! the orbit is with respect to laboratory (not element body) coordinates.
! That is, each time this routine is called it applies half the radiation kick for the entire element.
!
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
type (rad_int_ele_cache_struct), pointer :: ric
type (rad_map_struct) rad_mat

integer :: edge

real(rp) int_gx, int_gy, this_ran, mc2, int_g2, int_g3, ran6(6)
real(rp) gamma_0, dE_p, fact_d, fact_f, p_spin, spin_norm(3), norm, rel_p, c_radius
real(rp) damp_mat(6,6), stoc_mat(6,6), ref_orb(6)
real(rp), parameter :: rad_fluct_const = 55.0_rp * h_bar_planck * c_light / (24.0_rp * sqrt_3)
real(rp), parameter :: spin_const = 5.0_rp * sqrt_3 * h_bar_planck * c_light / 16.0_rp
real(rp), parameter :: damp_const = 2.0_rp / 3.0_rp
real(rp), parameter :: c1_spin = 2.0_rp / 9.0_rp, c2_spin = 8.0_rp / (5.0_rp * sqrt_3)

character(*), parameter :: r_name = 'track1_radiation'

!

if (.not. bmad_com%radiation_damping_on .and. .not. bmad_com%radiation_fluctuations_on) return
if (ele%value(l$) == 0) return
if (ele%tracking_method == taylor$ .and. ele%mat6_calc_method == taylor$) return
select case (ele%key)
case (drift$, taylor$, multipole$, ab_multipole$, mask$, marker$);  return
end select

! Use stochastic and damp mats

call this_rad_mat_setup(ele, edge)
ric => ele%rad_int_cache
if (edge == start_edge$) then
  rad_mat = ric%rm0
else
  rad_mat = ric%rm1
endif

if (bmad_com%radiation_damping_on) then
  orbit%vec = orbit%vec + synch_rad_com%scale * matmul(rad_mat%damp_mat, orbit%vec - rad_mat%ref_orb)
  if (.not. bmad_com%radiation_zero_average) orbit%vec = orbit%vec + synch_rad_com%scale * rad_mat%damp_vec
endif

if (bmad_com%radiation_fluctuations_on) then
  call ran_gauss (ran6)
  orbit%vec = orbit%vec + synch_rad_com%scale * matmul(rad_mat%stoc_mat, ran6)
endif

! Sokolov-Ternov Spin flip
! The equation is not correct

!if (bmad_com%spin_tracking_on .and. bmad_com%spin_sokolov_ternov_flipping_on) then
!  norm = norm2(orbit%spin)
!  if (norm /= 0) then
!    spin_norm = orbit%spin / norm
!    call ran_uniform (this_ran)
!    c_radius = classical_radius(orbit%species)
!    mc2 = mass_of(ele%ref_species)
!    gamma_0 = ele%value(e_tot$) / mc2
!    p_spin = (spin_const * c_radius * gamma_0**5 / (orbit%beta * mc2)) * &
!          (int_g3 - c1_spin * int_g3 * (spin_norm(3))**2 + c2_spin * int_g2 * dot_product([g_y, -g_x], spin_norm(1:2))) 
!    if (this_ran < p_spin) orbit%spin = -orbit%spin  ! spin flip
!  endif
!endif

!---------------------------------------------------------------------------
contains

subroutine this_rad_mat_setup (ele, edge)

use rad_6d_mod

type (coord_struct) orb1, orb2
type (ele_struct), target :: ele, runt
type (rad_int_ele_cache_struct), pointer :: ric
type (rad_map_struct) rad_mat
type (branch_struct), pointer :: branch
real(rp) tol, m_inv(6,6)
integer i, edge, info
logical err, rad_damp_on

!

if (.not. associated(ele%rad_int_cache)) allocate(ele%rad_int_cache)
ric => ele%rad_int_cache
if (.not. ric%stale) return
ric%stale = .false.

branch => pointer_to_branch(ele)

if (all(ric%rm0%ref_orb(2:4:2) == ric%rm1%ref_orb(2:4:2)) .and. &
        ele%key /= sbend$ .and. ele%key /= wiggler$ .and. ele%key /= undulator$) then
  ric%rm0 = rad_map_struct()
  ric%rm1 = rad_map_struct()
  ric%rm0%ref_orb = ele%map_ref_orb_in%vec
  ric%rm1%ref_orb = ele%map_ref_orb_out%vec
  return
endif

rad_damp_on = bmad_com%radiation_damping_on
bmad_com%radiation_damping_on = .false.

! Mats for first half of element

call create_element_slice (runt, ele, 0.5_rp*ele%value(l$), 0.0_rp, ele%branch%param, .true., .true., err, pointer_to_next_ele(ele, -1))
call make_mat6 (runt, branch%param, ele%map_ref_orb_in, orb1)
call tracking_rad_mat_setup (runt, 1e-4_rp, upstream_end$, ric%rm0)

! Mats for second half of element

call create_element_slice (runt, ele, 0.5_rp*ele%value(l$), 0.5_rp*ele%value(l$), ele%branch%param, .true., .true., err, runt)
call make_mat6 (runt, branch%param, orb1, orb2)
runt%map_ref_orb_out = ele%map_ref_orb_out  ! Important for if test above
call tracking_rad_mat_setup (runt, 1e-4_rp, downstream_end$, ric%rm1)

bmad_com%radiation_damping_on = rad_damp_on

end subroutine this_rad_mat_setup

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

end subroutine track1_radiation_center 

end module
