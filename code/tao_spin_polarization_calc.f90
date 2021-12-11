!+
! Subroutine tao_spin_polarization_calc (branch, tao_branch)
!
! Routine to calculate the spin equalibrium polarization in a ring along with the polarization rate and
! the depolarization rate due to emission of synchrotron radiation photons.
!
! From the Handbook of Accelerator Physics
!
! Input:
!   branch        -- branch_struct: Lattice branch to analyze.
!   tao_branch    -- tao_lattice_branch_struct: Contains %orbit
!
! Output:
!   tao_branch    -- tao_lattice_branch_struct: Calculated is:
!     %dn_dpz(:) 
!     %spin
!-

subroutine tao_spin_polarization_calc (branch, tao_branch)

use tao_data_and_eval_mod, dummy => tao_spin_polarization_calc
use radiation_mod
use ptc_interface_mod

implicit none

type (branch_struct), target :: branch
type (tao_lattice_branch_struct), target :: tao_branch
type (coord_struct), pointer :: orbit(:)
type (spin_orbit_map1_struct) :: q_1turn
type (spin_orbit_map1_struct), pointer :: q1
type (spin_orbit_map1_struct), target :: q_ele(branch%n_ele_track)
type (ele_struct), pointer :: ele

real(rp) n0(3), dn_dpz(3), integral_bdn_partial(3), partial(3,3)
real(rp) integral_bn, integral_bdn, integral_1ns, integral_dn2, integral_dn2_partial(3)
real(rp) int_gx, int_gy, int_g, int_g2, int_g3, b_vec(3), s_vec(3), del_p, cm_ratio, gamma, f
real(rp) old_int_g3, old_b_vec(3), old_dn_dpz(3), old_s_vec(3), old_n0(3), old_partial(3,3)
real(rp), parameter :: f_limit = 8 / (5 * sqrt(3.0_rp))
real(rp), parameter :: f_rate = 5 * sqrt(3.0_rp) * classical_radius_factor * h_bar_planck * c_light**2 / 8

integer ix1, ix2
integer i, j, k, kk, n, p, ie

logical valid_value, err

!

call tao_spin_tracking_turn_on()

orbit => tao_branch%orbit
if (.not. allocated(tao_branch%dn_dpz)) allocate (tao_branch%dn_dpz(0:branch%n_ele_track))
tao_branch%spin_valid = .true.

call spin_concat_linear_maps (q_1turn, branch, 0, branch%n_ele_track, q_ele, orbit)

tao_branch%spin%tune = 2.0_rp * atan2(norm2(q_1turn%spin_q(1:3,0)), q_1turn%spin_q(0,0))

! Loop over all elements.
! Assume that dn_dpz varies linearly within an element so dn_dpz varies quadratically.

integral_bn          = 0
integral_bdn         = 0
integral_1ns         = 0
integral_dn2         = 0
integral_bdn_partial = 0
integral_dn2_partial = 0

do ie = 0, branch%n_ele_track
  if (ie /= 0) q_1turn = q_ele(ie) * q_1turn * map1_inverse(q_ele(ie))
  
  dn_dpz = spin_dn_dpz_from_qmap(q_1turn%orb_mat, q_1turn%spin_q, partial, err)
  if (err) exit

  tao_branch%dn_dpz(ie)%vec = dn_dpz
  tao_branch%dn_dpz(ie)%partial = partial
  n0 = q_1turn%spin_q(1:3, 0)
  n0 = n0 / norm2(n0)

  del_p  = 1 + orbit(ie)%vec(6)
  s_vec(1:2) = [orbit(ie)%vec(2)/del_p, orbit(ie)%vec(4)/del_p]
  s_vec(3) = sqrt(1.0_rp - s_vec(1)**2 - s_vec(2)**2)

  ele => branch%ele(ie)
  call calc_radiation_tracking_integrals (ele, orbit(ie), branch%param, start_edge$, .true., int_gx, int_gy, int_g2, old_int_g3)
  old_b_vec = [int_gy, -int_gx, 0.0_rp]
  if (any(old_b_vec /= 0)) old_b_vec = old_b_vec / norm2(old_b_vec)

  call calc_radiation_tracking_integrals (ele, orbit(ie), branch%param, end_edge$, .true., int_gx, int_gy, int_g2, int_g3)
  b_vec = [int_gy, -int_gx, 0.0_rp]
  if (any(b_vec /= 0)) b_vec = b_vec / norm2(b_vec)

  if (int_g2 /= 0) then
    integral_bn   = integral_bn  + old_int_g3 * dot_product(old_b_vec, old_n0) + int_g3 * dot_product(b_vec, n0)
    integral_bdn  = integral_bdn + old_int_g3 * dot_product(old_b_vec, old_dn_dpz) + int_g3 * dot_product(b_vec, dn_dpz)
    integral_1ns  = integral_1ns + old_int_g3 * (1 - (2.0_rp/9.0_rp) * dot_product(old_n0, old_s_vec)**2) + &
                                   int_g3 * (1 - (2.0_rp/9.0_rp) * dot_product(n0, s_vec)**2)
    integral_dn2  = integral_dn2 + (11.0_rp/18.0_rp) * integ_dn2(old_int_g3, old_dn_dpz, int_g3, dn_dpz)

    do kk = 1, 3
      integral_bdn_partial(kk) = integral_bdn_partial(kk) + old_int_g3 * dot_product(old_b_vec, old_partial(kk,:)) + &
                                                            int_g3 * dot_product(b_vec, partial(kk,:))
      integral_dn2_partial(kk) = integral_dn2_partial(kk) + (11.0_rp/18.0_rp) * &
                                                      integ_dn2 (old_int_g3, old_partial(kk,:), int_g3, partial(kk,:))
    enddo
  endif

  old_dn_dpz = dn_dpz
  old_n0 = n0
  old_partial = partial
enddo

! Some toy lattices may not have any bends (EG: spin single resonance model lattice) or have lattice length zero.
! So test that integral_1ns is non-zero.

if (integral_1ns == 0 .or. err) then
  tao_branch%spin%pol_limit_st          = 0
  tao_branch%spin%pol_limit_dkm         = 0
  tao_branch%spin%pol_limit_dkm_partial = 0
  tao_branch%spin%pol_rate_bks          = 0
  tao_branch%spin%depol_rate            = 0
  tao_branch%spin%depol_rate_partial    = 0
else
  cm_ratio = charge_to_mass_of(branch%param%particle)
  call convert_pc_to ((1 + orbit(0)%vec(6)) * orbit(0)%p0c, branch%param%particle, gamma = gamma)
  f = f_rate * gamma**5 * cm_ratio**2 / branch%param%total_length
  tao_branch%spin%pol_limit_st          = f_limit * integral_bn / integral_1ns
  tao_branch%spin%pol_limit_dkm         = f_limit * (integral_bn - integral_bdn) / (integral_1ns + integral_dn2)
  tao_branch%spin%pol_limit_dkm_partial = f_limit * (integral_bn - integral_bdn_partial) / (integral_1ns + integral_dn2_partial)
  tao_branch%spin%pol_rate_bks          = f * integral_1ns
  tao_branch%spin%depol_rate            = f * integral_dn2
  tao_branch%spin%depol_rate_partial    = f * integral_dn2_partial
endif

!--------------------------------------
contains

function integ_dn2(old_g3int, old_dn, g3int, dn) result (integral)

real(rp) old_g3int, g3int, old_dn(3), dn(3), integral


integral = (old_g3int * dot_product(old_dn, old_dn) + sqrt(old_g3int*g3int) * dot_product(old_dn, dn) + &
                                                                         g3int * dot_product(dn, dn)) * (2.0_rp / 3.0_rp)

end function integ_dn2

end subroutine tao_spin_polarization_calc
