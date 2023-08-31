!+
! Subroutine tao_spin_polarization_calc (branch, tao_branch, excite_zero, ignore_kinetic)
!
! Routine to calculate the spin equalibrium polarization in a ring along with the polarization rate and
! the depolarization rate due to emission of synchrotron radiation photons.
!
! From the Handbook of Accelerator Physics
!
! Input:
!   branch            -- branch_struct: Lattice branch to analyze.
!   tao_branch        -- tao_lattice_branch_struct: Contains %orbit
!   excite_zero(3)    -- character(*), optional: See documentation on spin_concat_linear_maps.
!   ignore_kinetic    -- character(*), optional: List of elements to veto the kinetic polarization term.
!
! Output:
!   tao_branch    -- tao_lattice_branch_struct: Calculated is:
!     %spin%dn_dpz(:) 
!     %spin
!-

subroutine tao_spin_polarization_calc (branch, tao_branch, excite_zero, ignore_kinetic)

use tao_data_and_eval_mod, dummy => tao_spin_polarization_calc
use radiation_mod
use ptc_interface_mod

implicit none

type (branch_struct), target :: branch
type (tao_lattice_branch_struct), target :: tao_branch
type (coord_struct), pointer :: orbit(:)
type (spin_orbit_map1_struct) :: q_1turn
type (spin_orbit_map1_struct), pointer :: q1
type (spin_orbit_map1_struct), target :: q_ele(0:branch%n_ele_track)
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: eles(:)

real(rp) n0(3), dn_dpz(3), integral_bdn_partial(3), integral_bdn_partial2(3), partial(3,3), partial2(3,3)
real(rp) integral_bn, integral_bdn, integral_1ns, integral_dn2, integral_dn2_partial(3), integral_dn2_partial2(3)
real(rp) int_g(2), g, int_g2, int_g3, b_vec(3), s_vec(3), del_p, cm_ratio, gamma, f
real(rp) old_int_g3, old_b_vec(3), old_dn_dpz(3), old_s_vec(3), old_n0(3), old_partial(3,3), old_partial2(3,3)
real(rp) g_tol, g2_tol, g3_tol
real(rp), parameter :: f_limit = 8 / (5 * sqrt(3.0_rp))
real(rp), parameter :: f_rate = 5 * sqrt(3.0_rp) * classical_radius_factor * h_bar_planck * c_light**2 / 8

integer ix1, ix2, n_loc
integer i, j, k, kk, n, p, ie

logical valid_value, err
character(*), optional :: excite_zero(3), ignore_kinetic

!

g_tol  = 1e-4_rp * branch%param%g1_integral / branch%param%total_length
g2_tol = 1e-4_rp * branch%param%g2_integral / branch%param%total_length
g3_tol = 1e-4_rp * branch%param%g3_integral / branch%param%total_length

integral_bn           = 0
integral_bdn          = 0
integral_1ns          = 0
integral_dn2          = 0
integral_bdn_partial  = 0
integral_dn2_partial  = 0
integral_bdn_partial2 = 0
integral_dn2_partial2 = 0

!

call tao_spin_tracking_turn_on()

branch%ele%logic = .true.  ! False => Ignore kinetic term
if (present(ignore_kinetic)) then
  call lat_ele_locator (ignore_kinetic, branch%lat, eles, n_loc, err, ix_dflt_branch = branch%ix_branch)
  if (err) return
  do i = 1, n_loc
    call mark_logic_false(eles(i)%ele)
  enddo
endif

orbit => tao_branch%orbit
if (.not. allocated(tao_branch%spin_ele)) allocate (tao_branch%spin_ele(0:branch%n_ele_track))
tao_branch%spin_valid = .true.

call spin_concat_linear_maps (err, q_1turn, branch, 0, branch%n_ele_track, q_ele, orbit, excite_zero)
if (err) return

if (branch%param%geometry == closed$) then
  tao_branch%spin%tune = 2.0_rp * atan2(norm2(q_1turn%spin_q(1:3,0)), q_1turn%spin_q(0,0))
  n0       = 0
  dn_dpz   = 0
  partial  = 0
  partial2 = 0

else
  n0       = branch%lat%particle_start%spin
  dn_dpz   = branch%ele(0)%value(spin_dn_dpz_x$:spin_dn_dpz_z$)
  partial  = 0 ! Not sure this is computable.
  partial2 = 0
  tao_branch%spin_ele(0)%dn_dpz%vec      = dn_dpz
  tao_branch%spin_ele(0)%dn_dpz%partial  = partial
  tao_branch%spin_ele(0)%dn_dpz%partial2 = partial2
endif

! Loop over all elements.
! Assume that dn_dpz varies linearly within an element so dn_dpz varies quadratically.

do ie = 0, branch%n_ele_track

  old_dn_dpz   = dn_dpz
  old_n0       = n0
  old_partial  = partial
  old_partial2 = partial2

  if (branch%param%geometry == closed$) then
    q_1turn = q_ele(ie) * q_1turn * map1_inverse(q_ele(ie))
    dn_dpz = spin_dn_dpz_from_qmap(q_1turn%orb_mat, q_1turn%spin_q, partial, partial2, err)
    if (err) exit
    tao_branch%spin_ele(ie)%dn_dpz%vec      = dn_dpz
    tao_branch%spin_ele(ie)%dn_dpz%partial  = partial
    tao_branch%spin_ele(ie)%dn_dpz%partial2 = partial2
    n0 = q_1turn%spin_q(1:3, 0)
    n0 = n0 / norm2(n0)
  else
    if (ie == 0) cycle
    n0       = quat_rotate(q_ele(ie)%spin_q(:,0), n0)
    dn_dpz   = quat_rotate(q_ele(ie)%spin_q(:,0), tao_branch%spin_ele(ie-1)%dn_dpz%vec)
    partial  = 0
    partial2 = 0
    tao_branch%spin_ele(ie)%dn_dpz%vec      = dn_dpz
    tao_branch%spin_ele(ie)%dn_dpz%partial  = partial
    tao_branch%spin_ele(ie)%dn_dpz%partial2 = partial2
  endif

  del_p  = 1 + orbit(ie)%vec(6)
  s_vec(1:2) = [orbit(ie)%vec(2)/del_p, orbit(ie)%vec(4)/del_p]
  s_vec(3) = sqrt(1.0_rp - s_vec(1)**2 - s_vec(2)**2)

  ele => branch%ele(ie)
  call rad_g_integrals (ele, upstream$, orbit(ie-1), orbit(ie), int_g, int_g2, old_int_g3, g_tol, g2_tol, g3_tol)
  old_b_vec = [int_g(2), -int_g(1), 0.0_rp]
  if (any(old_b_vec /= 0)) old_b_vec = old_b_vec / norm2(old_b_vec)

  call rad_g_integrals (ele, downstream$, orbit(ie-1), orbit(ie), int_g, int_g2, int_g3, g_tol, g2_tol, g3_tol)
  b_vec = [int_g(2), -int_g(1), 0.0_rp]
  if (any(b_vec /= 0)) b_vec = b_vec / norm2(b_vec)

  if (int_g2 /= 0) then
    integral_bn   = integral_bn  + old_int_g3 * dot_product(old_b_vec, old_n0) + int_g3 * dot_product(b_vec, n0)
    if (ele%logic) integral_bdn = integral_bdn + &
                                   old_int_g3 * dot_product(old_b_vec, old_dn_dpz) + int_g3 * dot_product(b_vec, dn_dpz)
    
    integral_1ns  = integral_1ns + old_int_g3 * (1 - (2.0_rp/9.0_rp) * dot_product(old_n0, old_s_vec)**2) + &
                                   int_g3 * (1 - (2.0_rp/9.0_rp) * dot_product(n0, s_vec)**2)
    integral_dn2  = integral_dn2 + (11.0_rp/18.0_rp) * integ_dn2(old_int_g3, old_dn_dpz, int_g3, dn_dpz)

    do kk = 1, 3
      if (ele%logic) integral_bdn_partial(kk) = integral_bdn_partial(kk) + &
                      old_int_g3 * dot_product(old_b_vec, old_partial(kk,:)) + int_g3 * dot_product(b_vec, partial(kk,:))
      integral_dn2_partial(kk)  = integral_dn2_partial(kk) + (11.0_rp/18.0_rp) * &
                                                      integ_dn2 (old_int_g3, old_partial(kk,:), int_g3, partial(kk,:))
      integral_bdn_partial2(kk) = integral_bdn_partial2(kk) + old_int_g3 * dot_product(old_b_vec, old_partial2(kk,:)) + &
                                                            int_g3 * dot_product(b_vec, partial2(kk,:))
      integral_dn2_partial2(kk) = integral_dn2_partial2(kk) + (11.0_rp/18.0_rp) * &
                                                      integ_dn2 (old_int_g3, old_partial2(kk,:), int_g3, partial2(kk,:))
    enddo
  endif
enddo

! Some toy lattices may not have any bends (EG: spin single resonance model lattice) or have lattice length zero.
! So test that integral_1ns is non-zero.

if (integral_1ns == 0 .or. err) then
  tao_branch%spin%pol_limit_st          = 0
  tao_branch%spin%pol_limit_dk          = 0
  tao_branch%spin%pol_limit_dk_partial  = 0
  tao_branch%spin%pol_limit_dk_partial2 = 0
  tao_branch%spin%pol_rate_bks          = 0
  tao_branch%spin%depol_rate            = 0
  tao_branch%spin%depol_rate_partial    = 0
  tao_branch%spin%depol_rate_partial2   = 0
else
  cm_ratio = charge_to_mass_of(branch%param%particle)
  call convert_pc_to ((1 + orbit(0)%vec(6)) * orbit(0)%p0c, branch%param%particle, gamma = gamma)
  f = f_rate * gamma**5 * cm_ratio**2 / branch%param%total_length
  tao_branch%spin%pol_limit_st          = abs(f_limit * integral_bn / integral_1ns)
  tao_branch%spin%pol_limit_dk          = abs(f_limit * (integral_bn - integral_bdn) / (integral_1ns + integral_dn2))
  tao_branch%spin%pol_limit_dk_partial  = abs(f_limit * (integral_bn - integral_bdn_partial) / (integral_1ns + integral_dn2_partial))
  tao_branch%spin%pol_limit_dk_partial2 = abs(f_limit * (integral_bn - integral_bdn_partial2) / (integral_1ns + integral_dn2_partial2))
  tao_branch%spin%pol_rate_bks          = f * integral_1ns
  tao_branch%spin%depol_rate            = f * integral_dn2
  tao_branch%spin%depol_rate_partial    = f * integral_dn2_partial
  tao_branch%spin%depol_rate_partial2   = f * integral_dn2_partial2
  tao_branch%spin%integral_bn           = integral_bn
  tao_branch%spin%integral_bdn          = integral_bdn
  tao_branch%spin%integral_1ns          = integral_1ns
  tao_branch%spin%integral_dn2          = integral_dn2
endif

!--------------------------------------
contains

function integ_dn2(old_g3int, old_dn, g3int, dn) result (integral)
real(rp) old_g3int, g3int, old_dn(3), dn(3), integral
!
integral = (old_g3int + g3int) * dot_product(old_dn+dn, old_dn+dn) / 4
end function integ_dn2

!--------------------------------------
! contains

recursive subroutine mark_logic_false(ele)
type (ele_struct), target :: ele
type (ele_struct), pointer :: slave
integer ix
!
ele%logic = .false.
do ix = 1, ele%n_slave
  slave => pointer_to_slave(ele, ix)
  call mark_logic_false(slave)
enddo

end subroutine mark_logic_false

end subroutine tao_spin_polarization_calc
