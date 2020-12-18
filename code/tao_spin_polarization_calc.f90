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
use pointer_lattice, dummy2 => sqrt

implicit none

type (branch_struct), target :: branch
type (tao_lattice_branch_struct), target :: tao_branch
type (coord_struct), pointer :: orbit(:)
type (c_linear_map) :: q_1turn
type (c_linear_map), pointer :: q1
type (c_linear_map), target :: q_ele(branch%n_ele_track)
type (ele_struct), pointer :: ele
type (taylor_struct), pointer :: st

real(rp) vec0(6), mat6(6,6), n0(3), dn_dpz(3)
real(rp) integral_bn, integral_1minus, integral_dn_ddel
real(rp) int_gx, int_gy, int_g, int_g2, int_g3, b_vec(3), s_vec(3), del_p, cm_ratio, gamma, f
real(rp), parameter :: f_limit = 8 / (5 * sqrt(3.0_rp))
real(rp), parameter :: f_rate = 5 * sqrt(3.0_rp) * classical_radius_factor * h_bar_planck * c_light**2 / 8

integer ix1, ix2
integer i, j, k, kk, n, p, ie

logical valid_value, st_on, err

!

orbit => tao_branch%orbit
if (.not. allocated(tao_branch%dn_dpz)) allocate (tao_branch%dn_dpz(0:branch%n_ele_track))
tao_branch%spin_valid = .true.

q_1turn = 0
do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)
  if (.not. associated(ele%spin_taylor(0)%term)) then
    st_on = bmad_com%spin_tracking_on
    bmad_com%spin_tracking_on = .true.
    call ele_to_taylor(ele, branch%param, ele%map_ref_orb_in)
    bmad_com%spin_tracking_on = st_on
  endif

  q1 => q_ele(ie)
  q1%q = 0

  do i = 0, 3
    st => ele%spin_taylor(i)
    do k = 1, size(st%term)
      n = sum(st%term(k)%expn)
      select case (n)
      case (0)
        q1%q(i,0) = st%term(k)%coef
      case (1)
        do p = 1, 6
          if (st%term(k)%expn(p) == 0) cycle
          q1%q(i,p) = st%term(k)%coef
          exit
        enddo
      end select
    enddo
  enddo

  call taylor_to_mat6 (ele%taylor, ele%taylor%ref, vec0, mat6)
  q1%mat = mat6

  q_1turn = q1 * q_1turn
enddo

! Loop over all elements.

integral_bn      = 0 ! Integral of g^3 (b_hat dot (n - dn/ddelta))
integral_1minus  = 0 ! Integral of g^3 (1 - (2/9) (n dot s)^2)
integral_dn_ddel = 0 ! Integral of g^3 (11/18) dn_dpz

do ie = 0, branch%n_ele_track
  ele => branch%ele(ie)
  if (ie /= 0) q_1turn = q_ele(ie) * q_1turn * q_ele(ie)**(-1)
  
  dn_dpz = spin_dn_dpz_from_qmap(real(q_1turn%mat, rp), real(q_1turn%q, rp))
  tao_branch%dn_dpz(ie)%vec = dn_dpz

  n0 = q_1turn%q(1:3, 0)
  n0 = n0 / norm2(n0)

  del_p  = 1 + orbit(ie)%vec(6)
  s_vec(1:2) = [orbit(ie)%vec(2)/del_p, orbit(ie)%vec(4)/del_p]
  s_vec(3) = sqrt(1.0_rp - s_vec(1)**2 - s_vec(2)**2)

  ele => branch%ele(ie)
  call calc_radiation_tracking_integrals (ele, orbit(ie), branch%param, end_edge$, .true., int_gx, int_gy, int_g2, int_g3)
  if (int_g2 /= 0) then
    if (int_gx /= 0 .or. int_gy /= 0) then
      b_vec = [int_gy, -int_gx, 0.0_rp]
      b_vec = b_vec / norm2(b_vec)
      integral_bn = integral_bn + int_g3 * dot_product(b_vec, n0 - dn_dpz)
    endif
    integral_1minus  = integral_1minus  + int_g3 * (1 - (2.0_rp/9.0_rp) * dot_product(n0, s_vec)**2)
    integral_dn_ddel = integral_dn_ddel + int_g3 * (11.0_rp/18.0_rp) * dot_product(dn_dpz, dn_dpz)
  endif

  if (ie == branch%n_ele_track) cycle

  ele => branch%ele(ie+1)
  call calc_radiation_tracking_integrals (ele, orbit(ie), branch%param, start_edge$, .true., int_gx, int_gy, int_g2, int_g3)
  if (int_g2 /= 0) then
    if (int_gx /= 0 .or. int_gy /= 0) then
      b_vec = [int_gy, -int_gx, 0.0_rp]
      b_vec = b_vec / norm2(b_vec)
      integral_bn = integral_bn + int_g3 * dot_product(b_vec, n0 - dn_dpz)
    endif
    integral_1minus  = integral_1minus  + int_g3 * (1 - (2.0_rp/9.0_rp) * dot_product(n0, s_vec)**2)
    integral_dn_ddel = integral_dn_ddel + int_g3 * (11.0_rp/18.0_rp) * dot_product(dn_dpz, dn_dpz)
  endif

enddo

!

cm_ratio = charge_to_mass_of(branch%param%particle)
call convert_pc_to ((1 + orbit(0)%vec(6)) * orbit(0)%p0c, branch%param%particle, gamma = gamma)

f = f_rate * gamma**5 * cm_ratio**2 / branch%param%total_length

tao_branch%spin%pol_limit = f_limit * integral_bn / (integral_1minus + integral_dn_ddel)
tao_branch%spin%pol_rate = f * integral_1minus
tao_branch%spin%depol_rate = f * integral_dn_ddel

end subroutine tao_spin_polarization_calc
