!+
! Subroutine tao_spin_polarization_calc (branch, orbit, spin_pol)
!
! Routine to calculate the spin equalibrium polarization in a ring along with the polarization rate and
! the depolarization rate due to emission of synchrotron radiation photons.
!
! From the Handbook of Accelerator Physics
!
! Input:
!   branch        -- branch_struct: Lattice branch to analyze.
!   orbit(0:)     -- coord_struct: Reference orbit
!
! Output:
!   spin_pol      -- tao_spin_polarization_struct.
!-

subroutine tao_spin_polarization_calc (branch, orbit, spin_pol)

use tao_data_and_eval_mod, dummy => tao_spin_polarization_calc
use radiation_mod
use ptc_interface_mod
use pointer_lattice, dummy2 => sqrt
use eigen_mod
use f95_lapack

implicit none

type (branch_struct), target :: branch
type (coord_struct) :: orbit(0:)
type (tao_spin_polarization_struct) spin_pol
type (c_linear_map) :: q_1turn
type (c_linear_map), pointer :: q1
type (c_linear_map), target :: q_ele(branch%n_ele_track)
type (ele_struct), pointer :: ele
type (taylor_struct), pointer :: st

real(rp) vec0(6), mat6(6,6), n0(3), l0(3), m0(3), mat3(3,3)
real(rp) dn_ddelta(3), m_1turn(8,8)
real(rp) quat0(0:3), quat_lnm_to_xyz(0:3), q0_lnm(0:3), qq(0:3)
real(rp) integral_bn, integral_1minus, integral_dn_ddel
real(rp) int_gx, int_gy, int_g, int_g2, int_g3, b_vec(3), s_vec(3), del_p, cm_ratio, gamma, f
real(rp), parameter :: f_limit = 8 / (5 * sqrt(3.0_rp))
real(rp), parameter :: f_rate = 5 * sqrt(3.0_rp) * classical_radius_factor * h_bar_planck * c_light**2 / 8

complex(rp) eval(6), evec(6,6), dd(2,2), gv(2,1), w_vec(6,2), vv(6,6), aa(6,1)

integer ix1, ix2
integer i, j, k, kk, n, p, ie
integer pinfo, ipiv2(2), ipiv6(6)

logical st_on, err


!

spin_pol%valid_value = .true.

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
integral_dn_ddel = 0 ! Integral of g^3 (11/18) dn_ddelta

do ie = 0, branch%n_ele_track
  ele => branch%ele(ie)

  if (ie /= 0) q_1turn = q_ele(ie) * q_1turn * q_ele(ie)**(-1)

  ! Calculate n0, and dn_ddelta
  ! Construct coordinate systems (l0, n0, m0) and (l1, n1, m1)

  n0 = q_1turn%q(1:3, 0)
  n0 = n0 / norm2(n0)

  j = maxloc(abs(n0), 1)
  select case (j)
  case (1); l0 = [-n0(3), 0.0_rp, n0(1)]
  case (2); l0 = [n0(2), -n0(1), 0.0_rp]
  case (3); l0 = [0.0_rp, n0(3), -n0(2)]
  end select

  l0 = l0 / norm2(l0)
  m0 = cross_product(l0, n0)

  mat3(:,1) = l0
  mat3(:,2) = n0
  mat3(:,3) = m0
  quat_lnm_to_xyz = w_mat_to_quat(mat3)

  quat0 = q_1turn%q(:,0)
  q0_lnm = quat_mul(quat_mul(quat_inverse(quat_lnm_to_xyz), quat0), quat_lnm_to_xyz)
  mat3 = quat_to_w_mat(q0_lnm)

  M_1turn(1:6,1:6) = q_1turn%mat
  M_1turn(1:6,7:8) = 0
  M_1turn(7:8,7:8) = mat3(1:3:2,1:3:2)

  do p = 1, 6
    qq = q_1turn%q(:,p)
    qq = quat_mul(quat_mul(quat_inverse(quat_lnm_to_xyz), qq), quat_lnm_to_xyz)
    M_1turn(7,p) = 2 * (q0_lnm(1)*qq(2) + q0_lnm(2)*qq(1) - q0_lnm(0)*qq(3) - q0_lnm(3)*qq(0))
    M_1turn(8,p) = 2 * (q0_lnm(0)*qq(1) + q0_lnm(1)*qq(0) + q0_lnm(2)*qq(3) + q0_lnm(3)*qq(2))
  enddo

  call mat_eigen(m_1turn(1:6,1:6), eval, evec, err)

  do k = 1, 6
    dd = m_1turn(7:8,7:8)
    dd(1,1) = dd(1,1) - eval(k)
    dd(2,2) = dd(2,2) - eval(k)
    gv(:,1) = -matmul(m_1turn(7:8,1:6), evec(k,1:6))
    call zgesv_f95 (dd, gv, ipiv2, pinfo)
    w_vec(k,:) = gv(:,1)
  enddo

  vv = transpose(evec)
  aa(:,1) = [0, 0, 0, 0, 0, 1]
  call zgesv_f95(vv, aa, ipiv6, pinfo)
  dn_ddelta = [real(sum(w_vec(:,1)* aa(:,1)), rp), 0.0_rp, real(sum(w_vec(:,2)* aa(:,1)), rp)]
  dn_ddelta = rotate_vec_given_quat (quat_lnm_to_xyz, dn_ddelta)

  !

  del_p  = 1 + orbit(ie)%vec(6)
  s_vec(1:2) = [orbit(ie)%vec(2)/del_p, orbit(ie)%vec(4)/del_p]
  s_vec(3) = sqrt(1.0_rp - s_vec(1)**2 - s_vec(2)**2)


  ele => branch%ele(ie)
  call calc_radiation_tracking_integrals (ele, orbit(ie), branch%param, end_edge$, .true., int_gx, int_gy, int_g2, int_g3)
  if (int_g2 /= 0) then
    if (int_gx /= 0 .or. int_gy /= 0) then
      b_vec = [int_gy, -int_gx, 0.0_rp]
      b_vec = b_vec / norm2(b_vec)
      integral_bn = integral_bn + int_g3 * dot_product(b_vec, n0 - dn_ddelta)
    endif
    integral_1minus  = integral_1minus  + int_g3 * (1 - (2.0_rp/9.0_rp) * dot_product(n0, s_vec)**2)
    integral_dn_ddel = integral_dn_ddel + int_g3 * (11.0_rp/18.0_rp) * dot_product(dn_ddelta, dn_ddelta)
  endif

  if (ie == branch%n_ele_track) cycle

  ele => branch%ele(ie+1)
  call calc_radiation_tracking_integrals (ele, orbit(ie), branch%param, start_edge$, .true., int_gx, int_gy, int_g2, int_g3)
  if (int_g2 /= 0) then
    if (int_gx /= 0 .or. int_gy /= 0) then
      b_vec = [int_gy, -int_gx, 0.0_rp]
      b_vec = b_vec / norm2(b_vec)
      integral_bn = integral_bn + int_g3 * dot_product(b_vec, n0 - dn_ddelta)
    endif
    integral_1minus  = integral_1minus  + int_g3 * (1 - (2.0_rp/9.0_rp) * dot_product(n0, s_vec)**2)
    integral_dn_ddel = integral_dn_ddel + int_g3 * (11.0_rp/18.0_rp) * dot_product(dn_ddelta, dn_ddelta)
  endif

enddo

!

cm_ratio = charge_to_mass_of(branch%param%particle)
call convert_pc_to ((1 + orbit(0)%vec(6)) * orbit(0)%p0c, branch%param%particle, gamma = gamma)

f = f_rate * gamma**5 * cm_ratio**2 / branch%param%total_length

spin_pol%pol_limit = -f_limit * integral_bn / (integral_1minus + integral_dn_ddel)
spin_pol%pol_rate = f * (integral_1minus + integral_dn_ddel)
spin_pol%depol_rate = f * integral_dn_ddel

end subroutine tao_spin_polarization_calc
