!+
! Subroutine tao_spin_polarization_calc (branch, valid_value, why_invalid, pol_limit, pol_rate, depol_rate)
!
! Routine to calculate the spin equalibrium polarization in a ring along with the polarization rate and
! the depolarization rate due to emission of synchrotron radiation photons.
!
! From the Handbook of Accelerator Physics
!
! Input:
!   branch        -- branch_struct: Lattice branch to analyze.
!
! Output:
!   valid_value   -- Logical: Set false when there is a problem. Set true otherwise.
!   why_invalid   -- character(*), optional: Tells why datum value is invalid.
!   pol_limit     -- real(rp): Equalibrium Polarization calculated via the Derbenev–Kondratenko–Mane formula.
!
!-

subroutine tao_spin_polarization_calc (branch, orbit, valid_value, why_invalid, pol_limit, pol_rate, depol_rate)

use tao_data_and_eval_mod, dummy => tao_spin_polarization_calc
use radiation_mod
use ptc_interface_mod
use pointer_lattice, dummy2 => sqrt
use eigen_mod
use f95_lapack

implicit none

type (branch_struct), target :: branch
type (coord_struct) :: orbit(0:)
type (q_linear) :: q_1turn
type (q_linear), pointer :: q1
type (q_linear), target :: q_ele(branch%n_ele_track)
type (ele_struct), pointer :: ele
type (taylor_struct), pointer :: st

real(rp), optional :: pol_limit, pol_rate, depol_rate
real(rp) vec0(6), mat6(6,6), n0(3), l0(3), m0(3), mat3(3,3)
real(rp) dn_ddelta(3), m_1turn(8,8)
real(rp) quat0(0:3), quat_lnm_to_xyz(0:3), q0_lnm(0:3), qq(0:3)
real(rp) integral_bn, integral_1minus, integral_dn_ddel
real(rp) len2, g_x, g_y, g, g2, g3, b_vec(3), s_vec(3), del_p, mc2, q, gamma, f
real(rp), parameter :: f_limit = 8 / (5 * sqrt(3.0_rp))
real(rp), parameter :: f_rate = 5 * sqrt(3.0_rp) * classical_radius_factor * h_bar_planck * c_light**2 / 8


complex(rp) eval(6), evec(6,6), dd(2,2), gv(2,1), w_vec(6,2)

integer ix1, ix2
integer i, j, k, kk, n, p, ie
integer pinfo, ipiv(2), ipiv6(6)

logical valid_value
logical st_on, err

character(*) why_invalid

!

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

  m_1turn = 0
  M_1turn(7:8,7:8) = mat3(1:3:2,1:3:2)

  do p = 1, 6
    qq = q_1turn%q(:,p)
    qq = quat_mul(quat_mul(quat_inverse(quat_lnm_to_xyz), qq), quat_lnm_to_xyz)
    M_1turn(7,p) = 2 * (q0_lnm(1)*qq(2) + q0_lnm(2)*qq(1) - q0_lnm(0)*qq(3) - q0_lnm(3)*qq(0))
    M_1turn(8,p) = 2 * (q0_lnm(0)*qq(1) + q0_lnm(1)*qq(0) + q0_lnm(2)*qq(3) + q0_lnm(3)*qq(2))
  enddo

  call mat_eigen(m_1turn(1:6,1:6), eval, evec, err)

  do kk = 1, 3
    k = 2*kk - 1
    f = 2 * aimag(evec(k,1) * conjg(evec(k,2)) + evec(k,3) * conjg(evec(k,4)) + evec(k,5) * conjg(evec(k,6)))
    if (f < 0) then
      evec(k:k+1,:) = evec(k+1:k:-1, :)
      eval(k:k+1) = [eval(k+1), eval(k)]
    endif

    evec(k,:)   = evec(k,:) / sqrt(abs(f))
    evec(k+1,:) = evec(k+1,:) / sqrt(abs(f))
  enddo

  do k = 1, 6
    dd = m_1turn(7:8,7:8)
    dd(1,1) = dd(1,1) - eval(k)
    dd(2,2) = dd(2,2) - eval(k)
    gv(:,1) = -matmul(m_1turn(7:8,1:6), evec(k,1:6))
    call zgesv_f95 (dd, gv, ipiv, pinfo)
    w_vec(k,:) = gv(:,1)
  enddo

  dn_ddelta = 0
  do kk = 1, 3
    k = 2*kk - 1
    dn_ddelta(1:3:2) = dn_ddelta(1:3:2) - 2 * aimag(conjg(evec(k,5)) * w_vec(k,:))
  enddo

  dn_ddelta = rotate_vec_given_quat (quat_lnm_to_xyz, dn_ddelta)

  !

  del_p  = 1 + orbit(ie)%vec(6)
  s_vec(1:2) = [orbit(ie)%vec(2)/del_p, orbit(ie)%vec(4)/del_p]
  s_vec(3) = sqrt(1.0_rp - s_vec(1)**2 - s_vec(2)**2)


  ele => branch%ele(ie)
  call calc_radiation_tracking_g_factors (ele, orbit(ie), branch%param, end_edge$, len2, g_x, g_y, g2, g3)
  if (len2 /= 0) then
    b_vec = [g_y, -g_x, 0.0_rp]
    b_vec = b_vec / norm2(b_vec)

    integral_bn      = integral_bn      + len2 * g3 * dot_product(b_vec, n0 - dn_ddelta)
    integral_1minus  = integral_1minus  + len2 * g3 * (1 - (2.0_rp/9.0_rp) * dot_product(n0, s_vec)**2)
    integral_dn_ddel = integral_dn_ddel + len2 * g3 * (11.0_rp/18.0_rp) * dot_product(dn_ddelta, dn_ddelta)
  endif

  if (ie == branch%n_ele_track) cycle

  ele => branch%ele(ie+1)
  call calc_radiation_tracking_g_factors (ele, orbit(ie), branch%param, start_edge$, len2, g_x, g_y, g2, g3)
  if (len2 /= 0) then
    b_vec = [g_y, -g_x, 0.0_rp]
    b_vec = b_vec / norm2(b_vec)

    integral_bn      = integral_bn      + len2 * g3 * dot_product(b_vec, n0 - dn_ddelta)
    integral_1minus  = integral_1minus  + len2 * g3 * (1 - (2.0_rp/9.0_rp) * dot_product(n0, s_vec)**2)
    integral_dn_ddel = integral_dn_ddel + len2 * g3 * (11.0_rp/18.0_rp) * dot_product(dn_ddelta, dn_ddelta)
  endif

enddo

!

mc2 = mass_of(branch%param%particle)
q = charge_of(branch%param%particle)
call convert_pc_to ((1 + orbit(0)%vec(6)) * orbit(0)%p0c, branch%param%particle, gamma = gamma)

f = f_rate * gamma**5 * q**2 / (mc2**2 * branch%param%total_length)

pol_limit = -f_limit * integral_bn / (integral_1minus + integral_dn_ddel)
pol_rate = f * (integral_1minus + integral_dn_ddel)
depol_rate = f * integral_dn_ddel

end subroutine tao_spin_polarization_calc
