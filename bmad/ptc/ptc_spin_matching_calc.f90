!+
! Subroutine ptc_spin_matching_calc (branch, match_info)
!
! Routine to calculate the G-matrix and other stuff relavent to spin matching.
!
! Note: A call to lat_to_ptc_layout must be done before calling this routine.
!
! Input:
!   branch          -- branch_struct: Lattice branch to analyze.
!
! Output:
!   match_info(0:)  -- spin_matching_struct, allocatable: G-matrix and other stuff.
!                       The array will be allocated by this routine.
!-

subroutine ptc_spin_matching_calc (branch, match_info)

use ptc_layout_mod, dummy => ptc_spin_matching_calc
use pointer_lattice

implicit none

type (branch_struct), target :: branch
type (spin_matching_struct), allocatable, target :: match_info(:)
type (spin_matching_struct), pointer :: minfo, info0
type (ele_struct), pointer :: ele
type (internal_state) ptc_state
type (fibre), pointer :: ptc_fibre, fib_now, fib_next
type (layout), pointer :: ptc_layout
type (c_damap) cdamap, u, u_c, cdamap_1
type (probe_8) p8_a, p8_ele
type (probe) probe_orb
type (c_normal_form) cc_norm
type (real_8) r8
type (c_linear_map) q0, q2, n_axis, q_invar, q_nonlin, q_x, q_y, q_z, l_axis, m_axis, q_rot, q_lin, q_1turn
type (taylor_struct) spin_taylor, bmad_taylor(6)

type quat1_struct
  real(rp) q(0:3)
end type
type (quat1_struct) qr(6)

real(rp) spin_tune(2), phase(3), quat(0:3), mat3(3,3), vec0(6)
real(rp) quat0(0:3), quat_lnm_to_xyz(0:3), quat0_lnm_to_xyz(0:3), qq(0:3), q0_lnm(0:3)
real(rp) re_val(6), im_val(6), re_vec(6), im_vec(6)

integer ie, k, p, order, plane(3)

logical rf_on

!

rf_on = rf_is_on(branch)
if (rf_on) then
  ptc_state = ptc_private%base_state - NOCAVITY0 + SPIN0
else
  ptc_state = ptc_private%base_state + NOCAVITY0 + SPIN0
endif

!

if (allocated(match_info)) then
  if (ubound(match_info, 1) /= branch%n_ele_track) deallocate(match_info)
endif
if (.not. allocated(match_info)) allocate(match_info(0:branch%n_ele_track))

!

call init_all(ptc_state, 1, 0)   ! Only need first order map for this analysis

call alloc(cdamap, cdamap_1, u, u_c)
call alloc(p8_a)
call alloc(p8_ele)
call alloc(cc_norm)

q_x = 1
q_y = 2
q_z = 3
cdamap_1 = 1

ele => branch%ele(0)
minfo => match_info(0)

ptc_fibre => pointer_to_fibre(ele)
ptc_layout => ptc_fibre%parent_layout

!

if (branch%param%geometry == closed$) then
  call find_orbit_x (minfo%orb0, ptc_state, 1.e-8_rp, fibre1 = ptc_fibre) 

  probe_orb = minfo%orb0

  p8_a = probe_orb + cdamap_1

  call track_probe(p8_a, ptc_state, fibre1 = ptc_fibre)

  cdamap = p8_a

  call ptc_set_rf_state_for_c_normal(ptc_state%nocavity)
  call c_normal(cdamap, cc_norm, dospin = .true.)

  cc_norm%n = cc_norm%atot**(-1) * cdamap * cc_norm%atot
  p8_a = probe_orb + cc_norm%atot

  u = cc_norm%atot
  call c_fast_canonise(u, u_c, dospin = .true.)

  p8_a = probe_orb + u_c
  p8_ele = probe_orb + cdamap_1

else
  u = 1
  u%q = q_y
  u = u * u_c
  call c_fast_canonise (u, u_c, dospin = .true.)

  minfo%orb0 = 0
  probe_orb = minfo%orb0

  p8_a = probe_orb + u_c
  p8_ele = probe_orb + cdamap_1
endif

fib_now => ptc_fibre
q_invar = 1   ! Set %mat = unit matrix

!

do ie = 0, branch%n_ele_track
  if (ie /= 0) then
    fib_next => pointer_to_fibre(branch%ele(ie))
    p8_ele = probe_orb + cdamap_1

    call track_probe(p8_a, ptc_state, fibre1 = fib_now, fibre2 = fib_next)
    call track_probe(p8_ele,   ptc_state, fibre1 = fib_now, fibre2 = fib_next)

    fib_now => fib_next
  endif

  minfo => match_info(ie)
  minfo = spin_matching_struct()

  minfo%orb0 = p8_a%x

  u = p8_a
  q_invar = u

  q2 = q_invar * q_y * q_invar**(-1)
  minfo%dn_dpz = q2%q(1:3,6)

  cdamap = u * cc_norm%n * u**(-1)
  q_1turn = cdamap
  minfo%m_1turn(1:6,1:6) = q_1turn%mat

  bmad_taylor = p8_ele%x
  call taylor_to_mat6 (bmad_taylor, [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], vec0, minfo%m_ele(1:6,1:6))

  call c_fast_canonise (u, u_c, q_cs = q_lin, phase = phase, q_rot = q_rot, spin_tune = spin_tune, dospin = .true.)

  q0 = q_lin
  q0%q(0:3,1:6) = 0
  q_nonlin = q0**(-1) * q_lin

  minfo%alpha = q_nonlin%q(1,1:6)
  minfo%beta = q_nonlin%q(3,1:6)

  l_axis = q0 * q_x * q0**(-1)
  n_axis = q0 * q_y * q0**(-1)
  m_axis = q0 * q_z * q0**(-1)

  minfo%axis%l  = l_axis%q(1:3,0)
  minfo%axis%n0 = n_axis%q(1:3,0)
  minfo%axis%m  = m_axis%q(1:3,0)

  !

  if (ie /= 0) then
    mat3(:,1) = minfo%axis%l
    mat3(:,2) = minfo%axis%n0
    mat3(:,3) = minfo%axis%m
    quat_lnm_to_xyz = w_mat_to_quat(mat3)

    info0 => match_info(ie-1)
    mat3(:,1) = info0%axis%l
    mat3(:,2) = info0%axis%n0
    mat3(:,3) = info0%axis%m
    quat0_lnm_to_xyz = w_mat_to_quat(mat3)

    ! ele matrix

    do k = 0, 3
      spin_taylor = p8_ele%q%x(k)%t
      quat0(k) = taylor_coef(spin_taylor, [0,0,0,0,0,0])
      do p = 1, 6
        qr(p)%q(k) = taylor_coef(spin_taylor, taylor_expn([p]))
      enddo
    enddo

    q0_lnm = quat_mul(quat_inverse(quat_lnm_to_xyz), quat0, quat0_lnm_to_xyz)
    mat3 = quat_to_w_mat(q0_lnm)
    minfo%M_ele(7:8,7:8) = mat3(1:3:2,1:3:2)

    do p = 1, 6
      qq = quat_mul(quat_inverse(quat_lnm_to_xyz), qr(p)%q, quat0_lnm_to_xyz)
      minfo%M_ele(7,p) = 2 * (q0_lnm(1)*qq(2) + q0_lnm(2)*qq(1) - q0_lnm(0)*qq(3) - q0_lnm(3)*qq(0))
      minfo%M_ele(8,p) = 2 * (q0_lnm(0)*qq(1) + q0_lnm(1)*qq(0) + q0_lnm(2)*qq(3) + q0_lnm(3)*qq(2))
    enddo

    minfo%sq_ele = quat0

    ! 1-turn matrix

    quat0 = q_1turn%q(:,0)
    q0_lnm = quat_mul(quat_inverse(quat_lnm_to_xyz), quat0, quat_lnm_to_xyz)
    mat3 = quat_to_w_mat(q0_lnm)
    minfo%M_1turn(7:8,7:8) = mat3(1:3:2,1:3:2)

    do p = 1, 6
      qq = q_1turn%q(:,p)
      qq = quat_mul(quat_inverse(quat_lnm_to_xyz), qq, quat_lnm_to_xyz)
      minfo%M_1turn(7,p) = 2 * (q0_lnm(1)*qq(2) + q0_lnm(2)*qq(1) - q0_lnm(0)*qq(3) - q0_lnm(3)*qq(0))
      minfo%M_1turn(8,p) = 2 * (q0_lnm(0)*qq(1) + q0_lnm(1)*qq(0) + q0_lnm(2)*qq(3) + q0_lnm(3)*qq(2))
    enddo

    minfo%sq_1turn = quat0
  endif

  !

  probe_orb = p8_a
  p8_a = probe_orb + u_c
enddo

!

call kill (cdamap_1)
call kill (cdamap)
call kill (u)
call kill (u_c)
call kill (p8_a)
call kill (p8_ele)
call kill (cc_norm)

!!!call init (ptc_private%base_state, ptc_private%taylor_order_ptc, 0)
call init_all (ptc_private%base_state, ptc_private%taylor_order_ptc, 0)

end subroutine ptc_spin_matching_calc

