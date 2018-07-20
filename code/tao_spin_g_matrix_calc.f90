!+
! Subroutine tao_spin_g_matrix_calc (datum, u, ix_ref, ix_ele, spin_g, valid_value, why_invalid)
!
! Routine to calculate the spin g-matrix for a datum.
!
! Input:
!   datum         -- tao_data_struct: 
!   u             -- tao_universe_struct: Tao universe used for evaluation.
!   ix_ref        -- integer: Reference lattice element for the datum
!   ix_ele        -- integer: Lattice element to evaluate at.
!
! Output:
!   spin_g(2,6)   -- real(rp): G-matrix.
!   valid_value   -- logical: True if able to calculate the g-matrix 
!   why_invalid   -- Character(*), optional: Tells why datum value is invalid.
!-

subroutine tao_spin_g_matrix_calc (datum, u, ix_ref, ix_ele, spin_g, valid_value, why_invalid)

use tao_data_and_eval_mod, dummy => tao_spin_g_matrix_calc
use ptc_interface_mod
use pointer_lattice

implicit none

type (tao_data_struct) datum
type (tao_universe_struct), target :: u
type (tao_spin_map_struct), pointer :: sm
type (tao_spin_map_struct), allocatable ::sm_temp(:)
type (branch_struct), pointer :: branch
type (taylor_struct) orbit_taylor(6), spin_taylor(0:3)
type (c_damap) cd, cd2
type (taylor) ptc_taylor
type (real_8) r8(6)
type (q_linear) q_lin

integer ix_ref, ix_ele
integer ix_r

real(rp) spin_g(2,6)
real(rp) :: mat6(6,6), mat3(3,3), quat0(0:3), quat1(0:3), qq(0:3)
real(rp) :: v0(6) = 0, vec0(6), n0(3), n1(3), l0(3), m0(3), l1(3), m1(3)
real(rp) quat0_lnm_to_xyz(0:3), quat1_lnm_to_xyz(0:3)

integer i, j, n, p

logical valid_value

character(*) why_invalid

! Checks

valid_value = .false.
branch => u%model%lat%branch(datum%ix_branch)

if (ix_ref > -1 .and. all(datum%spin_n0 == 0)) then
  call tao_set_invalid (datum, 'SPIN_N0 NOT SET.', why_invalid, .true.)
  return
endif

! Has this been computed before? If so then use old calculation.

if (allocated(scratch%spin_map)) then
  do i = 1, size(scratch%spin_map)
    sm => scratch%spin_map(i)
    if (sm%ix_ele /= ix_ele .or. sm%ix_ref /= ix_ref .or. sm%ix_uni /= u%ix_uni) cycle
    if (ix_ref > -1 .and. all(sm%n0 /= datum%spin_n0)) cycle
    spin_g = sm%g_mat
    valid_value = .true.
    return
  enddo
endif

!----
! Calculate g-matrix...
! First concatenate the spin/orbital map

call alloc (cd)
call alloc (cd2)
call alloc (ptc_taylor)
call alloc (r8)

cd = 1

ix_r = ix_ref
if (ix_r < 0) ix_r = ix_ele

if (ix_r >= ix_ele) then
  call concat_this (ix_r, branch%n_ele_track)
  call concat_this (0, ix_ele)
else
  call concat_this (ix_r, ix_ele)
endif

q_lin = cd
quat0 = q_lin%q(:, 0)

!do i = 0, 3
!  ptc_taylor = cd%q%x(i)
!  spin_taylor(i) = ptc_taylor
!enddo

r8 = cd
orbit_taylor = r8

call taylor_to_mat6 (orbit_taylor, v0, vec0, mat6)

! If 1-turn then calculate n0. Otherwise take the value in the datum.

if (ix_r == ix_ele) then  ! 1-turn
  n0 = q_lin%q(1:3,0)   ! n0 is the rotation axis of the 0th order part of the map.
  n1 = n0
else
  n0 = datum%spin_n0 / norm2(datum%spin_n0)
  n1 = rotate_vec_given_quat(quat0, n0)
endif

! Construct coordinate systems (l0, n0, m0) and (l1, n1, m1)

j = maxloc(abs(n0), 1)
select case (j)
case (1); l0 = [-n0(3), 0.0_rp, n0(1)]
case (2); l0 = [n0(2), -n0(1), 0.0_rp]
case (3); l0 = [0.0_rp, n0(3), -n0(2)]
end select

l0 = l0 / norm2(l0)
m0 = cross_product(l0, n0)

if (ix_r == ix_ele) then
  l1 = l0
  m1 = m0
else
  l1 = rotate_vec_given_quat(quat0, l0)
  m1 = rotate_vec_given_quat(quat0, m0)
endif

mat3(:,1) = l0
mat3(:,2) = n0
mat3(:,3) = m0
quat0_lnm_to_xyz = w_mat_to_quat(mat3)

mat3(:,1) = l1
mat3(:,2) = n1
mat3(:,3) = m1
quat1_lnm_to_xyz = w_mat_to_quat(mat3)

! Calculate the g-matrix.

do p = 1, 6
  quat1 = q_lin%q(:, p)
  qq = quat_mul(quat_mul(quat_inverse(quat1_lnm_to_xyz), quat1), quat0_lnm_to_xyz)
  spin_g(1,p) = 2 * (quat0(1)*qq(2) + quat0(2)*qq(1) - quat0(0)*qq(3) - quat0(3)*qq(0))
  spin_g(2,p) = 2 * (quat0(0)*qq(1) + quat0(1)*qq(0) + quat0(2)*qq(3) + quat0(3)*qq(2))
enddo

! Cleanup

call kill (cd)
call kill (cd2)
call kill (ptc_taylor)
call kill (r8)

!---
! Save present results 

if (allocated(scratch%spin_map)) then
  n = size(scratch%spin_map)
  call move_alloc (scratch%spin_map, sm_temp)
  allocate (scratch%spin_map(n+1))
  scratch%spin_map(1:n) = sm_temp
else
  allocate (scratch%spin_map(1))
endif

sm => scratch%spin_map(size(scratch%spin_map))
sm%ix_ref = ix_ref
sm%ix_ele = ix_ele
sm%ix_uni = u%ix_uni
sm%g_mat = spin_g

valid_value = .true.

!--------------------------------------------
contains

subroutine concat_this (ix1, ix2)

integer ix1, ix2
integer i, j
!

do j = ix1+1, ix2
  call truncate_taylor_to_order (branch%ele(j)%spin_taylor, 1, spin_taylor)
  do i = 0, 3
    ptc_taylor = spin_taylor(i)
    cd2%q%x(i) = ptc_taylor
  enddo

  call truncate_taylor_to_order (branch%ele(j)%taylor, 1, orbit_taylor)
  r8 = orbit_taylor

  cd = cd2 * cd
enddo

end subroutine concat_this

end subroutine tao_spin_g_matrix_calc
