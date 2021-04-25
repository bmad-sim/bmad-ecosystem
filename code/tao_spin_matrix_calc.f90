!+
! Subroutine tao_spin_matrix_calc (datum, u, ix_ref, ix_ele, spin_map, valid_value, why_invalid)
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
!   spin_map      -- tao_spin_map_struct, pointer: G-matrix, etc.
!   valid_value   -- logical: True if able to calculate the g-matrix 
!   why_invalid   -- Character(*), optional: Tells why datum value is invalid.
!-

subroutine tao_spin_matrix_calc (datum, u, ix_ref, ix_ele, spin_map, valid_value, why_invalid)

use tao_data_and_eval_mod, dummy => tao_spin_matrix_calc
use ptc_interface_mod
use pointer_lattice

implicit none

type (tao_data_struct) datum
type (tao_universe_struct), target :: u
type (tao_spin_map_struct), pointer :: sm
type (tao_spin_map_struct), allocatable ::sm_temp(:)
type (tao_spin_map_struct), pointer :: spin_map
type (branch_struct), pointer :: branch
type (taylor_struct) orbit_taylor(6), spin_taylor(0:3)
type (c_linear_map) q_map

integer ix_ref, ix_ele
integer ix_r

real(rp) spin_g(2,6)
real(rp) :: mat3(3,3), quat0(0:3), quat1(0:3), qq(0:3), q0_lnm(0:3)
real(rp) :: n0(3), n1(3), l0(3), m0(3), l1(3), m1(3)
real(rp) quat0_lnm_to_xyz(0:3), quat1_lnm_to_xyz(0:3)

integer i, j, n, p

logical valid_value

character(*) why_invalid

! Checks

valid_value = .false.
branch => u%model%lat%branch(datum%ix_branch)

! Has this been computed before? If so then use old calculation.

if (allocated(scratch%spin_map)) then
  do i = 1, size(scratch%spin_map)
    sm => scratch%spin_map(i)
    if (sm%ix_ele /= ix_ele .or. sm%ix_ref /= ix_ref .or. sm%ix_uni /= u%ix_uni) cycle
    if (any(sm%axis_dat%n0 /= datum%spin_axis%n0) .or. any(sm%axis_dat%l /= datum%spin_axis%l) .or. &
                                                       any(sm%axis_dat%m /= datum%spin_axis%m)) cycle
    spin_map => sm
    valid_value = .true.
    return
  enddo
endif

!----
! Calculate g-matrix...

! Allocate space

if (allocated(scratch%spin_map)) then
  n = size(scratch%spin_map)
  call move_alloc (scratch%spin_map, sm_temp)
  allocate (scratch%spin_map(n+1))
  scratch%spin_map(1:n) = sm_temp
else
  allocate (scratch%spin_map(1))
endif

spin_map => scratch%spin_map(size(scratch%spin_map))
spin_map%mat8 = 0

! Concatenate the spin/orbital map

q_map = 0

ix_r = ix_ref
if (ix_r < 0) ix_r = ix_ele

if (ix_r >= ix_ele) then
  call concat_this (q_map, ix_r, branch%n_ele_track)
  call concat_this (q_map, 0, ix_ele)
else
  call concat_this (q_map, ix_r, ix_ele)
endif

spin_map%q_map%mat = q_map%mat
spin_map%q_map%q   = q_map%q

quat0 = q_map%q(:, 0)

! If 1-turn then calculate n0. Otherwise take the value in the datum.

if (any(datum%spin_axis%n0 /= 0)) then
  n0 = datum%spin_axis%n0 / norm2(datum%spin_axis%n0)
else  
  n0 = u%model%tao_branch(datum%ix_branch)%orbit(ix_r)%spin
endif

if (all(n0 == 0)) then
  call tao_set_invalid (datum, 'DATUM SPIN_AXIS%N0 NOT SET AND CANNOT BE COMPUTED.', why_invalid, .true.)
  return
endif


n1 = rotate_vec_given_quat(quat0, n0)

! Construct coordinate systems (l0, n0, m0) and (l1, n1, m1)

if (any(datum%spin_axis%l /= 0)) then
  l0 = datum%spin_axis%l - n0 * dot_product(datum%spin_axis%l, n0)  ! Make sure l0 is perpendicular to n0.
  l0 = l0 / norm2(l0)
  m0 = cross_product(l0, n0)
elseif (any(datum%spin_axis%m /= 0)) then
  m0 = datum%spin_axis%m - n0 * dot_product(datum%spin_axis%m, n0)  ! Make sure m0 is perpendicular to n0.
  m0 = m0 / norm2(m0)
  l0 = cross_product(n0, m0)
else
  j = maxloc(abs(n0), 1)
  select case (j)
  case (1); l0 = [-n0(3), 0.0_rp, n0(1)]
  case (2); l0 = [n0(2), -n0(1), 0.0_rp]
  case (3); l0 = [0.0_rp, n0(3), -n0(2)]
  end select
  l0 = l0 / norm2(l0)
  m0 = cross_product(l0, n0)
endif

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

! Calculate the 8x8 transfer matrix 

spin_map%mat8(1:6,1:6) = q_map%mat

q0_lnm = quat_mul(quat_mul(quat_inverse(quat1_lnm_to_xyz), quat0), quat0_lnm_to_xyz)
mat3 = quat_to_w_mat(q0_lnm)
spin_map%mat8(7:8,7:8) = mat3(1:3:2,1:3:2)

do p = 1, 6
  quat1 = q_map%q(:, p)
  qq = quat_mul(quat_mul(quat_inverse(quat1_lnm_to_xyz), quat1), quat0_lnm_to_xyz)
  spin_map%mat8(7,p) = 2 * (q0_lnm(1)*qq(2) + q0_lnm(2)*qq(1) - q0_lnm(0)*qq(3) - q0_lnm(3)*qq(0))
  spin_map%mat8(8,p) = 2 * (q0_lnm(0)*qq(1) + q0_lnm(1)*qq(0) + q0_lnm(2)*qq(3) + q0_lnm(3)*qq(2))
enddo

!

spin_map%ix_ref   = ix_ref
spin_map%ix_ele   = ix_ele
spin_map%ix_uni   = u%ix_uni
spin_map%axis_dat = datum%spin_axis
spin_map%axis0    = spin_axis_struct(l0, n0, m0)
spin_map%axis1    = spin_axis_struct(l1, n1, m1)

valid_value = .true.

!--------------------------------------------
contains

subroutine concat_this (q_map, ix1, ix2)

type (c_linear_map) q_map, q_ele
type (ele_struct), pointer :: ele
type (taylor_struct), pointer :: st
real(rp) vec0(6), mat6(6,6)
integer ix1, ix2
integer i, j, k, n, p
logical st_on

!

do j = ix1+1, ix2
  ele => branch%ele(j)
  if (.not. associated(ele%spin_taylor(0)%term)) then
    st_on = bmad_com%spin_tracking_on
    bmad_com%spin_tracking_on = .true.
    call ele_to_taylor(ele, branch%param, ele%map_ref_orb_in)
    bmad_com%spin_tracking_on = st_on
  endif

  q_ele%q = 0

  do i = 0, 3
    st => ele%spin_taylor(i)
    do k = 1, size(st%term)
      n = sum(st%term(k)%expn)
      select case (n)
      case (0)
        q_ele%q(i,0) = st%term(k)%coef
      case (1)
        do p = 1, 6
          if (st%term(k)%expn(p) == 0) cycle
          q_ele%q(i,p) = st%term(k)%coef
          exit
        enddo
      end select
    enddo
  enddo

  call taylor_to_mat6 (ele%taylor, ele%taylor%ref, vec0, mat6)
  q_ele%mat = mat6

  q_map = q_ele * q_map
enddo

end subroutine concat_this

end subroutine tao_spin_matrix_calc
