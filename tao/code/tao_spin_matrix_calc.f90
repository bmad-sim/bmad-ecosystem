!+
! Subroutine tao_spin_matrix_calc (datum, u, ele_ref, ele_eval, excite_zero)
!
! Routine to calculate the spin g-matrix for a datum.
!
! Input:
!   datum             -- tao_data_struct: 
!   u                 -- tao_universe_struct: Tao universe used for evaluation.
!   ele_ref           -- ele_struct: Reference lattice element for the datum
!   ele_eval          -- ele_struct: Lattice element to evaluate at.
!   excite_zero(3)    -- character(*): See documentation on spin_concat_linear_maps.
!
! Output:
!   datum%spin_map    -- tao_spin_map_struct, pointer: G-matrix, etc.
!-

subroutine tao_spin_matrix_calc (datum, u, ele_ref, ele_eval, excite_zero)

use tao_interface, dummy => tao_spin_matrix_calc

implicit none

type (tao_data_struct), target :: datum
type (tao_universe_struct), target :: u
type (ele_struct), pointer :: ele_ref, ele_eval, ele0, ele1
type (tao_spin_map_struct), pointer :: spin_map
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: tao_branch

integer ix_ref, ix_ele

real(rp) spin_g(2,6)
real(rp) :: mat3(3,3), quat0(0:3), quat1(0:3), qq(0:3), q0_lnm(0:3)
real(rp) :: n0(3), n1(3), l0(3), m0(3), l1(3), m1(3)
real(rp) quat0_lnm_to_xyz(0:3), quat1_lnm_to_xyz(0:3)

integer i, j, n, p
logical err

character(*), optional :: excite_zero(3)
character(*), parameter :: r_name = 'tao_spin_matrix_calc'

! Init

call tao_spin_tracking_turn_on()

branch => u%model%lat%branch(datum%ix_branch)
tao_branch => u%model%tao_branch(datum%ix_branch)

spin_map => datum%spin_map
spin_map%mat8 = 0
spin_map%valid = .false.

! Concatenate the spin/orbital map

ele1 => branch%ele(0)
if (associated(ele_eval)) ele1 => ele_eval
if (ele1%lord_status == super_lord$) ele1 => pointer_to_slave(ele1, ele1%n_slave)

ele0 => ele1
if (associated(ele_ref)) ele0 => ele_ref
if (ele0%lord_status == super_lord$) ele0 => pointer_to_slave(ele0, ele0%n_slave)

ix_ref = ele0%ix_ele
ix_ele = ele1%ix_ele

call spin_concat_linear_maps (err, spin_map%map1, branch, ix_ref, ix_ele, &
                                      orbit = tao_branch%orbit, excite_zero = excite_zero)
if (err) return

quat0 = spin_map%map1%spin_q(:, 0)

! If 1-turn then calculate n0. Otherwise take the value in the datum.

if (any(spin_map%axis_input%n0 /= 0)) then
  n0 = spin_map%axis_input%n0 / norm2(spin_map%axis_input%n0)
elseif (associated(ele_ref) .and. ix_ele == ix_ref) then
  n0 = real(spin_map%map1%spin_q(1:3,0), rp)
  n0 = n0 / norm2(n0)
else
  n0 = u%model%tao_branch(datum%ix_branch)%orbit(ix_ref)%spin
endif

if (all(n0 == 0)) then
  call tao_set_invalid (datum, 'DATUM SPIN_AXIS%N0 NOT SET AND CANNOT BE COMPUTED.', datum%why_invalid, .true.)
  return
endif

n1 = quat_rotate(quat0, n0)

! Construct coordinate systems (l0, n0, m0) and (l1, n1, m1)

if (any(spin_map%axis_input%l /= 0)) then
  l0 = spin_map%axis_input%l - n0 * dot_product(spin_map%axis_input%l, n0)  ! Make sure l0 is perpendicular to n0.
  l0 = l0 / norm2(l0)
  m0 = cross_product(l0, n0)
elseif (any(spin_map%axis_input%m /= 0)) then
  m0 = spin_map%axis_input%m - n0 * dot_product(spin_map%axis_input%m, n0)  ! Make sure m0 is perpendicular to n0.
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

if (ix_ref == ix_ele) then
  l1 = l0
  m1 = m0
else
  l1 = quat_rotate(quat0, l0)
  m1 = quat_rotate(quat0, m0)
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

spin_map%mat8(1:6,1:6) = spin_map%map1%orb_mat

q0_lnm = quat_mul(quat_inverse(quat1_lnm_to_xyz), quat0, quat0_lnm_to_xyz)
mat3 = quat_to_w_mat(q0_lnm)
spin_map%mat8(7:8,7:8) = mat3(1:3:2,1:3:2)

do p = 1, 6
  quat1 = spin_map%map1%spin_q(:, p)
  qq = quat_mul(quat_inverse(quat1_lnm_to_xyz), quat1, quat0_lnm_to_xyz)
  ! q0_lnm(1) & q0_lnm(3) should be 0 so could drop corresponding terms here.
  spin_map%mat8(7,p) = 2 * (q0_lnm(1)*qq(2) + q0_lnm(2)*qq(1) - q0_lnm(0)*qq(3) - q0_lnm(3)*qq(0))
  spin_map%mat8(8,p) = 2 * (q0_lnm(0)*qq(1) + q0_lnm(1)*qq(0) + q0_lnm(2)*qq(3) + q0_lnm(3)*qq(2))
enddo

!

if (associated(ele_ref)) then
  spin_map%ix_ref   = ele_ref%ix_ele
else
  spin_map%ix_ref   = -1
endif

spin_map%ix_ele   = ix_ele
spin_map%ix_uni   = u%ix_uni
spin_map%axis0    = spin_axis_struct(l0, n0, m0)
spin_map%axis1    = spin_axis_struct(l1, n1, m1)
spin_map%valid    = .true.

end subroutine tao_spin_matrix_calc
