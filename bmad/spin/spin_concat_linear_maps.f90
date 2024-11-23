!+
! Subroutine spin_concat_linear_maps (err_flag, map1, branch, n1, n2, q_ele, orbit, excite_zero)
!
! Routine to concatenate element spin/orbit maps in the range branch%ele(n1+1:n2)
! This routine will wrap around the ends of the lattice so n2 may be less than n1.
! In this case the range will be [n1+1:end] + [beginning:n2].
!
! If a Taylor map exists for an individual element, that map will be reused.
! If not, a new map will be made for the element. If a map is made, orbit(:) will be 
! used as the reference orbit. If it is not present, ele%map_ref_orb_in will be used.
!
! The excite_zero(:) argument is used to zero out linear ds_vec/dr_vec terms in the element's spin 
! transfer map where s_vec = spin vector and r_vec = orbital phase space vector. 
! This can be used to give insights into where spin matching is failing.
! The excite_zero(1) term is a list for nullifying of r_vec = x or px terms, excite_zero(2)
! is for r_vec = y or py terms, and excite_zero(3) is for r_vec = z or pz terms.
!
! Input:
!   branch            -- branch_struct: Lattice branch.
!   n1                -- integer: Starting element index. Start at element downstream end.
!   n2                -- integer: Ending element index. End at element downstream end
!   orbit(0:)         -- coord_struct, optional: Reference orbit used if maps must be created.
!   excite_zero(3)    -- character(*), optional: Three lists of elements where ds_vec/dr_vec terms are zeroed.
!
! Output:
!   err_flag          -- logical: Set True if there is an error. False otherwise.
!   map1              -- spin_orbit_map1_struct: Map with element spin/orbit maps concatenated.
!   map1_ele(:)       -- spin_orbit_map1_struct, optional: Individual spin/orbit maps.
!-

subroutine spin_concat_linear_maps (err_flag, map1, branch, n1, n2, map1_ele, orbit, excite_zero)

use ptc_interface_mod, dummy => spin_concat_linear_maps

implicit none

type (spin_orbit_map1_struct) map1
type (spin_orbit_map1_struct), optional :: map1_ele(0:)
type (branch_struct), target :: branch
type (coord_struct), optional :: orbit(0:)
type (ele_pointer_struct), allocatable :: eles(:)

integer n1, n2
integer i, j, n_loc

logical err_flag

character(*), optional :: excite_zero(3)

! Mark elements to be nullified

err_flag = .false.
branch%ele%ixx = 0
branch%lat%ele%ixx = 0   ! To take case of lord elements.

if (present(excite_zero)) then
  do i = 1, 3
    call lat_ele_locator (excite_zero(i), branch%lat, eles, n_loc, err_flag, .false., branch%ix_branch)
    if (err_flag) return
    do j = 1, n_loc
      call set_this_ele(eles(j)%ele, i)
    enddo
  enddo
endif

! Make the map

call map1_make_unit(map1)

if (n2 <= n1) then
  call concat_this_map(n1+1, branch%n_ele_track)
  call concat_this_map(1, n2)
else
  call concat_this_map(n1+1, n2)
  if (n1 == 0 .and. present(map1_ele)) call map1_make_unit(map1_ele(0))
endif

call spin_map1_normalize(map1%spin_q)

!------------------------------------------------------
contains

recursive subroutine set_this_ele(ele, ix_null)

type (ele_struct) :: ele
integer ix_null, k

!

ele%ixx = ibset(ele%ixx, ix_null)
do k = 1, ele%n_slave
  call set_this_ele(pointer_to_slave(ele, k), ix_null)
enddo

end subroutine set_this_ele
!------------------------------------------------------
! contains

subroutine concat_this_map(n1, n2)

type (ele_struct), pointer :: ele
type (ele_struct) ele2
type (taylor_struct), pointer :: st
type (spin_orbit_map1_struct) q1
type (coord_struct) :: map_ref_orb

real(rp) vec0(6), ref_orb(6), mat6(6,6), vec(0)
integer n1, n2
integer ie, i, k, n, p

!

do ie = n1, n2
  if (ie == 0) then
    if (present(map1_ele)) call map1_make_unit(map1_ele(0))
    cycle
  endif
  ele => branch%ele(ie)

  if (present(orbit)) then
    ref_orb = orbit(ie-1)%vec
    map_ref_orb = orbit(ie-1)
  else
    ref_orb = 0
    map_ref_orb = ele%map_ref_orb_in
  endif

  ! Spin map

  if (.not. associated(ele%spin_taylor(0)%term)) call ele_to_spin_taylor(ele, branch%param, map_ref_orb)
  q1%spin_q = spin_taylor_to_linear(ele%spin_taylor, .false., ref_orb - ele%spin_taylor_ref_orb_in, ele%is_on)

  ! Orbital map

  if (associated(ele%taylor(1)%term)) then
    call taylor_to_mat6 (ele%taylor, ref_orb, q1%vec0, q1%orb_mat)
  elseif (all(ele%map_ref_orb_in%vec == ref_orb)) then
    q1%orb_mat = ele%mat6
    q1%vec0 = ele%vec0
  else
    call transfer_ele(ele, ele2)
    map_ref_orb%vec = ref_orb
    call make_mat6(ele2, branch%param, map_ref_orb)
    q1%orb_mat = ele2%mat6
    q1%vec0 = ele2%vec0
  endif

  if (ele%ixx /= 0) then
    do i = 1, 3
      if (.not. btest(ele%ixx, i)) cycle
      k = 2*i-1
      q1%spin_q(:, k:k+1) = 0
    enddo
  endif

  if (present(map1_ele)) then
    map1_ele(ie) = q1
  endif

  map1 = q1 * map1
enddo

end subroutine

end subroutine
