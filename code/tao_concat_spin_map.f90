!+
! Subroutine tao_concat_spin_map (q_map, branch, n1, n2, q_ele)
!
! Routine to concatenate element spin/orbit maps in the range branch%ele(n1+1:n2)
!
! Input:
!   branch    -- branch_struct: Lattice branch.
!   n1        -- integer: Starting element index. Start at element downstream end.
!   n2        -- integer: Ending element index. End at element downstream end
!
! Output:
!   q_map     -- c_linear_map: Map with element spin/orbit maps concatenated.
!   q_ele(:)  -- c_linear_map, optional: Individual spin/orbit maps.
!-

subroutine tao_concat_spin_map (q_map, branch, n1, n2, q_ele)

use tao_interface, dummy => tao_concat_spin_map
use ptc_interface_mod
use pointer_lattice, only: c_linear_map, operator(*), assignment(=)

implicit none

type (c_linear_map) q_map
type (c_linear_map), optional :: q_ele(:)
type (branch_struct), target :: branch

integer n1, n2

!

q_map = 0

if (n2 <= n1) then
  call concat_this_map(n1+1, branch%n_ele_track)
  call concat_this_map(1, n2)
else
  call concat_this_map(n1+1, n2)
endif

!------------------------------------------------------
contains

subroutine concat_this_map(n1, n2)

type (ele_struct), pointer :: ele
type (taylor_struct), pointer :: st
type (c_linear_map) q1

real(rp) vec0(6), mat6(6,6)
integer n1, n2
integer ie, i, k, n, p
logical st_on

!

do ie = n1, n2
  if (ie == 0) cycle
  ele => branch%ele(ie)
  if (.not. associated(ele%spin_taylor(0)%term)) then
    st_on = bmad_com%spin_tracking_on
    bmad_com%spin_tracking_on = .true.
    call ele_to_taylor(ele, branch%param, ele%map_ref_orb_in)
    bmad_com%spin_tracking_on = st_on
  endif

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
  if (present(q_ele)) then
    q_ele(ie) = q1
  endif

  q_map = q1 * q_map
enddo

end subroutine

end subroutine
