!+
! Subroutine tracking_rad_map_setup (ele, tollerance, ref_edge, rad_map)
!
! Routine to setup the radiation damping and excitation matrices for an element.
! These matrices can then be used when tracking to simulate radiation effects.
!
! Restriction: ele must associated with a lattice (ele%branch must point to a branch).
!
! Input:
!   ele         -- ele_struct: Element to setup. Matrices will be with respect to the map reference orbit.
!   tollerance  -- real(rp): Tolerance used for the computation.
!   ref_edge    -- integer: Edge that the matrices are referenced to. upstream_end$ or downstream_end$.
!
! Output:
!   rad_map     -- rad_map_struct: Structure holding the matrices.
!-

subroutine tracking_rad_map_setup (ele, tollerance, ref_edge, rad_map)

use rad_6d_mod, dummy => tracking_rad_map_setup
use f95_lapack, only: dpotrf_f95

implicit none

type (ele_struct), target :: ele
type (rad_map_struct) rad_map
type (branch_struct), pointer :: branch
type (coord_struct) orb0, orb1

real(rp) tollerance
real(rp) tol, m_inv(6,6)

integer ref_edge
integer i, info

logical err

!

branch => pointer_to_branch(ele)
tol = tollerance / branch%param%total_length
orb0 = ele%map_ref_orb_in
orb1 = ele%map_ref_orb_out

rad_map = rad_map_struct()
if (orb0%vec(2) == orb1%vec(2) .and. orb0%vec(4) == orb1%vec(4) .and. ele%key /= sbend$) return

!

call rad1_damp_and_stoc_mats (ele, .true., orb0, orb1, rad_map, &
                                                  tol*branch%param%g2_integral, tol*branch%param%g3_integral)

select case (ref_edge)
case (upstream_end$)
  m_inv = mat_symp_conj(ele%mat6)
  rad_map%xfer_damp_mat = matmul(m_inv, rad_map%xfer_damp_mat)
  rad_map%xfer_damp_vec = matmul(m_inv, rad_map%xfer_damp_vec)
  rad_map%stoc_mat = matmul(matmul(m_inv, rad_map%stoc_mat), transpose(m_inv))
  rad_map%ref_orb = orb0%vec
case (downstream_end$)
  rad_map%ref_orb = orb1%vec
end select

! The Cholesky decomp can fail due to roundoff error in rad1_damp_and_stoc_mats. 
! This is especially true if there is little radiation. In this case just set the stoc mat to zero.

call dpotrf_f95 (rad_map%stoc_mat, 'L', info = info)
if (info /= 0) then
  rad_map%stoc_mat = 0  ! Cholesky failed
  return
endif

do i = 2, 6
  rad_map%stoc_mat(1:i-1, i) = 0
enddo

end subroutine
