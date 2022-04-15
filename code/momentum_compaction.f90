!+
! Function momentum_compaction(branch, mat6, branch_len) result (mom_comp)
!
! Routine to compute the momentum compaction factor.
!
! The momentum compaction factor can be calculated from either branch or mat6 and branch_len.
! Therefore only branch should be present or mat6 and branch_len should be present.
!
! Note:
!   slip_factor      = mom_comp - 1 / gamma^2
!   gamma_transition = sqrt(mom_comp)
!
! Input:
!   branch      -- branch_struct, optional: Lattice branch to calculate on.
!   mat6(6,6)   -- real(rp), optional: One turn matrix with RF off.
!   branch_len  -- real(rp), optional: Length of lattice branch.
!
! Output:
!   mom_comp    -- real(rp): Momentum compaction.
!-

function momentum_compaction(branch, mat6, branch_len) result (mom_comp)

use bmad_routine_interface, dummy => momentum_compaction

implicit none

type (branch_struct), optional, target :: branch
type (ele_struct), pointer :: ele
type (ele_struct) ele0
type (coord_struct) orbit

real(rp), optional :: mat6(6,6), branch_len
real(rp) mom_comp
real(rp) m(6,6), mat4(4,4), eta_vec(4), length
integer ie

! Form 1-turn matrix

if (present(branch)) then
  length = branch%param%total_length

  call mat_make_unit(m)
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%key == rfcavity$) then
      orbit = ele%map_ref_orb_in
      call track_a_drift(orbit, ele%value(l$), m, .true.)
    else
      m = matmul(ele%mat6, m)
    endif
  enddo

else
  length = branch_len
  m = mat6
endif

! Calc momentum compaction

call mat_make_unit(mat4)
mat4 = mat4 - m(1:4,1:4)
call mat_inverse (mat4, mat4)
eta_vec = matmul(mat4, m(1:4,6))

mom_comp = -(sum(eta_vec * m(5,1:4)) + m(5,6)) / length

end function
