!+
! Subroutine set_twiss(branch, twiss_ele, ix_ele, err_flag, print_err)
!
! Routine to set the beginning Twiss of a lattice branch so that the Twiss parameters at branch%ele(ix_ele) are
! the same as the Twiss parameters in twiss_ele.
!
! It is OK for twiss_ele to be an element in the branch (this routine makes a copy to be safe).
!
! Input:
!   branch      -- branch_struct: Branch to modify.
!   twiss_ele   -- ele_struct: Element with desired Twiss parameters.
!   ix_ele      -- integer: Match branch%ele(ix_ele) Twiss to twiss_ele.
!   err_flag    -- logical: Set True if there is an error. False otherwise.
!   print_err   -- logical, optional: Print an error message if there is an error? Default is True.
!-

subroutine set_twiss(branch, twiss_ele, ix_ele, err_flag, print_err)

use bmad, dummy => set_twiss

implicit none

type (branch_struct), target :: branch
type (ele_struct) twiss_ele, t_ele, ele0, ele1
type (ele_struct), pointer :: b_ele

real(rp) xmat(6,6), xvec(6)

integer ix_ele
integer ie

logical err_flag
logical, optional :: print_err

character(*), parameter :: r_name = 'set_twiss'

!

err_flag = .true.

if (branch%param%geometry == closed$) then
  if (logic_option(.true., print_err)) call out_io (s_error$, r_name, &
                                             'TWISS CANNOT BE SET IN A BRANCH WITH AN OPEN GEOMETRY.')
  return
endif

if (ix_ele < 0 .or. ix_ele > branch%n_ele_track) then
  if (logic_option(.true., print_err)) then
    call out_io (s_error$, r_name, 'IX_ELE ELEMENT INDEX (' // &
                                int_str(ix_ele) // ') OUT OF RANGE: [0, ' // int_str(branch%n_ele_track) // '].')
  endif
  return
endif

!

if (ix_ele == 0) then
  call transfer_twiss(twiss_ele, branch%ele(0))
  call twiss_propagate_all(branch%lat, branch%ix_branch, err_flag)
  return
endif

!

call transfer_twiss(twiss_ele, ele0)

call transfer_matrix_calc(branch%lat, ele1%mat6, ele1%vec0, 0, ix_ele, branch%ix_branch)
call mat_inverse(ele1%mat6, ele1%mat6)
ele1%key = quadrupole$
ele1%map_ref_orb_in = branch%ele(ix_ele)%map_ref_orb_out
ele1%map_ref_orb_out = branch%ele(1)%map_ref_orb_in
call twiss_propagate1(ele0, ele1, err_flag)
if (err_flag) then
  if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'CANNOT BACK PROPAGATE TWISS.')
  return
endif

call transfer_twiss(ele1, branch%ele(0))
call twiss_propagate_all(branch%lat, branch%ix_branch, err_flag)

end subroutine set_twiss
