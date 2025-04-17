!+
! Subroutine set_twiss(branch, twiss_ele, ix_ele, match_deta_ds, err_flag, print_err)
!
! For an open lattice branch, set the beginning Twiss, coupling, and dispersion of the branch so that the 
! Twiss parameters at the "target" element branch%ele(ix_ele) are the same as the Twiss parameters in twiss_ele.
!
! It is OK for twiss_ele to be an element in the branch (this routine makes a copy to be safe).
!
! Since the the reference orbit at twiss_ele may be different from the reference orbit at the target ele,
! it is not possible, in general, to match both etap and deta_ds. Which is matched is set by match_deta_ds.
!
! Additionally, It is not possible to simultaneously match both (x,y) deta_ds and (a,b) normal mode deta_ds.
! This routine matches (x,y) deta_ds only (when match_deta_ds = True).
!
! Input:
!   branch        -- branch_struct: Branch to modify.
!   twiss_ele     -- ele_struct: Element with desired Twiss parameters.
!   ix_ele        -- integer: Match branch%ele(ix_ele) Twiss to twiss_ele.
!   match_deta_ds -- logical: If True, match deta_ds. If False, match etap.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!   print_err     -- logical, optional: Print an error message if there is an error? Default is True.
!-

subroutine set_twiss(branch, twiss_ele, ix_ele, match_deta_ds, err_flag, print_err)

use bmad_routine_interface, dummy => set_twiss

implicit none

type (branch_struct), target :: branch
type (ele_struct) twiss_ele, ele0, ele1
type (ele_struct), pointer :: ele

real(rp) xmat(6,6), xvec(6), rel_p

integer ix_ele
integer ie

logical match_deta_ds, err_flag
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

! Note: twiss_propagate1 treats deta_ds as a dependent variable to etap

call transfer_twiss(twiss_ele, ele0)
ele0%map_ref_orb_out = twiss_ele%map_ref_orb_out

if (match_deta_ds) then
  ele => branch%ele(ix_ele)
  rel_p = 1.0_rp + ele%map_ref_orb_out%vec(6)
  ele0%x%etap = ele0%x%deta_ds * rel_p + ele%map_ref_orb_out%vec(2) / rel_p
  ele0%y%etap = ele0%y%deta_ds * rel_p + ele%map_ref_orb_out%vec(4) / rel_p
endif

do ie = ix_ele, 1, -1
  ele => branch%ele(ie)
  call mat_inverse(ele%mat6, ele1%mat6)
  ele1%key = ele%key
  ele1%map_ref_orb_in = ele%map_ref_orb_out
  ele1%map_ref_orb_out = ele%map_ref_orb_in
  call twiss_propagate1(ele0, ele1, err_flag)
  if (err_flag) then
    if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'CANNOT BACK PROPAGATE TWISS.')
    return
  endif
  call transfer_twiss(ele1, ele0)
enddo

call transfer_twiss(ele0, branch%ele(0))

call twiss_propagate_all(branch%lat, branch%ix_branch, err_flag)

end subroutine set_twiss
