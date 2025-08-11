module fixer_mod

use bmad_routine_interface

implicit none

contains

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!+
! Subroutine set_active_fixer(fixer, is_on)
!
! Set the acvitive fixer element.
! All other fixer/beginning_ele elements in the branch will be deactivated.
!
! If is_on is True (default), the fixer argument becomes the active fixer.
! If is_on is False, and fixer%is_on is also False, there is nothing to be done.
! If is_on is False, and fixer%is_on is True, turn this fixer off and turn on the beginning element.
!
! Input:
!   fixer         -- ele_struct: Fixer element to make active.
!   is_on         -- logical, optional: If True (default), make this fixer the active element. 
!                     If False, make the beginning element active.
!
! Output:
!   fixer         -- ele_struct: Element is now active.
!-

subroutine set_active_fixer(fixer, is_on)

type (ele_struct), target :: fixer
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch

integer ie
logical, optional :: is_on
logical on
character(*), parameter :: r_name = 'set_active_fixer'


!

if (fixer%key /= fixer$ .and. fixer%key /= beginning_ele$) then
  call out_io(s_error$, r_name, 'Element to set is not a fixer element nor a beginning element.')
  return
endif

branch => fixer%branch
on = logic_option(.true., is_on)

if (.not. is_on .and. .not. fixer%is_on) return

do ie = 0, branch%n_ele_track
  ele => branch%ele(ie)
  if (ele%key /= beginning_ele$ .and. ele%key /= fixer$) cycle
  ele%is_on = .false.
enddo

if (on) then
  fixer%is_on = .true.
  branch%ix_fixer = fixer%ix_ele
else
  branch%ele(0)%is_on = .true.
  branch%ix_fixer = 0
endif

end subroutine set_active_fixer

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

!+
! Subroutine fix(branch, to_active, who, err_flag)
!
! Input:
!   to_active   -- logical: If True, set active Twiss from stored. If False, set stored Twiss from active.
!   who         -- logical, optional: Who to set. Possibilities are:
!                   'all', ' ' (default and same as 'all'),
!                   'twiss', 'a-mode', 'b-mode', 'c-mode', 'x-twiss', 'y-twiss', 'cmat', 'dispersion', 'chromatic',
!                   'orbit', 'x-plane', 'y-plane', 'z-plane',
!                   and individual parameters like 'x', 'px', 'cmat_11', etc.
!-

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!+
! Function find_active_fixer(branch) result (fixer)
!
! Input:
!   branch
!
! Output:
!   fixer
!-

end module
