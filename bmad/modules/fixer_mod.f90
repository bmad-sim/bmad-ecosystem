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

if (.not. on .and. .not. fixer%is_on) return

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
! Function fix(branch, to_active) result (is_ok)
!
! Input:
!   to_active   -- logical: If True, set active Twiss from stored. If False, set stored Twiss from active.
!   who         -- logical, optional: Who to set. Possibilities are:
!                   'all', ' ' (default and same as 'all'),
!                   'twiss', 'a-mode', 'b-mode', 'c-mode', 'x-twiss', 'y-twiss', 'cmat', 'dispersion', 'chromatic',
!                   'orbit', 'x-plane', 'y-plane', 'z-plane',
!                   and individual parameters like 'x', 'px', 'cmat_11', etc.
!-

function fix(branch, to_active, who) result (is_ok)

type (branch_struct), target :: branch
type (ele_struct), pointer :: fixer

logical to_active, is_ok
character(*), optional :: who
character(40) whom

!

whom = ''
if (present(who)) whom = who

fixer => branch%ele(branch%ix_fixer)

if (to_active) then
  select case (whom)
  case ('all', '')
  case ('twiss')
  case ('a-mode')
  case ('b-mode')
  case ('c-mode')
  case ('x-twiss')
  case ('y-twiss')
  case ('cmat')
  case ('dispersion')
  case ('chromatic')
  case ('orbit')
  case ('x-plane')
  case ('y-plane')
  case ('z-plane')
  case ('x');               fixer%value(x_stored$)            = fixer%value(x_stored$)
  case ('px');              fixer%value(px_stored$)           = fixer%value(x_stored$)
  case ('y');               fixer%value(y_stored$)            = fixer%value(x_stored$)
  case ('py');              fixer%value(py_stored$)           = fixer%value(x_stored$)
  case ('z');               fixer%value(z_stored$)            = fixer%value(x_stored$)
  case ('pz');              fixer%value(pz_stored$)           = fixer%value(x_stored$)
  case ('beta_a');          fixer%value(beta_a_stored$)       = fixer%value(x_stored$)
  case ('alpha_a');         fixer%value(alpha_a_stored$)      = fixer%value(x_stored$)
  case ('beta_b');          fixer%value(beta_b_stored$)       = fixer%value(x_stored$)
  case ('alpha_b');         fixer%value(alpha_b_stored$)      = fixer%value(x_stored$)
  case ('phi_a');           fixer%value(phi_a_stored$)        = fixer%value(x_stored$)
  case ('phi_b');           fixer%value(phi_b_stored$)        = fixer%value(x_stored$)
  case ('mode_flip');       fixer%value(mode_flip_stored$)    = fixer%value(x_stored$)
  case ('eta_x');           fixer%value(eta_x_stored$)        = fixer%value(x_stored$)
  case ('etap_x');          fixer%value(etap_x_stored$)       = fixer%value(x_stored$)
  case ('eta_y');           fixer%value(eta_y_stored$)        = fixer%value(x_stored$)
  case ('etap_y');          fixer%value(etap_y_stored$)       = fixer%value(x_stored$)
  case ('cmat_11');         fixer%value(cmat_11_stored$)      = fixer%value(x_stored$)
  case ('cmat_12');         fixer%value(cmat_12_stored$)      = fixer%value(x_stored$)
  case ('cmat_21');         fixer%value(cmat_21_stored$)      = fixer%value(x_stored$)
  case ('cmat_22');         fixer%value(cmat_22_stored$)      = fixer%value(x_stored$)
  case ('dbeta_dpz_a');     fixer%value(dbeta_dpz_a_stored$)  = fixer%value(x_stored$)
  case ('dbeta_dpz_b');     fixer%value(dbeta_dpz_b_stored$)  = fixer%value(x_stored$)
  case ('dalpha_dpz_a');    fixer%value(dalpha_dpz_a_stored$) = fixer%value(x_stored$)
  case ('dalpha_dpz_b');    fixer%value(dalpha_dpz_b_stored$) = fixer%value(x_stored$)
  case ('deta_dpz_x');      fixer%value(deta_dpz_x_stored$)   = fixer%value(x_stored$)
  case ('deta_dpz_y');      fixer%value(deta_dpz_y_stored$)   = fixer%value(x_stored$)
  case ('detap_dpz_x');     fixer%value(detap_dpz_x_stored$)  = fixer%value(x_stored$)
  case ('detap_dpz_y');     fixer%value(detap_dpz_y_stored$)  = fixer%value(x_stored$)
  case default;             is_ok = .false.
end select

else

endif

end function fix

end module
