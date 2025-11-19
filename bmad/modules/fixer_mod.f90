module fixer_mod

use bmad_routine_interface

implicit none

contains

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!+
! Subroutine set_active_fixer(fixer, is_on, orbit)
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
!   orbit         -- coord_struct, optional: Load with stored fixer phase space and spin values.
!-

subroutine set_active_fixer(fixer, is_on, orbit)

type (ele_struct), target :: fixer
type (coord_struct), optional :: orbit
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch

integer ie
logical, optional :: is_on
logical on, is_ok
character(*), parameter :: r_name = 'set_active_fixer'

!

if (fixer%key /= fixer$ .and. fixer%key /= beginning_ele$) then
  call out_io(s_error$, r_name, 'Element to set is not a fixer element nor a beginning element.')
  return
endif

if (fixer%value(beta_a_stored$) <= 0 .or. fixer%value(beta_b_stored$) <= 0) then
  call out_io(s_error$, r_name, 'Fixer element does not have beta_a or beta_b set!. Active fixer not set.')
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

is_ok = transfer_fixer_params(branch%ele(branch%ix_fixer), .false., orbit)

end subroutine set_active_fixer

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

!+
! Function transfer_fixer_params(fixer, to_stored, orbit, who) result (is_ok)
!
! Set parameters of fixer.
!
! Input:
!   fixer       -- ele_struct: Fixer element to set.
!   to_stored   -- logical: If False, set real Twiss from stored. If True, set stored Twiss from real.
!   orbit       -- coord_struct, optional: Used for 'phase_space' transfers.
!   who         -- logical, optional: Who to set. Possibilities are:
!                   Groups:
!                     'all', ' ' (default and same as 'all') Note: This excludes all 'start' sets.,
!                     'twiss', 'a_twiss', 'b_twiss', 'cmat', 'x_dispersion', 'y_dispersion', 'dispersion', 'chromatic',
!                     'orbit', 'phase_space', 'spin', 'x_plane', 'y_plane', 'z_plane', 
!                     'start', 'start_spin', 'start_phase_space', 
!                   Individula Parameters:
!                     'x', 'px', 'cmat_11', etc.
!
! Output:
!   is_ok       -- logical
!-

function transfer_fixer_params(fixer, to_stored, orbit, who) result (is_ok)

type (ele_struct), target :: fixer
type (coord_struct), optional :: orbit
type (branch_struct), pointer :: branch

logical to_stored, is_ok
character(*), optional :: who
character(40) whom
character(*), parameter :: r_name = 'transfer_fixer_params'

!

whom = 'all'
if (present(who)) then
  if (who /= '') whom = who
endif

is_ok = .true.
branch => fixer%branch
call fix_this(branch, fixer, to_stored, is_ok, whom, orbit)

!----------------------------------------------------------------------------------------------------
contains

recursive subroutine fix_this(branch, fixer, to_stored, is_ok, whom, orbit)

type (branch_struct) branch
type (ele_struct) fixer
type (coord_struct), optional :: orbit

integer ix
logical to_stored, is_ok
character(*) whom
character(20) switch

!

call match_word(whom, [character(20):: &
                  'all', 'twiss', 'a_twiss', 'b_twiss', 'x_dispersion', 'y_dispersion', 'cmat', 'dispersion', &
                  'chromatic', 'orbit', 'spin', 'phase_space', 'x_plane', 'y_plane', 'z_plane', &
                  'start_orbit', 'start_spin', 'start_phase_space', 'start_x_plane',  'start_y_plane', 'start_z_plane', &
                  'spin_x', 'spin_y', 'spin_z', 'start_spin_x', 'start_spin_y', 'start_spin_z', &
                  'x', 'px', 'y', 'py', 'z', 'pz', 'start_x', 'start_px', 'start_y', 'start_py', 'start_z', 'start_pz', &
                  'beta_a', 'alpha_a', 'phi_a', 'dbeta_dpz_a', 'dalpha_dpz_a', &
                  'beta_b', 'alpha_b', 'phi_b', 'dbeta_dpz_b', 'dalpha_dpz_b', &
                  'eta_x', 'etap_x', 'deta_dpz_x', 'detap_dpz_x', &
                  'eta_y', 'etap_y', 'deta_dpz_y', 'detap_dpz_y', &
                  'mode_flip', 'cmat_11', 'cmat_12', 'cmat_21', 'cmat_22', &
                  'dcmat_dpz_11', 'dcmat_dpz_12', 'dcmat_dpz_21', 'dcmat_dpz_22'], &
                      ix, .false., .true., switch)

if (ix == 0) then
  is_ok = .false.
  call out_io (s_error$, r_name, 'Unknown fixer element parameter name: ' // whom)

elseif (ix < 0) then
  is_ok = .false.
  call out_io (s_error$, r_name, 'Fixer element parameter name has multiple matches: ' // whom)

else
  select case (switch)
  case ('all')
    call fix_this(branch, fixer, to_stored, is_ok, 'twiss', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'phase_space', orbit)

  case ('twiss')
    call fix_this(branch, fixer, to_stored, is_ok, 'a_twiss', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'b_twiss', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dispersion', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'cmat', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'chromatic', orbit)

  case ('a_twiss')
    call fix_this(branch, fixer, to_stored, is_ok, 'beta_a', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'alpha_a', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'phi_a', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dbeta_dpz_a', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dalpha_dpz_a', orbit)

  case ('b_twiss')
    call fix_this(branch, fixer, to_stored, is_ok, 'beta_b', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'alpha_b', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'phi_b', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dbeta_dpz_b', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dalpha_dpz_b', orbit)

  case ('dispersion')
    call fix_this(branch, fixer, to_stored, is_ok, 'x_dispersion', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'y_dispersion', orbit)

  case ('x_dispersion')
    call fix_this(branch, fixer, to_stored, is_ok, 'eta_x', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'etap_x', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'deta_dpz_x', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'detap_dpz_x', orbit)

  case ('y_dispersion')
    call fix_this(branch, fixer, to_stored, is_ok, 'eta_y', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'etap_y', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'deta_dpz_y', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'detap_dpz_y', orbit)

  case ('cmat')
    call fix_this(branch, fixer, to_stored, is_ok, 'cmat_11', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'cmat_12', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'cmat_21', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'cmat_22', orbit)

  case ('chromatic')
    call fix_this(branch, fixer, to_stored, is_ok, 'dbeta_dpz_a', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dalpha_dpz_a', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dbeta_dpz_b', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dalpha_dpz_b', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'deta_dpz_x', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'detap_dpz_x', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'deta_dpz_y', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'detap_dpz_y', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dcmat_dpz_11', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dcmat_dpz_12', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dcmat_dpz_21', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'dcmat_dpz_22', orbit)

  case ('phase_space')
    call fix_this(branch, fixer, to_stored, is_ok, 'spin', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'orbit', orbit)

  case ('orbit')
    call fix_this(branch, fixer, to_stored, is_ok, 'x_plane', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'y_plane', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'z_plane', orbit)

  case ('x_plane')
    call fix_this(branch, fixer, to_stored, is_ok, 'x', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'px', orbit)

  case ('y_plane')
    call fix_this(branch, fixer, to_stored, is_ok, 'y', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'py', orbit)

  case ('z_plane')
    call fix_this(branch, fixer, to_stored, is_ok, 'z', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'pz', orbit)

  case ('spin')
    call fix_this(branch, fixer, to_stored, is_ok, 'spin_x', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'spin_y', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'spin_z', orbit)

  case ('start_orbit')
    call fix_this(branch, fixer, to_stored, is_ok, 'start_x_plane', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'start_y_plane', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'start_z_plane', orbit)

  case ('start_spin')
    call fix_this(branch, fixer, to_stored, is_ok, 'start_spin_x', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'start_spin_y', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'start_spin_z', orbit)

  case ('start_phase_space')
    call fix_this(branch, fixer, to_stored, is_ok, 'start_spin', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'start_orbit', orbit)

  case ('start_x_plane')
    call fix_this(branch, fixer, to_stored, is_ok, 'start_x', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'start_px', orbit)

  case ('start_y_plane')
    call fix_this(branch, fixer, to_stored, is_ok, 'start_y', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'start_py', orbit)

  case ('start_z_plane')
    call fix_this(branch, fixer, to_stored, is_ok, 'start_z', orbit)
    call fix_this(branch, fixer, to_stored, is_ok, 'start_pz', orbit)

  case default
    call fix_this1(branch, fixer, to_stored, is_ok, whom, orbit)
    return

  end select
endif

end subroutine fix_this

!----------------------------------------------------------------------------------------------------
! contains

subroutine fix_this1(branch, fixer, to_stored, is_ok, whom, orbit)

type (branch_struct) branch
type (ele_struct) fixer
type (coord_struct), optional :: orbit
type (coord_struct) this_orb

real(rp) gc2
logical to_stored, is_ok
character(*) whom

!

if (present(orbit)) then
  this_orb = orbit
else
  this_orb = coord_struct()
endif

if (to_stored) then
  select case (whom)
  case ('spin_x');          fixer%value(spin_x_stored$)       = this_orb%spin(1)
  case ('spin_y');          fixer%value(spin_y_stored$)       = this_orb%spin(2)
  case ('spin_z');          fixer%value(spin_z_stored$)       = this_orb%spin(3)

  case ('start_spin_x');    fixer%value(spin_x_stored$)       = branch%particle_start%spin(1)
  case ('start_spin_y');    fixer%value(spin_y_stored$)       = branch%particle_start%spin(2)
  case ('start_spin_z');    fixer%value(spin_z_stored$)       = branch%particle_start%spin(3)

  case ('x');               fixer%value(x_stored$)            = this_orb%vec(1)
  case ('px');              fixer%value(px_stored$)           = this_orb%vec(2)
  case ('y');               fixer%value(y_stored$)            = this_orb%vec(3)
  case ('py');              fixer%value(py_stored$)           = this_orb%vec(4)
  case ('z');               fixer%value(z_stored$)            = this_orb%vec(5)
  case ('pz');              fixer%value(pz_stored$)           = this_orb%vec(6)

  case ('start_x');         fixer%value(x_stored$)            = branch%particle_start%vec(1)
  case ('start_px');        fixer%value(px_stored$)           = branch%particle_start%vec(2)
  case ('start_y');         fixer%value(y_stored$)            = branch%particle_start%vec(3)
  case ('start_py');        fixer%value(py_stored$)           = branch%particle_start%vec(4)
  case ('start_z');         fixer%value(z_stored$)            = branch%particle_start%vec(5)
  case ('start_pz');        fixer%value(pz_stored$)           = branch%particle_start%vec(6)

  case ('beta_a');          fixer%value(beta_a_stored$)       = fixer%a%beta
  case ('alpha_a');         fixer%value(alpha_a_stored$)      = fixer%a%alpha
  case ('phi_a');           fixer%value(phi_a_stored$)        = fixer%a%phi
  case ('dbeta_dpz_a');     fixer%value(dbeta_dpz_a_stored$)  = fixer%a%dbeta_dpz
  case ('dalpha_dpz_a');    fixer%value(dalpha_dpz_a_stored$) = fixer%a%dalpha_dpz

  case ('beta_b');          fixer%value(beta_b_stored$)       = fixer%b%beta
  case ('alpha_b');         fixer%value(alpha_b_stored$)      = fixer%b%alpha
  case ('phi_b');           fixer%value(phi_b_stored$)        = fixer%b%phi
  case ('dbeta_dpz_b');     fixer%value(dbeta_dpz_b_stored$)  = fixer%b%dbeta_dpz
  case ('dalpha_dpz_b');    fixer%value(dalpha_dpz_b_stored$) = fixer%b%dalpha_dpz

  case ('eta_x');           fixer%value(eta_x_stored$)        = fixer%x%eta
  case ('etap_x');          fixer%value(etap_x_stored$)       = fixer%x%etap
  case ('deta_dpz_x');      fixer%value(deta_dpz_x_stored$)   = fixer%x%deta_dpz
  case ('detap_dpz_x');     fixer%value(detap_dpz_x_stored$)  = fixer%x%detap_dpz

  case ('eta_y');           fixer%value(eta_y_stored$)        = fixer%y%eta
  case ('etap_y');          fixer%value(etap_y_stored$)       = fixer%y%etap
  case ('deta_dpz_y');      fixer%value(deta_dpz_y_stored$)   = fixer%y%deta_dpz
  case ('detap_dpz_y');     fixer%value(detap_dpz_y_stored$)  = fixer%y%detap_dpz

  case ('mode_flip');       fixer%value(mode_flip_stored$)    = int_logic(fixer%mode_flip)
  case ('cmat_11');         fixer%value(cmat_11_stored$)      = fixer%c_mat(1,1)
  case ('cmat_12');         fixer%value(cmat_12_stored$)      = fixer%c_mat(1,2)
  case ('cmat_21');         fixer%value(cmat_21_stored$)      = fixer%c_mat(2,1)
  case ('cmat_22');         fixer%value(cmat_22_stored$)      = fixer%c_mat(2,2)

  case ('dcmat_dpz_11');    fixer%value(dcmat_dpz_11_stored$) = fixer%dc_mat_dpz(1,1)
  case ('dcmat_dpz_12');    fixer%value(dcmat_dpz_12_stored$) = fixer%dc_mat_dpz(1,2)
  case ('dcmat_dpz_21');    fixer%value(dcmat_dpz_21_stored$) = fixer%dc_mat_dpz(2,1)
  case ('dcmat_dpz_22');    fixer%value(dcmat_dpz_22_stored$) = fixer%dc_mat_dpz(2,2)

  case default
    is_ok = .false.
    call out_io (s_error$, r_name, 'Fixer element parameter name not recognized: ' // whom)
  end select

else
  select case (whom)
  case ('spin_x');          this_orb%spin(1)                     = fixer%value(spin_x_stored$)
  case ('spin_y');          this_orb%spin(2)                     = fixer%value(spin_y_stored$)
  case ('spin_z');          this_orb%spin(3)                     = fixer%value(spin_z_stored$)

  case ('start_spin_x');    branch%particle_start%spin(1)     = fixer%value(spin_x_stored$)
  case ('start_spin_y');    branch%particle_start%spin(2)     = fixer%value(spin_y_stored$)
  case ('start_spin_z');    branch%particle_start%spin(3)     = fixer%value(spin_z_stored$)

  case ('x');               this_orb%vec(1)                      = fixer%value(x_stored$)
  case ('px');              this_orb%vec(2)                      = fixer%value(px_stored$)
  case ('y');               this_orb%vec(3)                      = fixer%value(y_stored$)
  case ('py');              this_orb%vec(4)                      = fixer%value(py_stored$)
  case ('z');               this_orb%vec(5)                      = fixer%value(z_stored$)
  case ('pz');              this_orb%vec(6)                      = fixer%value(pz_stored$)

  case ('start_x');         branch%particle_start%vec(1)      = fixer%value(x_stored$)
  case ('start_px');        branch%particle_start%vec(2)      = fixer%value(px_stored$)
  case ('start_y');         branch%particle_start%vec(3)      = fixer%value(y_stored$)
  case ('start_py');        branch%particle_start%vec(4)      = fixer%value(py_stored$)
  case ('start_z');         branch%particle_start%vec(5)      = fixer%value(z_stored$)
  case ('start_pz');        branch%particle_start%vec(6)      = fixer%value(pz_stored$)

  case ('beta_a');          fixer%a%beta                      = fixer%value(beta_a_stored$)
  case ('alpha_a');         fixer%a%alpha                     = fixer%value(alpha_a_stored$)
  case ('phi_a');           fixer%a%phi                       = fixer%value(phi_a_stored$)
  case ('dbeta_dpz_a');     fixer%a%dbeta_dpz                 = fixer%value(dbeta_dpz_a_stored$)
  case ('dalpha_dpz_a');    fixer%a%dalpha_dpz                = fixer%value(dalpha_dpz_a_stored$)

  case ('beta_b');          fixer%b%beta                      = fixer%value(beta_b_stored$)
  case ('alpha_b');         fixer%b%alpha                     = fixer%value(alpha_b_stored$)
  case ('phi_b');           fixer%b%phi                       = fixer%value(phi_b_stored$)
  case ('dbeta_dpz_b');     fixer%b%dbeta_dpz                 = fixer%value(dbeta_dpz_b_stored$)
  case ('dalpha_dpz_b');    fixer%b%dalpha_dpz                = fixer%value(dalpha_dpz_b_stored$)

  case ('eta_x');           fixer%x%eta                       = fixer%value(eta_x_stored$)
  case ('etap_x');          fixer%x%etap                      = fixer%value(etap_x_stored$)
  case ('deta_dpz_x');      fixer%x%deta_dpz                  = fixer%value(deta_dpz_x_stored$)
  case ('detap_dpz_x');     fixer%x%detap_dpz                 = fixer%value(detap_dpz_x_stored$)

  case ('eta_y');           fixer%y%eta                       = fixer%value(eta_y_stored$)
  case ('etap_y');          fixer%y%etap                      = fixer%value(etap_y_stored$)
  case ('deta_dpz_y');      fixer%y%deta_dpz                  = fixer%value(deta_dpz_y_stored$)
  case ('detap_dpz_y');     fixer%y%detap_dpz                 = fixer%value(detap_dpz_y_stored$)

  case ('mode_flip');       fixer%mode_flip                   = is_true(fixer%value(mode_flip_stored$))
  case ('cmat_11');         fixer%c_mat(1,1)                  = fixer%value(cmat_11_stored$)
  case ('cmat_12');         fixer%c_mat(1,2)                  = fixer%value(cmat_12_stored$)
  case ('cmat_21');         fixer%c_mat(2,1)                  = fixer%value(cmat_21_stored$)
  case ('cmat_22');         fixer%c_mat(2,2)                  = fixer%value(cmat_22_stored$)

  case ('dcmat_dpz_11');    fixer%dc_mat_dpz(1,1)             = fixer%value(dcmat_dpz_11_stored$)
  case ('dcmat_dpz_12');    fixer%dc_mat_dpz(1,2)             = fixer%value(dcmat_dpz_12_stored$)
  case ('dcmat_dpz_21');    fixer%dc_mat_dpz(2,1)             = fixer%value(dcmat_dpz_21_stored$)
  case ('dcmat_dpz_22');    fixer%dc_mat_dpz(2,2)             = fixer%value(dcmat_dpz_22_stored$)

  case default
    is_ok = .false.
    call out_io (s_error$, r_name, 'Fixer element parameter name not recognized: ' // whom)
  end select
endif

if (len(whom) > 4) then
  if (whom(1:4) == 'cmat') then
    gc2 = 1 - fixer%c_mat(1,1)*fixer%c_mat(2,2) + fixer%c_mat(1,2)*fixer%c_mat(2,1)
    if (gc2 > 0) fixer%gamma_c = sqrt(gc2)
  endif
endif

end subroutine fix_this1

end function transfer_fixer_params

end module
