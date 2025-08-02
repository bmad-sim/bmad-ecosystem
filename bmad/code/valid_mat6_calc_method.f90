!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function valid_mat6_calc_method (ele, species, mat6_calc_method) result (is_valid)
!
! Routine to return whether a given mat6_calc method is valid for a given element.
!
! Input:
!   ele              -- ele_struct: Lattice element.
!   species          -- Type of particle being tracked. electron$, etc. or not_set$
!   mat6_calc_method -- integer: bmad_standard$, etc.
!
! Output:
!   is_valid  -- logical: True if a valid method. False otherwise.
!-

function valid_mat6_calc_method (ele, species, mat6_calc_method) result (is_valid)

use bmad_interface, dummy => valid_mat6_calc_method

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: field_ele
integer mat6_calc_method, species
logical is_valid

!

is_valid = .false.

! If tracking photons...

if (species == photon$) then
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, tracking$, custom$)
    is_valid = .true.
  end select

  return
endif

! Save some writting since fixed_step_* versions are valid when non fixed_step_* versions are valid.

select case (ele%key)

case (ab_multipole$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (ac_kicker$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, tracking$, custom$)
    is_valid = .true.
  end select

case (beambeam$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, tracking$, custom$)
    is_valid = .true.
  end select

case (beginning_ele$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$)
    is_valid = .true.
  end select

case (converter$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, custom$)
    is_valid = .true.
  end select

case (crab_cavity$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, tracking$, custom$, taylor$)
    is_valid = .true.
  end select

case (fork$, photon_fork$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, tracking$, custom$)
    is_valid = .true.
  end select

case (crystal$, diffraction_plate$, mirror$, multilayer_mirror$, capillary$)
  if (species == not_set$) then
    select case (mat6_calc_method)
    case (auto$, bmad_standard$, tracking$, custom$)
      is_valid = .true.
    end select
  endif

case (custom$)
  select case (mat6_calc_method)
  case (auto$, tracking$, custom$)
    is_valid = .true.
  end select

case (drift$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, mad$, tracking$, custom$)
    is_valid = .true.
  end select

case (e_gun$)
  select case (mat6_calc_method)
  case (auto$, tracking$, custom$)
    is_valid = .true.
  end select

case (ecollimator$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (elseparator$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, mad$, tracking$, custom$)
    is_valid = .true.
  end select

case (em_field$)
  select case (mat6_calc_method)
  case (auto$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (fiducial$, floor_shift$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, tracking$, custom$)
    is_valid = .true.
  end select

case (gkicker$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select


case (group$, overlay$, girder$, ramper$)
  ! No valid methods

case (hkicker$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (hybrid$)
  select case (mat6_calc_method)
  case (auto$, taylor$)
    is_valid = .true.
  end select

case (instrument$, pipe$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (kicker$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (lcavity$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (mask$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, tracking$, custom$)
    is_valid = .true.
  end select

case (marker$, fixer$, detector$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (match$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, tracking$, custom$)
    is_valid = .true.
  end select

case (monitor$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (multipole$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (null_ele$)
  ! Nothing to do

case (octupole$, thick_multipole$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (patch$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, tracking$, custom$)
    is_valid = .true.
  end select

case (pickup$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, tracking$, custom$)
    is_valid = .true.
  end select

case (photon_init$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, tracking$, custom$)
    is_valid = .true.
  end select

case (quadrupole$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, mad$, tracking$, custom$)
    is_valid = .true.
  end select

case (rcollimator$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (rf_bend$)
  select case (mat6_calc_method)
  case (auto$, tracking$, custom$)
    is_valid = .true.
  end select

case (sad_mult$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, tracking$, custom$, symp_lie_ptc$, taylor$)
    is_valid = .true.
  end select

case (sample$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, tracking$, custom$)
    is_valid = .true.
  end select

case (sbend$, rbend$, sextupole$, rfcavity$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, mad$, tracking$, custom$)
    is_valid = .true.
  end select

case (solenoid$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, mad$, tracking$, custom$)
    is_valid = .true.
  end select

case (sol_quad$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, tracking$, custom$)
    is_valid = .true.
  end select

case (foil$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, custom$)
    is_valid = .true.
  end select

case (taylor$)
  select case (mat6_calc_method)
  case (auto$, taylor$, custom$, symp_lie_ptc$)
    is_valid = .true.
  end select

case (vkicker$)
  select case (mat6_calc_method)
  case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, tracking$, custom$)
    is_valid = .true.
  end select

case (wiggler$, undulator$)
  ! %field_calc = int_garbage during parsing. Must accept any possible mat6_calc_method in this case.
  field_ele => pointer_to_field_ele(ele, 1)
  select case (field_ele%field_calc)
  case (int_garbage$)
    select case (mat6_calc_method)
    case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, tracking$, custom$)
      is_valid = .true.
    end select
  case (fieldmap$)       ! Is map type
    select case (mat6_calc_method)
    case (auto$, symp_lie_ptc$, taylor$, symp_lie_bmad$, tracking$, custom$)
      is_valid = .true.
    end select
  case (planar_model$)   ! Is periodic type
    select case (mat6_calc_method)
    case (auto$, bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, tracking$, custom$)
      is_valid = .true.
    end select
  case (helical_model$)  ! Is periodic type
    select case (mat6_calc_method)
    case (auto$, bmad_standard$, symp_lie_bmad$, tracking$, custom$)
      is_valid = .true.
    end select
  end select

case default
  call err_exit

end select

end function valid_mat6_calc_method 

