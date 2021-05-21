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

use bmad_struct

implicit none

type (ele_struct) ele
integer mat6_calc_method, species
logical is_valid

!

is_valid = .false.

! If tracking photons...

if (species == photon$) then
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

  return
endif

! Save some writting since fixed_step_* versions are valid when non fixed_step_* versions are valid.

select case (ele%key)

case (ab_multipole$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (ac_kicker$)
  select case (mat6_calc_method)
  case (static$, tracking$, custom$)
    is_valid = .true.
  end select

case (beambeam$)
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (beginning_ele$)
  select case (mat6_calc_method)
  case (bmad_standard$)
    is_valid = .true.
  end select

case (crab_cavity$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, static$, tracking$, custom$, taylor$)
    is_valid = .true.
  end select

case (fork$, photon_fork$)
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (crystal$, diffraction_plate$, mirror$, multilayer_mirror$, capillary$)
  if (species == not_set$) then
    select case (mat6_calc_method)
    case (bmad_standard$, static$, tracking$, custom$)
      is_valid = .true.
    end select
  endif

case (custom$)
  select case (mat6_calc_method)
  case (static$, tracking$, custom$)
    is_valid = .true.
  end select

case (drift$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, mad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (e_gun$)
  select case (mat6_calc_method)
  case (static$, tracking$, custom$)
    is_valid = .true.
  end select

case (ecollimator$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (elseparator$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, mad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (em_field$)
  select case (mat6_calc_method)
  case (symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (fiducial$, floor_shift$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (group$, overlay$, girder$, ramper$)
  ! No valid methods

case (hkicker$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (hybrid$)
  select case (mat6_calc_method)
  case (taylor$, static$)
    is_valid = .true.
  end select

case (instrument$, pipe$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (kicker$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (lcavity$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (mask$)
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (marker$, detector$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (match$)
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (monitor$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (multipole$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (octupole$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (patch$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (photon_init$)
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (quadrupole$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, mad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (rcollimator$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (sad_mult$)
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$, symp_lie_ptc$, taylor$)
    is_valid = .true.
  end select

case (sample$)
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (sbend$, rbend$, sextupole$, rfcavity$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, mad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (solenoid$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, mad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (sol_quad$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (taylor$)
  select case (mat6_calc_method)
  case (taylor$, static$, custom$, symp_lie_ptc$)
    is_valid = .true.
  end select

case (vkicker$)
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (wiggler$, undulator$)
  ! %field_calc = int_garbage during parsing. Must accept any possible mat6_calc_method in this case.
  if (ele%field_calc == fieldmap$ .or. ele%field_calc == int_garbage$) then
    select case (mat6_calc_method)
    case (symp_lie_ptc$, taylor$, symp_lie_bmad$, static$, tracking$, custom$)
      is_valid = .true.
    end select
  elseif (ele%field_calc == planar_model$ .or. ele%field_calc == int_garbage$) then
    select case (mat6_calc_method)
    case (bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, static$, tracking$, custom$)
      is_valid = .true.
    end select
  elseif (ele%field_calc == helical_model$ .or. ele%field_calc == int_garbage$) then
    select case (mat6_calc_method)
    case (bmad_standard$, symp_lie_bmad$, static$, tracking$, custom$)
      is_valid = .true.
    end select
  endif

case default
  call err_exit

end select

end function valid_mat6_calc_method 

