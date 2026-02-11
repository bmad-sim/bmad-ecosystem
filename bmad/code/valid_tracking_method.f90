!+
! Function valid_tracking_method (ele, species, tracking_method) result (is_valid)
!
! Routine to return whether a given tracking method is valid for a given element.
!
! Input:
!   ele             -- ele_struct: Lattice element.
!   species         -- Type of particle being tracked. electron$, etc. or not_set$
!   tracking_method -- integer: bmad_standard$, etc.
!
! Output:
!   is_valid  -- logical: True if a valid method. False otherwise.
!-

function valid_tracking_method (ele, species, tracking_method) result (is_valid)

use bmad_interface, dummy => valid_tracking_method

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: field_ele
integer tracking_method, species, method
logical is_valid

!

is_valid = .false.

! If tracking photons...

if (species == photon$) then
  select case (tracking_method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

  return
endif

! Save some writting since fixed_step_* versions are valid when non fixed_step_* versions are valid.
! Exception: patch elements.

method = tracking_method
if (ele%key /= patch$) then
  if (method == fixed_step_runge_kutta$) method = runge_kutta$
  if (method == fixed_step_time_runge_kutta$) method = time_runge_kutta$
endif

!

select case (ele%key)

case (ab_multipole$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, linear$, taylor$, custom$)
    is_valid = .true.
  end select

case (ac_kicker$)
  select case (method)
  case (bmad_standard$, runge_kutta$, time_runge_kutta$, linear$, custom$, symp_lie_ptc$, taylor$)
    is_valid = (method /= symp_lie_ptc$ .and. method /= taylor$)
    if (.not. is_valid .and. allocated(ele%ac_kick%frequency)) then
      if (size(ele%ac_kick%frequency) == 1) is_valid = .true.
    endif
  end select

case (beambeam$)
  select case (method)
  case (bmad_standard$, linear$, custom$, symp_lie_ptc$, taylor$)
    is_valid = .true.
  end select

case (beginning_ele$)
  select case (method)
  case (bmad_standard$)
    is_valid = .true.
  end select

case (converter$)
  select case (method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (crab_cavity$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, linear$, custom$, runge_kutta$, time_runge_kutta$)
    is_valid = .true.
  end select

case (crystal$, mirror$, multilayer_mirror$, capillary$)
  if (species == not_set$) then
    select case (method)
    case (bmad_standard$, custom$)
      is_valid = .true.
    end select
  endif

case (custom$)
  select case (method)
  case (custom$)
    is_valid = .true.
  end select

case (diffraction_plate$, mask$)
  select case (method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (drift$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (e_gun$)
  select case (method)
  case (runge_kutta$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (ecollimator$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (elseparator$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (em_field$)
  select case (method)
  case (symp_lie_ptc$, taylor$, runge_kutta$, time_runge_kutta$, custom$)
    is_valid = .true.
  case (symp_lie_bmad$)
    if (associated(ele%cartesian_map) .and. ele%field_calc == fieldmap$) is_valid = .true.
  end select

case (floor_shift$, fiducial$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, custom$)
    is_valid = .true.
  end select

case (fork$, photon_fork$)
  select case (method)
  case (bmad_standard$, linear$, custom$)
    is_valid = .true.
  end select

case (gkicker$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, linear$, taylor$, mad$, custom$)
    is_valid = .true.
  end select

case (group$, overlay$, girder$, ramper$)
  ! No valid methods

case (hkicker$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (hybrid$)
  select case (method)
  case (linear$, taylor$)
    is_valid = .true.
  end select

case (instrument$, pipe$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (kicker$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (lcavity$)
  select case (method)
  case (symp_lie_ptc$, runge_kutta$, time_runge_kutta$)
    if (ele%value(l_active$) /= 0) is_valid = .true.
  case (bmad_standard$, linear$, taylor$, custom$)
    is_valid = .true.
  end select

case (marker$, detector$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, linear$, taylor$, custom$)
    is_valid = .true.
  end select

case (fixer$)
  select case (method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (match$)
  select case (method)
  case (bmad_standard$, taylor$, custom$)
    is_valid = .true.
  end select

case (monitor$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (multipole$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, linear$, taylor$, custom$)
    is_valid = .true.
  end select

case (null_ele$)
  ! Nothing to do

case (octupole$, thick_multipole$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (patch$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, custom$, runge_kutta$, time_runge_kutta$)
    is_valid = .true.
    if ((method == runge_kutta$ .or. method == time_runge_kutta$) .and. &
                      cos(ele%value(x_pitch$)) * cos(ele%value(y_pitch$)) < 0) is_valid = .false.
  end select

case (photon_init$)
  select case (method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (pickup$)
  select case (method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (quadrupole$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, symp_lie_bmad$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (rcollimator$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (rf_bend$)
  select case (method)
  case (runge_kutta$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (rfcavity$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (sad_mult$)
  select case (method)
  case (bmad_standard$, custom$, symp_lie_ptc$, linear$, taylor$, runge_kutta$)
    is_valid = .true.
  end select

case (sample$)
  select case (method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (sbend$, rbend$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, linear$, taylor$, mad$, custom$, runge_kutta$, time_runge_kutta$)
    is_valid = .true.
  end select

case (sextupole$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (solenoid$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, symp_lie_bmad$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (sol_quad$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, symp_lie_bmad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (foil$)
  select case (method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (taylor$)
  select case (method)
  case (taylor$, linear$, custom$, symp_lie_ptc$)
    is_valid = .true.
  end select

case (vkicker$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (wiggler$, undulator$)
  ! %field_calc = int_garbage during parsing. Must accept any possible tracking_method in this case.
  field_ele => pointer_to_field_ele(ele, 1)
  select case (field_ele%field_calc)
  case (int_garbage$)
    select case (method)
    case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, symp_lie_bmad$, time_runge_kutta$, custom$)
      is_valid = .true.
    end select
  case (fieldmap$)      ! Is map type
    select case (method)
    case (symp_lie_ptc$, runge_kutta$, linear$, taylor$, symp_lie_bmad$, time_runge_kutta$, custom$)
      is_valid = .true.
    end select
  case (planar_model$)  ! Is periodic type
    select case (method)
    case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, taylor$, symp_lie_bmad$, time_runge_kutta$, custom$)
      is_valid = .true.
    end select
  case (helical_model$)  ! Is periodic type
    select case (method)
    case (bmad_standard$, runge_kutta$, linear$, symp_lie_bmad$, time_runge_kutta$, custom$)
      is_valid = .true.
    end select
  end select

case default
  call err_exit   ! Should not be here

end select

end function valid_tracking_method 

