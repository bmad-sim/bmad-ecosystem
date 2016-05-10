module attribute_mod

use bmad_struct
use basic_bmad_interface
use multipole_mod
use lat_ele_loc_mod

type(ele_struct), private, pointer, save :: ele0 ! For Error message purposes
character(40), private, save :: attrib_name0     ! For Error message purposes

private check_this_attribute_free, print_error

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free
!
! Overloaded function for:
!   Function attribute_free1 (ix_ele, attrib_name, lat,
!                                err_print_flag, except_overlay) result (free)
!   Function attribute_free2 (ele, attrib_name, 
!                                err_print_flag, except_overlay) result (free)
!   Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, 
!                                err_print_flag, except_overlay) result (free)
!
! Routine to check if an attribute is free to vary.
!
! Attributes that cannot be changed directly include super_slave attributes (since
! these attributes are controlled by their super_lords) and attributes that
! are controlled by an overlay.
!
! Also dependent variables such as the angle of a bend cannot be 
!   freely variable.
!
! Modules needed:
!   use bmad
!
! Input:
!   ix_ele          -- Integer: Index of element in element array.
!   ix_branch       -- Integer: Branch index of element. 
!   ele             -- Ele_struct: Element containing the attribute
!   attrib_name     -- Character(*): Name of the attribute. Assumed upper case.
!   lat             -- lat_struct: Lattice structure.
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message if attribute is not free.
!   except_overlay  -- Logical, optional: If present and True then an attribute that
!                       is controlled by an overlay will be treated as free. 
!                       This is used by, for example, the create_overlay routine.
!
! Output:
!   free   -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!-

interface attribute_free
  module procedure attribute_free1
  module procedure attribute_free2
  module procedure attribute_free3
end interface

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function valid_tracking_method (ele, species, tracking_method, num_valid) result (is_valid)
!
! Routine to return whether a given tracking method is valid for a given element.
!
! Module needed:
!   use track1_mod
!
! Input:
!   ele             -- ele_struct: Lattice element.
!   species         -- Type of particle being tracked. electron$, etc. or not_set$
!   tracking_method -- integer: bmad_standard$, etc.
!
! Output:
!   num_valid -- integer, optional: Number of valid methods. 
!   is_valid  -- logical: True if a valid method. False otherwise.
!-

function valid_tracking_method (ele, species, tracking_method, num_valid) result (is_valid)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord, field_ele
integer tracking_method, species
integer, optional :: num_valid
logical is_valid

!

is_valid = .false.
if (present(num_valid)) num_valid = 0

! If tracking photons...

if (species == photon$) then

  select case (ele%key)
  case (crystal$, mirror$, multilayer_mirror$, drift$, fork$, photon_fork$, capillary$)
    if (present(num_valid)) num_valid = 2
    select case (tracking_method)
    case (bmad_standard$, custom$)
      is_valid = .true.
    end select

  case default
    ! Nothing is valid
  end select

  return

endif

!

select case (ele%key)

case (ab_multipole$)
  if (present(num_valid)) num_valid = 6
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, linear$, symp_map$, taylor$, custom$)
    is_valid = .true.
  end select

case (beambeam$)
  if (present(num_valid)) num_valid = 3
  select case (tracking_method)
  case (bmad_standard$, linear$, custom$)
    is_valid = .true.
  end select

case (fork$, photon_fork$)
  if (present(num_valid)) num_valid = 3
  select case (tracking_method)
  case (bmad_standard$, linear$, custom$)
    is_valid = .true.
  end select
  
case (crystal$, mirror$, multilayer_mirror$, capillary$)
  if (species == not_set$) then
    if (present(num_valid)) num_valid = 2
    select case (tracking_method)
    case (bmad_standard$, custom$)
      is_valid = .true.
    end select
  endif

case (diffraction_plate$, mask$)
  if (present(num_valid)) num_valid = 2
  select case (tracking_method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select
  

case (drift$)
  if (present(num_valid)) num_valid = 10
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (e_gun$)
  if (present(num_valid)) num_valid = 2
  select case (tracking_method)
  case (time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (ecollimator$)
  if (present(num_valid)) num_valid = 9
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

! There is no bmad_standard field calculation.

case (elseparator$)
  if (present(num_valid)) num_valid = 10
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (em_field$)
  if (present(num_valid)) num_valid = 4
  select case (tracking_method)
  case (runge_kutta$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (floor_shift$)
  if (present(num_valid)) num_valid = 2
  select case (tracking_method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (hkicker$)
  if (present(num_valid)) num_valid = 9
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (hybrid$)
  if (present(num_valid)) num_valid = 2
  select case (tracking_method)
  case (linear$, taylor$)
    is_valid = .true.
  end select

case (instrument$, pipe$)
  if (present(num_valid)) num_valid = 9
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (kicker$)
  if (present(num_valid)) num_valid = 9
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (lcavity$)
  if (present(num_valid)) num_valid = 9
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (marker$, detector$)
  if (present(num_valid)) num_valid = 6
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, linear$, symp_map$, taylor$, custom$)
    is_valid = .true.
  end select

case (match$)
  if (present(num_valid)) num_valid = 3
  select case (tracking_method)
  case (bmad_standard$, taylor$, custom$)
    is_valid = .true.
  end select

case (monitor$)
  if (present(num_valid)) num_valid = 9
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (multipole$)
  if (present(num_valid)) num_valid = 6
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, linear$, symp_map$, taylor$, custom$)
    is_valid = .true.
  end select

case (octupole$)
  if (present(num_valid)) num_valid = 9
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (patch$)
  if (present(num_valid)) num_valid = 5
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, custom$, runge_kutta$)
    is_valid = .true.
  end select

case (photon_init$)
  if (present(num_valid)) num_valid = 2
  select case (tracking_method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (quadrupole$)
  if (present(num_valid)) num_valid = 11
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, symp_lie_bmad$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (rcollimator$)
  if (present(num_valid)) num_valid = 9
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (rfcavity$)
  if (present(num_valid)) num_valid = 10
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (sad_mult$)
  if (present(num_valid)) num_valid = 6
  select case (tracking_method)
  case (bmad_standard$, custom$, symp_lie_ptc$, linear$, symp_map$, taylor$)
    is_valid = .true.
  end select

case (sample$)
  if (present(num_valid)) num_valid = 2
  select case (tracking_method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (sbend$, rbend$)
  if (present(num_valid)) num_valid = 8
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, linear$, symp_map$, taylor$, mad$, custom$, runge_kutta$, time_runge_kutta$)
    is_valid = .true.
  end select

case (sextupole$)
  if (present(num_valid)) num_valid = 10
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (solenoid$)
  if (present(num_valid)) num_valid = 11
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, symp_lie_bmad$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (sol_quad$)
  if (present(num_valid)) num_valid = 11
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, symp_lie_bmad$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (taylor$)
  if (present(num_valid)) num_valid = 4
  select case (tracking_method)
  case (taylor$, linear$, custom$, symp_lie_ptc$)
    is_valid = .true.
  end select

case (vkicker$)
  if (present(num_valid)) num_valid = 9
  select case (tracking_method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (wiggler$, undulator$)
  if (ele%sub_key == map_type$) then
    if (present(num_valid)) num_valid = 9
    select case (tracking_method)
    case (symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, symp_lie_bmad$, boris$, time_runge_kutta$, custom$)
      is_valid = .true.
    end select
  elseif (ele%sub_key == periodic_type$) then
    if (present(num_valid)) num_valid = 9
    select case (tracking_method)
    case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, symp_lie_bmad$, time_runge_kutta$, custom$)
      is_valid = .true.
    end select
  endif

case default
  call err_exit   ! Should not be here

end select

end function valid_tracking_method 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function valid_mat6_calc_method (ele, species, mat6_calc_method, num_valid) result (is_valid)
!
! Routine to return whether a given mat6_calc method is valid for a given element.
!
! Module needed:
!   use track1_mod
!
! Input:
!   ele              -- ele_struct: Lattice element.
!   species          -- Type of particle being tracked. electron$, etc. or not_set$
!   mat6_calc_method -- integer: bmad_standard$, etc.
!
! Output:
!   num_valid -- integer, optional: Number of valid methods.
!   is_valid  -- logical: True if a valid method. False otherwise.
!-

function valid_mat6_calc_method (ele, species, mat6_calc_method, num_valid) result (is_valid)

implicit none

type (ele_struct) ele
integer mat6_calc_method, species
integer, optional :: num_valid
logical is_valid

!

is_valid = .false.
if (present(num_valid)) num_valid = 0

! If tracking photons...

if (species == photon$) then

  select case (ele%key)
  case (crystal$, mirror$, multilayer_mirror$, drift$, fork$, photon_fork$, capillary$)
    if (present(num_valid)) num_valid = 4
    select case (mat6_calc_method)
    case (bmad_standard$, static$, tracking$, custom$)
      is_valid = .true.
    end select

  case default
    ! Nothing is valid
  end select

  return

endif

!

select case (ele%key)

case (ab_multipole$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (beambeam$)
  if (present(num_valid)) num_valid = 4
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (bend_sol_quad$)
  if (present(num_valid)) num_valid = 4
  select case (mat6_calc_method)
  case (symp_lie_bmad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (fork$, photon_fork$)
  if (present(num_valid)) num_valid = 4
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (crystal$, mirror$, multilayer_mirror$, capillary$)
  if (species == not_set$) then
    select case (mat6_calc_method)
    case (bmad_standard$, static$, tracking$, custom$)
      is_valid = .true.
    end select
  endif

case (diffraction_plate$, mask$)
  if (present(num_valid)) num_valid = 4
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (drift$)
  if (present(num_valid)) num_valid = 7
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, mad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (e_gun$)
  if (present(num_valid)) num_valid = 3
  select case (mat6_calc_method)
  case (static$, tracking$, custom$)
    is_valid = .true.
  end select

case (ecollimator$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (elseparator$)
  if (present(num_valid)) num_valid = 7
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, mad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (em_field$)
  if (present(num_valid)) num_valid = 3
  select case (mat6_calc_method)
  case (static$, tracking$, custom$)
    is_valid = .true.
  end select

case (floor_shift$)
  if (present(num_valid)) num_valid = 4
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (hkicker$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (hybrid$)
  if (present(num_valid)) num_valid = 2
  select case (mat6_calc_method)
  case (taylor$, static$)
    is_valid = .true.
  end select

case (instrument$, pipe$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (kicker$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (lcavity$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (marker$, detector$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (match$)
  if (present(num_valid)) num_valid = 4
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (monitor$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (multipole$)
  if (present(num_valid)) num_valid = 5
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (octupole$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (patch$)
  if (present(num_valid)) num_valid = 5
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (photon_init$)
  if (present(num_valid)) num_valid = 4
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (quadrupole$)
  if (present(num_valid)) num_valid = 8
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, mad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (rcollimator$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (sad_mult$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$, symp_lie_ptc$, taylor$)
    is_valid = .true.
  end select

case (sample$)
  if (present(num_valid)) num_valid = 4
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (sbend$, rbend$, sextupole$, rfcavity$)
  if (present(num_valid)) num_valid = 7
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, mad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (solenoid$)
  if (present(num_valid)) num_valid = 8
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, mad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (sol_quad$)
  if (present(num_valid)) num_valid = 7
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (taylor$)
  if (present(num_valid)) num_valid = 4
  select case (mat6_calc_method)
  case (taylor$, static$, custom$, symp_lie_ptc$)
    is_valid = .true.
  end select

case (vkicker$)
  if (present(num_valid)) num_valid = 6
  select case (mat6_calc_method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (wiggler$, undulator$)
  if (ele%sub_key == map_type$) then
    if (present(num_valid)) num_valid = 6
    select case (mat6_calc_method)
    case (symp_lie_ptc$, taylor$, symp_lie_bmad$, static$, tracking$, custom$)
      is_valid = .true.
    end select
  elseif (ele%sub_key == periodic_type$) then
    if (present(num_valid)) num_valid = 7
    select case (mat6_calc_method)
    case (bmad_standard$, symp_lie_ptc$, taylor$, symp_lie_bmad$, static$, tracking$, custom$)
      is_valid = .true.
    end select
  endif

case default
  call err_exit

end select

end function valid_mat6_calc_method 

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function valid_spin_tracking_method (ele, mat6_calc_method, num_valid) result (is_valid)
!
! Routine to return whether a given spin_tracking method is valid for a given element.
!
! Module needed:
!   use track1_mod
!
! Input:
!   ele              -- ele_struct: Lattice element.
!   spin_tracking_method -- integer: bmad_standard$, etc.
!
! Output:
!   num_valid -- integer, optional: Number of valid methods.
!   is_valid  -- logical: True if a valid method. False otherwise.
!-

function valid_spin_tracking_method (ele, spin_tracking_method, num_valid) result (is_valid)

implicit none

type (ele_struct) ele
integer spin_tracking_method
integer, optional :: num_valid
logical is_valid

!


select case (ele%key)
case (ab_multipole$)
  if (present(num_valid)) num_valid = 2
  select case (spin_tracking_method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (capillary$, crystal$, mirror$, multilayer_mirror$, taylor$)
  if (present(num_valid)) num_valid = 0
  is_valid = .false.

case (custom$)
  if (present(num_valid)) num_valid = 2
  select case (spin_tracking_method)
  case (tracking$, custom$)
    is_valid = .true.
  end select

case (hybrid$)
  if (present(num_valid)) num_valid = 0
  is_valid = .false.

case (sad_mult$, patch$)
  if (present(num_valid)) num_valid = 1
  select case (spin_tracking_method)
  case (custom$, symp_lie_ptc$)
    is_valid = .true.
  end select

case default
  if (present(num_valid)) num_valid = 3
  select case (spin_tracking_method)
  case (bmad_standard$, custom$, symp_lie_ptc$, tracking$)
    is_valid = .true.
  end select
end select

end function valid_spin_tracking_method 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine pointer_to_indexed_attribute (ele, ix_attrib, do_allocation,
!                                                a_ptr, err_flag, err_print_flag)
!
! Returns a pointer to an attribute of an element ele with attribute index ix_attrib.
! 
! Use of this routine is restricted to attributes that have an index like k1$, tracking_method$, etc.
! A more general routine is pointer_to_attribute.
! Alternatively, consider the routine pointers_to_attribute.
! Note: Use attribute_free to see if the attribute may be varied independently.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele             -- Ele_struct: After this routine finishes A_ptr 
!                        will point to a variable within this element.
!   ix_attrib       -- Integer, Attribute index.
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message on error.
!
! Output:
!   a_ptr      -- all_pointer_struct: Pointer to the attribute. 
!     %r           -- pointer to real attribute. Nullified if error or attribute is not real.               
!     %i           -- pointer to integer attribute. Nullified if error or attribute is not integer.
!     %l           -- pointer to logical attribute. Nullified if error or attribute is not logical.               
!   err_flag   -- Logical: Set True if attribtute not found. False otherwise.
!-

subroutine pointer_to_indexed_attribute (ele, ix_attrib, do_allocation, a_ptr, err_flag, err_print_flag)

implicit none

type (ele_struct), target :: ele
type (all_pointer_struct) :: a_ptr

integer :: ix_attrib
integer ix, iy

character(40) a_name
character(*), parameter :: r_name = 'pointer_to_indexed_attribute'

logical err_flag, do_allocation, do_print
logical, optional :: err_print_flag

! Init

err_flag = .true.
nullify (a_ptr%r, a_ptr%i, a_ptr%l)
do_print = logic_option (.true., err_print_flag)

! overlay or group

if (ele%key == overlay$ .or. ele%key == group$) then
  if (is_attribute(ix_attrib, control_var$)) then
    ix = ix_attrib - var_offset$ 
    if (ix > size(ele%control_var)) return
    a_ptr%r => ele%control_var(ix)%value
    err_flag = .false.
    return
  endif

  if (is_attribute(ix_attrib, old_control_var$)) then
    ix = ix_attrib - old_control_var_offset$ 
    if (ix > size(ele%control_var)) return
    a_ptr%r => ele%control_var(ix)%old_value
    err_flag = .false.
    return
  endif
endif

! Multipole or curvature attribute
! Note that elements that have surface curvature automatically have ele%photon allocated.

if (ix_attrib >= a0$ .and. ix_attrib <= b21$) then   
  a_name = attribute_name(ele, ix_attrib)

  ! Curvature
  if (a_name(1:4) == 'CURV') then
    read (a_name(12:12), *) ix
    read (a_name(15:15), *) iy
    if (ix > ubound(ele%photon%surface%curvature_xy, 1) .or. iy > ubound(ele%photon%surface%curvature_xy, 2)) return
    a_ptr%r => ele%photon%surface%curvature_xy(ix,iy)

  ! Multipole
  else
    if (.not. associated(ele%a_pole)) then
      if (do_allocation) then
        call multipole_init (ele)
      else
        if (do_print) call out_io (s_error$, r_name, 'MULTIPOLE NOT ALLOCATED FOR ELEMENT: ' // ele%name)
        return
      endif
    endif

    if (ix_attrib >= b0$) then
      a_ptr%r => ele%b_pole(ix_attrib-b0$)
    else
      a_ptr%r => ele%a_pole(ix_attrib-a0$)
    endif
  endif

! Electric Multipole 

if (ix_attrib >= a0_elec$ .and. ix_attrib <= b21_elec$) then   
  a_name = attribute_name(ele, ix_attrib)

  if (.not. associated(ele%a_pole_elec)) then
    if (do_allocation) then
      call multipole_init (ele, .true.)
    else
      if (do_print) call out_io (s_error$, r_name, 'MULTIPOLE NOT ALLOCATED FOR ELEMENT: ' // ele%name)
      return
    endif
  endif

  if (ix_attrib >= b0_elec$) then
    a_ptr%r => ele%b_pole_elec(ix_attrib-b0_elec$)
  else
    a_ptr%r => ele%a_pole_elec(ix_attrib-a0_elec$)
  endif
endif

! If none of the above

else
  select case (ix_attrib)
  case (is_on$);                          a_ptr%l => ele%is_on
  case (symplectify$);                    a_ptr%l => ele%symplectify
  case (mat6_calc_method$);               a_ptr%i => ele%mat6_calc_method
  case (tracking_method$);                a_ptr%i => ele%tracking_method

  case default
    ! Out of bounds
    if (ix_attrib < 1 .or. ix_attrib > num_ele_attrib$) then
      if (do_print) call out_io (s_error$, r_name, 'INVALID ATTRIBUTE INDEX: \i0\ ', 'FOR THIS ELEMENT: ' // ele%name, &
          i_array = [ix_attrib])
      return
    endif
    a_ptr%r => ele%value(ix_attrib)
  end select
endif

err_flag = .false.

end subroutine pointer_to_indexed_attribute 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag, except_overlay) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag, except_overlay) result (free)

implicit none

type (lat_struct) :: lat

integer ix_ele

character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay

!

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)
free = .false.

call check_this_attribute_free (lat%ele(ix_ele), attrib_name, lat, &
                                             do_print, do_except_overlay, free, 0)

end function attribute_free1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free2 (ele, attrib_name, err_print_flag, except_overlay) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free2 (ele, attrib_name, err_print_flag, except_overlay) result (free)

implicit none

type (lat_struct), target :: lat
type (ele_struct) ele

character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay

character(16) :: r_name = 'attribute_free'

! Elements not assocaited with a lattice are considered free.

if (.not. associated(ele%branch)) then
  free = .true.
  return
endif

! init & check

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)
free = .false.

call check_this_attribute_free (ele, attrib_name, ele%branch%lat, do_print, do_except_overlay, free, 0)

end function attribute_free2

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, 
!                                                 err_print_flag, except_overlay) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, &
                                                  err_print_flag, except_overlay) result (free)

implicit none

type (lat_struct) :: lat

integer ix_ele, ix_branch
character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay

!

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)
free = .false.

call check_this_attribute_free (lat%branch(ix_branch)%ele(ix_ele), attrib_name, &
                                             lat, do_print, do_except_overlay, free, 0)

end function attribute_free3

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

recursive subroutine check_this_attribute_free (ele, attrib_name, lat, &
                            do_print, do_except_overlay, free, ix_recursion, ix_lord)

implicit none

type (ele_struct), target :: ele
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_p, lord
type (branch_struct), pointer :: branch
type (ele_attribute_struct) attrib_info
type (control_struct), pointer :: control

integer ix_branch, ix_recursion, i, ir, ix_attrib, ix, ic
integer, optional :: ix_lord

character(*) attrib_name
character(40) a_name

logical free, do_print, do_except_overlay

! If this is first time then set pointers used in error message printing.

branch => lat%branch(ele%ix_branch)

if (ix_recursion == 0) then
  ele0 => branch%ele(ele%ix_ele)
  attrib_name0 = attrib_name
endif

! Init & check that the name corresponds to an attribute

ix_attrib = attribute_index(ele, attrib_name)
attrib_info = attribute_info(ele, ix_attrib)

a_name = attribute_name (ele, ix_attrib)

if (attrib_info%type == does_not_exist$ .or. attrib_info%type == private$) then
  if (do_print) call print_error (ele, ix_attrib, &
          'THIS NAME DOES NOT CORRESPOND TO A VALID ATTRIBUTE.')
  return
endif

! only one particular attribute of an overlay lord is allowed to be adjusted

if (ele%key == overlay$ .or. ele%key == group$) then
  if (all(attrib_name /= ele%control_var%name)) then
    if (do_print) call print_error (ele, ix_attrib, &
           'FOR THIS OVERLAY ELEMENT THE ATTRIBUTE TO VARY IS: ' // ele%component_name)
  else
    free = .true.
  endif
  return
endif

! Here if checking something that is not an overlay or group lord... 

if (attrib_info%type == dependent$) then
  if (do_print) call print_error (ele, ix_attrib, 'THIS ATTRIBUTE CANNOT BE VARIED.')
  return
endif

! csr_calc_on, etc. are always free.
! x_offset_tot, etc are never free.

select case (a_name)
case ('NUM_STEPS')
  return

case ('FIELD_SCALE', 'PHI0_FIELDMAP')
  free = .true.   ! This may not be true with autoscaling
  return

case ('E_TOT', 'P0C')
  if (ele%key == beginning_ele$) then
    free = .true.
    return
  endif

  if (ele%lord_status /= multipass_lord$) return
  if (ele%field_master) return
  if (ele%value(n_ref_pass$) /= 0) return
  select case (ele%key)
  case (quadrupole$, sextupole$, octupole$, solenoid$, sol_quad$, sbend$, &
        hkicker$, vkicker$, kicker$, elseparator$, bend_sol_quad$, lcavity$, rfcavity$)
    free = .true.
  end select
  return

end select

! if the attribute is controled by an overlay lord then it cannot be varied.
! Exception: Multiple overlays can control the same attribute.

if (.not. do_except_overlay) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i, control)
    if (present(ix_lord)) then
      if (ix_lord == lord%ix_ele) cycle
      if (lat%ele(ix_lord)%key == overlay$) cycle
    endif
    if (lord%key == overlay$) then
      if (control%ix_attrib == ix_attrib) then 
        if (do_print) call print_error (ele, ix_attrib, 'IT IS CONTROLLED BY THE OVERLAY: ' // lord%name)
        return
      endif
    endif
  enddo
endif

! Super_slaves attributes cannot be varied except for L.

if (ele%slave_status == super_slave$ .and. a_name /= 'L') then
  if (do_print) call print_error (ele, ix_attrib, 'THIS ELEMENT IS A SUPER_SLAVE.')
  return
endif

! Check for a multipass_slave.
! Exception: phi0_multipass can be varied for lcavity and rfcavity slaves, etc.

if (ele%slave_status == multipass_slave$) then
  free = .true.
  select case (ele%key)
  case (lcavity$, rfcavity$) 
    if (ix_attrib == phi0_multipass$) return
  case (patch$)
    lord => pointer_to_lord(ele, 1)
  end select

  free = .false.
  if (do_print) call print_error (ele, ix_attrib, 'THIS ELEMENT IS A MULTIPASS_SLAVE.')
  return
endif

! check if it is a dependent variable.

free = .true.
if (attrib_info%type == not_a_child$) return

select case (ele%key)
case (sbend$)
  if (any(ix_attrib == [angle$, l_chord$, rho$])) free = .false.
case (rfcavity$)
  if (ix_attrib == rf_frequency$ .and. ele%field_master) free = .false.
  if (ix_attrib == harmon$ .and. .not. ele%field_master) free = .false.
  if (ix_attrib == gradient$) free = .false.
case (lcavity$)
  if (ix_attrib == voltage$ .and. ele%value(l$) /= 0) free = .false.
  if (ix_attrib == gradient$ .and. ele%value(l$) == 0) free = .false.
case (elseparator$)
  if (ix_attrib == voltage$) free = .false.
end select

if (ele%key == sbend$ .and. ele%lord_status == multipass_lord$ .and. &
    ele%value(n_ref_pass$) == 0 .and. ix_attrib == p0c$) free = .true.

if (.not.free) then
  if (do_print) call print_error (ele, ix_attrib, 'THE ATTRIBUTE IS A DEPENDENT VARIABLE.')
  return
endif

! field_master on means that the b_field and bn_gradient values control
! the strength.

if (ele%field_master) then
  select case (ele%key)
  case (quadrupole$)
    if (ix_attrib == k1$) free = .false.
  case (sextupole$)
    if (ix_attrib == k2$) free = .false.
  case (octupole$)
    if (ix_attrib == k3$) free = .false.
  case (solenoid$)
    if (ix_attrib == ks$) free = .false.
  case (sol_quad$)
    if (ix_attrib == ks$) free = .false.
    if (ix_attrib == k1$) free = .false.
    if (ix_attrib == k2$) free = .false.
  case (sbend$)
    if (ix_attrib == g$) free = .false.
    if (ix_attrib == g_err$) free = .false.
  case (hkicker$, vkicker$)
    if (ix_attrib == kick$) free = .false.
  end select

  if (has_hkick_attributes(ele%key)) then
    if (ix_attrib == hkick$) free = .false.
    if (ix_attrib == vkick$) free = .false.
  endif

else
  select case (ele%key)
  case (elseparator$)
    if (ix_attrib == e_field$) free = .false.
  case (quadrupole$)
    if (ix_attrib == b1_gradient$) free = .false.
  case (sextupole$)
    if (ix_attrib == b2_gradient$) free = .false.
  case (octupole$)
    if (ix_attrib == b3_gradient$) free = .false.
  case (solenoid$)
    if (ix_attrib == bs_field$) free = .false.
  case (sol_quad$)
    if (ix_attrib == bs_field$) free = .false.
    if (ix_attrib == b1_gradient$) free = .false.
  case (sbend$)
    if (ix_attrib == b_field$) free = .false.
    if (ix_attrib == b_field_err$) free = .false.
    if (ix_attrib == b1_gradient$) free = .false.
    if (ix_attrib == b2_gradient$) free = .false.
  case (hkicker$, vkicker$)
    if (ix_attrib == bl_kick$) free = .false.
  end select

  if (has_hkick_attributes(ele%key)) then
    if (ix_attrib == bl_hkick$) free = .false.
    if (ix_attrib == bl_vkick$) free = .false.
  endif

endif

if (.not. free) then
  if (do_print) call print_error (ele, ix_attrib, &
       "THE ATTRIBUTE IS A DEPENDENT VARIABLE SINCE", &
       "THE ELEMENT'S FIELD_MASTER IS " // on_off_logic (ele%field_master))
  return
endif

end subroutine check_this_attribute_free

!-------------------------------------------------------

subroutine print_error (ele, ix_attrib, l1, l2)

type (ele_struct) ele

integer ix_attrib, nl

character(*) l1
character(*), optional :: l2
character(100) li(8)
character (20) :: r_name = 'attribute_free'

!

nl = 0

nl=nl+1; li(nl) =   'THE ATTRIBUTE: ' // attrib_name0
nl=nl+1; li(nl) =   'OF THE ELEMENT: ' // ele0%name

if (ele%ix_branch == 0 .and. ele%ix_ele == ele0%ix_ele) then
  nl=nl+1; li(nl) = 'IS NOT FREE TO VARY SINCE:'
else 
  nl=nl+1; li(nl) = 'IS NOT FREE TO VARY SINCE IT IS TRYING TO CONTROL:'
  nl=nl+1; li(nl) = 'THE ATTRIBUTE: ' // attrib_name0
  nl=nl+1; li(nl) = 'OF THE ELEMENT: ' // ele%name
  nl=nl+1; li(nl) = 'AND THIS IS NOT FREE TO VARY SINCE'
endif

nl=nl+1; li(nl) = l1
if (present(l2)) then
  nl=nl+1; li(nl) = l2
endif

call out_io (s_error$, r_name, li(1:nl))   

end subroutine print_error

end module
