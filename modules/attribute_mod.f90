module attribute_mod

use basic_attribute_mod
use multipole_mod
use lat_ele_loc_mod

private check_this_attribute_free

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free
!
! Overloaded function for:
!   Function attribute_free1 (ix_ele, attrib_name, lat,
!                                err_print_flag, except_overlay, ignore_field_master) result (free)
!   Function attribute_free2 (ele, attrib_name, 
!                                err_print_flag, except_overlay, ignore_field_master) result (free)
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
!   ix_ele                  -- integer: Index of element in element array.
!   ix_branch               -- integer: Branch index of element. 
!   ele                     -- ele_struct: Element containing the attribute
!   attrib_name             -- character(*): Name of the attribute. Assumed upper case.
!   lat                     -- lat_struct: Lattice structure.
!   err_print_flag          -- logical, optional: If present and False then suppress
!                               printing of an error message if attribute is not free.
!   except_overlay          -- logical, optional: If present and True then an attribute that
!                               is controlled by an overlay will be treated as free. 
!                               This is used by, for example, the create_overlay routine.
!   ignore_fielde_master    -- logical, optional: If present and True then do not mark as not free 
!                               attributes that are slaved due to the setting of ele%field_master.
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
! Function valid_tracking_method (ele, species, tracking_method) result (is_valid)
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
!   is_valid  -- logical: True if a valid method. False otherwise.
!-

function valid_tracking_method (ele, species, tracking_method) result (is_valid)

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord, field_ele
integer tracking_method, species, method
logical is_valid

!

is_valid = .false.

! If tracking photons...

if (species == photon$) then

  select case (ele%key)
  case (crystal$, mirror$, multilayer_mirror$, drift$, fork$, photon_fork$, capillary$)
    select case (tracking_method)
    case (bmad_standard$, custom$)
      is_valid = .true.
    end select

  case default
    ! Nothing is valid
  end select

  return

endif

! Save some writting since fixed_step_* versions are valid when non fixed_step_* versions are valid.

method = tracking_method
if (method == fixed_step_runge_kutta$) method = runge_kutta$
if (method == fixed_step_time_runge_kutta$) method = time_runge_kutta$

select case (ele%key)

case (ab_multipole$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, linear$, symp_map$, taylor$, custom$)
    is_valid = .true.
  end select

case (ac_kicker$)
  select case (method)
  case (bmad_standard$, runge_kutta$, time_runge_kutta$, boris$, linear$, custom$)
    is_valid = .true.
  end select
  

case (beambeam$)
  select case (method)
  case (bmad_standard$, linear$, custom$)
    is_valid = .true.
  end select

case (bend_sol_quad$)
  select case (method)
  case (symp_lie_bmad$, custom$)
    is_valid = .true.
  end select

case (crystal$, mirror$, multilayer_mirror$, capillary$, fiducial$)
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
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (e_gun$)
  select case (method)
  case (time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (ecollimator$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (elseparator$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (em_field$)
  select case (method)
  case (runge_kutta$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (floor_shift$)
  select case (method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (fork$, photon_fork$)
  select case (method)
  case (bmad_standard$, linear$, custom$)
    is_valid = .true.
  end select
  
case (hkicker$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (hybrid$)
  select case (method)
  case (linear$, taylor$)
    is_valid = .true.
  end select

case (instrument$, pipe$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (kicker$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (lcavity$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (marker$, detector$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, linear$, symp_map$, taylor$, custom$)
    is_valid = .true.
  end select

case (match$)
  select case (method)
  case (bmad_standard$, taylor$, custom$)
    is_valid = .true.
  end select

case (monitor$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (multipole$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, linear$, symp_map$, taylor$, custom$)
    is_valid = .true.
  end select

case (octupole$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (patch$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, taylor$, custom$, runge_kutta$, time_runge_kutta$)
    is_valid = .true.
  end select

case (photon_init$)
  select case (method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (quadrupole$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, symp_lie_bmad$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (rcollimator$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (rfcavity$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (sad_mult$)
  select case (method)
  case (bmad_standard$, custom$, symp_lie_ptc$, linear$, symp_map$, taylor$)
    is_valid = .true.
  end select

case (sample$)
  select case (method)
  case (bmad_standard$, custom$)
    is_valid = .true.
  end select

case (sbend$, rbend$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, linear$, symp_map$, taylor$, mad$, custom$, runge_kutta$, time_runge_kutta$)
    is_valid = .true.
  end select

case (sextupole$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (solenoid$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, symp_lie_bmad$, boris$, mad$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (sol_quad$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, symp_lie_bmad$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (taylor$)
  select case (method)
  case (taylor$, linear$, custom$, symp_lie_ptc$)
    is_valid = .true.
  end select

case (vkicker$)
  select case (method)
  case (bmad_standard$, symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, boris$, time_runge_kutta$, custom$)
    is_valid = .true.
  end select

case (wiggler$, undulator$)
  if (ele%sub_key == map_type$) then
    select case (method)
    case (symp_lie_ptc$, runge_kutta$, linear$, symp_map$, taylor$, symp_lie_bmad$, boris$, time_runge_kutta$, custom$)
      is_valid = .true.
    end select
  elseif (ele%sub_key == periodic_type$) then
    select case (method)
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
! Function valid_mat6_calc_method (ele, species, mat6_calc_method) result (is_valid)
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
!   is_valid  -- logical: True if a valid method. False otherwise.
!-

function valid_mat6_calc_method (ele, species, mat6_calc_method) result (is_valid)

implicit none

type (ele_struct) ele
integer mat6_calc_method, species
logical is_valid

!

is_valid = .false.

! If tracking photons...

if (species == photon$) then
  select case (ele%key)
  case (crystal$, mirror$, multilayer_mirror$, drift$, fork$, photon_fork$, capillary$)
    select case (mat6_calc_method)
    case (bmad_standard$, static$, tracking$, custom$)
      is_valid = .true.
    end select

  case default
    ! Nothing is valid
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

case (bend_sol_quad$)
  select case (mat6_calc_method)
  case (symp_lie_bmad$, static$, tracking$, custom$)
    is_valid = .true.
  end select

case (fork$, photon_fork$)
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

case (custom$)
  select case (mat6_calc_method)
  case (static$, tracking$, custom$)
    is_valid = .true.
  end select

case (diffraction_plate$, mask$, fiducial$, floor_shift$)
  select case (mat6_calc_method)
  case (bmad_standard$, static$, tracking$, custom$)
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
  case (static$, tracking$, custom$)
    is_valid = .true.
  end select

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
  if (ele%sub_key == map_type$) then
    select case (mat6_calc_method)
    case (symp_lie_ptc$, taylor$, symp_lie_bmad$, static$, tracking$, custom$)
      is_valid = .true.
    end select
  elseif (ele%sub_key == periodic_type$) then
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
! Function valid_spin_tracking_method (ele, spin_tracking_method) result (is_valid)
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
!   is_valid  -- logical: True if a valid method. False otherwise.
!-

function valid_spin_tracking_method (ele, spin_tracking_method) result (is_valid)

implicit none

type (ele_struct) ele
integer spin_tracking_method
logical is_valid

! 

select case (ele%key)

case (ab_multipole$)
  select case (spin_tracking_method)
  case (tracking$, custom$)
    is_valid = .true.
  end select

case (ac_kicker$)
  select case (spin_tracking_method)
  case (tracking$, custom$)
    is_valid = .true.
  end select

case (capillary$, crystal$, mirror$, multilayer_mirror$, taylor$)
  is_valid = .false.

case (custom$)
  select case (spin_tracking_method)
  case (tracking$, custom$)
    is_valid = .true.
  end select

case (hybrid$)
  is_valid = .false.

case (sad_mult$, patch$)
  select case (spin_tracking_method)
  case (custom$, symp_lie_ptc$)
    is_valid = .true.
  end select

case default
  select case (spin_tracking_method)
  case (custom$, symp_lie_ptc$, tracking$)
    is_valid = .true.
  end select
end select

end function valid_spin_tracking_method 

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine pointer_to_indexed_attribute (ele, ix_attrib, do_allocation, a_ptr, err_flag, err_print_flag)
!
! DEPRECATED ROUTINE! DO NOT USE!
! Consider instead: pointer_to_attribute.
!
! Returns a pointer to an attribute of an element ele with attribute index ix_attrib.
!
! Use of this routine is restricted to attributes that have an index like k1$, tracking_method$, etc.
!
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
integer i, j, n, ix, iy, expn(6)

character(40) a_name
character(10) str
character(*), parameter :: r_name = 'pointer_to_indexed_attribute'

logical err_flag, do_allocation, do_print
logical, optional :: err_print_flag

! Init

err_flag = .true.
nullify (a_ptr%r, a_ptr%i, a_ptr%l)
do_print = logic_option (.true., err_print_flag)

! Taylor

if (ele%key == taylor$ .and. ix_attrib > taylor_offset$) then
  write (str, '(i0)') ix_attrib - taylor_offset$
  n = index('123456', str(1:1))
  if (.not. associated(ele%taylor(1)%term)) then
    if (.not. do_allocation) return
    do i = 1, 6
      call init_taylor_series(ele%taylor(i), 0)
    enddo
  endif

  expn = 0
  do i = 2, len_trim(str)
    j = index('123456', str(i:i))
    expn(j) = expn(j) + 1
  enddo

  i = taylor_term_index(ele%taylor(n), expn, do_allocation)
  if (i /= 0) then
    a_ptr%r => ele%taylor(n)%term(i)%coef
    err_flag = .false.
  endif
  return
endif

! Overlay or Group

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
        call multipole_init (ele, magnetic$)
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
      call multipole_init (ele, electric$)
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
! Function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag, except_overlay, ignore_field_master) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free1 (ix_ele, attrib_name, lat, err_print_flag, except_overlay, ignore_field_master) result (free)

implicit none

type (lat_struct) :: lat

integer ix_ele

character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay, ignore_field_master

!

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)

call check_this_attribute_free (lat%ele(ix_ele), attrib_name, lat, do_print, do_except_overlay, ignore_field_master, free)

end function attribute_free1

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free2 (ele, attrib_name, err_print_flag, except_overlay, ignore_field_master) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free2 (ele, attrib_name, err_print_flag, except_overlay, ignore_field_master) result (free)

implicit none

type (lat_struct), target :: lat
type (ele_struct) ele

character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay, ignore_field_master

character(16) :: r_name = 'attribute_free'

! Elements not assocaited with a lattice are considered free.

if (.not. associated(ele%branch)) then
  free = .true.
  return
endif

! init & check

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)

call check_this_attribute_free (ele, attrib_name, ele%branch%lat, do_print, do_except_overlay, ignore_field_master, free)

end function attribute_free2

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, 
!                                         err_print_flag, except_overlay, ignore_field_master) result (free)
!
! This function overloaded by attribute_free. See attribute_free for more details.
!-

function attribute_free3 (ix_ele, ix_branch, attrib_name, lat, &
                                            err_print_flag, except_overlay, ignore_field_master) result (free)

implicit none

type (lat_struct) :: lat

integer ix_ele, ix_branch
character(*) attrib_name

logical free, do_print, do_except_overlay
logical, optional :: err_print_flag, except_overlay, ignore_field_master

!

do_print = logic_option (.true., err_print_flag)
do_except_overlay = logic_option(.false., except_overlay)

call check_this_attribute_free (lat%branch(ix_branch)%ele(ix_ele), attrib_name, &
                                             lat, do_print, do_except_overlay, ignore_field_master, free)

end function attribute_free3

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine check_this_attribute_free (ele, attrib_name, lat, do_print, do_except_overlay, &
                                                                      ignore_field_master, free, ix_lord)

implicit none

type (ele_struct), target :: ele
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele_p, lord
type (branch_struct), pointer :: branch
type (ele_attribute_struct) attrib_info
type (control_struct), pointer :: control
type (all_pointer_struct) a_ptr

integer ix_branch, i, ir, ix_attrib, ix, ic
integer, optional :: ix_lord

character(*) attrib_name
character(40) a_name

logical free, do_print, do_except_overlay, ignore_field_master, err_flag

!

free = .true.

branch => lat%branch(ele%ix_branch)

! Init & check that the name corresponds to an attribute

ix_attrib = attribute_index(ele, attrib_name)
attrib_info = attribute_info(ele, ix_attrib)

a_name = attribute_name (ele, ix_attrib)

if (attrib_info%type == private$) then
  call it_is_not_free (ele, ix_attrib, 'THIS ATTRIBUTE IS PRIVATE.')
  return
endif

if (attrib_info%type == does_not_exist$) then
  ! Calculated quantities like the Twiss function do not have an entry in the attribute table but
  ! pointer_to_attribute will return a pointer.
  ! Note: Something like beginning element beta_a do have an entry in the attribute table. 
  select case (attrib_name)
  case ('ALPHA_A', 'ALPHA_B', 'BETA_A', 'BETA_B', 'PHI_A', 'PHI_B', 'DPHI_A', 'DPHI_B', &
        'ETA_A', 'ETAP_A', 'ETA_B', 'ETAP_B')
    call it_is_not_free (ele, ix_attrib, 'THIS ATTRIBUTE IS A COMPUTED PARAMETER.')
  case default
    ! Something like 'cartesian_map(1)%field_scale' does not have an attribute index
    call pointer_to_attribute (ele, attrib_name, .true., a_ptr, err_flag, .false.)
    if (.not. err_flag) return
    call it_is_not_free (ele, ix_attrib, 'THIS NAME DOES NOT CORRESPOND TO A VALID ATTRIBUTE.')
  end select
  return
endif

! only one particular attribute of an overlay lord is allowed to be adjusted

if (ele%key == overlay$ .or. ele%key == group$) then
  if (all(attrib_name /= ele%control_var%name)) then
    call it_is_not_free (ele, ix_attrib, 'IT IS NOT A VALID CONTROL VARIABLE')
  endif
  return
endif

! Here if checking something that is not an overlay or group lord... 

if (attrib_info%type == dependent$) then
  call it_is_not_free (ele, ix_attrib, 'THIS IS A DEPENDENT ATTRIBUTE.')
  return
endif

! If the attribute is controled by an overlay lord then it cannot be varied.
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
        call it_is_not_free (ele, ix_attrib, 'IT IS CONTROLLED BY THE OVERLAY: ' // lord%name)
        return
      endif
    endif

  enddo
endif

! With a few exceptions, super_slaves attributes cannot be varied.

if (ele%slave_status == super_slave$) then
  select case (a_name)
  case ('L', 'CSR_CALC_ON')
  case default
    call it_is_not_free (ele, ix_attrib, 'THIS ELEMENT IS A SUPER_SLAVE.', &
                                         '[ATTRIBUTES OF SUPER_SLAVE ELEMENTS ARE GENERALLY NOT FREE TO VARY.]')
  end select
  return
endif

! Check for a multipass_slave.
! Exception: phi0_multipass can be varied for lcavity and rfcavity slaves, etc.

if (ele%slave_status == multipass_slave$) then
  if (a_name == 'CSR_CALC_ON') return
  select case (ele%key)
  case (lcavity$, rfcavity$) 
    if (ix_attrib == phi0_multipass$) return
  case (patch$)
    lord => pointer_to_lord(ele, 1)
  end select

  call it_is_not_free (ele, ix_attrib, 'THIS ELEMENT IS A MULTIPASS_SLAVE.', &
                                       '[ATTRIBUTES OF MULTIPASS_SLAVE ELEMENTS ARE GENERALLY NOT FREE TO VARY.]')

  return
endif

! 

select case (a_name)
case ('NUM_STEPS')
  call it_is_not_free (ele, ix_attrib, 'THIS ATTRIBUTE CANNOT BE VARIED.')
  return

case ('FIELD_SCALE', 'PHI0_FIELDMAP', 'CSR_CALC_ON')
  return

case ('E_TOT', 'P0C')
  if (ele%key == beginning_ele$) return

  if (.not. logic_option(.false., ignore_field_master) .and. ele%lord_status == multipass_lord$ .and. &
                                                     .not. ele%field_master .and. ele%value(n_ref_pass$) == 0) then
    select case (ele%key)
    case (quadrupole$, sextupole$, octupole$, solenoid$, sol_quad$, sbend$, &
          hkicker$, vkicker$, kicker$, ac_kicker$, elseparator$, bend_sol_quad$, lcavity$, rfcavity$)
      return  ! Is free
    end select
  endif

  call it_is_not_free (ele, ix_attrib, 'THIS IS A DEPENDENT ATTRIBUTE.')
  return
end select

! check if it is a dependent variable.

if (attrib_info%type == is_free$) return

select case (ele%key)
case (sbend$)
  if (any(ix_attrib == [angle$, l_chord$, rho$])) free = .false.
case (rfcavity$)
  if (.not. logic_option(.false., ignore_field_master)) then
    if (ix_attrib == rf_frequency$ .and. ele%field_master) free = .false.
    if (ix_attrib == harmon$ .and. .not. ele%field_master) free = .false.
  endif
  if (ix_attrib == gradient$) free = .false.
case (lcavity$)
  if (ix_attrib == voltage$ .and. ele%value(l$) /= 0) free = .false.
  if (ix_attrib == gradient$ .and. ele%value(l$) == 0) free = .false.
case (elseparator$)
  if (ix_attrib == voltage$) free = .false.
end select

if (ele%key == sbend$ .and. ele%lord_status == multipass_lord$ .and. &
    ele%value(n_ref_pass$) == 0 .and. ix_attrib == p0c$) free = .true.

if (.not. free) then
  call it_is_not_free (ele, ix_attrib, 'THIS IS A DEPENDENT ATTRIBUTE.')
  return
endif

! field_master on means that the b_field and bn_gradient values control the strength.

if (.not. logic_option(.false., ignore_field_master)) then
  free = field_attribute_free (ele, a_name)
  if (.not. free) then
    call it_is_not_free (ele, ix_attrib, &
         "THIS IS A DEPENDENT ATTRIBUTE SINCE", &
         "THE ELEMENT'S FIELD_MASTER IS SET TO: " // on_off_logic (ele%field_master))
    return
  endif
endif

!-------------------------------------------------------
contains

subroutine it_is_not_free (ele, ix_attrib, l1, l2)

type (ele_struct) ele

integer ix_attrib, nl

character(*) l1
character(*), optional :: l2
character(100) li(8)
character(*), parameter :: r_name = 'attribute_free'

!

free = .false.

if (.not. do_print) return

nl = 0

nl=nl+1; li(nl) =   'THE ATTRIBUTE: ' // attrib_name
nl=nl+1; li(nl) =   'OF THE ELEMENT: ' // ele%name
nl=nl+1; li(nl) = 'IS NOT FREE TO VARY SINCE:'

nl=nl+1; li(nl) = l1
if (present(l2)) then
  nl=nl+1; li(nl) = l2
endif

call out_io (s_error$, r_name, li(1:nl))   

end subroutine it_is_not_free

end subroutine check_this_attribute_free

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function field_attribute_free (ele, attrib_name) result (free)
!
! Routine to check if a field attribute is free to vary.
!
! Field attributes are either normalized (EG K2 of a sextupole) or unnormalized (EG B2_GRADIENT of a sextupole).
! Whether normalized or unnormalized attributes are free to vary will depend on the setting  of ele%field_master.
!
! Generally, this routine should not be called directly. Use the routine attribute_free instead.
!
! Input:
!   ele           -- ele_struct: Element containing the attribute
!   attrib_name   -- character(*): Name of the field attribute. Assumed upper case.
!
! Output:
!   free          -- logical: Is the attribute free to vary? 
!                     If the attribute is not recognized, free = True will be returned.
!-

function field_attribute_free (ele, attrib_name) result (free)

implicit none

type (ele_struct) ele
character(*) attrib_name
logical free

!

free = .true.

if (ele%field_master) then
  select case (ele%key)
  case (quadrupole$)
    if (attrib_name == 'K1') free = .false.
  case (sextupole$)
    if (attrib_name == 'K2') free = .false.
  case (octupole$)
    if (attrib_name == 'K3') free = .false.
  case (solenoid$)
    if (attrib_name == 'KS') free = .false.
  case (sol_quad$)
    if (attrib_name == 'KS') free = .false.
    if (attrib_name == 'K1') free = .false.
    if (attrib_name == 'K2') free = .false.
  case (sbend$)
    if (attrib_name == 'G') free = .false.
    if (attrib_name == 'G_ERR') free = .false.
  case (hkicker$, vkicker$)
    if (attrib_name == 'KICK') free = .false.
  end select

  if (has_hkick_attributes(ele%key)) then
    if (attrib_name == 'HKICK') free = .false.
    if (attrib_name == 'VKICK') free = .false.
  endif

else
  select case (ele%key)
  case (elseparator$)
    if (attrib_name == 'E_FIELD') free = .false.
  case (quadrupole$)
    if (attrib_name == 'B1_GRADIENT') free = .false.
  case (sextupole$)
    if (attrib_name == 'B2_GRADIENT') free = .false.
  case (octupole$)
    if (attrib_name == 'B3_GRADIENT') free = .false.
  case (solenoid$)
    if (attrib_name == 'BS_FIELD') free = .false.
  case (sol_quad$)
    if (attrib_name == 'BS_FIELD') free = .false.
    if (attrib_name == 'B1_GRADIENT') free = .false.
  case (sbend$)
    if (attrib_name == 'B_FIELD') free = .false.
    if (attrib_name == 'B_FIELD_ERR') free = .false.
    if (attrib_name == 'B1_GRADIENT') free = .false.
    if (attrib_name == 'B2_GRADIENT') free = .false.
  case (hkicker$, vkicker$)
    if (attrib_name == 'BL_KICK') free = .false.
  end select

  if (has_hkick_attributes(ele%key)) then
    if (attrib_name == 'BL_HKICK') free = .false.
    if (attrib_name == 'BL_VKICK') free = .false.
  endif
endif

end function field_attribute_free

end module
