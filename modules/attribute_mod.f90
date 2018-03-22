module attribute_mod

use basic_attribute_mod
use multipole_mod
use lat_ele_loc_mod

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

end module
