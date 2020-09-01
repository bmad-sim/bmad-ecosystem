!+
! Subroutine tao_evaluate_element_parameters (err, param_name, values, print_err, dflt_source, dflt_component, dflt_uni)
!
! Routine to evaluate a lattice element parameter of the form 
!     <universe>@ele::{<class>}::<ele_name_or_num>[<parameter>]{|<component>}
! or to evaluate at the middle of the element
!     <universe>@ele_mid::{<class>}::<ele_name_or_num>[<parameter>]{|<component>}
! Note: size(values) can be zero without an error
! 
! Input:
!   param_name      -- character(*): parameter name.
!   print_err       -- logical: Print error message? 
!   dflt_source     -- character(*): Default source
!   dflt_component  -- character(*), optional: Default component
!   dflt_uni        -- integer, optional :: Default universe to use.
!
! Output:
!   err       -- Logical: True if there is an error in syntax. False otherwise
!   values(:) -- Real(rp), allocatable: Array of datum valuse.
!-

subroutine tao_evaluate_element_parameters (err, param_name, values, print_err, dflt_source, dflt_component, dflt_uni)

use tao_interface, except_dummy => tao_evaluate_element_parameters

implicit none

type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: ele0
type (ele_struct) ele3
type (tao_lattice_struct), pointer :: tao_lat
type (coord_struct) orb
type (coord_struct), pointer :: orb0
type (branch_struct), pointer :: branch
type (all_pointer_struct) a_ptr
type (tao_lattice_branch_struct), pointer :: tao_branch

character(*) param_name
character(*) dflt_source
character(*), optional :: dflt_component
character(60) name, class_ele, parameter, component, why_invalid
character(*), parameter :: r_name = 'tao_evaluate_element_parameters'

real(rp), allocatable :: values(:)
real(rp) :: real_val

integer, optional :: dflt_uni
integer i, j, ix, num, ixe, ix1, ios, n_tot, ix_start, ib, where

integer, parameter :: beginning$ = 1, middle$ = 2, end$ = 3

logical err, valid
logical :: print_err

!

call tao_pick_universe (param_name, name, scratch%this_u, err, dflt_uni = dflt_uni)
if (err) return

err = .true.

if (name(1:5) == 'ele::') then
  name = name(6:)  ! Strip off 'ele::'
  where = end$
elseif (name(1:9) == 'ele_mid::') then   
  name = name(10:)  ! Strip off 'ele_mid::'
  where = middle$
elseif (name(1:11) == 'ele_begin::') then   
  name = name(12:)  ! Strip off 'ele_begin::'
  where = beginning$
elseif (dflt_source /= 'element') then
  return
endif

! Get component

ix = index(name, '|')
if (ix == 0) then
  component = 'model'
  if (present(dflt_component)) then
    if (dflt_component /= '') component = dflt_component
  endif
else
  component = name(ix+1:)
  name = name(1:ix-1)
endif

! Get class:name

ix1 = index(name, '[');  
if (ix1 == 0) return

ix1 = index(name, '[');  if (ix1 == 0) return
class_ele = name(1:ix1-1)
name = name(ix1+1:)
if (class_ele(1:2) == '::') class_ele = class_ele(3:)
ix1 = index(name, ']');  if (ix1 == 0) return
parameter = name(1:ix1-1)

! "Intrinsic" element parameter values are not affected by evaluation in the middle.
! It is easier to list what is not intrinsic.

if (where /= end$) then
  select case (parameter)
  ! These are non-intrinsic
  case ('x_position', 'y_position', 'z_position', 'theta_position', 'phi_position', 'psi_position', &
        'beta_a', 'beta_b', 'alpha_a', 'alpha_b', 'gamma_a', 'gamma_b', 'phi_a', 'phi_b', &
        'eta_a', 'eta_b', 'eta_x', 'eta_y', 'eta_z', 'etap_a', 'etap_b', 'etap_x', 'etap_y', 'etap_z', &
        'cmat_11', 'cmat_12', 'cmat_21', 'cmat_22', &
        'orbit_x', 'orbit.x', 'orbit_px', 'orbit.px', 'orbit_y', 'orbit.y', 'orbit_py', 'orbit.py', &
        'orbit_z', 'orbit.z', 'orbit_pz', 'orbit.pz', 'spin.x', 'spin_x', 'spin.y', 'spin_y', &
        'spin.z', 'spin_z', 'intensity', 'intensity_x', 'intensity.x', 'intensity_y', 'intensity.y', &
        'phase_x', 'phase.x', 'phase_y', 'phase.y', 't', 'time', 'beta', 'energy', 'pc', 's')

  case default
    where = end$
  end select
endif

! Evaluate

n_tot = 0
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. scratch%this_u(i)) cycle
  u => s%u(i)
  call tao_locate_elements (class_ele, u%ix_uni, scratch%eles, err)
  if (err) return
  call re_allocate (values, n_tot + size(scratch%eles))

  do j = 1, size(scratch%eles)
    ixe = scratch%eles(j)%ele%ix_ele

    if (parameter == 'index') then
      values(n_tot+j) = ixe
      cycle
    endif

    ib = scratch%eles(j)%ele%ix_branch
    select case (component)
    case ('model');     tao_lat => u%model
    case ('base');      tao_lat => u%base
    case ('design');    tao_lat => u%design
    case default
      call out_io (s_error$, r_name, 'BAD DATUM COMPONENT FOR: ' // param_name)
      return
    end select

      lat        => tao_lat%lat
      tao_branch => tao_lat%tao_branch(ib)
      branch     => lat%branch(ib)

    if (where /= end$ .and. ixe /= 0) then
      ! Need to find element just before the element under consideration. 
      ! This is complicated if the element under consideration is a lord.
      ele0 => branch%ele(ixe)
      do 
        if (ele0%ix_ele <= branch%lat%branch(ele0%ix_branch)%n_ele_track) exit
        ele0 => pointer_to_slave(ele0, 1)
      enddo
      ele0 => pointer_to_next_ele(ele0, -1)
      orb0 => tao_lat%tao_branch(ele0%ix_branch)%orbit(ele0%ix_ele)

      if (where == middle$) then
        select case (parameter)
        case ('x_position', 'y_position', 'z_position', 'theta_position', 'phi_position', 'psi_position')
          call twiss_and_track_intra_ele (branch%ele(ixe), lat%param, 0.0_rp, branch%ele(ixe)%value(l$)/2, &
                                                         .true., .false., orb0, orb, ele0, ele3, compute_floor_coords = .true.)
          err = .true. ! To trigger call to pointer_to_attribute
        case default
          call twiss_and_track_intra_ele (branch%ele(ixe), lat%param, 0.0_rp, branch%ele(ixe)%value(l$)/2, &
                                                                        .true., .false., orb0, orb, ele0, ele3, err)
          call tao_orbit_value (parameter, orb, values(n_tot+j), err)
        end select

        if (err) then
          call pointer_to_attribute (ele3, parameter, .true., a_ptr, err, print_err)
          if (err) return
          values(n_tot+j) = value_of_all_ptr(a_ptr)
        endif

      else
        call tao_orbit_value (parameter, orb0, values(n_tot+j), err)
        if (err) then
          call pointer_to_attribute (ele0, parameter, .true., a_ptr, err, print_err)
          if (err) return
          values(n_tot+j) = value_of_all_ptr(a_ptr)
        endif
      endif

    else
      select case (parameter)
      case ('spin_dn_dpz.x')
        if (.not. tao_branch%spin_valid) call tao_spin_polarization_calc(branch, tao_branch); err = .not. tao_branch%spin_valid
        values(n_tot+j) = tao_branch%dn_dpz(ixe)%vec(1)
      case ('spin_dn_dpz.y')
        if (.not. tao_branch%spin_valid) call tao_spin_polarization_calc(branch, tao_branch); err = .not. tao_branch%spin_valid
        values(n_tot+j) = tao_branch%dn_dpz(ixe)%vec(2)
      case ('spin_dn_dpz.z')
        if (.not. tao_branch%spin_valid) call tao_spin_polarization_calc(branch, tao_branch) ;err = .not. tao_branch%spin_valid
        values(n_tot+j) = tao_branch%dn_dpz(ixe)%vec(3)
      case ('spin_dn_dpz.amp')
        if (.not. tao_branch%spin_valid) call tao_spin_polarization_calc(branch, tao_branch); err = .not. tao_branch%spin_valid
        values(n_tot+j) = norm2(tao_branch%dn_dpz(ixe)%vec)
      case default
        call tao_orbit_value (parameter, tao_branch%orbit(ixe), values(n_tot+j), err)
        if (err) then
          call pointer_to_attribute (branch%ele(ixe), parameter, .true., a_ptr, err, print_err)
          if (err) return
          values(n_tot+j) = value_of_all_ptr(a_ptr)
        endif
      end select
    endif

  enddo

  n_tot = n_tot + size(scratch%eles)
enddo

err = .false.

end subroutine tao_evaluate_element_parameters
