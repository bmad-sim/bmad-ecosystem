!+
! Subroutine tao_evaluate_element_parameters (err, param_name, values, print_err, dflt_ele,
!                                               dflt_source, dflt_component, dflt_uni, dflt_eval_point, info)
!
! Routine to evaluate a lattice element parameter of the form 
!     <universe>@ele::{<ele_class>}::<ele_name_or_num>[<parameter>]{|<component>}
! or to evaluate at the middle of the element
!     <universe>@ele_mid::{<ele_class>}::<ele_name_or_num>[<parameter>]{|<component>}
! Note: size(values) can be zero without an error.
! 
! Input:
!   param_name      -- character(*): parameter name.
!   print_err       -- logical: Print error message? 
!   dflt_ele        -- ele_struct, pointer, optional: Default element if not specified by param_name.
!   dflt_source     -- character(*): Default source
!   dflt_component  -- character(*), optional: Default component
!   dflt_uni        -- integer, optional: Default universe to use.
!   dflt_eval_point -- integer, optional: Evaluation point: anchor_beginning$, anchor_center$, or anchor_end$.
!
! Output:
!   err       -- logical: True if there is an error in syntax. False otherwise
!   values(:) -- real(rp), allocatable: Array of datum values.
!   info(:)   -- tao_expression_info_struct), allocatable, optional:  
!-

subroutine tao_evaluate_element_parameters (err, param_name, values, print_err, dflt_ele, &
                                             dflt_source, dflt_component, dflt_uni, dflt_eval_point, info)

use tao_interface, except_dummy => tao_evaluate_element_parameters

implicit none

type (tao_universe_struct), pointer :: u
type (ele_struct), pointer, optional :: dflt_ele
type (tao_expression_info_struct), allocatable, optional :: info(:)

character(*) param_name
character(*) dflt_source
character(*), optional :: dflt_component
character(60) name, class_ele, parameter, component, why_invalid
character(*), parameter :: r_name = 'tao_evaluate_element_parameters'

real(rp), allocatable :: values(:)
real(rp) :: real_val

integer, optional :: dflt_uni, dflt_eval_point
integer i, j, ix, num, ix1, ios, n_tot, ix_start, where

logical err, valid, use_dflt_ele
logical :: print_err
logical, allocatable :: this_u(:)

!

call tao_pick_universe (param_name, name, this_u, err, dflt_uni = dflt_uni)
if (err) return

err = .true.
where = integer_option(anchor_end$, dflt_eval_point)
use_dflt_ele = .false.

if (name(1:5) == 'ele::') then
  name = name(6:)  ! Strip off 'ele::'
elseif (name(1:9) == 'ele_mid::') then   
  name = name(10:)  ! Strip off 'ele_mid::'
  where = anchor_center$
elseif (name(1:11) == 'ele_begin::') then   
  name = name(12:)  ! Strip off 'ele_begin::'
  where = anchor_beginning$
elseif (present(dflt_ele)) then
  use_dflt_ele = .true.
elseif (dflt_source /= 'ele') then
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

! Get class::name

if (use_dflt_ele) then
  parameter = name
else
  ix1 = index(name, '[');  if (ix1 == 0) return
  class_ele = name(1:ix1-1)
  name = name(ix1+1:)
  if (class_ele(1:2) == '::') class_ele = class_ele(3:)
  ix1 = index(name, ']');  if (ix1 == 0) return
  parameter = name(1:ix1-1)
endif

! "Intrinsic" element parameter values are not affected by evaluation in the middle.
! It is easier to list what is not intrinsic.

if (where /= anchor_end$) then
  select case (parameter)
  ! These are non-intrinsic
  case ('x_position', 'y_position', 'z_position', 'theta_position', 'phi_position', 'psi_position', &
        'beta_a', 'beta_b', 'alpha_a', 'alpha_b', 'gamma_a', 'gamma_b', 'phi_a', 'phi_b', &
        'eta_a', 'eta_b', 'eta_x', 'eta_y', 'eta_z', 'etap_a', 'etap_b', 'etap_x', 'etap_y', 'etap_z', &
        'deta_a_ds', 'deta_b_ds', 'deta_x_ds', 'deta_y_ds', 'deta_z_ds', &
        'cmat_11', 'cmat_12', 'cmat_21', 'cmat_22', 'cmat.11', 'cmat.12', 'cmat.21', 'cmat.22', &
        'cbar_11', 'cbar_12', 'cbar_21', 'cbar_22', 'cbar.11', 'cbar.12', 'cbar.21', 'cbar.22', &
        'orbit_x', 'orbit.x', 'orbit_px', 'orbit.px', 'orbit_y', 'orbit.y', 'orbit_py', 'orbit.py', &
        'orbit_z', 'orbit.z', 'orbit_pz', 'orbit.pz', 'spin.x', 'spin_x', 'spin.y', 'spin_y', &
        'spin.z', 'spin_z', 'intensity', 'intensity_x', 'intensity.x', 'intensity_y', 'intensity.y', &
        'phase_x', 'phase.x', 'phase_y', 'phase.y', 't', 'time', 'beta', 'energy', 'pc', 's')

  case default
    where = anchor_end$
  end select
endif

! Evaluate

n_tot = 0

if (use_dflt_ele) then
  call re_allocate (values, 1)
  call evaluate_this_parameter(dflt_ele, parameter, component, where, s%u(dflt_uni), n_tot, values, err)
  if (present(info)) then
    call tao_re_allocate_expression_info(info, 1)
    info(1)%ele => dflt_ele
  endif

else
  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. this_u(i)) cycle
    u => s%u(i)
    call tao_locate_elements (class_ele, u%ix_uni, scratch%eles, err)
    if (err) return
    call re_allocate (values, n_tot + size(scratch%eles))
    if (present(info)) call tao_re_allocate_expression_info(info, n_tot+size(scratch%eles))

    do j = 1, size(scratch%eles)
      call evaluate_this_parameter(scratch%eles(j)%ele, parameter, component, where, u, n_tot, values, err)
      if (present(info)) info(n_tot)%ele => scratch%eles(j)%ele
      if (err) return
    enddo
  enddo
endif

!-------------------------------------------------------------------------------
contains

subroutine evaluate_this_parameter(ele, parameter, component, where, u, n_tot, values, err)

type (ele_struct) ele
type (tao_universe_struct) :: u
type (tao_lattice_struct), pointer :: tao_lat
type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (tao_lattice_branch_struct), pointer :: tao_branch
type (ele_struct), pointer :: ele0
type (ele_struct) ele3
type (coord_struct) orb
type (coord_struct), pointer :: orb0
type (all_pointer_struct) a_ptr

real(rp), allocatable :: values(:)
integer n_tot, where, ib, ixe
character(*) parameter, component
logical err

! Note: ele may be from the wrong component lattice so do not use directly.

n_tot = n_tot + 1
ixe = ele%ix_ele
err = .true.

if (parameter == 'index') then
  values(n_tot) = ixe
  err = .false.
  return
endif

ib = ele%ix_branch
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

if (where /= anchor_end$ .and. ixe /= 0) then
  ! Need to find element just before the element under consideration. 
  ! This is complicated if the element under consideration is a lord.
  ele0 => branch%ele(ixe)
  do 
    if (ele0%ix_ele <= branch%lat%branch(ele0%ix_branch)%n_ele_track) exit
    ele0 => pointer_to_slave(ele0, 1)
  enddo
  ele0 => pointer_to_next_ele(ele0, -1)
  orb0 => tao_lat%tao_branch(ele0%ix_branch)%orbit(ele0%ix_ele)

  if (where == anchor_center$) then
    select case (parameter)
    case ('x_position', 'y_position', 'z_position', 'theta_position', 'phi_position', 'psi_position')
      call twiss_and_track_intra_ele (branch%ele(ixe), lat%param, 0.0_rp, branch%ele(ixe)%value(l$)/2, &
                                                     .true., .false., orb0, orb, ele0, ele3, compute_floor_coords = .true.)
      err = .true. ! To trigger call to pointer_to_attribute
    case default
      call twiss_and_track_intra_ele (branch%ele(ixe), lat%param, 0.0_rp, branch%ele(ixe)%value(l$)/2, &
                                                                    .true., .false., orb0, orb, ele0, ele3, err)
      values(n_tot) = tao_param_value_at_s (parameter, ele3, ele3, orb, err)
    end select

    if (err) then
      call pointer_to_attribute (ele3, parameter, .true., a_ptr, err, print_err)
      if (err) return
      values(n_tot) = value_of_all_ptr(a_ptr)
    endif

  else
    values(n_tot) = tao_param_value_at_s (parameter, ele0, ele0, orb0, err)
    if (err) then
      call pointer_to_attribute (ele0, parameter, .true., a_ptr, err, print_err)
      if (err) return
      values(n_tot) = value_of_all_ptr(a_ptr)
    endif
  endif

else
  if (parameter(1:12) == 'spin_dn_dpz.') then
    if (.not. allocated(tao_branch%spin_ele)) call tao_spin_polarization_calc(branch, tao_branch)
    if (.not. tao_branch%spin_ele(ixe)%valid) call tao_spin_polarization_calc(branch, tao_branch)
    err = (.not. tao_branch%spin_ele(ixe)%valid)
    if (err) return

    select case (parameter)
    case ('spin_dn_dpz.x')
      values(n_tot) = tao_branch%spin_ele(ixe)%dn_dpz%vec(1)
    case ('spin_dn_dpz.y')
      values(n_tot) = tao_branch%spin_ele(ixe)%dn_dpz%vec(2)
    case ('spin_dn_dpz.z')
      values(n_tot) = tao_branch%spin_ele(ixe)%dn_dpz%vec(3)
    case ('spin_dn_dpz.amp')
      values(n_tot) = norm2(tao_branch%spin_ele(ixe)%dn_dpz%vec)
    case default
      err = .true.
      return  
    end select

  else
    values(n_tot) = tao_param_value_at_s (parameter, branch%ele(ixe), branch%ele(ixe), tao_branch%orbit(ixe), err, print_err = print_err)
    if (err) then
      call pointer_to_attribute (branch%ele(ixe), parameter, .true., a_ptr, err, print_err)
      if (err) return
      values(n_tot) = value_of_all_ptr(a_ptr)
    endif
  endif
endif

err = .false.

end subroutine evaluate_this_parameter

end subroutine tao_evaluate_element_parameters
