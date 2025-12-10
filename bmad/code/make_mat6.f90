!+
! Subroutine make_mat6 (ele, param, start_orb, end_orb, err_flag)
!
! Subroutine to make the 6x6 1st order transfer matrix for an element 
! along with the 0th order transfer vector. Also optionally track the particle.
!
! Note: Radiation fluctuations (but not damping) is turned off for the calculation.
!
! Input:
!   ele       -- Ele_struct: Element holding the transfer matrix.
!   param     -- lat_param_struct: Lattice global parameters.
!   start_orb -- Coord_struct, optional: Reference coordinates at the beginning of element. 
!                  If not present, default is to use the zero orbit.
!
! Output:
!   ele      -- Ele_struct: Element
!     %mat6    -- Real(rp): 1st order 6x6 transfer matrix.
!     %vec0    -- Real(rp): 0th order transfer vector.
!   end_orb  -- Coord_struct, optional: Reference coordinates at the end of element.
!   err_flag -- Logical, optional: Set True if there is an error. False otherwise.
!-

recursive subroutine make_mat6 (ele, param, start_orb, end_orb, err_flag)

use bookkeeper_mod, only: attribute_bookkeeper
use mad_mod, only: make_mat6_mad
use high_energy_space_charge_mod, except_dummy => make_mat6
use equality_mod, only: operator(==)

implicit none

type (ele_struct), target :: ele
type (coord_struct), optional :: start_orb, end_orb
type (lat_param_struct)  param
type (coord_struct) a_start_orb, a_end_orb, b_start_orb

real(rp), parameter :: zero_vec(6) = 0
integer mat6_calc_method, species

logical, optional :: err_flag
logical rad_fluct_save, err, finished

character(*), parameter :: r_name = 'make_mat6'

!--------------------------------------------------------
! The beginning element is handled specially.
! Also see twiss_propagate1.

if (ele%key == beginning_ele$) then
  return
endif

! X-ray mat6 calc not implemented.
! Since X-ray tracking is highly nonlinear (think cyrstal diffraction and apertures, etc.), a linear
! transfer matrix is generally useful.

if (ele%ref_species == photon$) return

! Some init.
! If start_orb is in its not_set state (can happen if a particle is lost in 
! tracking and ele is downstream from the loss point), init the orbit to zero.

if (present(err_flag)) err_flag = .true.

if (.not. present(start_orb)) then
  call init_coord (a_start_orb, zero_vec, ele, upstream_end$, default_tracking_species(param))
else if (start_orb%state == not_set$ .or. significant_difference(start_orb%p0c, ele%value(p0c_start$), rel_tol = small_rel_change$)) then
  call init_coord(a_start_orb, start_orb, ele, upstream_end$, default_tracking_species(param))
else
  a_start_orb = start_orb
endif

if (a_start_orb%direction == -1) then  ! Can only happen if start_orb is present
  call out_io (s_fatal$, r_name, 'TRANSFER MATRIX CALCULATION NOT ALLOWED WITH BACKWARD TRACKING.')
  if (global_com%exit_on_error) call err_exit
  return
endif

! To keep the code simple, only use forward tracking.

a_start_orb%time_dir = 1

! init

if (bmad_com%auto_bookkeeper) call attribute_bookkeeper (ele)

mat6_calc_method = ele%mat6_calc_method
if (.not. ele%is_on) mat6_calc_method = bmad_standard$
if (mat6_calc_method == auto$) then
  select case (ele%tracking_method)
  case (bmad_standard$, linear$);   mat6_calc_method = bmad_standard$
  case (custom$);                   mat6_calc_method = custom$
  case (mad$);                      mat6_calc_method = mad$
  case (symp_lie_bmad$);            mat6_calc_method = symp_lie_bmad$
  case (symp_lie_ptc$);             mat6_calc_method = symp_lie_ptc$
  case (taylor$);                   mat6_calc_method = taylor$
  case (runge_kutta$, fixed_step_runge_kutta$, time_runge_kutta$, fixed_step_time_runge_kutta$)
    mat6_calc_method = tracking$
  case default
    call out_io (s_fatal$, r_name, 'UNKNOWN TRACKING_METHOD: \i0\ ', ele%tracking_method)
    if (global_com%exit_on_error) call err_exit
    return
  end select

  if (ele%key == foil$) mat6_calc_method = tracking$
endif

ele%map_ref_orb_in = a_start_orb

rad_fluct_save = bmad_com%radiation_fluctuations_on
bmad_com%radiation_fluctuations_on = .false.

! Compute matrix
! Matrix must be made around the zero orbit for linear tracking.

err = .false.
if (ele%tracking_method == linear$) then
  b_start_orb = a_start_orb
  a_start_orb%vec = 0  
endif

select case (mat6_calc_method)

case (custom$)
  if (.not. associated(make_mat6_custom_ptr)) then
    call out_io (s_error$, r_name, 'CUSTOM MAT6_CALC_METHOD OR CUSTOM ELEMENT NEEDS MAKE_MAT6_CUSTOM_PTR SET IN THIS PROGRAM!', &
                                   'NEEDED FOR ELEMENT: ' // ele%name)
    a_end_orb%state = lost$
  endif

  call make_mat6_custom_ptr (ele, param, a_start_orb, a_end_orb, err)

case (taylor$)
  call make_mat6_taylor (ele, a_start_orb, a_end_orb, err)

case (bmad_standard$)
  if (a_start_orb%species == photon$) then
    call make_mat6_bmad_photon (ele, param, a_start_orb, a_end_orb, err)
  else
    call make_mat6_bmad (ele, param, a_start_orb, a_end_orb, err)
  endif

case (symp_lie_ptc$)
  call make_mat6_symp_lie_ptc (ele, a_start_orb, a_end_orb)

case (symp_lie_bmad$)
  a_end_orb = a_start_orb
  call symp_lie_bmad (ele, param, a_end_orb, mat6 = ele%mat6, make_matrix = .true.)

case (tracking$)
  call make_mat6_tracking (ele, param, a_start_orb, a_end_orb, err)

case (mad$)
  call make_mat6_mad (ele, param, a_start_orb, a_end_orb)

! Static is used with hybrid elements since, in this case, the transfer map is not recomputable.

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN MAT6_CALC_METHOD: ' // mat6_calc_method_name(ele%mat6_calc_method))
  if (global_com%exit_on_error) call err_exit
  return
end select

if (err) then
  if (present(end_orb)) end_orb = a_end_orb
  return
endif

! Add space charge effects

if (bmad_com%high_energy_space_charge_on) call make_mat6_high_energy_space_charge (ele, param)

! symplectify if wanted

if (ele%symplectify) call mat_symplectify (ele%mat6, ele%mat6, ele%value(p0c$)/ele%value(p0c_start$))

! If the tracking_method is not consistant with the mat6_calc_method then need to track.

if (present(end_orb)) then
  if (ele%tracking_method == linear$) a_start_orb = b_start_orb
  if (.not. ele%is_on .or. mat6_calc_method == tracking$ .or. mat6_calc_method == ele%tracking_method) then
    end_orb = a_end_orb
  else
    call track1 (a_start_orb, ele, param, end_orb)
  endif

  if (end_orb%state /= alive$) then
    end_orb%location = inside$
  elseif (end_orb%direction == 1) then
    end_orb%location = downstream_end$
  else
    end_orb%location = upstream_end$
  endif
endif


! Finish up mat6 calc

ele%map_ref_orb_out = a_end_orb
if (ele%bookkeeping_state%mat6 == stale$) ele%bookkeeping_state%mat6 = ok$

! Spin

!if (bmad_com%spin_tracking_on .and. a_start_orb%species /= photon$) then
!  if (.not. associated(ele%spin_taylor(0)%term)) then
!    call ele_to_spin_taylor(ele, param, a_start_orb)
!  endif
!
!  if (associated(ele%spin_taylor(0)%term)) then
!    ele%spin_q = spin_taylor_to_linear (ele%spin_taylor, .true., a_start_orb%vec-ele%spin_taylor_ref_orb_in, ele%is_on)
!  else
!    call out_io (s_error$, r_name, 'CANNOT CONSTRUCT SPIN TAYLOR MAP FOR ELEMENT: ' // ele_full_name(ele))
!  endif
!endif

! And end

bmad_com%radiation_fluctuations_on = rad_fluct_save
if (present(err_flag)) err_flag = .false.

end subroutine

