!+
! Subroutine make_mat6 (ele, param, start_orb, end_orb, end_in, err_flag)
!
! Subroutine to make the 6x6 1st order transfer matrix for an element 
! along with the 0th order transfer vector.
!
! Note: end_orb is used as input (when end_in = True) when the orbit through the element is known.
! This makes the calculation faster. The default is to treat end_orb as output.
!
! Note: Radiation fluctuations (but not damping) is turned off for the calculation.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele     -- Ele_struct: Element holding the transfer matrix.
!   param   -- lat_param_struct: Lattice global parameters.
!   start   -- Coord_struct, optional: Coordinates at the beginning of element. 
!                If not present then default is start = 0.
!   end_orb -- Coord_struct, optional: Coordinates at the end of element.
!                end_orb is an input only if end_in is set to True.
!   end_in  -- Logical, optional: If present and True then the end coords
!                will be taken as input. not output as normal. 
!                If end_orb is not present then end_in will be ignored.
!
!
! Output:
!   ele      -- Ele_struct: Element
!     %mat6    -- Real(rp): 1st order 6x6 transfer matrix.
!     %vec0    -- Real(rp): 0th order transfer vector.
!   end_orb  -- Coord_struct, optional: If end_in is not set to True:
!                 Coordinates at the end of element.
!   err_flag -- Logical, optional: Set True if there is an error. False otherwise.
!-

recursive subroutine make_mat6 (ele, param, start_orb, end_orb, end_in, err_flag)

use symp_lie_mod, only: symp_lie_bmad
use bookkeeper_mod, only: attribute_bookkeeper
use mad_mod, only: make_mat6_mad
use space_charge_mod, except_dummy => make_mat6
use equality_mod, only: operator(==)

implicit none

type (ele_struct), target :: ele
type (coord_struct), optional :: start_orb, end_orb
type (lat_param_struct)  param
type (coord_struct) a_start_orb, a_end_orb

real(rp), parameter :: zero_vec(6) = 0
integer mat6_calc_method, species

logical, optional :: end_in, err_flag
logical end_input, rad_fluct_save, err, finished

character(16), parameter :: r_name = 'make_mat6'

!--------------------------------------------------------
! Some init.
! If start_orb is in its not_set state (can happen if a particle is lost in 
! tracking and ele is downstream from the loss point), init the orbit to zero.

if (present(err_flag)) err_flag = .true.
if (ele%bookkeeping_state%mat6 == stale$) ele%bookkeeping_state%mat6 = ok$

if (present(start_orb)) then
  if (start_orb%state == not_set$) then
    call init_coord(a_start_orb, zero_vec, ele, upstream_end$, param%particle)
  else
    a_start_orb = start_orb
  endif
else
  call init_coord (a_start_orb, ele = ele, element_end = upstream_end$, particle = default_tracking_species(param))
endif

end_input = (logic_option (.false., end_in) .and. present(end_orb))
if (end_input) then
  if (end_orb%state == not_set$) then
    end_input = .false.
  else
    a_end_orb = end_orb
  endif
endif

! init

if (bmad_com%auto_bookkeeper) call attribute_bookkeeper (ele, param)

mat6_calc_method = ele%mat6_calc_method
if (.not. ele%is_on) mat6_calc_method = bmad_standard$

!

rad_fluct_save = bmad_com%radiation_fluctuations_on
bmad_com%radiation_fluctuations_on = .false.

!

err = .false.

select case (mat6_calc_method)

case (custom$)
  call make_mat6_custom (ele, param, a_start_orb, a_end_orb, err)

case (taylor$)
  call make_mat6_taylor (ele, param, a_start_orb)
  if (.not. end_input) call track1_taylor (a_start_orb, ele, param, a_end_orb)

case (bmad_standard$)
  if (a_start_orb%species == photon$) then
    call make_mat6_bmad_photon (ele, param, a_start_orb, a_end_orb, end_input, err)
  else
    call make_mat6_bmad (ele, param, a_start_orb, a_end_orb, end_input, err)
  endif

case (symp_lie_ptc$)
  call make_mat6_symp_lie_ptc (ele, param, a_start_orb, a_end_orb)

case (symp_lie_bmad$)
  call symp_lie_bmad (ele, param, a_start_orb, a_end_orb, .true.)

case (tracking$)
  call make_mat6_tracking (ele, param, a_start_orb, a_end_orb)

case (mad$)
  call make_mat6_mad (ele, param, a_start_orb, a_end_orb)

! Static is used with hybrid elements since, in this case, the transfer map is not recomputable.

case (static$)
  if (present(err_flag)) err_flag = .false.
  if (present(end_orb) .and. .not. logic_option (.false., end_in)) call track1 (a_end_orb, ele, param, end_orb)
  return

case default
  call out_io (s_fatal$, r_name, 'UNKNOWN MAT6_CALC_METHOD: ' // mat6_calc_method_name(ele%mat6_calc_method))
  if (global_com%exit_on_error) call err_exit
  return
end select

if (err) then
  if (present(end_orb) .and. .not. logic_option (.false., end_in)) end_orb = a_end_orb
  return
endif

! Add space charge effects

if (bmad_com%space_charge_on) call make_mat6_ultra_rel_space_charge (ele, param)

! symplectify if wanted

if (ele%symplectify) call mat_symplectify (ele%mat6, ele%mat6, ele%value(p0c$)/ele%value(p0c_start$))

! Finish up

if (any(ele%map_ref_orb_in%vec /= a_start_orb%vec)) then
  if (associated(ele%rad_int_cache)) ele%rad_int_cache%stale = .true.
endif

ele%map_ref_orb_in = a_start_orb
ele%map_ref_orb_out = a_end_orb
if (present(end_orb) .and. .not. logic_option (.false., end_in)) end_orb = a_end_orb

bmad_com%radiation_fluctuations_on = rad_fluct_save

if (present(err_flag)) err_flag = .false.

end subroutine

