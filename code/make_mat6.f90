!+
! Subroutine make_mat6 (ele, param, start_orb, end_orb, end_in, err)
!
! Subroutine to make the 6x6 1st order transfer matrix for an element 
! along with the 0th order transfer vector.
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
!                end is an input only if end_in is set to True.
!   end_in  -- Logical, optional: If present and True then the end coords
!                will be taken as input. not output as normal. 
!                If end_orb is not present then end_in will be ignored.
!
!
! Output:
!   ele     -- Ele_struct: Element
!     %mat6   -- Real(rp): 1st order 6x6 transfer matrix.
!     %vec0   -- Real(rp): 0th order transfer vector.
!   end_orb -- Coord_struct, optional: Coordinates at the end of element.
!                end is an output if end_in is not set to True.
!   param   -- lat_param_struct:
!     %lost   -- Since make_mat6 may do tracking %lost may be set to True if
!                  tracking was unsuccessful. %lost set to False otherwise.
!   err     -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine make_mat6 (ele, param, start_orb, end_orb, end_in, err)

use symp_lie_mod, only: symp_lie_bmad
use bookkeeper_mod, only: attribute_bookkeeper
use mad_mod, only: make_mat6_mad
use space_charge_mod, except_dummy => make_mat6
use equality_mod

implicit none

type (ele_struct), target :: ele
type (coord_struct), optional :: start_orb, end_orb
type (lat_param_struct)  param
type (coord_struct) a_start_orb, a_end_orb

integer mat6_calc_method

logical, optional :: end_in, err
logical end_input, rad_fluct_save

!--------------------------------------------------------
! Some Init

mat6_calc_method = ele%mat6_calc_method

if (present(start_orb)) then
  a_start_orb = start_orb
else
  call init_coord (a_start_orb)
endif

end_input = (logic_option (.false., end_in) .and. present(end_orb))
if (end_input) a_end_orb = end_orb

! custom calc 

if (mat6_calc_method == custom$) then
  call make_mat6_custom (ele, param, a_start_orb, a_end_orb)
  return
endif

! init

param%lost = .false.
if (bmad_com%auto_bookkeeper) call attribute_bookkeeper (ele, param)

if (.not. ele%is_on) mat6_calc_method = bmad_standard$

!

if (present(err)) err = .false.

rad_fluct_save = bmad_com%radiation_fluctuations_on
bmad_com%radiation_fluctuations_on = .false.

!

select case (mat6_calc_method)

case (taylor$)
  call make_mat6_taylor (ele, param, a_start_orb)
  if (.not. end_input) call track1_taylor (a_start_orb, ele, param, a_end_orb)

case (bmad_standard$)
  call make_mat6_bmad (ele, param, a_start_orb, a_end_orb, end_in, err)

case (symp_lie_ptc$)
  call make_mat6_symp_lie_ptc (ele, param, a_start_orb)
  if (.not. end_input) call track1_taylor (a_start_orb, ele, param, a_end_orb)

case (symp_lie_bmad$)
  call symp_lie_bmad (ele, param, a_start_orb, a_end_orb, .true.)

case (tracking$)
  call make_mat6_tracking (ele, param, a_start_orb, a_end_orb)

case (mad$)
  call make_mat6_mad (ele, param, a_start_orb, a_end_orb)

case (no_method$)
  return

case default
  print *, 'ERROR IN MAKE_MAT6: UNKNOWN MAT6_CALC_METHOD: ', &
                                  calc_method_name(ele%mat6_calc_method)
  call err_exit
end select

! Add space charge effects

if (bmad_com%space_charge_on) call make_mat6_ultra_rel_space_charge (ele, param)

! symplectify if wanted

if (ele%symplectify) call mat_symplectify (ele%mat6, ele%mat6)

! Finish up

if (.not. ele%map_ref_orb_in == a_start_orb) then
  ele%map_ref_orb_in = a_start_orb
  if (associated(ele%const)) ele%const = -1
endif

ele%map_ref_orb_out = a_end_orb
if (present(end_orb) .and. .not. end_input) end_orb = a_end_orb

bmad_com%radiation_fluctuations_on = rad_fluct_save


end subroutine

