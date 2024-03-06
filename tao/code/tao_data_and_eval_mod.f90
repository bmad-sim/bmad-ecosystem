module tao_data_and_eval_mod

use tao_interface
use bmad_interface

implicit none

! Used for parsing expressions

integer, parameter :: var_num$ = 101, lat_num$ = 102, data_num$ = 103, ele_num$ = 104

private tao_scratch_values_calc, tao_eval_floor_orbit, tao_ele_geometry_with_misalignments

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_lat_or_beam_data (err, data_name, values, print_err, default_source, default_source,  
!               dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_component, dflt_uni, dflt_eval_point, dflt_s_offset)
!
! Routine to evaluate data with a lat or beam source of the form:
!     <universe>@lat::<data_type>[<ix_ele_start>&<ix_ele>]|<component>
!
! Input:
!   data_name      -- character(*): data name.
!   print_err      -- logical: Print error message?
!   dflt_source    -- character(*): If not blank: Default source: 'lat' or 'beam'.
!   dflt_ele_ref   -- ele_struct, pointer, optional: Default reference element.
!   dflt_ele_start -- ele_struct, pointer, optional: Default start element.
!   dflt_ele       -- ele_struct, pointer, optional: Default element to evaluate at.
!   dflt_component -- character(*), optional: Default component: 'model' (default), 'base', or 'design'.
!   dflt_uni       -- integer, optional: Default universe to use.
!   dflt_eval_point -- integer, optional: Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.
!   dflt_s_offset   -- real(rp), optional: Default offset of eval_point. Default = 0.
!
! Output:
!   err       -- Logical: True if there is an error. False otherwise.
!   values(:) -- Real(rp), allocatable: Array of datum valuse.
!-

subroutine tao_evaluate_lat_or_beam_data (err, data_name, values, print_err, default_source, &
                dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_component, dflt_uni, dflt_eval_point, dflt_s_offset)

type (tao_data_struct) datum
type (tao_universe_struct), pointer :: u
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer, optional :: dflt_ele_ref, dflt_ele_start, dflt_ele
type (ele_struct), pointer :: this_ele

character(*) data_name
character(*) default_source
character(*), optional :: dflt_component
character(100) name, ele_name, component, offset_str
character(*), parameter :: r_name = 'tao_evaluate_lat_or_beam_data'

real(rp), allocatable :: values(:), off_val(:)
real(rp), optional :: dflt_s_offset
real(rp) s_offset

integer, optional :: dflt_uni, dflt_eval_point
integer j, num, ix, ix1, ios, n_tot, n_loc, iu

logical err, valid, err_flag
logical print_err, use_dflt_ele, has_assoc_ele
logical, allocatable :: this_u(:)

!

datum%exists = .true.
datum%ele_start_name = ''
datum%ele_ref_name = ''
datum%ix_ele_start = -1
datum%ix_ele_ref = -1

call tao_pick_universe (data_name, name, this_u, err, dflt_uni = dflt_uni)
if (err) return

err = .true.
if (name(1:5) == 'lat::') then
  datum%data_source = 'lat'
  name = name(6:)  ! Strip off 'lat:'
elseif (name(1:6) == 'beam::') then
  datum%data_source = 'beam'
  name = name(7:)  ! Strip off 'beam:'
elseif (default_source /= '') then
  datum%data_source = default_source
else
  if (print_err) call out_io (s_error$, r_name, 'DATUM NOT "LAT::" OR "BEAM::"' // data_name)
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

! Get data type

ele_name = ''
offset_str = ''
s_offset = real_option(0.0_rp, dflt_s_offset)
use_dflt_ele = .true.
has_assoc_ele = .true.
ix1 = index(name, '[')

if (ix1 == 0) then
  datum%data_type = name
  if (present(dflt_ele)) then
    if (associated(dflt_ele)) ele_name = dflt_ele%name
  endif

  has_assoc_ele = (tao_datum_has_associated_ele(name) == yes$)
  if (ele_name == '' .and. has_assoc_ele) then
    if (print_err) call out_io (s_error$, r_name, 'NO "[" FOUND IN: ' // data_name)
    return
  endif

else
  datum%data_type = name(1:ix1-1)
  name = name(ix1+1:)
  ix1 = str_first_in_set(name, ']', .true.)
  if (ix1 == 0) then
    if (print_err) call out_io (s_error$, r_name, 'NO "]" FOUND IN: ' // data_name)
    return
  endif
  name(ix1:ix1) = ''
  if (name(ix1+1:) /= '') then
    if (print_err) call out_io (s_error$, r_name, 'MANGLED CONSTRUCT: ' // data_name)
    return
  endif
  use_dflt_ele = .false.

  ! Get ele_ref & ele -> s_offset

  ix = index(name, '->')
  if (ix /= 0) then
    offset_str = name(ix+2:)
    name = name(:ix-1)
  endif

  ix = index(name, '&')
  if (ix /= 0) then
    if (is_integer(name)) then
      read (name(:ix-1), *, iostat = ios) datum%ix_ele_ref
      if (ios /= 0) then
        if (print_err) call out_io (s_error$, r_name, 'BAD ELE_REF: ' // data_name)
        return
      endif
    else
      datum%ele_ref_name = name(:ix-1)
    endif
    ele_name = name(ix+1:)
  else
    ele_name = name
  endif
endif

! Evaluate

n_tot = 0
do iu = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. this_u(iu)) cycle
  u => s%u(iu)

  if (has_assoc_ele) then
    if (use_dflt_ele) then
      n_loc = 1
      if (present(dflt_ele_start)) then
        if (associated(dflt_ele_start)) then
          n_loc = dflt_ele%ix_ele - dflt_ele_start%ix_ele + 1
          if (n_loc < 0) n_loc = n_loc + dflt_ele%branch%n_ele_track + 1
        endif
      endif

      if (present(dflt_ele_ref)) then
        if (associated(dflt_ele_ref)) then
          datum%ix_ele_ref = dflt_ele_ref%ix_ele
        endif
      endif

    else
      if (datum%ele_ref_name /= '') then
        call lat_ele_locator (datum%ele_ref_name, u%model%lat, eles, n_loc, err_flag)
        if (err_flag) return
        if (n_loc /= 1) then
          if (print_err) call out_io (s_error$, r_name, &
                          'MULTIPLE ELEMENTS MATCH REFERENCE NAME: ' // datum%ele_ref_name)
          return
        endif
        datum%ix_ele_ref = eles(1)%ele%ix_ele
      endif

      ! If ele_name = "*" throw out group, overlay, and girder elements
      if (ele_name == '*') ele_name = trim(ele_name) // ' ~group::* ~overlay::* ~girder::*'
      call lat_ele_locator (ele_name, u%model%lat, eles, n_loc, err_flag)
      if (err_flag) return
    endif

  else
    n_loc = 1
  endif

  call re_allocate (values, n_tot + n_loc)

  !

  do j = 1, n_loc
    if (has_assoc_ele) then
      if (use_dflt_ele) then
        this_ele => dflt_ele
        if (present(dflt_ele_start)) then
          if (associated(dflt_ele_start)) then
            datum%ix_ele = dflt_ele_start%ix_ele + j - 1
            if (datum%ix_ele > dflt_ele%branch%n_ele_track) datum%ix_ele = datum%ix_ele - dflt_ele%branch%n_ele_track + 1
          endif
        endif

      else
        this_ele => eles(j)%ele
        datum%ix_ele = eles(j)%ele%ix_ele
        datum%ix_branch = eles(j)%ele%ix_branch
      endif

      datum%ele_name  = this_ele%name
      datum%ix_branch = this_ele%ix_branch
      datum%ix_ele    = this_ele%ix_ele

      ! Offset_str may be something like "L/2" where L is the element length.
      if (offset_str /= '') then
        call tao_evaluate_expression(offset_str, 1, .false., off_val, err_flag, .true., &
                                             dflt_source = 'ele', dflt_ele = this_ele, dflt_uni = iu)
        if (err_flag) then
          if (print_err) call out_io (s_error$, r_name, 'BAD S_OFFSET: ' // data_name)
          return
        endif
        s_offset = off_val(1)
      endif

      datum%eval_point = integer_option(anchor_end$, dflt_eval_point)
      datum%s_offset = s_offset
    endif

    select case (component)
    case ('model')   
      call tao_evaluate_a_datum (datum, u, u%model, values(n_tot+j), valid)
    case ('base')  
      call tao_evaluate_a_datum (datum, u, u%base, values(n_tot+j), valid)
    case ('design')  
      call tao_evaluate_a_datum (datum, u, u%design, values(n_tot+j), valid)
    case default
      if (print_err) call out_io (s_error$, r_name, 'BAD DATUM COMPONENT FOR: ' // data_name)
      return
    end select

    if (valid) err = .false.
  enddo

  n_tot = n_tot + n_loc
enddo

if (n_tot == 0) then
  if (print_err) call out_io (s_error$, r_name, 'ELEMENT NOT FOUND: ' // ele_name)
  return
endif

end subroutine tao_evaluate_lat_or_beam_data

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value)
!
! Buffer routine for to_phase_and_coupling_reading.
!
! Input:
!   ele -- Ele_struct: The monitor.
!
! Output:
!   bpm_data     -- Bpm_phase_coupling_struct: Monitor values
!   valid_value  -- Logical: Valid data value?
!-

subroutine tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value, why_invalid, datum)

use measurement_mod, only: to_phase_and_coupling_reading, ele_is_monitor

type (ele_struct) ele
type (bpm_phase_coupling_struct) bpm_data
type (bpm_phase_coupling_struct), save :: old_bpm_data
type (tao_data_struct) datum

integer, save :: ix_ele_old = -1

logical valid_value
logical, save :: err

character(*) why_invalid

!

valid_value = .false.

if (ix_ele_old /= ele%ix_ele) then
  call to_phase_and_coupling_reading (ele, s%com%add_measurement_noise, old_bpm_data, err)
  if (err) then
    if (ele%a%beta == 0) then;               call tao_set_invalid (datum, 'UNSTABLE LATTICE', why_invalid)
    elseif (.not. ele%is_on) then;           call tao_set_invalid (datum, 'ELEMENT IS TURNED OFF.', why_invalid)
    elseif (.not. ele_is_monitor(ele)) then; call tao_set_invalid (datum, &
              'ELEMENT TYPE IS NOT SUITABLE FOR BEAM MONITORING: ' // key_name(ele%key), why_invalid, .true.)
    endif
    return
  endif
  ix_ele_old = ele%ix_ele
endif

bpm_data = old_bpm_data
valid_value = .true.

end subroutine tao_to_phase_and_coupling_reading

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_get_data (data_value, data_weight, data_meas_value, dat_ix_dModel)
!
! Subroutine to get the values of the data used in optimization and put them
! in an array. The data is ordered starting with the first universe
!
! Input:
! 
! Output:
!   data_value(:)       -- Real(Rp), allocatable, optional: Data model values.
!   data_weight(:)      -- Real(Rp), allocatable, optional: Data weights in the merit function.
!   data_meas_value(:)  -- Real(Rp), allocatable, optional: Data values when the data was taken.
!   data_ix_dModel(:)   -- Integer, allocatable, optional: Data ix_dModel indices
!-

subroutine tao_get_data (data_value, data_weight, data_meas_value, data_ix_dModel)

real(rp), allocatable, optional :: data_value(:), data_meas_value(:), data_weight(:)
integer, allocatable, optional :: data_ix_dModel(:)

integer i, j, iu
integer n_data

!
  
n_data = 0
do iu = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(iu)%is_on) cycle
  n_data  = n_data + count(s%u(iu)%data(:)%useit_opt)
enddo
if (present(data_value))      call re_allocate (data_value, n_data)
if (present(data_meas_value)) call re_allocate (data_meas_value, n_data)
if (present(data_weight))     call re_allocate (data_weight, n_data)
if (present(data_ix_dModel))  call re_allocate (data_ix_dModel, n_data)

j = 0
do iu = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(iu)%is_on) cycle
  do i = 1, size(s%u(iu)%data)
    if (.not. s%u(iu)%data(i)%useit_opt) cycle
    j = j + 1
    if (present(data_value))        data_value(j)      = s%u(iu)%data(i)%model_value
    if (present(data_weight))       data_weight(j)     = s%u(iu)%data(i)%weight
    if (present(data_meas_value))   data_meas_value(j) = s%u(iu)%data(i)%meas_value
    if (present(data_ix_dModel))    data_ix_dModel(j)  = s%u(iu)%data(i)%ix_dModel
  enddo
enddo

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_data_coupling_init (u)
!
! Routine to initialize the coupling structure for a lattice branch.
! This routine is called by tao_lattic_calc and is not meant for general use.
!
! Input:
!   branch -- branch_struct: New lattice branch.
!-

subroutine tao_data_coupling_init (branch)

type (branch_struct) branch
integer m

! 

m = branch%n_ele_max
if (.not. allocated(scratch%cc)) allocate (scratch%cc(0:m))
if (ubound(scratch%cc, 1) < m) then
  deallocate(scratch%cc)
  allocate(scratch%cc(0:m))
endif

scratch%cc%coupling_calc_done = .false.
scratch%cc%amp_calc_done = .false.

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, why_invalid)
!
! Subroutine to put the proper data in the specified datum
!
! Input:
!   datum          -- Tao_data_struct: What type of datum
!   u              -- Tao_universe_struct: Which universe to use.
!   tao_lat        -- Tao_lattice_struct: Lattice to use.
!
! Output:
!   datum          -- Tao_data_struct: 
!     %ix_ele_merit   -- For max/min type constraints: Place where value is max/min. 
!   datum_value   -- Real(rp): Value of the datum.
!   valid_value   -- Logical: Set false when there is a problem. Set true otherwise.
!   why_invalid   -- Character(*), optional: Tells why datum value is invalid.
!-

recursive subroutine tao_evaluate_a_datum (datum, u, tao_lat, datum_value, valid_value, why_invalid)

use pointer_lattice, only: operator(.sub.)
use ptc_interface_mod, only: taylor_inverse
use ptc_layout_mod, only: normal_form_rd_terms
use measurement_mod, only: to_orbit_reading, to_eta_reading, ele_is_monitor
use expression_mod, only: numeric$

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (tao_lattice_branch_struct), pointer :: tao_branch
type (tao_building_wall_section_struct), pointer :: section
type (tao_building_wall_point_struct) :: pt
type (tao_data_struct), pointer :: dp
type (tao_data_array_struct), allocatable :: d_array(:)
type (tao_lattice_struct), target :: tao_lat
type (tao_expression_info_struct), allocatable :: info(:)
type (lat_struct), pointer :: lat
type (tao_dynamic_aperture_struct), pointer :: da
type (aperture_scan_struct), pointer :: scan
type (normal_modes_struct) mode
type (ele_struct), pointer :: ele, ele_start, ele_ref, ele2
type (ele_pointer_struct), allocatable :: eles(:)
type (coord_struct), pointer :: orb0, orbit(:), orb
type (coord_struct) :: orb_at_s, orb1
type (bpm_phase_coupling_struct) bpm_data
type (taylor_struct), save :: taylor_save(6), taylor(6) ! Saved taylor map
type (floor_position_struct) floor
type (branch_struct), pointer :: branch
type (bunch_params_struct), pointer :: bunch_params(:)
type (bmad_normal_form_struct), pointer :: bmad_nf
type (ptc_normal_form_struct), pointer :: ptc_nf
type (taylor_struct), pointer :: taylor_ptr
type (complex_taylor_struct), pointer :: complex_taylor_ptr
type (all_pointer_struct) a_ptr
type (rad_int_branch_struct), pointer :: rad_int_branch
type (c_taylor), pointer :: phase_map
type (twiss_struct), pointer :: z0, z1, z2

real(rp) datum_value, mat6(6,6), vec0(6), angle, px, py, vec2(2)
real(rp) eta_vec(4), v_mat(4,4), v_inv_mat(4,4), a_vec(4), mc2, charge
real(rp) beta_gamma, one_pz, xi_sum, xi_diff, w0_mat(3,3), w_mat(3,3), vec3(3), value, s_len, n0(3)
real(rp) dz, dx, cos_theta, sin_theta, zz_pt, xx_pt, zz0_pt, xx0_pt, dE, s_offset
real(rp) zz_center, xx_center, xx_wall, phase, amp, dalpha, dbeta, aa, bb, g2
real(rp) xx_a, xx_b, dxx1, dzz1, drad, ang_a, ang_b, ang_c, dphi, amp_a, amp_b
real(rp), allocatable :: value_vec(:)
real(rp), allocatable :: expression_value_vec(:)
real(rp) theta, phi, psi

complex(rp) eval(6), evec(6,6), n_eigen(6,3)
complex(rp) temp_cplx

integer i, j, jj, k, m, n, k_old, ix, ie, is, iz, ix_ele, ix_start, ix_ref, ie0, ie1
integer n_size, ix0, which, expnt(6), n_track, n_max, ix_branch, expo(6), n_da

character(*), optional :: why_invalid
character(6) expn_str
character(16) constraint
character(20) :: r_name = 'tao_evaluate_a_datum'
character(40) head_data_type, sub_data_type, data_source, name, dflt_dat_index
character(100) data_type, str
character(:), allocatable :: e_str

logical found, valid_value, err, taylor_is_complex, use_real_part, term_found, ok
logical particle_lost, exterminate, printit
logical, allocatable :: good(:)

! If does not exist

valid_value = .false.
datum%why_invalid = ''

if (.not. datum%exists) then
  datum_value = real_garbage$
  call tao_set_invalid(datum, 'Datum does not exist.', why_invalid)
  return
endif

! To save time, don't evaluate if unnecessary when the running an optimizer.
! Exception: When there are datums that use expressions, things are 
!   complicated so don't try to save time in this case.

if (s%com%optimizer_running .and. .not. datum%useit_opt .and. .not. s%com%have_datums_using_expressions) then
  datum_value = 0
  return
endif

! See if there is a hook for this datum

if (associated(tao_hook_evaluate_a_datum_ptr)) then
  call tao_hook_evaluate_a_datum_ptr (found, datum, u, tao_lat, datum_value, valid_value, why_invalid)
  if (found) return
endif

! Set ix_ele, ix_ref, and ix_start. 
! Note: To check that a start element was set, need to look at datum%ele_start_name, not ix_start.

data_source = datum%data_source
data_type = datum%data_type  ! Needed since %data_type is a var length str so evaluating something like %data_type(1:10) is problematical
head_data_type = datum%data_type
lat => tao_lat%lat

if (head_data_type == 'null') then
  datum_value = 0
  call tao_set_invalid (datum, 'Datum data_type is set to "null".', why_invalid)
  valid_value = .false.
  return
endif

ele => tao_pointer_to_datum_ele (lat, datum%ele_name, datum%ix_ele, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_ele = -1
ix_branch = datum%ix_branch
if (associated(ele)) ix_ele = tao_tracking_ele_index(ele, datum, ix_branch)

ele_ref => tao_pointer_to_datum_ele (lat, datum%ele_ref_name, datum%ix_ele_ref, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_ref = -1
if (associated(ele_ref)) ix_ref = tao_tracking_ele_index(ele_ref, datum)

ele_start => tao_pointer_to_datum_ele (lat, datum%ele_start_name, datum%ix_ele_start, datum, valid_value, why_invalid)
if (.not. valid_value) return
ix_start = ix_ele
if (associated(ele_start)) ix_start = tao_tracking_ele_index(ele_start, datum)

! Some inits

valid_value = .false.

datum_value = 0           ! default
datum%ix_ele_merit = -1   ! default

branch => lat%branch(ix_branch)
tao_branch => tao_lat%tao_branch(ix_branch)
orbit => tao_branch%orbit
bunch_params => tao_branch%bunch_params

n_track = branch%n_ele_track
n_max   = branch%n_ele_max

call re_allocate2 (value_vec, 0, n_track, .false.) ! Make sure is large enough if used.
call re_allocate2 (good,      0, n_track, .false.) ! Make sure is large enough if used.

ix = index(head_data_type, '.')
if (head_data_type(1:11) == 'expression:') then
  head_data_type = 'expression:'
elseif (ix /= 0) then
  sub_data_type  = head_data_type(ix+1:)
  head_data_type = head_data_type(1:ix) 
endif

select case (data_source)
case ('lat', 'beam')
  ! Valid data source
case default
  if ( head_data_type /= 'expression:') then
    call tao_set_invalid (datum, 'UNKNOWN DATA_SOURCE: ' // data_source, why_invalid, .true.)
    return
  endif
end select

if (index(head_data_type, 'stable') == 0 .and. head_data_type /= 'expression:') then
  if (associated(ele)) then
    if (orbit(ix_ele)%state /= alive$) then
      call tao_set_invalid (datum, 'UNSTABLE ORBIT AT EVALUATION POINT', why_invalid)
      return
    endif
  endif

  if (associated(ele_ref)) then
    if (orbit(ix_ref)%state /= alive$) then
      call tao_set_invalid (datum, 'UNSTABLE ORBIT AT REFERENCE EVALUATION POINT', why_invalid)
      return
    endif
  endif

  if (associated(ele_start)) then
    if (orbit(ix_start)%state /= alive$) then
      call tao_set_invalid (datum, 'UNSTABLE ORBIT AT START EVALUATION POINT', why_invalid)
      return
    endif
  endif
endif

! ele_ref must not be specified for some data types. Check this.

select case (head_data_type)
case ('wall.')
  if (datum%ele_ref_name /= '') then
    call tao_set_invalid (datum, 'SPECIFYING ELE_REF NOT VALID', why_invalid, .true.)
    return
 endif
end select


! ele_start must not be specified for some data types. Check this.

if (data_type(1:11) == 'periodic.tt' .or. data_type == 'sigma.pz') then
  if (datum%ele_start_name /= '') then
    call tao_set_invalid (datum, 'SPECIFYING ELE_START NOT VALID', why_invalid, .true.)
    return
  endif
endif

if (data_type(1:11) == 'periodic.tt') then
  if (datum%ele_ref_name /= '') then
    call tao_set_invalid (datum, 'SPECIFYING ELE_REF NOT VALID', why_invalid, .true.)
    return
  endif
endif

if (tao_branch%track_state /= moving_forward$ .and. ix_ele >= tao_branch%track_state) then
  if ((data_source == 'beam' .and. head_data_type /= 'n_particle_loss') .or. &
                         head_data_type(1:4) == 'bpm_' .or. head_data_type == 'orbit.') then
    call tao_set_invalid (datum, 'CANNOT EVALUATE DUE TO PARTICLE LOSS.', why_invalid)
    return
  endif
endif

if (data_source == 'beam' .and. .not. s%com%have_tracked_beam) then
  call tao_set_invalid (datum, 'DATA_SOURCE FOR DATUM SET TO "beam". BUT NO BEAM TRACKING HAS BEEN DONE!', why_invalid, err_level = s_warn$)
  return
endif

!-------------------------------------------------------------
! Case where evaluation point not at the end of the element.

if (head_data_type /= 'expression:' .and. (datum%s_offset /= 0 .or. datum%eval_point /= anchor_end$)) then
  if (data_source /= 'lat') then
    call tao_set_invalid (datum, 'CANNOT USE A BEAM DATA_SOURCE WITH A FINITE S_OFFSET OR EVAL_POINT = CENTER.', why_invalid, .true.)
    return
  endif

  if (datum%ele_start_name /= '') then
    call out_io (s_warn$, r_name, 'If there is an evaluation range (that is, ele_start is set), s_offset and', & 
                                  ' eval_point are ignored. For datum:' // tao_datum_name(datum))
  else
    datum_value = tao_evaluate_datum_at_s (datum, tao_lat, ele, ele_ref, valid_value, str, exterminate)
    if (.not. valid_value) call tao_set_invalid (datum, str, why_invalid, exterminate)
    return
  endif
endif

!---------------------------------------------------

select case (head_data_type)

!-----------

case ('alpha.')

  select case (data_type)

  case ('alpha.a')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%a%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%a%norm_emit /= 0)
      if (.not. valid_value) call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif
    
  case ('alpha.b')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%b%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%b%norm_emit /= 0)
      if (.not. valid_value) call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif

  case ('alpha.z')
    if (data_source == 'lat') goto 9001  ! Error message and return
    call tao_load_this_datum (bunch_params(:)%z%alpha, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID.', why_invalid, .true.)
    return

  end select

!-----------

case ('apparent_emit.', 'norm_apparent_emit.')

  select case (data_type)

  case ('apparent_emit.x', 'norm_apparent_emit.x')
    do i = ix_start, ix_ele
      if (data_source == 'lat') then
        value_vec(i) = tao_lat_emit_calc (x_plane$, apparent_emit$, branch%ele(i), tao_branch%modes_6d)
      else
        value_vec(i) = tao_beam_emit_calc (x_plane$, apparent_emit$, branch%ele(i), bunch_params(i))
      endif
    enddo

    if (ix_ref > -1) then
      if (data_source == 'lat') then
        value_vec(ix_ref) = tao_lat_emit_calc (x_plane$, apparent_emit$, branch%ele(ix_ref), tao_branch%modes_6d)
      else
        value_vec(ix_ref) = tao_beam_emit_calc (x_plane$, apparent_emit$, branch%ele(i), bunch_params(ix_ref))
      endif
    endif

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

    if (data_type == 'norm_apparent_emit.x') then
      beta_gamma = ele%value(p0c$) /  mass_of(branch%param%particle)
      datum_value = datum_value * beta_gamma
    endif

  case ('apparent_emit.y', 'norm_apparent_emit.y')
    do i = ix_start, ix_ele
      if (data_source == 'lat') then
        value_vec(i) = tao_lat_emit_calc (y_plane$, apparent_emit$, branch%ele(i), tao_branch%modes_6d)
      else
        value_vec(i) = tao_beam_emit_calc (y_plane$, apparent_emit$, branch%ele(i), bunch_params(i))
      endif
    enddo

    if (ix_ref > -1) then
      if (data_source == 'lat') then
        value_vec(ix_ref) = tao_lat_emit_calc (y_plane$, apparent_emit$, branch%ele(ix_ref), tao_branch%modes_6d)
        call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      else
        value_vec(ix_ref) = tao_beam_emit_calc (y_plane$, apparent_emit$, branch%ele(i), bunch_params(ix_ref))
        call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      endif
    endif


    if (data_type == 'norm_apparent_emit.y') then
      beta_gamma = ele%value(p0c$) /  mass_of(branch%param%particle)
      datum_value = datum_value * beta_gamma
    endif


  case default
    call tao_set_invalid (datum, 'UNKNOWN DATUM TYPE: ' // data_type, why_invalid, .true.)
    return

  end select

!-----------

case ('beta.')

  select case (data_type)

  case ('beta.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      if (bunch_params(ix_ele)%x%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    else
      call tao_set_invalid (datum, 'DATA_SOURCE: ' // trim(data_source) // ' INVALID WITH: beta.x DATA_TYPE', why_invalid, .true.)
    endif
    
  case ('beta.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      if (bunch_params(ix_ele)%y%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    else
      call tao_set_invalid (datum, 'DATA_SOURCE: ' // trim(data_source) // ' INVALID WITH: beta.y DATA_TYPE', why_invalid, .true.)
    endif

  case ('beta.z')
    if (data_source == 'lat') goto 9001  ! Error message and return
    call tao_load_this_datum (bunch_params(:)%z%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    if (bunch_params(ix_ele)%z%norm_emit == 0) then
      valid_value = .false.
      call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif

  case ('beta.a')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%a%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      if (bunch_params(ix_ele)%a%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    endif
    
  case ('beta.b')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%b%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      if (bunch_params(ix_ele)%b%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    endif

  case ('beta.c')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%z%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%c%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      if (bunch_params(ix_ele)%a%norm_emit == 0) then
        valid_value = .false.
        call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
      endif
    endif
    
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('bpm_orbit.')

  select case (data_type)
  case ('bpm_orbit.x')
    which = x_plane$
  case ('bpm_orbit.y')
    which = y_plane$
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

  if (data_source == 'beam') goto 9000  ! Set error message and return

  valid_value = .true.
  particle_lost = .false.

  do i = ix_start, ix_ele
    if (i /= ix_ele .and. .not. ele_is_monitor(branch%ele(i), .false.)) then
      value_vec(i) = 0
      cycle
    endif
    call to_orbit_reading (orbit(i), branch%ele(i), which, s%com%add_measurement_noise, value_vec(i), err)
    particle_lost = particle_lost .or. (tao_branch%track_state /= moving_forward$ .and. i > tao_branch%track_state)
    valid_value = valid_value .and. .not. err
  enddo

  if (ix_ref > -1) then
    call to_orbit_reading (orbit(ix_ref), branch%ele(ix_ref), which, s%com%add_measurement_noise, value_vec(ix_ref), err)
    particle_lost = particle_lost .or. (tao_branch%track_state /= moving_forward$ .and. ix_ref > tao_branch%track_state)
    valid_value = valid_value .and. .not. err
  endif

  if (particle_lost) then
    call tao_set_invalid (datum, 'CANNOT EVALUATE DUE TO PARTICLE LOSS', why_invalid)
    return
  elseif (.not. valid_value) then
    call tao_set_invalid (datum, 'NO VALID MONITOR ELEMENT', why_invalid, .true.)
    return
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)


!-----------

case ('bpm_eta.')

  select case (data_type)
  case ('bpm_eta.x')
    which = x_plane$
  case ('bpm_eta.y')
    which = y_plane$
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

  if (data_source == 'beam') goto 9000  ! Set error message and return
  vec2 = [ele%x%eta, ele%y%eta]
  call to_eta_reading (vec2, ele, which, s%com%add_measurement_noise, datum_value, err)
  valid_value = .not. err

!-----------

case ('bpm_phase.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value, why_invalid, datum); if (.not. valid_value) return

  select case (data_type)
  case ('bpm_phase.a')
    datum_value = bpm_data%phi_a
  case ('bpm_phase.b')
    datum_value = bpm_data%phi_b
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    valid_value = .false.
    return
  end select

!-----------

case ('bpm_k.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value, why_invalid, datum); if (.not. valid_value) return

  select case (data_type)
  case ('bpm_k.22a')
    datum_value = bpm_data%k_22a
  case ('bpm_k.12a')
    datum_value = bpm_data%k_12a
  case ('bpm_k.11b')
    datum_value = bpm_data%k_11b
  case ('bpm_k.12b')
    datum_value = bpm_data%k_12b
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    valid_value = .false.
    return
  end select

!-----------

case ('bpm_cbar.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  call tao_to_phase_and_coupling_reading (ele, bpm_data, valid_value, why_invalid, datum); if (.not. valid_value) return

  select case (data_type)
  case ('bpm_cbar.22a')
    datum_value = bpm_data%cbar22_a
  case ('bpm_cbar.12a')
    datum_value = bpm_data%cbar12_a
  case ('bpm_cbar.11b')
    datum_value = bpm_data%cbar11_b
  case ('bpm_cbar.12b')
    datum_value = bpm_data%cbar12_b
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    valid_value = .false.
    return
  end select

!-----------

case ('bunch_charge.')

  call tao_load_this_datum (bunch_params(:)%charge_live, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  if (data_type == 'bunch_charge.live_relative') then
    charge = bunch_params(ele%ix_ele)%charge_tot
    if (charge == 0) then
      call tao_set_invalid (datum, 'BUNCH HAS NO CHARGE FOR EVALUATING A DATUM OF TYPE "bunch_charge_live.percent', why_invalid)
      valid_value = .false.
      return
    endif
    datum_value = datum_value / charge

  elseif (data_type /= 'bunch_charge.live') then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    valid_value = .false.
    return
  endif

!-----------

case ('bunch_max.', 'bunch_min.')
  if (data_source /= 'beam') goto 9000  ! Set error message and return
  select case (data_type(11:))
  case ('x');  i = 1
  case ('px'); i = 2
  case ('y');  i = 3
  case ('py'); i = 4
  case ('z');  i = 5
  case ('pz'); i = 6
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select
  
  select case (data_type(1:10))
  case ('bunch_max.')
    call tao_load_this_datum (bunch_params(:)%rel_max(i), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
  case ('bunch_min.')
    call tao_load_this_datum (bunch_params(:)%rel_min(i), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  end select

!-----------

case ('c_mat.', 'cmat.')

  if (data_type(1:5) == 'c_mat') then
    data_type = 'cmat' // data_type(6:)
    call out_io (s_warn$, r_name, 'Note: "c_mat" data type is now called "cmat"')
  endif

  select case (data_type)

  case ('cmat.11')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (branch%ele(:)%c_mat(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cmat.12')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (branch%ele(:)%c_mat(1,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cmat.21')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (branch%ele(:)%c_mat(2,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cmat.22')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (branch%ele(:)%c_mat(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('cbar.')

  if (ele%a%beta == 0) then ! Can happen if the lattice is unstable.
    call tao_set_invalid (datum, 'UNSTABLE LATTICE', why_invalid)
    return
  endif

  select case (data_type)

  case ('cbar.11')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%cbar(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cbar.12')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%cbar(1,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cbar.21')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%cbar(2,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('cbar.22')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%cbar(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('chrom.')
  
  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (branch%param%geometry == open$) then
    select case (data_type)
    case ('chrom.dtune.a', 'chrom.a', 'chrom.dtune.b', 'chrom.b')
      call tao_set_invalid (datum, 'Cannot calc ' // trim(data_type) // ' with an open geometry.', why_invalid)
      return
    end select
  endif

  if (.not. tao_lat%chrom_calc_ok) then
    call tao_set_invalid (datum, 'Chrom calc failed.', why_invalid)
    return
  elseif (.not. allocated(tao_lat%low_E_lat%branch)) then
    if (branch%param%unstable_factor == 0) then
      call tao_set_invalid (datum, 'Chrom bookkeeping problem. Please contact DCS.', why_invalid)
    else
      call tao_set_invalid (datum, 'Unstable lattice.', why_invalid)
    endif
    return
  endif

  !----

  dE = 2 * s%global%delta_e_chrom  ! Actually this is the change in pz

  select case (data_type)

  case ('chrom.dtune.a', 'chrom.a')
    if (data_type == 'chrom.dtune.a') call out_io (s_warn$, r_name, '"chrom.dtune.a" IS DEPRECATED. PLEASE CHANGE TO "chrom.a".')
    datum_value = tao_branch%a%chrom
    valid_value = .true.

  case ('chrom.dtune.b', 'chrom.b')
    if (data_type == 'chrom.dtune.b') call out_io (s_warn$, r_name, '"chrom.dtune.b" IS DEPRECATED. PLEASE CHANGE TO "chrom.b".')
    datum_value = tao_branch%b%chrom
    valid_value = .true.

  case ('chrom.dbeta.a')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_lat%high_E_lat%branch(ix_branch)%ele(i)%a%beta - tao_lat%low_E_lat%branch(ix_branch)%ele(i)%a%beta) / (tao_lat%lat%ele(i)%a%beta * dE)
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.dbeta.b')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_lat%high_E_lat%branch(ix_branch)%ele(i)%b%beta - tao_lat%low_E_lat%branch(ix_branch)%ele(i)%b%beta) / (tao_lat%lat%ele(i)%b%beta * dE)
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  
  case ('chrom.dphi.a')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_lat%high_E_lat%branch(ix_branch)%ele(i)%a%phi - tao_lat%low_E_lat%branch(ix_branch)%ele(i)%a%phi)/ dE
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.dphi.b')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_lat%high_E_lat%branch(ix_branch)%ele(i)%b%phi - tao_lat%low_E_lat%branch(ix_branch)%ele(i)%b%phi)/ dE
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.deta.x')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_lat%high_E_lat%branch(ix_branch)%ele(i)%x%eta - tao_lat%low_E_lat%branch(ix_branch)%ele(i)%x%eta)/ dE
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.deta.y')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_lat%high_E_lat%branch(ix_branch)%ele(i)%y%eta - tao_lat%low_E_lat%branch(ix_branch)%ele(i)%y%eta)/ dE
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.detap.x')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_lat%high_E_lat%branch(ix_branch)%ele(i)%x%etap - tao_lat%low_E_lat%branch(ix_branch)%ele(i)%x%etap)/ dE
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.detap.y')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = (tao_lat%high_E_lat%branch(ix_branch)%ele(i)%y%etap - tao_lat%low_E_lat%branch(ix_branch)%ele(i)%y%etap)/ dE
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('chrom.w.a', 'chrom.w.b')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        if (data_type == 'chrom.w.a') then
          z2 => tao_lat%high_E_lat%branch(ix_branch)%ele(i)%a
          z1 => tao_lat%low_E_lat%branch(ix_branch)%ele(i)%a
          z0 => branch%ele(i)%a
        else
          z2 => tao_lat%high_E_lat%branch(ix_branch)%ele(i)%b
          z1 => tao_lat%low_E_lat%branch(ix_branch)%ele(i)%b
          z0 => branch%ele(i)%b
        endif
        dalpha = (z2%alpha - z1%alpha) / dE
        dbeta  = (z2%beta - z1%beta) / dE
        aa = dalpha - z0%alpha * dbeta / z0%beta
        bb = dbeta / z0%beta
        value_vec(i) = sqrt(aa**2 + bb**2)
      end do
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif      

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('chrom_ptc.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  ptc_nf => tao_branch%ptc_normal_form

  if (.not. ptc_nf%valid_map) then
    if (.not. u%calc%one_turn_map) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE ONE_TURN_MAP_CALC IS NOT SET TO TRUE.', why_invalid)
    elseif (branch%param%geometry /= closed$) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE LATTICE GEOMETRY IS NOT CLOSED.', why_invalid)
    else
      call tao_set_invalid (datum, '?????', why_invalid)
    endif
    return
  endif

  select case (data_type(1:12))
  case ('chrom_ptc.a.')
    phase_map => ptc_nf%phase(1)
  case ('chrom_ptc.b.')
    phase_map => ptc_nf%phase(2)
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

  if (.not. is_integer(data_type(13:), n)) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  expo = 0
  expo(6) = n 

  datum_value = real(phase_map .sub. expo)
  valid_value = .true.

!-----------

case ('curly_h.')

  if (data_source == 'beam')  goto 9000  ! Set error message and return

  select case (data_type)
  case ('curly_h.a')
    do i = ix_start, ix_ele
      ele => branch%ele(i)
      value_vec(i) = ele%a%gamma * ele%a%eta**2 + 2 * ele%a%alpha * ele%a%eta * ele%a%etap + ele%a%beta * ele%a%etap**2
    enddo
    if (ix_ref > -1) then
      ele => branch%ele(ix_ref)
      value_vec(ix_ref) = ele%a%gamma * ele%a%eta**2 + 2 * ele%a%alpha * ele%a%eta * ele%a%etap + ele%a%beta * ele%a%etap**2
    endif
    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('curly_h.b')
    do i = ix_start, ix_ele
      ele => branch%ele(i)
      value_vec(i) = ele%b%gamma * ele%b%eta**2 + 2 * ele%b%alpha * ele%b%eta * ele%b%etap + ele%b%beta * ele%b%etap**2
    enddo
    if (ix_ref > -1) then
      ele => branch%ele(ix_ref)
      value_vec(ix_ref) = ele%b%gamma * ele%b%eta**2 + 2 * ele%b%alpha * ele%b%eta * ele%b%etap + ele%b%beta * ele%b%etap**2
    endif
    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  end select

!-----------

case ('damp.')

  if (data_source == 'beam') goto 9000  ! Set error message and return

  select case (data_type)

  case ('damp.j_a')
    datum_value = tao_branch%modes_6d%a%j_damp
    valid_value = .true.

  case ('damp.j_b')
    datum_value = tao_branch%modes_6d%b%j_damp
    valid_value = .true.

  case ('damp.j_z')
    datum_value = tao_branch%modes_6d%z%j_damp
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('deta_ds.')

  select case (data_type)

  case ('deta_ds.a')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%a%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('deta_ds.b')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%b%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('deta_ds.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%x%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('deta_ds.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%y%deta_ds, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('dpx_dx') 
  if (data_source == 'lat') then
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) = tao_branch%lat_sigma(ix_ref)%mat(1,2) / tao_branch%lat_sigma(ix_ref)%mat(1,1)
      value_vec(ix_ele) = tao_branch%lat_sigma(ix_ele)%mat(1,2) / tao_branch%lat_sigma(ix_ele)%mat(1,1)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (tao_branch%lat_sigma%mat(1,2) / tao_branch%lat_sigma%mat(1,1), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif
  else
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) = bunch_params(ix_ref)%sigma(1,2) / bunch_params(ix_ref)%sigma(1,1)
      value_vec(ix_ele) = bunch_params(ix_ele)%sigma(1,2) / bunch_params(ix_ele)%sigma(1,1)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    else
      call tao_load_this_datum (bunch_params%sigma(1,2) / bunch_params%sigma(1,1), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif
  endif

case ('dpy_dy') 
  if (data_source == 'lat') then
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  tao_branch%lat_sigma(ix_ref)%mat(3,4) / tao_branch%lat_sigma(ix_ref)%mat(3,3)
      value_vec(ix_ele) = tao_branch%lat_sigma(ix_ele)%mat(3,4) / tao_branch%lat_sigma(ix_ele)%mat(3,3)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (tao_branch%lat_sigma%mat(3,4) / tao_branch%lat_sigma%mat(3,3), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  else
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  bunch_params(ix_ref)%sigma(3,4) / bunch_params(ix_ref)%sigma(3,3)
      value_vec(ix_ele) = bunch_params(ix_ele)%sigma(3,4) / bunch_params(ix_ele)%sigma(3,3)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    else
      call tao_load_this_datum (bunch_params%sigma(3,4) / bunch_params%sigma(3,3), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif
  endif

case ('dpz_dz') 
  if (data_source == 'lat') then
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  tao_branch%lat_sigma(ix_ref)%mat(5,6) / tao_branch%lat_sigma(ix_ref)%mat(5,5)
      value_vec(ix_ele) = tao_branch%lat_sigma(ix_ele)%mat(5,6) / tao_branch%lat_sigma(ix_ele)%mat(5,5)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (tao_branch%lat_sigma%mat(5,6) / tao_branch%lat_sigma%mat(5,5), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif
  else
    if (ix_start == ix_ele) then
      if (ix_ref > -1) value_vec(ix_ref) =  bunch_params(ix_ref)%sigma(5,6) / bunch_params(ix_ref)%sigma(5,5)
      value_vec(ix_ele) = bunch_params(ix_ele)%sigma(5,6) / bunch_params(ix_ele)%sigma(5,5)
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    else
      call tao_load_this_datum (bunch_params%sigma(5,6) / bunch_params%sigma(5,5), &
                          ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif
  endif

!-----------

case ('dynamic_aperture.')
  da => u%dynamic_aperture
  if (allocated(da%scan)) then
    if (.not. allocated(da%scan(1)%point)) then
      call tao_set_invalid (datum, 'DYNAMIC APERTURE NOT CALCULATED', why_invalid)
      return
    endif
  else
    call tao_set_invalid (datum, 'DYNAMIC APERTURE NOT CALCULATED', why_invalid)
    return
  endif

  if (da%a_emit <= 0 .or. da%b_emit <= 0) then
    call tao_set_invalid (datum, 'A_EMIT OR B_EMIT NOT SET IN TAO_DYNAMIC_APERTURE STRUCTURE.', why_invalid)
    return
  endif

  n_da = size(da%scan)

  if (.not. is_integer(sub_data_type, n)) then
    call tao_set_invalid (datum, 'MALFORMED DATA_TYPE: ' // quote(data_type) // '. ' // quote(sub_data_type) // ' IS NOT AN INTEGER.', why_invalid, .true.)
    return
  endif

  if (n < 1 .or. n > n_da) then
    call tao_set_invalid (datum, 'SCAN INDEX OUT OF RANGE FOR DATA_TYPE: ' // quote(data_type), why_invalid, .true.)
    return
  endif

  scan => da%scan(n)
  if (da%param%start_ele == '') then
    ele => lat%ele(0)
  else
    call lat_ele_locator (da%param%start_ele, lat, eles, n)
    ele => eles(1)%ele
  endif

  datum_value = 1d100   ! Something large
  do j = 1, size(scan%point)
    orb1 = scan%ref_orb
    orb1%vec(1:4) = [scan%point(j)%x, 0.0_rp, scan%point(j)%y, 0.0_rp]
    call orbit_amplitude_calc (ele, orb1, amp_a, amp_b)
    amp = sqrt(2 * (amp_a / da%a_emit + amp_b / da%b_emit))
    datum_value = min(datum_value, amp)
  enddo

  valid_value = .true.

!-----------

case ('e_tot_ref')
  if (data_source == 'beam') goto 9000  ! Set error message and return
  call tao_load_this_datum (branch%ele(:)%value(e_tot$), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('element_attrib.')

  name = upcase(data_type(16:))
  value_vec = 0
  good = .false.

  do i = ix_start, ix_ele
    call pointer_to_attribute (branch%ele(i), name, .false., a_ptr, err, .false.)
    if (.not. associated (a_ptr%r)) cycle
    value_vec(i) = a_ptr%r
    good(i) = .true.
  enddo

  if (all(.not. good(ix_start:ix_ele))) then
    call tao_set_invalid (datum, 'CANNOT EVALUATE DATUM WITH DATA_TYPE = "' // trim(data_type) // '" AT ASSOCIATED ELEMENT', why_invalid, .true.)
    return
  endif

  if (ix_ref > -1) then
    call pointer_to_attribute (ele_ref, name, .false., a_ptr, err, .false.)
    if (associated (a_ptr%r)) then
      value_vec(ix_ref) = a_ptr%r
      good(ix_ref) = .true.
    else
      if (.not. good(ix_ref)) then
        call tao_set_invalid (datum, 'CANNOT EVALUATE DATUM WITH DATA_TYPE = "' // trim(data_type) // '" AT REFERENCE ELEMENT', why_invalid, .true.)
        return
      endif
    endif
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, good)

!-----------

case ('emit.', 'norm_emit.')

  if (associated(ele)) then
    beta_gamma = ele%value(p0c$) / mass_of(branch%param%particle)
  else
    beta_gamma = 0
  endif

  select case (data_type)

  case ('emit.x', 'norm_emit.x')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = tao_lat_emit_calc (x_plane$, projected_emit$, branch%ele(i), tao_branch%modes_6d)
      enddo
      if (ix_ref > -1) then
        value_vec(ix_ref) = tao_lat_emit_calc (x_plane$, projected_emit$, branch%ele(ix_ref), tao_branch%modes_6d)
      endif
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.x') datum_value = datum_value * beta_gamma

  case ('emit.y', 'norm_emit.y')
    if (data_source == 'lat') then
      do i = ix_start, ix_ele
        value_vec(i) = tao_lat_emit_calc (y_plane$, projected_emit$, branch%ele(i), tao_branch%modes_6d)
      enddo
      if (ix_ref > -1) then
        value_vec(ix_ref) = tao_lat_emit_calc (y_plane$, projected_emit$, branch%ele(ix_ref), tao_branch%modes_6d)
      endif
      call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.y') datum_value = datum_value * beta_gamma

  case ('emit.z', 'norm_emit.z')
    if (data_source == 'lat') then
      goto 9001   ! Error message and return
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%z%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.y') datum_value = datum_value * beta_gamma

  case ('emit.a', 'norm_emit.a')
    if (data_source == 'lat') then
      if (lat%param%geometry == open$ .and. ix_ele > -1) then
        if (.not. allocated(tao_lat%rad_int%branch)) then
          call out_io (s_error$, r_name, 'tao_lat%rad_int%branch not allocated')
          return
        endif
        rad_int_branch => tao_lat%rad_int%branch(ix_branch)
        call tao_load_this_datum (rad_int_branch%ele%lin_norm_emit_a, &
                                ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
        datum_value = datum_value / beta_gamma
      else
        datum_value = tao_branch%modes_6d%a%emittance
        valid_value = .true.
      endif
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.a') datum_value = datum_value * beta_gamma
    
  case ('emit.b', 'norm_emit.b')
    if (data_source == 'lat') then
      if (lat%param%geometry == open$ .and. ix_ele > -1) then
        if (.not. allocated(tao_lat%rad_int%branch)) then
          call out_io (s_error$, r_name, 'tao_lat%rad_int%branch not allocated')
          return
        endif
        rad_int_branch => tao_lat%rad_int%branch(ix_branch)
        call tao_load_this_datum (rad_int_branch%ele%lin_norm_emit_b, &
                                ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
        datum_value = datum_value / beta_gamma
      else
        datum_value = tao_branch%modes_6d%b%emittance
        valid_value = .true.
      endif
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.b') datum_value = datum_value * beta_gamma

  case ('emit.c', 'norm_emit.c')
    if (data_source == 'lat') then
      goto 9001     ! Error message and return
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%c%emit, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

    if (data_type == 'norm_emit.y') datum_value = datum_value * beta_gamma

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('eta.')

  select case (data_type)

  case ('eta.a')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%a%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.b')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%b%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%x%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%y%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('eta.z')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%z%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%z%eta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('etap.')

  select case (data_type)

  case ('etap.a')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%a%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('etap.b')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%b%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('etap.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%x%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%x%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('etap.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%y%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = .true.
    else
      call tao_load_this_datum (branch%ele(:)%y%etap, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('expression:', 'expression.')

  write (dflt_dat_index, '(i0)') datum%ix_d1
  e_str = datum%data_type(12:)
  do
    ix = index(e_str, 'ele::#[')
    if (ix == 0) exit
    if (ix_ele == -1) then
      call tao_set_invalid (datum, 'NO ASSOCIATED ELEMENT' // e_str, why_invalid, .true.)
      return
    endif
    e_str = e_str(1:ix+4) // trim(ele_loc_name(ele)) // e_str(ix+6:)
  enddo

  printit = (s%com%n_err_messages_printed < s%global%datum_err_messages_max) 
  call tao_evaluate_expression (e_str, 0, .false., expression_value_vec, err, printit, info, &
                  datum%stack, tao_lat%name, datum%data_source, ele_ref, ele_start, ele, &
                  dflt_dat_index, u%ix_uni, datum%eval_point, datum%s_offset, datum = datum)
  if (err) then
    call tao_set_invalid (datum, 'CANNOT EVALUATE EXPRESSION: ' // e_str, why_invalid)
    return
  endif

  select case (datum%merit_type)
  case ('min');      datum_value = minval(expression_value_vec)
  case ('max');      datum_value = maxval(expression_value_vec)
  case ('abs_min');  datum_value = minval(abs(expression_value_vec))
  case ('abs_max');  datum_value = maxval(abs(expression_value_vec))
  case ('max-min');  datum_value = maxval(expression_value_vec) - minval(expression_value_vec)

  case ('integral', 'average', 'rms')
    s_offset = 0
    if (.not. associated(ele_start)) then
      call tao_set_invalid (datum, 'ELE_START NOT SET. THIS IS NEEDED WHEN MERIT_TYPE IS SET TO "integral", "average", OR "rms".', why_invalid)
      return
    endif

    do i = 1, size(info)
      j = i + ele_start%ix_ele - 1
      if (j > branch%n_ele_track) then
        j = j - branch%n_ele_track - 1
        s_offset = branch%param%total_length
      endif
      info(i)%s = tao_datum_s_position(datum, branch%ele(j)) + s_offset
    enddo

    if (j /= ele%ix_ele) then
      call out_io (s_error$, r_name, 'BOOKKEEPING ERROR IN EVALUATING INTEGRAL/AVERAGE OF EXPRESSION.', &
                                     'PLEASE REPORT!')
      call tao_set_invalid (datum, 'CANNOT EVALUATE EXPRESSION: ' // datum%data_type, why_invalid)
    endif

    datum_value = tao_datum_integrate(datum, branch, info(:)%s, expression_value_vec, valid_value, why_invalid)
    return

  case ('target')
    if (size(expression_value_vec) /= 1) then
      call tao_set_invalid (datum, 'MERIT_TYPE IS SET TO "TARGET" BUT DATUM DOES NOT EVALUATE TO A SINGLE NUMBER!', why_invalid, .true.)
      return
    endif
    datum_value = expression_value_vec(1)
  case default
    call out_io (s_error$, r_name, &
                'SINCE THIS DATUM: ' // tao_datum_name(datum), &
                'SPECIFIES A RANGE OF ELEMENTS, THEN THIS MERIT_TYPE: ' // datum%merit_type, &
                'IS NOT VALID. VALID MERIT_TYPES ARE MIN, MAX, ABS_MIN, AND ABS_MAX.')
    call tao_set_invalid (datum, 'MERIT_TYPE: ' // quote(datum%merit_type) // ' IS NOT VALID WHEN THERE IS AN EVALUATION RANGE', why_invalid, .true.)
    return
  end select

  ! Make sure that any datums used in the expression have already been evaluated.
  do i = 1, size(datum%stack)
    if (datum%stack(i)%type /= numeric$) cycle
    call tao_find_data (err, datum%stack(i)%name, d_array = d_array, print_err = .false.)
    if (err .or. size(d_array) == 0) cycle  ! Err -> This is not associated then not a datum.
    dp => d_array(1)%d
    if (dp%d1%d2%ix_universe < u%ix_uni) cycle ! OK
    if (dp%d1%d2%ix_universe == u%ix_uni .and. dp%ix_data < datum%ix_data) cycle
    call out_io (s_error$, r_name, 'DATUM: ' // tao_datum_name(datum), &
                    'WHICH IS OF TYPE EXPRESSION:' // datum%data_type, &
                    'THE EXPRESSION HAS A COMPONENT: ' // datum%stack(i)%name, &
                    'AND THIS COMPONENT IS EVALUATED AFTER THE EXPRESSION!', &
                    'TO FIX: MOVE THE EXPRESSION DATUM TO BE AFTER THE COMPONENT DATUM IN THE FILE THAT DEFINES THE DATA.')
    return
  enddo
  valid_value = .true.

!-----------

case ('floor.')

  select case (data_type)

  case ('floor.x')
    call tao_load_this_datum (branch%ele(:)%floor%r(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.y')
    call tao_load_this_datum (branch%ele(:)%floor%r(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.z')
    call tao_load_this_datum (branch%ele(:)%floor%r(3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.theta')
    call tao_load_this_datum (branch%ele(:)%floor%theta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.phi')
    call tao_load_this_datum (branch%ele(:)%floor%phi, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('floor.psi')
    call tao_load_this_datum (branch%ele(:)%floor%psi, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('floor_actual.')

  value_vec(ix_ele) = tao_ele_geometry_with_misalignments(datum, ele, valid_value, why_invalid)
  if (.not. valid_value) return
  if (associated(ele_ref)) value_vec(ix_ref) = tao_ele_geometry_with_misalignments(datum, ele_ref, valid_value, why_invalid)
  if (.not. valid_value) return

  if (associated(ele_start)) then
    do ie = ix_start, ix_ele - 1
      value_vec(ie) = tao_ele_geometry_with_misalignments(datum, branch%ele(ie), valid_value, why_invalid)
      if (.not. valid_value) return
    enddo
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('floor_orbit.')

  value_vec(ix_ele) = tao_eval_floor_orbit (datum, ele, orbit(ix_ele), bunch_params(ix_ele), valid_value, why_invalid)
  if (.not. valid_value) return
  if (associated(ele_ref)) value_vec(ix_ref) = tao_eval_floor_orbit (datum, ele_ref, orbit(ix_ref), bunch_params(ix_ref), valid_value, why_invalid)
  if (.not. valid_value) return

  if (associated(ele_start)) then
    do ie = ix_start, ix_ele - 1
      value_vec(ie) = tao_eval_floor_orbit (datum, branch%ele(ie), orbit(ie), bunch_params(ie), valid_value, why_invalid)
      if (.not. valid_value) return
    enddo
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('gamma.')

  select case (data_type)

  case ('gamma.a')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%a%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%a%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%a%norm_emit /= 0)
      call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif
    
  case ('gamma.b')
    if (data_source == 'lat') then
      call tao_load_this_datum (branch%ele(:)%b%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    elseif (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%b%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      valid_value = valid_value .and. (bunch_params(ix_ele)%b%norm_emit /= 0)
      call tao_set_invalid (datum, 'CANNOT EVALUATE SINCE THE EMITTANCE IS ZERO!', why_invalid)
    endif

  case ('gamma.z')
    if (data_source == 'lat') goto 9001  ! Error message and return
    call tao_load_this_datum (bunch_params(:)%z%gamma, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('k.')

  select case (data_type)

  case ('k.11b')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%k_11a, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  case ('k.12a')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%k_12a, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  case ('k.12b')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%k_12b, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  case ('k.22a')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%k_22b, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('momentum')
  if (data_source == 'beam') goto 9000  ! Set error message and return
  call tao_load_this_datum (branch%ele(0:n_track)%value(p0c$) * (1+orbit(0:n_track)%vec(6)), &
                            ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('momentum_compaction')
  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (ix_ref < 0) then
    ix_ref = 0
    ele_ref => lat%branch(branch%ix_branch)%ele(ix_ref)
  endif

  g2 = (mass_of(ele_ref%ref_species) / ele_ref%value(E_tot$))**2   ! 1/gamma^2

  orb0 => orbit(ix_ref)
  call make_v_mats (ele_ref, v_mat, v_inv_mat)
  eta_vec = [ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap]
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb0%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb0%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb0%vec(4) / one_pz

  call transfer_matrix_calc (lat, mat6, vec0, ix_ref, ix_start, branch%ix_branch)

  do i = ix_start, ix_ele
    s_len = branch%ele(i)%s - branch%ele(ix_ref)%s
    if (s_len == 0) then
      value_vec(i) = 0
    else
      value_vec(i) = g2 - (sum(mat6(5,1:4) * eta_vec) + mat6(5,6)) / s_len
    endif
    if (i /= ix_ele) mat6 = matmul(branch%ele(i+1)%mat6, mat6)
  enddo
  call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('momentum_compaction_ptc.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  ptc_nf => tao_branch%ptc_normal_form

  if (.not. ptc_nf%valid_map) then
    if (.not. u%calc%one_turn_map) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE ONE_TURN_MAP_CALC IS NOT SET TO TRUE.', why_invalid)
    elseif (branch%param%geometry /= closed$) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE LATTICE GEOMETRY IS NOT CLOSED.', why_invalid)
    else
      call tao_set_invalid (datum, '?????', why_invalid)
    endif
    return
  endif

  if (.not. is_integer(data_type(25:), n)) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  expo = 0
  expo(6) = n 

  datum_value = real(ptc_nf%path_length .sub. expo) / branch%param%total_length
  valid_value = .true.

!-----------

case ('n_particle_loss')
  if (data_source /= 'beam') goto 9001  ! Set error message and return
  if (ix_ele < 0) ix_ele = branch%n_ele_track
  datum_value = sum(bunch_params(ix_ref+1:ix_ele)%n_particle_lost_in_ele)
  valid_value = .true.

!--------------

case ('normal.')

  ! Fetches normal_form components.
  if (data_source == 'beam') goto 9000  ! Set error message and return

  ptc_nf => tao_branch%ptc_normal_form
  bmad_nf => tao_branch%bmad_normal_form

  ! Do nothing it the map wasn't made
  if (.not. ptc_nf%valid_map) then
    call tao_set_invalid (datum, 'PTC one-turn map not calculated.', why_invalid)
    return
  endif


  if (.not. associated(bmad_nf%ele_origin)) then
    ! Get resonant driving terms
    call normal_form_rd_terms(ptc_nf%one_turn_map, bmad_nf, rf_is_on(branch))
    bmad_nf%ele_origin => ptc_nf%ele_origin
  endif

  ! Expect: taylor.#.######
  ! Example: normal.dhdj.2.000001 is the b-mode chromaticity
  !          head   sub
  ! Get position of first number. 
  iz = index(sub_data_type, '.') + 1
  
  if(sub_data_type(1:2) == 'h.') then
    valid_value = .true.
    term_found = .false.
    do i=1, size(bmad_nf%h)
      if(sub_data_type(3:8) == bmad_nf%h(i)%id) then
        temp_cplx = bmad_nf%h(i)%c_val
        term_found = .true.
      endif
    enddo
    if(term_found) then
      select case (sub_data_type(10:10))
      case('r')
        datum_value = real(temp_cplx)
      case('i')
        datum_value = aimag(temp_cplx)
      case('a')
        datum_value = abs(temp_cplx)
      case default
        call tao_set_invalid (datum, 'Data_type not ending in .r, .i, or .a.', why_invalid, .true.)
        valid_value = .false.
        return
      end select
    else
      call tao_set_invalid (datum, 'Data_type not found in normal_form_struct', why_invalid, .true.)
      valid_value = .false.
      return
    endif
  else
    i = tao_read_phase_space_index (sub_data_type, iz, .false.)
    if (i == 0) then
      call tao_set_invalid (datum, 'Bad phase space index.', why_invalid, .true.)
      return
    endif
    ! Point to taylor
    taylor_is_complex = .false.
    if (sub_data_type(1:5) == 'dhdj.') then
      taylor_ptr => bmad_nf%dhdj(i)
    else if (sub_data_type(1:2) == 'A.') then
      taylor_ptr => bmad_nf%A(i)
    else if (sub_data_type(1:6) == 'A_inv.') then
      taylor_ptr => bmad_nf%A_inv(i)
    else if (sub_data_type(1:2) == 'M.') then
      taylor_ptr => bmad_nf%M(i)
    else if (sub_data_type(1:4) == 'ReF.') then
      taylor_is_complex = .true.
      use_real_part = .true.
      complex_taylor_ptr => bmad_nf%f(i)
    else if (sub_data_type(1:4) == 'ImF.') then
      taylor_is_complex = .true.
      use_real_part = .false.
      complex_taylor_ptr => bmad_nf%F(i)
    else if (sub_data_type(1:4) == 'ReL.') then
      taylor_is_complex = .true.
      use_real_part = .true.
      complex_taylor_ptr => bmad_nf%L(i)
    else if (sub_data_type(1:4) == 'ImL.') then
      taylor_is_complex = .true.
      use_real_part = .false.
      complex_taylor_ptr => bmad_nf%L(i)
    endif
   
    ! Check for second dot
    if (sub_data_type(iz+1:iz+1) /= '.') then
     call tao_set_invalid (datum, 'Missing dot "." in data_type', why_invalid, .true.)
     call out_io (s_error$, r_name, 'data_type: '//trim(data_type) )
     call out_io (s_error$, r_name, 'expect dot: ', sub_data_type(1:iz)//'.######' )
    endif
   
    ! Get exponent
    expn_str = sub_data_type(iz+2:iz+7)
    expnt = 0
    do j = 1, 6
      if (expn_str(j:j) == ' ') exit
      expnt(j) = index('0123456789', expn_str(j:j)) - 1
    enddo
    
    ! Coefficient
    if (taylor_is_complex) then
      if (use_real_part) then
        datum_value = real(complex_taylor_coef(complex_taylor_ptr, expnt))
      else
        datum_value = aimag(complex_taylor_coef(complex_taylor_ptr, expnt))
      endif
    else
      datum_value = taylor_coef(taylor_ptr, expnt)
    endif
    valid_value = .true.  
  endif

!-----------

case ('orbit.')

  if (tao_branch%track_state /= moving_forward$ .and. ix_ele > tao_branch%track_state) then
    valid_value = .false.
    call tao_set_invalid (datum, 'Particle lost.', why_invalid)
    return
  endif

  select case (data_type)

  case ('orbit.energy', 'orbit.e_tot', 'orbit.kinetic')  ! orbit.e_tot is old style
    if (ix_ref > -1) then
      if (data_source == 'beam') then
        orb => bunch_params(ix_ref)%centroid
      else
        orb => orbit(ix_ref)
      endif
      if (orb%state == not_set$) goto 7000  ! Set error message and return
      value_vec(ix_ref) = (1 + orb%vec(6)) * orb%p0c / orb%beta
    endif

    do i = ix_start, ix_ele
      if (data_source == 'beam') then
        orb => bunch_params(i)%centroid
      else
        orb => orbit(i)
      endif
      if (orb%state == not_set$) goto 7000  ! Set error message and return
      call convert_pc_to ((1 + orb%vec(6))*orb%p0c, orb%species, e_tot = value_vec(i))
    enddo

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

    if (data_type == 'orbit.kinetic') datum_value = datum_value - mass_of(orb%species)

  case ('orbit.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.z')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.px')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.py')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.pz')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('orbit.amp_a')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%amp_a, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)

  case ('orbit.amp_b')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%amp_b, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)

  case ('orbit.norm_amp_a')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%amp_na, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)

  case ('orbit.norm_amp_b')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    call tao_load_this_datum (scratch%cc%amp_nb, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('pc')
  if (ix_ref > -1) then
    if (data_source == 'beam') then
      value_vec(ix_ref) = (1 + bunch_params(ix_ref)%centroid%vec(6)) * bunch_params(ix_ref)%centroid%p0c
    else
      value_vec(ix_ref) = (1 + orbit(ix_ref)%vec(6)) * orbit(ix_ref)%p0c
    endif
  endif

  if (data_source == 'beam') then
    do i = ix_start, ix_ele
      value_vec(i) = (1 + bunch_params(i)%centroid%vec(6)) * bunch_params(ix_ref)%centroid%p0c
    enddo
  else
    do i = ix_start, ix_ele
      value_vec(i) = (1 + orbit(i)%vec(6)) * orbit(ix_ref)%p0c
    enddo
  endif

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('periodic.')

  ix = index(data_type(10:), '.') + 9
  select case (data_type(1:ix))

  case ('periodic.tt.')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    if (lat%param%geometry /= closed$ .and. .not. associated(ele_ref)) then
      call tao_set_invalid (datum, 'LATTICE MUST BE CIRCULAR FOR A DATUM LIKE: ' // data_type, why_invalid)
      call err_exit
    endif

    call transfer_map_calc (lat, taylor, err, ix_ele, ix_ele, orbit(ix_ele), branch%ix_branch, &
                                                       one_turn = .true., concat_if_possible = s%global%concatenate_maps)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT INVERT MAP', why_invalid)
      return
    endif

    do i = 1, 4
      call add_taylor_term (taylor(i), -1.0_rp, taylor_expn([i]))
    enddo
    call taylor_inverse (taylor, taylor, err)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT INVERT MAP', why_invalid)
      return
    endif

    expnt = 0
    i = tao_read_phase_space_index (data_type, 13, .false.)
    if (i == 0) then
      call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
      return
    endif

    do j = 14, 24
      if (data_type(j:j) == ' ') exit
      k = tao_read_phase_space_index (data_type, j, .false.)
      if (k == 0) then
        call tao_set_invalid (datum, 'BAD DATA_TYPE = "' // trim(data_type), why_invalid, .true.)
        return
      endif
      expnt(k) = expnt(k) + 1
    enddo

    datum_value = taylor_coef (taylor(i), expnt)
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('phase.', 'phase_frac.')

  select case (data_type)

  case ('phase.a', 'phase_frac.a')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    if (ix_ref < 0) then
      datum_value = ele%a%phi
    else
      datum_value = ele%a%phi - ele_ref%a%phi
      if (ix_ref > ix_ele .and. branch%param%geometry == closed$) then
        dphi = branch%ele(n_track)%a%phi - branch%ele(0)%a%phi
        if (2*datum_value < -dphi) datum_value = datum_value + dphi
      endif
    if (data_type == 'phase_frac.a') datum_value = modulo2(datum_value, pi)
    endif
    valid_value = .true.

  case ('phase.b', 'phase_frac.b')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    if (ix_ref < 0) then
      datum_value = ele%b%phi
    else
      datum_value = ele%b%phi - ele_ref%b%phi
      if (ix_ref > ix_ele .and. branch%param%geometry == closed$) then
        dphi = branch%ele(n_track)%b%phi - branch%ele(0)%b%phi
        if (2*datum_value < -dphi) datum_value = datum_value + dphi
      endif
    endif
    if (data_type == 'phase_frac.b') datum_value = modulo2(datum_value, pi)
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('phase_frac_diff')
  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (ix_ref < 0) then
    px = ele%a%phi 
    py = ele%b%phi 
  else
    px = ele%a%phi - ele_ref%a%phi
    py = ele%b%phi - ele_ref%b%phi
    if (ix_ref > ix_ele .and. branch%param%geometry == closed$) then
      dphi = branch%ele(n_track)%a%phi - branch%ele(0)%a%phi
      if (2*px < -dphi) px = px + dphi
      dphi = branch%ele(n_track)%b%phi - branch%ele(0)%b%phi
      if (2*px < -dphi) py = py + dphi
    endif
  endif

  datum_value = modulo2 (px - py, pi)
  valid_value = .true.

!-----------

case ('photon.')

  select case (data_type)

  case ('photon.intensity_x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%centroid%field(1)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%field(1)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('photon.intensity_y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%centroid%field(2)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%field(2)**2, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('photon.intensity')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params(:)%centroid%field(1)**2+bunch_params(:)%centroid%field(2)**2, &
                                                                      ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%field(1)**2 + orbit(:)%field(2)**2, ele_ref, ele_start, ele, &
                                                                           datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('photon.phase_x')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (orbit(:)%phase(1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('photon.phase_y')
    if (data_source == 'beam') goto 9000  ! Set error message and return
    call tao_load_this_datum (orbit(:)%phase(2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  end select

!-----------

case ('ping_a.')
  select case (data_type)
  case ('ping_a.amp_x')
    datum_value = ele%gamma_c * sqrt(ele%a%beta)
    if (associated(ele_ref)) datum_value = datum_value - ele%gamma_c * sqrt(ele_ref%a%beta)
    valid_value = .true.

  case ('ping_a.phase_x')
    datum_value = ele%a%phi
    if (associated(ele_ref)) datum_value = datum_value - ele_ref%a%phi
    valid_value = .true.

  case ('ping_a.amp_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = sqrt(ele%b%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(2,2)**2))
    if (associated(ele_ref)) then
      datum_value = datum_value - sqrt(ele_ref%b%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(2,2)**2))
    endif
    valid_value = .true.

  case ('ping_a.phase_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = ele%a%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), - scratch%cc(ix_ele)%cbar(2,2))
    if (associated(ele_ref)) then
      datum_value = datum_value - ele_ref%a%phi - atan2(scratch%cc(ix_ref)%cbar(1,2), - scratch%cc(ix_ref)%cbar(2,2))
    endif
    valid_value = .true.

  case ('ping_a.amp_sin_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    amp = sqrt(ele%b%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(2,2)**2))
    phase = ele%a%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), -scratch%cc(ix_ele)%cbar(2,2))
    datum_value = amp * sin(phase)
    if (associated(ele_ref)) then
      amp = sqrt(ele_ref%b%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(2,2)**2))
      phase = ele_ref%a%phi + atan2(scratch%cc(ix_ref)%cbar(1,2), -scratch%cc(ix_ref)%cbar(2,2))
      datum_value = datum_value - amp * sin(phase)
    endif
    valid_value = .true.

  case ('ping_a.amp_cos_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    amp = sqrt(ele%b%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(2,2)**2))
    phase = ele%a%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), -scratch%cc(ix_ele)%cbar(2,2))
    datum_value = amp * cos(phase)
    if (associated(ele_ref)) then
      amp = sqrt(ele_ref%b%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(2,2)**2))
      phase = ele_ref%a%phi + atan2(scratch%cc(ix_ref)%cbar(1,2), -scratch%cc(ix_ref)%cbar(2,2))
      datum_value = datum_value - amp * cos(phase)
    endif
    valid_value = .true.

  case ('ping_a.amp_sin_rel_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = -sqrt(ele%b%beta) * scratch%cc(ix_ele)%cbar(1,2)

    if (associated(ele_ref)) then
      datum_value = datum_value + sqrt(ele_ref%b%beta) * scratch%cc(ix_ref)%cbar(1,2)
    endif
    valid_value = .true.

  case ('ping_a.amp_cos_rel_y')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = -sqrt(ele%b%beta) * scratch%cc(ix_ele)%cbar(2,2)

    if (associated(ele_ref)) then
      datum_value = datum_value + sqrt(ele_ref%b%beta) * scratch%cc(ix_ref)%cbar(2,2)
    endif
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  end select

!-----------

case ('ping_b.')
  select case (data_type)
  case ('ping_b.amp_y')
    datum_value = ele%gamma_c * sqrt(ele%b%beta)
    if (associated(ele_ref)) datum_value = datum_value - ele%gamma_c * sqrt(ele_ref%b%beta)
    valid_value = .true.

  case ('ping_b.phase_y')
    datum_value = ele%b%phi
    if (associated(ele_ref)) datum_value = datum_value - ele_ref%b%phi
    valid_value = .true.

  case ('ping_b.amp_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = sqrt(ele%a%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(1,1)**2))
    if (associated(ele_ref)) then
      datum_value = datum_value - sqrt(ele_ref%a%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(1,1)**2))
    endif
    valid_value = .true.

  case ('ping_b.phase_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = ele%b%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), scratch%cc(ix_ele)%cbar(1,1))
    if (associated(ele_ref)) then
      datum_value = datum_value - ele_ref%b%phi - atan2(scratch%cc(ix_ref)%cbar(1,2), scratch%cc(ix_ref)%cbar(1,1))
    endif
    valid_value = .true.

  case ('ping_b.amp_sin_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    amp = sqrt(ele%a%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(1,1)**2))
    phase = ele%b%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), scratch%cc(ix_ele)%cbar(1,1))
    datum_value = amp * sin(phase)
    if (associated(ele_ref)) then
      amp = sqrt(ele_ref%a%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(1,1)**2))
      phase = ele_ref%b%phi + atan2(scratch%cc(ix_ref)%cbar(1,2), scratch%cc(ix_ref)%cbar(1,1))
      datum_value = datum_value - amp * sin(phase)
    endif
    valid_value = .true.

  case ('ping_b.amp_cos_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    amp = sqrt(ele%a%beta * (scratch%cc(ix_ele)%cbar(1,2)**2 + scratch%cc(ix_ele)%cbar(1,1)**2))
    phase = ele%b%phi + atan2(scratch%cc(ix_ele)%cbar(1,2), scratch%cc(ix_ele)%cbar(1,1))
    datum_value = amp * cos(phase)
    if (associated(ele_ref)) then
      amp = sqrt(ele_ref%a%beta * (scratch%cc(ix_ref)%cbar(1,2)**2 + scratch%cc(ix_ref)%cbar(1,1)**2))
      phase = ele_ref%b%phi + atan2(scratch%cc(ix_ref)%cbar(1,2), scratch%cc(ix_ref)%cbar(1,1))
      datum_value = datum_value - amp * cos(phase)
    endif
    valid_value = .true.

  case ('ping_b.amp_sin_rel_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = -sqrt(ele%a%beta) * scratch%cc(ix_ele)%cbar(1,2)
    if (associated(ele_ref)) then
      datum_value = datum_value + sqrt(ele_ref%a%beta) * scratch%cc(ix_ref)%cbar(1,2)
    endif
    valid_value = .true.

  case ('ping_b.amp_cos_rel_x')
    call tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)
    datum_value = sqrt(ele%a%beta) * scratch%cc(ix_ele)%cbar(1,1)
    if (associated(ele_ref)) then
      datum_value = datum_value - sqrt(ele_ref%a%beta) * scratch%cc(ix_ref)%cbar(1,1)
    endif
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  end select

!-----------

case ('r.')
  if (data_source == 'beam') goto 9000  ! Set error message and return

  i = tao_read_phase_space_index (data_type, 3, .false.)
  j = tao_read_phase_space_index (data_type, 4, .false.)
  if (i == 0 .or. j == 0 .or. len_trim(data_type) /= 4) then
    call tao_set_invalid (datum, 'BAD DATA_TYPE = "' // trim(data_type), why_invalid, .true.)
    return
  endif

  if (ix_ref < 0) ix_ref = 0
  call transfer_matrix_calc (lat, mat6, vec0, ix_ref, ix_start, branch%ix_branch)

  k = ix_start
  do 
    value_vec(k) = mat6(i, j)
    if (k == ix_ele) exit
    k = k + 1
    if (k > n_track) k = 0
    mat6 = matmul(branch%ele(k)%mat6, mat6)
  enddo

  call tao_load_this_datum (value_vec, NULL(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('r56_compaction')
  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (ix_ref < 0) then
    ix_ref = 0
    ele_ref => lat%branch(branch%ix_branch)%ele(ix_ref)
  endif

  orb0 => orbit(ix_ref)
  call make_v_mats (ele_ref, v_mat, v_inv_mat)
  eta_vec = [ele_ref%a%eta, ele_ref%a%etap, ele_ref%b%eta, ele_ref%b%etap]
  eta_vec = matmul (v_mat, eta_vec)
  one_pz = 1 + orb0%vec(6)
  eta_vec(2) = eta_vec(2) * one_pz + orb0%vec(2) / one_pz
  eta_vec(4) = eta_vec(4) * one_pz + orb0%vec(4) / one_pz

  call transfer_matrix_calc (lat, mat6, vec0, ix_ref, ix_start, branch%ix_branch)

  do i = ix_start, ix_ele
    value_vec(i) = sum(mat6(5,1:4) * eta_vec) + mat6(5,6)
    if (i /= ix_ele) mat6 = matmul(branch%ele(i+1)%mat6, mat6)
  enddo
  call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

!-----------

case ('rad_int.')

  if (ix_ref > -1 .or. ix_ele > -1) then
    if (ix_ele < 0) ix_ele = branch%n_ele_track
    if (ix_ref < 0) ix_ref = 0
  endif

  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (.not. allocated(tao_lat%rad_int%branch)) then
    call out_io (s_error$, r_name, 'tao_lat%rad_int%branch not allocated')
    return
  endif
  rad_int_branch => tao_lat%rad_int%branch(ix_branch)

  select case (data_type)
  case ('rad_int.i0')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%i0)
    else
      datum_value = tao_branch%modes_ri%synch_int(0)
    endif

  case ('rad_int.i1')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%i1)
    else
      datum_value = tao_branch%modes_ri%synch_int(1)
    endif

  case ('rad_int.i2')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%i2)
    else
      datum_value = tao_branch%modes_ri%synch_int(2)
    endif

  case ('rad_int.i2_e4')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%lin_i2_e4)
    else
      datum_value = tao_branch%modes_ri%lin%i2_e4
    endif

  case ('rad_int.i3')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%i3)
    else
      datum_value = tao_branch%modes_ri%synch_int(3)
    endif

  case ('rad_int.i3_e7')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%lin_i3_e7)
    else
      datum_value = tao_branch%modes_ri%lin%i3_e7
    endif

  case ('rad_int.i4a')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%i4a)
    else
      datum_value = tao_branch%modes_ri%a%synch_int(4)
    endif

  case ('rad_int.i4b')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%i4b)
    else
      datum_value = tao_branch%modes_ri%b%synch_int(4)
    endif

  case ('rad_int.i4z')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%i4z)
    else
      datum_value = tao_branch%modes_ri%z%synch_int(4)
    endif

  case ('rad_int.i5a')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%i5a)
    else
      datum_value = tao_branch%modes_ri%a%synch_int(5)
    endif

  case ('rad_int.i5a_e6')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%lin_i5a_e6)
    else
      datum_value = tao_branch%modes_ri%lin%i5a_e6
    endif

  case ('rad_int.i5b')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%i5b)
    else
      datum_value = tao_branch%modes_ri%b%synch_int(5)
    endif

  case ('rad_int.i5b_e6')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%lin_i5b_e6)
    else
      datum_value = tao_branch%modes_ri%lin%i5b_e6
    endif

  case ('rad_int.i6b')
    if (ix_ele > -1) then
      datum_value = sum(rad_int_branch%ele(ix_ref:ix_ele)%i6b)
    else
      datum_value = tao_branch%modes_ri%b%synch_int(6)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

  valid_value = .true.

!-----------

case ('rad_int1.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  if (ix_ele < 0) return
  rad_int_branch => tao_lat%rad_int%branch(ix_branch)
  if (.not. allocated(rad_int_branch%ele)) return

  select case (data_type)
  case ('rad_int1.i0')
    datum_value = rad_int_branch%ele(ix_ele)%i0
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%i0

  case ('rad_int1.i1')
    datum_value = rad_int_branch%ele(ix_ele)%i1
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%i1

  case ('rad_int1.i2')
    datum_value = rad_int_branch%ele(ix_ele)%i2
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%i2

  case ('rad_int1.i2_e4')
    datum_value = rad_int_branch%ele(ix_ele)%lin_i2_e4
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%lin_i2_e4

  case ('rad_int1.i3')
    datum_value = rad_int_branch%ele(ix_ele)%i3
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%i3

  case ('rad_int1.i3_e7')
    datum_value = rad_int_branch%ele(ix_ele)%lin_i3_e7
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%lin_i3_e7

  case ('rad_int1.i4a')
    datum_value = rad_int_branch%ele(ix_ele)%i4a
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%i4a

  case ('rad_int1.i5a')
    datum_value = rad_int_branch%ele(ix_ele)%i5a
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%i5a

  case ('rad_int1.i5a_e6')
    datum_value = rad_int_branch%ele(ix_ele)%lin_i5a_e6
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%lin_i5a_e6

  case ('rad_int1.i4b')
    datum_value = rad_int_branch%ele(ix_ele)%i4b
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%i4b

  case ('rad_int1.i5b')
    datum_value = rad_int_branch%ele(ix_ele)%i5b
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%i5b

  case ('rad_int1.i5b_e6')
    datum_value = rad_int_branch%ele(ix_ele)%lin_i5b_e6
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%lin_i5b_e6

  case ('rad_int1.i6b')
    datum_value = rad_int_branch%ele(ix_ele)%i6b
    if (ix_ref > -1) datum_value = datum_value - rad_int_branch%ele(ix_ref)%i6b

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

  valid_value = .true.

!-----------

case ('ref_time')
    call tao_load_this_datum (branch%ele(:)%ref_time, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    
!-----------

case ('rel_floor.')

  select case (data_type)

  case ('rel_floor.x', 'rel_floor.y', 'rel_floor.z')

    if (ix_ref < 0) ele_ref => lat%branch(branch%ix_branch)%ele(0)

    call floor_angles_to_w_mat (-ele_ref%floor%theta, -ele_ref%floor%phi, -ele_ref%floor%psi, w0_mat)

    i = ix_start
    do 
      ele2 => branch%ele(i)
      vec3 = ele2%floor%r - ele_ref%floor%r
      vec3 = matmul (w0_mat, vec3)
      select case (data_type)
      case ('rel_floor.x')
        value_vec(i) = vec3(1)
      case ('rel_floor.y')
        value_vec(i) = vec3(2)
      case ('rel_floor.z')
        value_vec(i) = vec3(3)
      end select
      if (i == ix_ele) exit
      i = i + 1
      if (i > n_track) i = 0
    enddo

    call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case ('rel_floor.theta', 'rel_floor.phi', 'rel_floor.psi')

    if (ix_ref < 0) ele_ref => lat%branch(branch%ix_branch)%ele(0)

    call floor_angles_to_w_mat (-ele_ref%floor%theta, -ele_ref%floor%phi, -ele_ref%floor%psi, w0_mat)

    i = ix_start
    do 
      ele2 => branch%ele(i)
      call floor_angles_to_w_mat (ele2%floor%theta, ele2%floor%phi, ele2%floor%psi, w_mat)
      w_mat = matmul (w0_mat, w_mat)
      call floor_w_mat_to_angles (w_mat, theta, phi, psi)

      select case (data_type)
      case ('rel_floor.theta')
        value_vec(i) = theta
      case ('rel_floor.phi')
        value_vec(i) = phi
      case ('rel_floor.psi')
        value_vec(i) = psi
      end select
      if (i == ix_ele) exit
      i = i + 1
      if (i > n_track) i = 0
    enddo

    call tao_load_this_datum (value_vec, null(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('s_position') 
  if (data_source == 'beam') goto 9000  ! Set error message and return
  if (ix_ref >= 0) then
    datum_value = ele%s - ele_ref%s
  else
    datum_value = ele%s 
  endif
  valid_value = .true.

!-----------

case ('sigma.')

  ! Looks for numbers: e.g. sigma.13
  i = index('123456', data_type(7:7))
  j = index('123456', data_type(8:8))
  if (i > 0 .and. j > 0) then
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(i,j), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(i,j), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif
    return
  endif

  select case (data_type)

  case ('sigma.x')
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(1,1), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif

  case ('sigma.px')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(2,2), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif

  case ('sigma.y')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(3,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(3,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif
    
  case ('sigma.py')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(4,4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(4,4), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif
    
  case ('sigma.z')
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(5,5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
      datum_value = sqrt(datum_value)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(5,5), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif

  case ('sigma.pz')  
    if (data_source == 'lat') then
      if (lat%param%geometry == closed$) then
        call tao_load_this_datum (tao_branch%lat_sigma%mat(6,6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
        datum_value = sqrt(datum_value)
      else
        if (ix_ele == -1) ix_ele = branch%n_ele_track
        rad_int_branch => tao_lat%rad_int%branch(ix_branch)
        datum_value = rad_int_branch%ele(ix_ele)%lin_sig_E / ele%value(E_tot$)
        if (ix_ref > 0) datum_value = datum_value - rad_int_branch%ele(ix_ref)%lin_sig_E / ele_ref%value(E_tot$)
        valid_value = .true.
      endif
    else
      call tao_load_this_datum (bunch_params(:)%sigma(6,6), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
      datum_value = sqrt(datum_value)
    endif
    
  case ('sigma.xy')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(1,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(1,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

  case ('sigma.Lxy')  
    if (data_source == 'lat') then
      call tao_load_this_datum (tao_branch%lat_sigma%mat(1,4) - tao_branch%lat_sigma%mat(2,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    else
      call tao_load_this_datum (bunch_params(:)%sigma(1,4) - bunch_params(:)%sigma(2,3), ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%twiss_valid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('slip_factor_ptc.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  ptc_nf => tao_branch%ptc_normal_form

  if (.not. ptc_nf%valid_map) then
    if (.not. u%calc%one_turn_map) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE ONE_TURN_MAP_CALC IS NOT SET TO TRUE.', why_invalid)
    elseif (branch%param%geometry /= closed$) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE LATTICE GEOMETRY IS NOT CLOSED.', why_invalid)
    else
      call tao_set_invalid (datum, '?????', why_invalid)
    endif
    return
  endif

  if (.not. is_integer(data_type(17:), n)) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  expo = 0
  expo(6) = n 

  datum_value = -real(ptc_nf%phase(3) .sub. expo) / branch%param%total_length
  valid_value = .true.

!-----------

case ('spin.')

  if (.not. bmad_com%spin_tracking_on) then
    call tao_set_invalid (datum, 'NO SPIN TRACKING WHEN BMAD_COM%SPIN_TRACKING_ON = FALSE!', why_invalid)
    return
  endif

  select case (data_type)

  case ('spin.x', 'spin.y', 'spin.z', 'spin.amp')
    do i = ix_start, ix_ele
      if (data_source == 'beam') then
        vec3 = bunch_params(i)%centroid%spin
      else
        vec3 = orbit(i)%spin
      endif

      select case (data_type)
      case ('spin.x');    value_vec(i) = vec3(1)
      case ('spin.y');    value_vec(i) = vec3(2)
      case ('spin.z');    value_vec(i) = vec3(3)
      case ('spin.amp');  value_vec(i) = norm2(vec3)
      end select
    enddo

    if (ix_ref > -1) then
      if (data_source == 'beam') then
        vec3 = bunch_params(ix_ref)%centroid%spin
      else
        vec3 = orbit(ix_ref)%spin
      endif

      select case (data_type)
      case ('spin.x');    value_vec(ix_ref) = vec3(1)
      case ('spin.y');    value_vec(ix_ref) = vec3(2)
      case ('spin.z');    value_vec(ix_ref) = vec3(3)
      case ('spin.amp');  value_vec(ix_ref) = norm2(vec3)
      end select
    endif

    call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    
  case ('spin.depolarization_rate', 'spin.polarization_rate', 'spin.polarization_limit')
    if (.not. tao_branch%spin%valid) call tao_spin_polarization_calc(branch, tao_branch, err_flag = err)
    valid_value = tao_branch%spin%valid

    if (.not. valid_value) then
      call tao_set_invalid (datum, 'ERROR IN SPIN POLARIZAITON CALC.', why_invalid, .false.)
      return
    endif
    select case (data_type)
    case ('spin.depolarization_rate')
      datum_value = tao_branch%spin%depol_rate
    case ('spin.polarization_rate')
      datum_value = tao_branch%spin%pol_rate_bks
    case ('spin.polarization_limit')
      datum_value = tao_branch%spin%pol_limit_dk
    end select

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('spin_dn_dpz.')

  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (.not. bmad_com%spin_tracking_on) then
    call tao_set_invalid (datum, 'NO SPIN TRACKING WHEN BMAD_COM%SPIN_TRACKING_ON = FALSE!', why_invalid)
    return
  endif

  call tao_spin_polarization_calc(branch, tao_branch)

  valid_value = (tao_branch%spin_ele(ix_start)%valid .and. tao_branch%spin_ele(ix_ele)%valid)
  if (ix_ref > -1) valid_value = (valid_value .and. tao_branch%spin_ele(ix_ref)%valid)

  if (.not. valid_value) then
     call tao_set_invalid (datum, 'ERROR IN SPIN POLARIZAITON CALC.', why_invalid, .false.)
     return
  endif

  select case (data_type)
  case ('spin_dn_dpz.x')
    do i = ix_start, ix_ele
      value_vec(i) = tao_branch%spin_ele(i)%dn_dpz%vec(1)
    enddo
    value_vec(ix_ele) = tao_branch%spin_ele(ix_ele)%dn_dpz%vec(1)
  case ('spin_dn_dpz.y')
    do i = ix_start, ix_ele
      value_vec(i) = tao_branch%spin_ele(i)%dn_dpz%vec(2)
    enddo
    value_vec(ix_ele) = tao_branch%spin_ele(ix_ele)%dn_dpz%vec(2)
  case ('spin_dn_dpz.z')
    do i = ix_start, ix_ele
      value_vec(i) = tao_branch%spin_ele(i)%dn_dpz%vec(3)
    enddo
    value_vec(ix_ele) = tao_branch%spin_ele(ix_ele)%dn_dpz%vec(3)
  case ('spin_dn_dpz.amp')
    do i = ix_start, ix_ele
      value_vec(i) = norm2(tao_branch%spin_ele(i)%dn_dpz%vec)
    enddo
    value_vec(ix_ele) = norm2(tao_branch%spin_ele(ix_ele)%dn_dpz%vec)
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

  call tao_load_this_datum (value_vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)


!-----------

case ('spin_g_matrix.')

  call tao_spin_matrix_calc (datum, u, ele_ref, ele)
  valid_value = datum%spin_map%valid

  if (.not. valid_value) then
    call tao_set_invalid (datum, datum%why_invalid)
    return
  endif

  select case (data_type)
  case ('spin_g_matrix.11');  datum_value = datum%spin_map%mat8(7,1)
  case ('spin_g_matrix.12');  datum_value = datum%spin_map%mat8(7,2)
  case ('spin_g_matrix.13');  datum_value = datum%spin_map%mat8(7,3)
  case ('spin_g_matrix.14');  datum_value = datum%spin_map%mat8(7,4)
  case ('spin_g_matrix.15');  datum_value = datum%spin_map%mat8(7,5)
  case ('spin_g_matrix.16');  datum_value = datum%spin_map%mat8(7,6)
  case ('spin_g_matrix.21');  datum_value = datum%spin_map%mat8(8,1)
  case ('spin_g_matrix.22');  datum_value = datum%spin_map%mat8(8,2)
  case ('spin_g_matrix.23');  datum_value = datum%spin_map%mat8(8,3)
  case ('spin_g_matrix.24');  datum_value = datum%spin_map%mat8(8,4)
  case ('spin_g_matrix.25');  datum_value = datum%spin_map%mat8(8,5)
  case ('spin_g_matrix.26');  datum_value = datum%spin_map%mat8(8,6)

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    valid_value = .false.
    return
  end select

!-----------

case ('spin_res.')

  if (data_source == 'beam') goto 9000  ! Set error message and return

  if (.not. bmad_com%spin_tracking_on) then
    call out_io (s_info$, r_name, 'Note: Turning on spin tracking (setting: bmad_com%spin_tracking_on = T)')
    bmad_com%spin_tracking_on = .true.
  endif

  call tao_spin_matrix_calc (datum, u, ele, ele)

  call spin_mat_to_eigen (datum%spin_map%map1%orb_mat, datum%spin_map%map1%spin_q, eval, evec, n0, n_eigen, err)
  if (err) then
    call tao_set_invalid (datum, 'ERROR CONVERTING SPIN/ORBIT 1-TURN MATRIX TO EIGEN VALUES.', why_invalid)
    return
  endif

  j = index('abc', data_type(10:10))
  if (j == 0) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  call spin_quat_resonance_strengths(evec(2*j-1,:), datum%spin_map%map1%spin_q, xi_sum, xi_diff)

  select case (data_type(11:))
  case ('.sum');   datum_value = xi_sum
  case ('.diff');  datum_value = xi_diff
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

  valid_value = .true.

!-----------

case ('spin_tune_ptc.')

  if (data_source == 'beam') goto 9000  ! Set error message and return
  ptc_nf => tao_branch%ptc_normal_form

  if (.not. ptc_nf%valid_map) then
    if (.not. u%calc%one_turn_map) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE ONE_TURN_MAP_CALC IS NOT SET TO TRUE.', why_invalid)
    elseif (branch%param%geometry /= closed$) then
      call tao_set_invalid (datum, 'MAP IS NOT BEING CALCULATED SINCE LATTICE GEOMETRY IS NOT CLOSED.', why_invalid)
    else
      call tao_set_invalid (datum, '?????', why_invalid)
    endif
    return
  endif

  if (.not. is_integer(data_type(15:), n)) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  expo = 0
  expo(6) = n 

  datum_value = real(ptc_nf%spin_tune .sub. expo)
  valid_value = .true.

!-----------

case ('srdt.')
  select case(sub_data_type(1:6))
  case('h00111');  temp_cplx = tao_branch%srdt%h00111
  case('h00201');  temp_cplx = tao_branch%srdt%h00201
  case('h00220');  temp_cplx = tao_branch%srdt%h00220
  case('h00310');  temp_cplx = tao_branch%srdt%h00310
  case('h00400');  temp_cplx = tao_branch%srdt%h00400
  case('h10002');  temp_cplx = tao_branch%srdt%h10002
  case('h10020');  temp_cplx = tao_branch%srdt%h10020
  case('h10110');  temp_cplx = tao_branch%srdt%h10110
  case('h10200');  temp_cplx = tao_branch%srdt%h10200
  case('h11001');  temp_cplx = tao_branch%srdt%h11001
  case('h11110');  temp_cplx = tao_branch%srdt%h11110
  case('h11200');  temp_cplx = tao_branch%srdt%h11200
  case('h20001');  temp_cplx = tao_branch%srdt%h20001
  case('h20020');  temp_cplx = tao_branch%srdt%h20020
  case('h20110');  temp_cplx = tao_branch%srdt%h20110
  case('h20200');  temp_cplx = tao_branch%srdt%h20200
  case('h21000');  temp_cplx = tao_branch%srdt%h21000
  case('h22000');  temp_cplx = tao_branch%srdt%h22000
  case('h30000');  temp_cplx = tao_branch%srdt%h30000
  case('h31000');  temp_cplx = tao_branch%srdt%h31000
  case('h40000');  temp_cplx = tao_branch%srdt%h40000
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    term_found = .false.
    valid_value = .false.
    return
  end select

  valid_value = .true.
  select case (sub_data_type(8:8))
  case('r')
    datum_value = real(temp_cplx)
  case('i')
    datum_value = aimag(temp_cplx)
  case('a')
    datum_value = abs(temp_cplx)
  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID (DATA_TYPE NOT ENDING IN .r, .i, or .a).', why_invalid, .true.)
    valid_value = .false.
  end select

!-----------

case ('time')
  if (data_source == 'beam') then
    call tao_load_this_datum (bunch_params%centroid%t, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
  else
    call tao_load_this_datum (orbit(:)%t, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  endif

!-----------

case ('tune.')

  if (data_source == 'beam') goto 9000  ! Set error message and return

  select case (data_type)
  case ('tune.a')
    datum_value = branch%ele(branch%n_ele_track)%a%phi
    valid_value = .true.

  case ('tune.b')
    datum_value = branch%ele(branch%n_ele_track)%b%phi
    valid_value = .true.

  case ('tune.z')
    call calc_z_tune (branch)
    datum_value = -branch%z%tune
    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  end select

!-----------

case ('t.', 'tt.')
  if (branch%ix_branch /= 0) then
    call out_io (s_fatal$, r_name, 'TRANSFER MAP CALC NOT YET MODIFIED FOR BRANCHES.')
    return
  endif
  if (data_source == 'beam') goto 9000  ! Set error message and return

  expnt = 0
  i = tao_read_phase_space_index (sub_data_type, 1, .false.)
  do j = 2, 20
    if (sub_data_type(j:j) == ' ') exit
    k = tao_read_phase_space_index (sub_data_type, j, .false.); if (k == 0) exit
    expnt(k) = expnt(k) + 1
  enddo

  if (i == 0 .or. k == 0) then
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return
  endif

  if (ix_ref < 0) ix_ref = 0

  ! Computation if there is no range

  if (ix_start == ix_ele) then
    if (s%com%ix_ref_taylor /= ix_ref .or. s%com%ix_ele_taylor /= ix_ele) then
      ix0 = s%com%ix_ele_taylor
      if (s%com%ix_ref_taylor == ix_ref .and. ix_ele > ix0) then
        call transfer_map_calc (lat, taylor_save, err, ix0, ix_ele, orbit(ix0), &
                                                  unit_start = .false., concat_if_possible = s%global%concatenate_maps)
      else
        call transfer_map_calc (lat, taylor_save, err, ix_ref, ix_ele, orbit(ix_ref), concat_if_possible = s%global%concatenate_maps)
      endif

      if (err) then
        call tao_set_invalid (datum, 'MAP TERM OVERFLOW', why_invalid)
        return
      endif

      s%com%ix_ref_taylor = ix_ref
      s%com%ix_ele_taylor = ix_ele
    endif
    datum_value = taylor_coef (taylor_save(i), expnt)
    valid_value = .true.

  ! Here if there is a range.
  else
    k = ix_start
    call transfer_map_calc (lat, taylor, err, ix_ref, k, orbit(ix_ref), concat_if_possible = s%global%concatenate_maps)
    if (err) then
      call tao_set_invalid (datum, 'MAP TERM OVERFLOW', why_invalid)
      return
    endif

    do
      value_vec(k) = taylor_coef (taylor(i), expnt)
      if (k == ix_ele) exit
      k_old = k
      k = k + 1
      if (k > branch%n_ele_track) k = 0
      call transfer_map_calc (lat, taylor, err, k_old, k, unit_start = .false., concat_if_possible = s%global%concatenate_maps)
      if (err) then
        call tao_set_invalid (datum, 'MAP TERM OVERFLOW', why_invalid)
        return
      endif
    enddo
    call tao_load_this_datum (value_vec, NULL(), ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
  endif

!-----------

case ('unstable.')

  select case (data_type)

  case ('unstable.eigen', 'unstable.eigen.a', 'unstable.eigen.b', 'unstable.eigen.c')
    call transfer_matrix_calc (lat, mat6, vec0, 0, branch%n_ele_track, branch%ix_branch, one_turn = .true.)
    call mat_eigen (mat6, eval, evec, err)
    if (err) then
      call tao_set_invalid (datum, 'CANNOT COMPUTE EIGENVALUES FOR TRANSFER MATRIX', why_invalid)
      return
    endif

    select case (data_type)
    case ('unstable.eigen');    datum_value = maxval(abs(eval))
    case ('unstable.eigen.a');  datum_value = max(abs(eval(1)), abs(eval(2)))
    case ('unstable.eigen.b');  datum_value = max(abs(eval(3)), abs(eval(4)))
    case ('unstable.eigen.c');  datum_value = max(abs(eval(5)), abs(eval(6)))
    end select
    valid_value = .true.
    

  case ('unstable.orbit')

    if (data_source == 'beam') then
      ie0 = u%model_branch(ix_branch)%beam%ix_track_start

      if (datum%ele_name == '') then
        ie1 = u%model_branch(ix_branch)%beam%ix_track_end
      else
        ie1 = ix_ele
      endif

      if (ie0 == not_set$) then
        call tao_set_invalid (datum, 'NO TRACKING DONE IN BRANCH', why_invalid)
        return
      endif

      if (ie1 > ie0) then
        n = ie1 - ie0
      else
        n = branch%n_ele_track - ie0 + ie1
      endif

      datum_value = 0
      do j = 1, branch%n_ele_track
        jj = j + ie0
        if (jj > branch%n_ele_track) jj = jj - branch%n_ele_track
        datum_value = datum_value + (n - j + 1) * bunch_params(jj)%n_particle_lost_in_ele
        if (jj == ie1) exit
      enddo
      datum_value = datum_value / bunch_params(ie0)%n_particle_tot
      datum%ix_ele_merit = -1

    elseif (lat%param%geometry == open$) then
      if (datum%ele_name == '') ix_ele = branch%n_ele_track
      iz = tao_branch%track_state
      if (iz /= moving_forward$ .and. iz <= ix_ele) then
        datum_value = 1 + ix_ele - iz
        if (orbit(iz)%s < branch%ele(iz)%s) then
          orb => orbit(iz-1)
          datum%ix_ele_merit = iz - 1
          if (branch%ele(iz)%value(L$) /= 0) then
            ! Add s_rel/L 
            datum_value = datum_value + (branch%ele(iz)%s - orbit(iz)%s)/branch%ele(iz)%value(L$)
          endif
        else
          datum_value = datum_value - 0.5
          orb => orbit(iz)
          datum%ix_ele_merit = iz
        endif

        datum_value = datum_value + 0.5 * tanh(lat%param%unstable_factor)
      endif

    else   ! closed geometry
      if (tao_branch%track_state == moving_forward$) then
        datum_value = 0
      else
        datum_value = 1
      endif
      datum%ix_ele_merit = 0
    endif

    valid_value = .true.

  case ('unstable.ring', 'unstable.lattice')
    if (data_source == 'beam') goto 9000  ! Set error message and return

    if (data_type == 'unstable.ring') then
      call out_io (s_error$, r_name, '"unstable.ring" has been replaced by "unstable.lattice". Please change this in your input file.')
    endif

    if (lat%param%geometry == closed$ .and. tao_branch%track_state /= moving_forward$) then
      datum_value = 1
    else
      datum_value = lat%param%unstable_factor
      ! unstable_penalty is needed since at the metastable borderline the growth rate is zero.
      if (.not. lat%param%stable) datum_value = datum_value + s%global%unstable_penalty
    endif

    valid_value = .true.

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
    return

  end select

!-----------

case ('velocity', 'velocity.')

  if (tao_branch%track_state /= moving_forward$ .and. ix_ele > tao_branch%track_state) then
    valid_value = .false.
    call tao_set_invalid (datum, 'Particle lost.', why_invalid)
    return
  endif

  select case (data_type)

  case ('velocity')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%beta, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('velocity.x')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(2)*(1+bunch_params%centroid%vec(6))*bunch_params%centroid%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(2)*(1+orbit(:)%vec(6))*orbit(:)%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('velocity.y')
    if (data_source == 'beam') then
      call tao_load_this_datum (bunch_params%centroid%vec(4)*(1+bunch_params%centroid%vec(6))*bunch_params%centroid%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (orbit(:)%vec(4)*(1+orbit(:)%vec(6))*orbit(:)%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case ('velocity.z')
    if (data_source == 'beam') then
      call tao_load_this_datum (sqrt(1 - (bunch_params%centroid%vec(2)*(1+bunch_params%centroid%vec(6)))**2 - (bunch_params%centroid%vec(4)*(1+bunch_params%centroid%vec(6)))**2)*(1+bunch_params%centroid%vec(6))*bunch_params%centroid%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, bunch_params%n_particle_live > 0)
    else
      call tao_load_this_datum (sqrt(1 - (orbit(:)%vec(2)*(1+orbit(:)%vec(6)))**2 - (orbit(:)%vec(4)*(1+orbit(:)%vec(6)))**2)*(1+orbit(:)%vec(6))*orbit(:)%beta, &
                        ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid)
    endif

  case default
    call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  end select

!-----------

case ('wall.')
  if (data_source == 'beam') goto 9000  ! Set error message and return

  constraint = data_type(6:)
  zz0_pt = 0
  datum_value = 1e10   ! Something large
  found = .false.

  if (.not. allocated(s%building_wall%section)) then
    valid_value = .false.
    call tao_set_invalid (datum, 'No building wall sections defined.', why_invalid)
    return
  endif

  do i = 1, size(s%building_wall%section)
    section => s%building_wall%section(i)
    if (section%constraint /= constraint) cycle
    do ie = ix_start, ix_ele
      ele => branch%ele(ie)

      ! (zz, xx) is the local reference coordinate system at the element

      do is = 1, size(section%point)
        pt = tao_oreint_building_wall_pt(section%point(is))
        dz = pt%z - ele%floor%r(3); dx = pt%x - ele%floor%r(1)
        cos_theta = cos(ele%floor%theta); sin_theta = sin(ele%floor%theta)
        zz_pt =  dz * cos_theta + dx * sin_theta
        xx_pt = -dz * sin_theta + dx * cos_theta

        if (is == 1) then
          zz0_pt = zz_pt
          xx0_pt = xx_pt
          cycle
        endif


        if (pt%radius == 0 .or. zz_pt == 0 .or. zz0_pt == 0) then
          ! The perpendicular to the machine line intersects the segment if
          ! zz_pt and zz0_pt have a different signs.
          if (zz_pt * zz0_pt > 0) cycle
          xx_wall = (xx0_pt * zz_pt - xx_pt * zz0_pt) / (zz_pt - zz0_pt)

        else  ! Circular arc
          dz = pt%z_center - ele%floor%r(3); dx = pt%x_center - ele%floor%r(1)
          zz_center =  dz * cos_theta + dx * sin_theta
          xx_center = -dz * sin_theta + dx * cos_theta
          drad = pt%radius**2 -zz_center**2
          if (drad <= 0) cycle
          drad = sqrt(drad)
          ! There are two possible points at (0, xx_a) and (0, xx_b) where the perpendicular 
          ! to the machine line intersects the arc.
          xx_a = xx_center - drad
          xx_b = xx_center + drad
          dxx1 = xx_pt - xx0_pt
          dzz1 = zz_pt - zz0_pt
          ang_a = dzz1 * (xx_a - xx0_pt) + dxx1 * zz0_pt
          ang_b = dzz1 * (xx_b - xx0_pt) + dxx1 * zz0_pt
          ang_c = dzz1 * (xx_center - xx0_pt) - dxx1 * (zz_center - zz0_pt)
          ! A point is within the circular arc if it and the center point are on opposite
          ! sides of the chord from (zz0_pt, xx0_pt) to (zz_pt, xx_pt). 
          ! This assumes the arc is less than 180^deg.
          ! It should not be that both intersection points are within the arc.
          if (ang_a * ang_c < 0) then
            xx_wall = xx_a
          elseif (ang_b * ang_c < 0) then
            xx_wall = xx_b
          else
            cycle
          endif
        endif

        if (data_type =='wall.right_side')  xx_wall = -xx_wall
        datum_value = min(datum_value, xx_wall)
        valid_value = .true.

        zz0_pt = zz_pt
        xx0_pt = xx_pt
      enddo
    enddo

  enddo

  if (.not. valid_value) call tao_set_invalid (datum, 'No wall section found in the transverse plane of the evaluation point.', why_invalid)

!-----------

case ('wire.')  
  if (data_source == 'lat') goto 9001  ! Error message and return
  read (data_type(6:), '(a)') angle
  datum_value = tao_do_wire_scan (ele, angle, u%model_branch(branch%ix_branch)%ele(ix_ele)%beam)
  valid_value = .true.
  
case default
  call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(data_type) // '" IS NOT VALID', why_invalid, .true.)
  return
end select

!-----------------------------------------------------------------------
! End stuff

if (datum%ix_ele_merit > -1) then
  datum%s = tao_datum_s_position(datum, branch%ele(datum%ix_ele_merit))
elseif (associated(ele)) then
  datum%s = tao_datum_s_position(datum, ele)
else
  datum%s = real_garbage$
endif

if (valid_value) datum%err_message_printed = .false.  ! Reset

return

!----------------------------------------------------------------------

7000 continue
call tao_set_invalid (datum, 'PARTICLE SPECIES TYPE NOT SET ??!! PLEASE SEEK HELP!', why_invalid)
return

9000 continue
call tao_set_invalid (datum, 'DATA_SOURCE = "beam" NOT VALID FOR THIS DATA_TYPE: ' // datum%data_type, why_invalid, .true.)
return

9001 continue
call tao_set_invalid (datum, 'DATA_SOURCE = "lat" NOT VALID FOR THIS DATA_TYPE: ' // datum%data_type, why_invalid, .true.)
return

9100 continue
call tao_set_invalid (datum, 'DATA_TYPE: ' // quote(datum%data_type) // ' NOT APPLICABLE TO A LATTICE BRANCH WITH AN OPEN GEOMETRY.', why_invalid, .true.)
return

9101 continue
call tao_set_invalid (datum, 'DATA_TYPE: ' // quote(datum%data_type) // ' NOT APPLICABLE TO A LATTICE BRANCH WITH A CLOSED GEOMETRY.', why_invalid, .true.)
return

end subroutine tao_evaluate_a_datum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine tao_load_this_datum (vec, ele_ref, ele_start, ele, datum_value, valid_value, datum, branch, why_invalid, good)

type (tao_data_struct) datum
type (branch_struct) branch
type (ele_struct), pointer :: ele_ref, ele_start, ele

real(rp), target :: vec(0:)
real(rp) datum_value, ref_value, l_sum
real(rp), pointer :: vec_ptr(:)
real(rp), allocatable :: s_pos(:), value(:)

character(*), parameter :: r_name = 'tao_load_this_datum'
character(*), optional :: why_invalid

integer ix_m, i, n_track, ix_m2, ix_ref, ix_start, ix_ele, n

logical valid_value
logical, optional :: good(0:)

!

valid_value = .true.
datum_value = 0

n_track = branch%n_ele_track
ix_start = -1; ix_ref = -1; ix_ele = -1
if (associated(ele)) ix_ele = tao_tracking_ele_index(ele, datum)
if (associated(ele_ref)) ix_ref = tao_tracking_ele_index(ele_ref, datum)
if (associated(ele_start)) ix_start = tao_tracking_ele_index(ele_start, datum)

if (ix_ele < 0) then
  datum%exists = .false.
  valid_value = .false.
  if (.not. associated(ele)) then
    call tao_set_invalid (datum, 'NO ASSOCIATED LATTICE ELEMENT TO EVALUATE AT!', why_invalid, .true.)
  else
    call tao_set_invalid (datum, 'BAD LATTICE ELEMENT NAME: ' // datum%ele_name, why_invalid, .true.)
  endif
  return
endif

if (ix_ele < ix_start .and. branch%param%geometry == open$) then
  if (datum%useit_opt) call out_io (s_error$, r_name, &
                'ERROR: ELEMENTS ARE REVERSED FOR: ' // tao_datum_name(datum), &
                'STARTING ELEMENT: ' // ele_start%name, &
                'IS AFTER ENDING ELEMENT: ' // ele%name)
  call tao_set_invalid (datum, 'ELEMENTS ARE REVERSED FOR: ' // tao_datum_name(datum), why_invalid, .true.)
  valid_value = .false.
  return
endif

! Set up refrence value

if (ix_ref > -1) then
  ref_value = vec(ix_ref)
  if (present(good)) then
    if (.not. good(ix_ref)) then
      call tao_set_invalid (datum, 'CANNOT EVALUATE DATUM AT REFERENCE ELEMENT.', why_invalid)
      valid_value = .false.
      return
    endif
  endif
else
  ref_value = 0
endif

 
! If ele_start does not exist

if (datum%ele_start_name == '' .or. ix_start == ix_ele .or. datum%merit_type == 'target') then
  datum_value = vec(ix_ele) - ref_value
  if (datum%merit_type(1:4) == 'abs_') datum_value = abs(vec(ele%ix_ele))
  if (present(good)) valid_value = good(ix_ele)
  if (.not. valid_value) call tao_set_invalid (datum, 'CANNOT EVALUATE DATUM.', why_invalid)
  return
endif

! Set up the vector of values with the reference subtracted off

if (ref_value == 0) then
  vec_ptr => vec
else
  allocate(vec_ptr(0:ubound(vec,1)))
  vec_ptr = vec - ref_value
endif

!------------------------
! If there is a range

if (ix_ele < ix_start) then   ! wrap around
  if (present(good)) then
    if (.not. all(good(0:ix_ele)) .or. .not. all(good(ix_start:n_track))) then
      call tao_set_invalid (datum, 'CANNOT EVALUATE OVER EVALUATION RANGE.', why_invalid)
      valid_value = .false.
      return
    endif
  endif

  select case (datum%merit_type)
  case ('min')
    ix_m = minloc (vec_ptr(0:ix_ele), 1) - 1
    ix_m2 = minloc (vec_ptr(ix_start:n_track), 1) + ix_start - 1
    if (vec_ptr(ix_m2) < vec_ptr(ix_m2)) ix_m = ix_m2
    datum_value = vec_ptr(ix_m)

  case ('max')
    ix_m = maxloc (vec_ptr(0:ix_ele), 1) - 1
    ix_m2 = maxloc (vec_ptr(ix_start:n_track), 1) + ix_start - 1
    if (vec_ptr(ix_m2) > vec_ptr(ix_m2)) ix_m = ix_m2
    datum_value = vec_ptr(ix_m)

  case ('abs_min')
    ix_m = minloc (abs(vec_ptr(0:ix_ele)), 1) - 1
    ix_m2 = minloc (abs(vec_ptr(ix_start:n_track)), 1) + ix_start - 1
    if (abs(vec_ptr(ix_m2)) < abs(vec_ptr(ix_m2))) ix_m = ix_m2
    datum_value = abs(vec_ptr(ix_m))

  case ('abs_max')
    ix_m = maxloc (abs(vec_ptr(0:ix_ele)), 1) - 1
    ix_m2 = maxloc (abs(vec_ptr(ix_start:n_track)), 1) + ix_start - 1
    if (abs(vec_ptr(ix_m2)) > abs(vec_ptr(ix_m2))) ix_m = ix_m2
    datum_value = abs(vec_ptr(ix_m))

  case ('max-min')
    ix_m  = maxloc (vec_ptr(0:ix_ele), 1) - 1
    ix_m2 = maxloc (vec_ptr(ix_start:n_track), 1) + ix_start - 1
    if (vec_ptr(ix_m2) > vec_ptr(ix_m2)) ix_m = ix_m2
    datum_value = vec_ptr(ix_m)

    ix_m  = minloc (vec_ptr(0:ix_ele), 1) - 1
    ix_m2 = minloc (vec_ptr(ix_start:n_track), 1) + ix_start - 1
    if (vec_ptr(ix_m2) < vec_ptr(ix_m2)) ix_m = ix_m2
    datum_value = datum_value - vec_ptr(ix_m)

  case ('average', 'integral', 'rms')
    n = ix_ele + n_track - ix_start + 2
    allocate(s_pos(n), value(n))

    n = 0
    do i = ix_start, n_track
      n = n + 1
      s_pos(n) = tao_datum_s_position(datum, branch%ele(i))
      value(n) = vec_ptr(i)
    enddo

    do i = 0, ix_ele
      n = n + 1
      s_pos(n) = tao_datum_s_position(datum, branch%ele(i)) + branch%param%total_length
      value(n) = vec_ptr(i)
    enddo

    datum_value = tao_datum_integrate(datum, branch, s_pos, vec_ptr, valid_value, why_invalid)

  case default
    call tao_set_invalid (datum, 'BAD MERIT_TYPE WHEN THERE IS A RANGE OF ELEMENTS: ' // datum%merit_type, why_invalid, .true.)
    valid_value = .false.
    return
  end select

  if (.not. valid_value) call tao_set_invalid (datum, 'INVALID DATA IN RANGE FROM ELE_START TO ELE_REF', why_invalid)

! no wrap case

else

  if (present(good)) then
    if (.not. all(good(ix_start:ix_ele))) then
      call tao_set_invalid (datum, 'CANNOT EVALUATE OVER EVALUATION RANGE.', why_invalid)
      valid_value = .false.
      return
    endif
  endif

  select case (datum%merit_type)
  case ('min')
    ix_m = minloc (vec_ptr(ix_start:ix_ele), 1) + ix_start - 1
    datum_value = vec_ptr(ix_m)

  case ('max')
    ix_m = maxloc (vec_ptr(ix_start:ix_ele), 1) + ix_start - 1
    datum_value = vec_ptr(ix_m)

  case ('abs_min')
    ix_m = minloc (abs(vec_ptr(ix_start:ix_ele)), 1) + ix_start - 1
    datum_value = abs(vec_ptr(ix_m))

  case ('abs_max')
    ix_m = maxloc (abs(vec_ptr(ix_start:ix_ele)), 1) + ix_start - 1
    datum_value = abs(vec_ptr(ix_m))

  case ('max-min')
    ix_m = maxloc (vec_ptr(ix_start:ix_ele), 1) + ix_start - 1
    datum_value = vec_ptr(ix_m)

    ix_m = minloc (vec_ptr(ix_start:ix_ele), 1) + ix_start - 1
    datum_value = datum_value - vec_ptr(ix_m)

  case ('average', 'integral', 'rms')
    n = ix_ele - ix_start + 1
    allocate(s_pos(n), value(n))

    n = 0
    do i = ix_start, ix_ele
      n = n + 1
      s_pos(n) = tao_datum_s_position(datum, branch%ele(i))
      value(n) = vec_ptr(i)
    enddo

    datum_value = tao_datum_integrate(datum, branch, s_pos, value, valid_value, why_invalid)
    ix_m = -1

  case default
    call tao_set_invalid (datum, 'BAD MERIT_TYPE WHEN THERE IS A RANGE OF ELEMENTS: ' // datum%merit_type, why_invalid, .true.)
    valid_value = .false.
    return
  end select

  if (.not. valid_value) call tao_set_invalid (datum, 'INVALID DATA IN RANGE FROM ELE_START TO ELE_REF', why_invalid)
endif

datum%ix_ele_merit = ix_m
if (ref_value /= 0) deallocate (vec_ptr)

end subroutine tao_load_this_datum

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_datum_s_position (datum, ele) result (s_pos)
!
! Routine to calculate the longitudinal position associated with a datum.
!
! Input:
!   datum     -- tao_data_struct: Datum under conideration.
!   ele       -- ele_struct: Associated lattice element.
!
! Output
!   s_pos     -- real(rp): Associated longitudinal position.
!-

function tao_datum_s_position (datum, ele) result (s_pos)

type (tao_data_struct) datum
type (ele_struct) ele

real(rp) s_pos

!

select case (datum%eval_point)
case (anchor_beginning$); s_pos = ele%s_start
case (anchor_center$);    s_pos = 0.5_rp * (ele%s_start + ele%s)
case (anchor_end$);       s_pos = ele%s
end select

end function tao_datum_s_position

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_datum_integrate (datum, branch, s_pos, values, valid_value, why_invalid) result (result)
!
! Routine to calculate the integral, rms, or average of an array of values associated with a datum.
!
! Input:
!   datum         -- tao_data_struct: Datum under consideration.
!   branch        -- branch_struct: Associated lattice branch.
!   s_pos(:)      -- real(rp): Array of s-positions of the values.
!   values(:)     -- real(rp): Array of values.
!
! Output:
!   valid_value   -- logical: Set false if, for example, all s_pos(:) are the same.
!   why_invalid   -- character(*): Information string if there is a problem.
!   result      -- real(rp): Integral, rms, or average depending upon datum%merit_type.
!-

function tao_datum_integrate (datum, branch, s_pos, values, valid_value, why_invalid) result (result)

type (tao_data_struct) datum
type (branch_struct) branch

real(rp) values(:), s_pos(:)
real(rp) result, integ, integ2, ds1, ds2
integer i, n_pt

logical valid_value
character(*) why_invalid

!

valid_value = .false.

if (size(values) /= size(s_pos)) then
  call tao_set_invalid (datum, 'INTERNAL ERROR. PLEASE REPORT!', why_invalid)
  return
endif

n_pt = size(values)
if (n_pt < 2) then
  call tao_set_invalid (datum, 'NUMBER OF POINTS TO INTEGRATE/AVERAGE OVER LESS THAN 2!')
  return
endif

ds1 = s_pos(2) - s_pos(1)
ds2 = s_pos(n_pt) - s_pos(n_pt-1)

integ = 0.5_rp * (values(1) * ds1 + values(n_pt) * ds2)
integ2 = 0.5_rp * (values(1)**2 * ds1 + values(n_pt)**2 * ds2)

do i = 2, n_pt-1
  ds2 = s_pos(i+1) - s_pos(i)
  integ = integ + 0.5_rp * values(i) * (ds1 + ds2)
  integ2 = integ2 + 0.5_rp * values(i)**2 * (ds1 + ds2)
  ds1 = ds2
enddo

select case (datum%merit_type)
case ('integral')
  result = integ

case ('average', 'rms')
  ds2 = s_pos(n_pt) - s_pos(1)
  if (ds2 == 0) then
    call tao_set_invalid (datum, 'INTERVAL TO AVERAGE OVER HAS ZERO LENGTH!')
    return
  endif

  if (datum%merit_type == 'average') then
    result = integ / ds2
  else
    result = sqrt(max(0.0_rp, integ2/ds2 - (integ / ds2)**2))
  endif

case default
  call err_exit
end select

valid_value = .true.

end function tao_datum_integrate

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_tracking_ele_index(ele, datum, ix_branch) result (ix_ele)
!
! Routine to return the index in the tracking part of a lattice that corresponds to ele.
!
! Input:
!   ele     -- ele_struct: Lattice element.
!   datum   -- tao_data_struct: Datum
!
! Output:
!   ix_branch -- integer, optional: Lattice branch associated with element
!   ix_ele    -- integer: Element index associated with ele.
!-

function tao_tracking_ele_index(ele, datum, ix_branch) result (ix_ele)

type (ele_struct), pointer :: ele, ele2
type (tao_data_struct) datum
integer ix_ele
integer, optional :: ix_branch

!

select case (ele%lord_status)
case (super_lord$, overlay_lord$, group_lord$)
  ele2 => pointer_to_slave(ele, ele%n_slave)

case default
  ele2 => ele
end select

ix_ele = ele2%ix_ele
if (present(ix_branch)) ix_branch = ele2%ix_branch

end function tao_tracking_ele_index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine integrate_min (ix_start, ix_ele, datum_value, ix_m, branch, vec, datum)

type (branch_struct) branch
type (tao_data_struct) datum

real(rp) vec(0:), val0, dv0, dv1, ds
real(rp) datum_value

integer ix_m, i, ix_start, ix_ele

!

val0 = datum%meas_value

do i = ix_start, ix_ele
  if (ix_m < 0) ix_m = i
  if (vec(i) < vec(ix_m)) ix_m = i
  if (i == ix_start) cycle
  dv0 = val0 - vec(i-1) 
  dv1 = val0 - vec(i)
  ds = branch%ele(i)%s - branch%ele(i)%s_start
  if (dv0 < 0 .and. dv1 < 0) cycle
  if (dv0 > 0 .and. dv1 > 0) then
    datum_value = datum_value + 0.5 * ds * (dv0 + dv1)
  elseif (dv0 > 0) then
    datum_value = datum_value + 0.5 * ds * dv0**2 / (dv0 - dv1)
  else
    datum_value = datum_value + 0.5 * ds * dv1**2 / (dv1 - dv0)
  endif
enddo

end subroutine integrate_min

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine integrate_max (ix_start, ix_ele, datum_value, ix_m, branch, vec, datum)

type (branch_struct) branch
type (tao_data_struct) datum

real(rp) vec(0:), val0, dv0, dv1, ds
real(rp) datum_value

integer ix_m, i, ix_start, ix_ele

!

val0 = datum%meas_value

do i = ix_start, ix_ele
  if (ix_m < 0) ix_m = i
  if (vec(i) > vec(ix_m)) ix_m = i
  if (i == ix_start) cycle
  dv0 = vec(i-1) - val0
  dv1 = vec(i) - val0
  ds = branch%ele(i)%s - branch%ele(i)%s_start
  if (dv0 < 0 .and. dv1 < 0) cycle
  if (dv0 > 0 .and. dv1 > 0) then
    datum_value = datum_value + 0.5 * ds * (dv0 + dv1)
  elseif (dv0 > 0) then
    datum_value = datum_value + 0.5 * ds * dv0**2 / (dv0 - dv1)
  else
    datum_value = datum_value + 0.5 * ds * dv1**2 / (dv1 - dv0)
  endif
enddo

end subroutine integrate_max

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine tao_scratch_values_calc (ele_ref, ele_start, ele, datum, branch, orbit)

type (ele_struct), pointer :: ele_ref, ele_start, ele
type (tao_beam_shake_struct), pointer :: cc_p
type (tao_data_struct) datum
type (branch_struct) branch
type (coord_struct) orbit(0:)

integer i, ip, ix_ele, ix_start

!

ip = index(datum%data_type, '.')
if (ip == 0) return

call this_scratch_calc(ele%ix_ele, datum%data_type(1:ip), datum, branch, orbit)
if (associated(ele_ref)) call this_scratch_calc(ele_ref%ix_ele, datum%data_type(1:ip), datum, branch, orbit)

ix_ele = tao_tracking_ele_index(ele, datum)

if (associated(ele_start)) then
  ix_start = tao_tracking_ele_index(ele_start, datum)
  if (ix_start <= ix_ele) then
    do i = ix_start, ix_ele
      call this_scratch_calc(i, datum%data_type(1:ip), datum, branch, orbit)
    enddo
  else
    do i = ix_start, branch%n_ele_track
      call this_scratch_calc(i, datum%data_type(1:ip), datum, branch, orbit)
    enddo
    do i = 1, ix_ele
      call this_scratch_calc(i, datum%data_type(1:ip), datum, branch, orbit)
    enddo
  endif
endif

!--------------------------------
contains

subroutine this_scratch_calc(ix_ele, data_class, datum, branch, orbit)

type (tao_data_struct) datum
type (branch_struct), target :: branch
type (coord_struct) orbit(0:)
type (ele_struct), pointer :: ele
type (tao_beam_shake_struct), pointer :: cc_p

real(rp) f, f1, f2
integer ix_ele
character(*) data_class

cc_p => scratch%cc(ix_ele)
ele => branch%ele(ix_ele)

select case (data_class)

case ('k.', 'cbar.', 'ping_a.', 'ping_b.')
  if (cc_p%coupling_calc_done) return
  cc_p%coupling_calc_done = .true.

  call c_to_cbar (ele, cc_p%cbar)
  f = sqrt(ele%a%beta/ele%b%beta) 
  f1 = f / ele%gamma_c
  f2 = 1 / (f * ele%gamma_c)

  cc_p%k_11a = cc_p%cbar(1,1) * f1
  cc_p%k_12a = cc_p%cbar(1,2) * f2
  cc_p%k_12b = cc_p%cbar(1,2) * f1
  cc_p%k_22b = cc_p%cbar(2,2) * f2

! Amplitude calc
! 'orbit.amp_a', 'orbit.amp_b', 'orbit.norm_amp_a', 'orbit.norm_amp_b'

case ('orbit.')
  if (index(datum%data_type, 'amp_') == 0) return
  if (cc_p%amp_calc_done) return
  cc_p%amp_calc_done = .true.

  call orbit_amplitude_calc (ele, orbit(ix_ele), cc_p%amp_a, cc_p%amp_b, cc_p%amp_na, cc_p%amp_nb)

end select

end subroutine this_scratch_calc

end subroutine tao_scratch_values_calc

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Subroutine tao_do_wire_scan (ele, wire_params, theta, beam) result (moment)
!
! Returns the beam's second moment using the wire along the specified angle.
! Keep in mind that the actual correlation axis is 90 degrees off of the 
! wire angle
!
! This simulates a fast wire scanner that performs the scan over only one
! bunch. Obviously, this isn't realistic. Any dynamic effects will not be
! accounted for!
!
! Input:
!  ele         -- Element_struct: 
!    %value(noise$) -- relative wire resolution RMS 
!    %value(tilt$)  -- wire angle error in radians rms.
!  theta       -- Real(rp): wire angle wrt x axis (in degrees)
!  beam        -- Beam_struct: contains the beam distribution
!
! Output:
!   moment  -- Real(rp): second moment along axis specified by angle.
!-

function tao_do_wire_scan (ele, theta, beam) result (moment)

use random_mod

type (ele_struct) ele
type (beam_struct) beam

real(rp), allocatable :: dist(:)
real(rp) theta, theta_rad, moment, ran_num(2)
real(rp) avg

!

call ran_gauss (ran_num)

! angle in radians and correlation angle is 90 off from wire angle
theta_rad = ele%value(tilt_tot$)+ (theta - 90) * (2.0*pi / 360.0)

if (.not. allocated (dist)) then
  allocate (dist(size(beam%bunch(1)%particle)))
elseif (size(dist) /= size(beam%bunch(1)%particle)) then
  deallocate (dist)
  allocate (dist(size(beam%bunch(1)%particle)))
endif

! Rotating the wire scanner is equivalent to rotating the beam by -theta

dist =  beam%bunch(1)%particle%vec(1) * cos(-theta_rad ) &
        + beam%bunch(1)%particle%vec(3) * sin(-theta_rad)
  
avg = sum (dist, mask = (beam%bunch(1)%particle%state == alive$)) &
          / count (beam%bunch(1)%particle%state == alive$)
        
moment = (1 + ele%value(noise$)*ran_num(2)) * sum ((dist-avg)*(dist-avg), &
                 mask = (beam%bunch(1)%particle%state == alive$)) &
          / count (beam%bunch(1)%particle%state == alive$)

end function tao_do_wire_scan

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+         
! Function tao_pointer_to_datum_ele (lat, ix_ele, datum, valid, why_invalid) result (ele)
! 
! Routine to see if an element index corresponds to an element with a definite 
! location such as an overlay or multipass element.
!
! If the element is a super_lord then the super_slave element at the exit end
! of the lord will be returned. Otherwise ix_loc will be set to ix_ele.
!
! Input:
!   lat    -- Lat_struct: Lattice
!   ix_ele -- Integer: Index of element.
!   datum  -- Tao_data_struct: Used for error messages and gives branch index.
!
! Output:
!   ele          -- Ele_struct, pointer :: Pointer to the element. Set to NULL if not valid 
!                     or no associated element.
!   valid        -- Logical: Set False if element does not have a definite location.
!                     Set True otherwise
!   why_invalid  -- Character(*), optional: Tells why datum value is invalid.
!-

function tao_pointer_to_datum_ele (lat, ele_name, ix_ele, datum, valid, why_invalid) result (ele)

type (lat_struct) lat
type (tao_data_struct) datum
type (ele_struct), pointer :: ele

integer ix_ele, ixc, n_track, n_max

logical valid

character(*) ele_name
character(*), optional :: why_invalid
character(*), parameter :: r_name = 'tao_pointer_to_datum_ele'

! If ele_name is blank but ix_ele is not -1 then this datum was constructed internally
! by Tao (as opposed being constructed from tao init file info).
! If this is the case, assume that everything is OK and just return a pointer to the element.

valid = .true.
nullify (ele)

if (ele_name == '') then
  if (ix_ele /= -1) ele => pointer_to_ele (lat, ix_ele, datum%ix_branch)
  return
endif

! Here if ele_name is not blank.
! Do some checking...

n_track = lat%branch(datum%ix_branch)%n_ele_track
n_max   = lat%branch(datum%ix_branch)%n_ele_max

if (ix_ele < 0 .or. ix_ele > n_max) then
  call out_io (s_error$, r_name, 'ELEMENT INDEX OUT OF RANGE! \i5\ ', ix_ele)
  call tao_set_invalid (datum, 'ELEMENT INDEX OUT OF RANGE FOR: ' // tao_datum_name(datum), why_invalid, .true.)
  valid = .false.
  return
endif

ele => pointer_to_ele (lat, ix_ele, datum%ix_branch)

if (ele%lord_status == multipass_lord$ .or. ele%lord_status == girder_lord$) then
  call out_io (s_error$, r_name, &
            'ELEMENT: ' // trim(ele%name) // &
            '    WHICH IS A: ' // control_name(ele%lord_status), &
            'CANNOT BE USED IN DEFINING A DATUM SINCE IT DOES NOT HAVE ', &
            '   A DEFINITE LOCATION IN THE LATTICE.', &
            'FOR DATUM: ' // tao_datum_name(datum) )
  call tao_set_invalid (datum, 'NO DEFINITE LOCATION IN LATTICE FOR: ' // tao_datum_name(datum), why_invalid, .true.)
  return
endif

end function tao_pointer_to_datum_ele

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_to_real (expression, value, err_flag)
!
! Mathematically evaluates an expression.
!
! Input:
!   expression    -- character(*): arithmetic expression
!
! Output:
!   value        -- real(rp): Value of arithmetic expression.
!   err_flag     -- Logical: TRUE on error.
!-

subroutine tao_to_real (expression, value, err_flag)

character(*) :: expression

real(rp) value
real(rp), allocatable :: vec(:)

logical err_flag

!

call tao_evaluate_expression (expression, 1, .false., vec, err_flag)
if (err_flag) return
value = vec(1)

end subroutine tao_to_real

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_expression (expression, n_size, use_good_user, value, err_flag, print_err, &
!                   info, stack, dflt_component, dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, &
!                   dflt_dat_or_var_index, dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)
!
! Mathematically evaluates a character expression.
!
! Input:
!   expression      -- character(*): Arithmetic expression.
!   n_size          -- integer: Size of the value array. If the expression evaluates to a
!                       a scalar, each value in the value array will get this value.
!                       If n_size = 0, the natural size is determined by the expression itself.
!   use_good_user   -- logical: Use the good_user logical in evaluating good(:)
!   print_err       -- logical, optional: If False then supress evaluation error messages.
!                       This does not affect syntax error messages. Default is True.
!   dflt_component  -- character(*), optional: Component to use if not specified in the expression. 
!                        'model' (default), 'base', or 'design'.
!   dflt_source     -- character(*), optional: Default source ('lat', 'data', etc.). Default is ''.
!   dflt_ele_ref    -- ele_struct, pointer, optional: Default reference element.
!   dflt_ele_start  -- ele_struct, pointer, optional: Default start element for ranges.
!   dflt_ele        -- ele_struct, pointer, optional: Default element to evaluate at.
!   dflt_dat_or_var_index -- character(*), optional: Default datum or variable index to use.
!   dflt_uni        -- integer, optional: Default universe to use. If 0 or not present, use viewed universe.
!   dflt_eval_point -- integer, optional: Default eval_point. anchor_end$ (default), anchor_center$, or anchor_beginning$.
!   dflt_s_offset   -- real(rp), optional: Default offset of eval_point. Default = 0.
!   dflt_orbit      -- coord_struct, optional: Default orbit to evaluate at.
!   datum           -- tao_data_struct, optional: If present, check to see that the expression does not depend upon
!                       a datum that will be evaluated after this datum. If so, this is an error.
!
! Output:
!   value(:)  -- Real(rp), allocatable: Value of arithmetic expression.
!   err_flag  -- Logical: True on an error. EG: Invalid expression.
!                  A divide by zero is not an error but good(:) will be set to False.
!   info(:)    -- tao_expression_info_struct, allocatable, optional: Is the value valid?, etc.
!                  Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                  orbit.x[23]|good_user is False.
!   stack(:)  -- Tao_eval_stack1_struct, allocatable, optional: Evaluation stack for the
!                  expression. This is useful to save if the same expression is
!                  to be evaluated repeatedly. 
!                  With this, tao_evaluate_stack can be called directly.
!-

recursive &
subroutine tao_evaluate_expression (expression, n_size, use_good_user, value, err_flag, print_err, &
                      info, stack, dflt_component, dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, &
                      dflt_dat_or_var_index, dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)

use random_mod
use expression_mod

type expression_func_struct
  character(12) :: name = ''      ! Name of function
  integer :: n_arg_target = 0     ! Number of arguments the function should have. -1 => 0 or 1 arg
  integer :: n_arg_count = 0      ! Number of arguments found.
end type

type (tao_eval_stack1_struct), allocatable :: stk(:)
type (tao_eval_stack1_struct), allocatable, optional :: stack(:)
type (ele_struct), optional, pointer :: dflt_ele_ref, dflt_ele_start, dflt_ele
type (coord_struct), optional :: dflt_orbit
type (tao_expression_info_struct), allocatable, optional :: info(:)
type (tao_data_struct), optional :: datum
type (expression_func_struct) func(0:20)

integer, optional :: dflt_uni, dflt_eval_point
integer, allocatable :: op(:)
integer i_lev, i_op, i, ios, n, n_size, ix0, ix1, ix2, ix3, ix4, n_func
integer ix_word, i_delim, i2, ix, ix_word2

real(rp), allocatable :: value(:)
real(rp), optional :: dflt_s_offset

character(*) :: expression
character(*), optional :: dflt_component, dflt_source
character(*), optional :: dflt_dat_or_var_index

character(len(expression)+20) :: phrase, word, word2
character(1) delim, cc
character(80) default_source
character(40) saved_prefix
character(*), parameter :: r_name = "tao_evaluate_expression"

logical delim_found, do_combine, use_good_user
logical err_flag, err, wild, printit, found, species_here
logical, optional :: print_err

! Don't destroy the input expression

err_flag = .true.
saved_prefix = ''
printit = logic_option(.true., print_err)
default_source = ''
if (present(dflt_source)) default_source = dflt_source

phrase = expression
if (len(phrase) > 11) then
  if (phrase(1:11) == 'expression:') phrase = phrase(12:)
endif

! if phrase is blank then return 0.0

call string_trim (phrase, phrase, ios)
if (ios == 0) then
  call out_io (s_warn$, r_name, "Expression is blank")
  call re_allocate (value, max(1, n_size))
  value = 0.0
  return
endif
 
! General idea: Create a reverse polish stack that represents the expression.
! Reverse polish means that the operand goes last so that 2 * 3 is written 
! on the stack as: [2, 3, *]

! The stack is called: stk
! Since operations move towards the end of the stack we need a separate
! stack called op which keeps track of what operations have not yet
! been put on stk.

! init

n_func = 0
i_lev = 0
i_op = 0

allocate (stk(20), op(20))

! parsing loop to build up the stack.

parsing_loop: do

  ! get a word

  call word_read (phrase, '+-*/()^,}[ ', word, ix_word, delim, delim_found, phrase)

  ! Args are counted counted at the beginning of the function and at each comma.

  if (n_func > 0 .and. (ix_word /= 0 .or. delim /= ')')) then
    if (func(n_func)%n_arg_count == 0) func(n_func)%n_arg_count = 1
  endif

  !--------------------------
  ! Preliminary: If we have split up something that should have not been split
  ! then put it back together again...

  ! just make sure we are not chopping a number in two, e.g. "3.5e-7" should not
  ! get split at the "-" even though "-" is a delimiter. Also "Q01+4[k1]" should not be split.

  do_combine = (delim == '-' .or. delim == '+') 
  if (do_combine .and. ix_word == 0) do_combine = .false.

  if (do_combine) then
    found = .false.   ! Found "NNN[" like construct where NNN is an integer?
    ix = index(phrase, '[')
    if (ix /= 0) then
      if (is_integer(phrase(:ix-1))) found = .true.
    endif

    if (.not. found) then ! Test if a number
      cc = upcase(word(ix_word:ix_word))
      if (cc == 'E' .or. cc == 'D') then
        do i = 1, ix_word-1
          if (index('.0123456789', word(i:i)) == 0) do_combine = .false.
        enddo
      else
        do_combine = .false.
      endif
    endif
  endif

  ! If still DO_COMBINE = True then we need to unsplit

  if (do_combine) then
    word = trim(word) // delim
    call word_read (phrase, '+-*/()^,}[', word2, ix_word2, delim, delim_found, phrase)
    word = trim(word) // word2
    ix_word = len_trim(word)
  endif

  ! Something like "lcav[lr(2).freq]" or "[2,4]@orbit.x[1,4] or "[ele::q20w[hkick], ele::q20w[hkick]]"
  ! will get split on the "["

  do
    if (delim /= '[') exit

    call word_read (phrase, ']', word2, ix_word2, delim, delim_found, phrase, ignore_interior = .true.)
    if (.not. delim_found) then
      call out_io (s_error$, r_name, "NO MATCHING ']' FOR OPENING '[':" // expression)
      return
    endif
    word = trim(word) // '[' // trim(word2) // ']'
    ix_word = len_trim(word)
    if (phrase == ' ') then  
      delim_found = .false.
      delim = ' '
    elseif (phrase(1:1) == ' ') then  
      call string_trim (phrase, phrase, ix)
      if (index('+-*/()^,}[', phrase(1:1)) == 0) then
        delim = ' '
      else
        delim = phrase(1:1)
        phrase = phrase(2:)
      endif
    else          ! even more...
      call word_read (phrase, '[+-*/()^,}', word2, ix_word2, delim, delim_found, phrase)
      word = trim(word) // trim(word2)       
      ix_word = len_trim(word)
    endif
  enddo

  ! If delim = "*" then see if this is being used as a wildcard
  ! Examples: "[*]|", "*.*|", "*.x|", "*@orbit.x|", "*@*|", "orbit.*[3]|", "ele::q*1[beta_a]", "var::*d|model"
  ! If so, we have split in the wrong place and we need to correct this. 
  ! Something like "3*[1,2]" or "3*.42" does not get split.

  do
    ix0 = index(word, '::')
    ix1 = index(phrase, '[')
    ix2 = index(phrase, ']')
    ix3 = index(phrase, '|')
    ix4 = index(word, '|')

    if (delim /= '*' .or. (phrase(1:1) == '[' .and. ((ix0 == 0) .eqv. (ix4 == 0)))) exit

    ! If in "[...*...]" construct is wild
    wild = .false.
    if (ix2 /= 0 .and. (ix1 == 0 .or. ix1 > ix2)) wild = .true.
    if (ix3 /= 0 .and. (ix1 == 0 .or. ix1 > ix3)) wild = .true.

    if (.not. wild) then
      select case (phrase(1:1))
      case ( ']', '[', '|', '@')
        wild = .true.
      case ('.')
        wild = (index('0123456789', phrase(2:2)) == 0) ! Wild if not a number
      case default
        ! If in "::xxx*yyy[" construct where each x and y is not one of ":", " ", etc.
        found = .false.
        if (ix0 /= 0 .and. ix1 /= 0) then
          do ix = ix0+2, ix_word
            if (index(': ]|', word(ix:ix)) /= 0) found = .true.
          enddo
          do ix = 1, ix1-1
            if (index(': ]|', phrase(ix:ix)) /= 0) found = .true.
          enddo
          if (.not. found) wild = .true.
        endif
      end select
    endif

    if (.not. wild) exit

    word = word(:ix_word) // '*'
    call word_read (phrase, '+-*/()^,}', word2, ix_word2, delim, delim_found, phrase, .true.)
    word = trim(word) // trim(word2)       
    ix_word = len_trim(word)
  enddo

  !---------------------------
  ! Now see what we got...

  ! For a word ending in '|' then must be a construct like 'orbit.x|-model'.
  ! So store the 'orbit.x|' prefix

  if (ix_word /= 0) then
    if (word(ix_word:ix_word) == '|') then
      saved_prefix = word
      word = ''
      ix_word = 0
    endif
  endif

  ! if the word is a datum without an "|", and dflt_component is present, 
  ! Then use the dflt_component.

  if (present(dflt_component) .and. index(word, '|') == 0) then
    call tao_find_data (err, word, print_err = .false.)
    if (.not. err) then
      phrase = trim(word) // '|' // trim(dflt_component) // delim // trim(phrase)
      cycle   ! Try again
    endif
  endif

  ! For a "(" delim we must have a function

  if (delim == '(') then

    if (ix_word /= 0) then
      word2 = word
      call downcase_string (word2)
      n_func = n_func + 1
      func(n_func) = expression_func_struct(word2, 1, 0)
      select case (word2)
      case ('cot');             call push_op_stack (op, i_op, cot$)
      case ('csc');             call push_op_stack (op, i_op, csc$)
      case ('sec');             call push_op_stack (op, i_op, sec$)
      case ('sin');             call push_op_stack (op, i_op, sin$)
      case ('sinc');            call push_op_stack (op, i_op, sinc$)
      case ('cos');             call push_op_stack (op, i_op, cos$)
      case ('tan');             call push_op_stack (op, i_op, tan$)
      case ('asin');            call push_op_stack (op, i_op, asin$)
      case ('acos');            call push_op_stack (op, i_op, acos$)
      case ('atan');            call push_op_stack (op, i_op, atan$)
      case ('atan2')
        call push_op_stack (op, i_op, atan2$)
        func(n_func)%n_arg_target = 2
      case ('modulo')
        call push_op_stack (op, i_op, modulo$)
        func(n_func)%n_arg_target = 2
      case ('sinh');            call push_op_stack (op, i_op, sinh$)
      case ('cosh');            call push_op_stack (op, i_op, cosh$)
      case ('tanh');            call push_op_stack (op, i_op, tanh$)
      case ('coth');            call push_op_stack (op, i_op, coth$)
      case ('asinh');           call push_op_stack (op, i_op, asinh$)
      case ('acosh');           call push_op_stack (op, i_op, acosh$)
      case ('atanh');           call push_op_stack (op, i_op, atanh$)
      case ('acoth');           call push_op_stack (op, i_op, acoth$)
      case ('abs');             call push_op_stack (op, i_op, abs$)
      case ('min');             call push_op_stack (op, i_op, min$)
      case ('max');             call push_op_stack (op, i_op, max$)
      case ('rms');             call push_op_stack (op, i_op, rms$)
      case ('average', 'mean'); call push_op_stack (op, i_op, average$)
      case ('sum');             call push_op_stack (op, i_op, sum$)
      case ('sqrt');            call push_op_stack (op, i_op, sqrt$)
      case ('log');             call push_op_stack (op, i_op, log$)
      case ('exp');             call push_op_stack (op, i_op, exp$)
      case ('factorial');       call push_op_stack (op, i_op, factorial$)
      case ('ran')         
        call push_op_stack (op, i_op, ran$)
        func(n_func)%n_arg_target = 0
      case ('ran_gauss')
        call push_op_stack (op, i_op, ran_gauss$)
        func(n_func)%n_arg_target = -1      ! 0 or 1 args
      case ('int');             call push_op_stack (op, i_op, int$)
      case ('sign');            call push_op_stack (op, i_op, sign$)
      case ('nint');            call push_op_stack (op, i_op, nint$)
      case ('floor');           call push_op_stack (op, i_op, floor$)
      case ('ceiling');         call push_op_stack (op, i_op, ceiling$)
      case ('mass_of');         call push_op_stack (op, i_op, mass_of$)
      case ('charge_of');       call push_op_stack (op, i_op, charge_of$)
      case ('anomalous_moment_of'); call push_op_stack (op, i_op, anomalous_moment_of$)
      case ('species');         call push_op_stack (op, i_op, species$)
      case ('antiparticle');    call push_op_stack (op, i_op, antiparticle$)
      case default
        call out_io (s_error$, r_name, 'UNEXPECTED CHARACTERS (BAD FUNCTION NAME?) BEFORE "(": ', 'IN EXPRESSION: ' // expression)
        return
      end select

      call push_op_stack (op, i_op, l_func_parens$)

    else
      call push_op_stack (op, i_op, l_parens$)
    endif

    cycle parsing_loop

  ! for a unary "-"

  elseif (delim == '-' .and. ix_word == 0) then
    call push_op_stack (op, i_op, unary_minus$)
    cycle parsing_loop

  ! for a unary "+"

  elseif (delim == '+' .and. ix_word == 0) then
    call push_op_stack (op, i_op, unary_plus$)
    cycle parsing_loop

  ! for a ")" delim

  elseif (delim == ')') then
    if (ix_word == 0) then
      if (n_func == 0 .or. (func(n_func)%n_arg_target /= 0 .and. func(n_func)%n_arg_target /= -1)) then
        if (printit) call out_io (s_error$, r_name, 'CONSTANT OR VARIABLE MISSING BEFORE ")"', &
                                                    'IN EXPRESSION: ' // expression)
        return
      endif

    else
      species_here = .false.
      if (i_op > 1) then
        select case(op(i_op-1))   ! op(i_op) will be l_func_parens$
        case (mass_of$, charge_of$, anomalous_moment_of$, antiparticle$, species$);  species_here = .true.
        end select
      endif

      if (species_here) then
        call push_stack (stk, i_lev, species_const$)
        stk(i_lev)%name = word
      else
        call push_stack (stk, i_lev, numeric$)
        call tao_param_value_routine (word, use_good_user, saved_prefix, stk(i_lev), err, printit, &
               dflt_component, default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, &
               dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)
        if (err) then
          if (printit) call out_io (s_error$, r_name, &
                          'ERROR IN EVALUATING EXPRESSION: ' // expression, &
                          'CANNOT EVALUATE: ' // word)
          return
        endif
      endif
    endif

    do
      do i = i_op, 1, -1       ! release pending ops
        if (op(i) == l_parens$ .or. op(i) == l_func_parens$) exit            ! break do loop
        call push_stack (stk, i_lev, op(i))
      enddo

      if (i == 0) then
        call out_io (s_error$, r_name, 'UNMATCHED ")" IN EXPRESSION: ' // expression)
        return
      endif

      i_op = i - 1

      if (op(i) == l_func_parens$) then
        if (func(n_func)%n_arg_target == -1) then
          if (func(n_func)%n_arg_count /= 0 .and. func(n_func)%n_arg_count /= 1) then
            if (printit) call out_io(s_error$, r_name, &
                            'FUNCTION: ' // trim(func(n_func)%name) // ' DOES NOT HAVE 0 OR 1 ARGUMENTS.', &
                            'IN EXPRESSION: ' // expression)
            return
          endif
          call push_stack (stk, i_lev, arg_count$)
          call re_allocate (stk(i_lev)%value, 1)
          stk(i_lev)%value(1) = func(n_func)%n_arg_count

        else
          if (func(n_func)%n_arg_count /= func(n_func)%n_arg_target) then
            if (printit) call out_io(s_error$, r_name, &
                          'FUNCTION: ' // trim(func(n_func)%name) // ' DOES NOT HAVE THE CORRECT NUMBER OF ARGUMENTS.', &
                          'IN EXPRESSION: ' // expression)
            return
          endif
        endif

        n_func = n_func - 1
      endif

      call word_read (phrase, '+-*/()^,}', word, ix_word, delim, delim_found, phrase)
      if (ix_word /= 0) then
        if (printit) call out_io (s_error$, r_name, 'UNEXPECTED CHARACTERS AFTER ")" IN EXPRESSION: ' // expression)
        return
      endif

      if (delim /= ')') exit    ! if no more ')' then no need to release more
    enddo


    if (delim == '(') then
      if (printit) call out_io (s_error$, r_name, '")(" CONSTRUCT DOES NOT MAKE SENSE IN EXPRESSION: ' // expression)
      return
    endif

  ! For binary "+-/*^" delims

  else
    if (ix_word == 0) then
      call out_io (s_error$, r_name, 'CONSTANT OR VARIABLE MISSING IN EXPRESSION: ' // expression)
      return
    endif
    call push_stack (stk, i_lev, numeric$)
    call tao_param_value_routine (word, use_good_user, saved_prefix, stk(i_lev), err, printit, &
            dflt_component, default_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, &
            dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit, datum)
    if (err) then
      if (printit) call out_io (s_error$, r_name, &
                        'ERROR IN EXPRESSION: ' // expression, &
                        'CANNOT EVALUATE: ' // word)
      return
    endif
  endif

  ! If we are here then we have an operation that is waiting to be identified

  select case (delim)
  case ('+')
    i_delim = plus$
  case ('-')
    i_delim = minus$
  case ('*')
    i_delim = times$
  case ('/')
    i_delim = divide$
  case (')')
    i_delim = r_parens$
  case ('^')
    i_delim = power$
  case (',')
    i_delim = comma$
    func(n_func)%n_arg_count = func(n_func)%n_arg_count + 1
  case ('}')
    i_delim = no_delim$
    call out_io (s_error$, r_name, &
                      'DELIMITOR FOUND OUT OF PLACE: ' // delim, &
                      'IN EXPRESSION: ' // expression)
    return
  case default
    if (delim_found) then
      if (delim == ' ') then
        if (printit) call out_io (s_error$, r_name, 'MALFORMED EXPRESSION: ' // expression)
        return
      endif

      call out_io (s_error$, r_name, 'INTERNAL ERROR')
      call err_exit
    endif
    i_delim = no_delim$
  end select

  ! Now see if there are operations on the OP stack that need to be transferred
  ! to the STK stack

  do i = i_op, 1, -1
    if (expression_eval_level(op(i)) < expression_eval_level(i_delim)) exit

    if (op(i) == l_parens$) then
      if (printit) call out_io (s_error$, r_name, 'UNMATCHED "(" IN EXPRESSION: ' // expression)
      return
    endif

    if (op(i) == l_func_parens$) then
      if (i_delim /= comma$) then
        if (printit) call out_io (s_error$, r_name, 'UNMATCHED "("', 'IN EXPRESSION: ' // expression)
        return
      endif
      i_op = i
      cycle parsing_loop
    endif

    call push_stack (stk, i_lev, op(i))
  enddo

  ! put the pending operation on the OP stack

  i_op = i
  select case (i_delim)
  case (no_delim$); exit parsing_loop
  case (comma$)
    if (printit) call out_io (s_error$, r_name, 'COMMA AT END OF EXPRESSION IS OUT OF place: ' // expression, &
                                   '(NEEDS "[...]" BRACKETS IF AN ARRAY.)')
    return
  case default; call push_op_stack (op, i_op, i_delim)
  end select

enddo parsing_loop

!------------------------------------------------------------------
! Some error checks

if (i_op /= 0) then
  call out_io (s_error$, r_name, 'UNMATCHED "(" IN EXPRESSION: ' // expression)
  return
endif

if (i_lev == 0) then
  call out_io (s_error$, r_name, 'NO VALUE FOUND IN EXPRESSION: ' // expression)
  return
endif

if (phrase /= '') then
  call out_io (s_error$, r_name, 'EXTRA STUFF AFTER EXPRESSION: ' // phrase)
  return
endif

call tao_evaluate_stack (stk(1:i_lev), n_size, use_good_user, value, err_flag, printit, expression, info)

! If the stack argument is present then copy stk to stack

if (present(stack)) then
  if (allocated(stack)) deallocate(stack)
  allocate (stack(i_lev))
  do i = 1, i_lev
    if (allocated (stk(i)%value)) then
      n = size(stk(i)%value)
      allocate (stack(i)%value(n), stack(i)%info(n))
      if (allocated (stack(i)%value_ptr)) allocate (stack(i)%value_ptr(n))
    endif
    stack(i) = stk(i)
  enddo
endif

!-------------------------------------------------------------------------
! The op_stack is for operators and functions.

contains

subroutine push_op_stack (op_stack, i_lev, this_type)

integer, allocatable :: op_stack(:)
integer i_lev, this_type

character(*), parameter :: r_name = "push_op_stack"

!

i_lev = i_lev + 1
if (i_lev > size(op_stack)) call re_allocate(op_stack, 2*i_lev)
op_stack(i_lev) = this_type

end subroutine push_op_stack

!-------------------------------------------------------------------------
! contains

subroutine push_stack (stack, i_lev, this_type)

type (tao_eval_stack1_struct), allocatable :: stack(:), tmp_stk(:)
integer i_lev, this_type

character(*), parameter :: r_name = "push_stack"

!

i_lev = i_lev + 1

if (i_lev > size(stack)) then
  call move_alloc(stack, tmp_stk)
  allocate(stack(2*i_lev))
  stack(1:i_lev-1) = tmp_stk
endif

stack(i_lev)%type = this_type
stack(i_lev)%name = expression_op_name(this_type)
stack(i_lev)%scale = 1

end subroutine push_stack
                       
end subroutine tao_evaluate_expression

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

recursive &
subroutine tao_param_value_routine (str, use_good_user, saved_prefix, stack, err_flag, print_err, dflt_component, &
                    dflt_source, dflt_ele_ref, dflt_ele_start, dflt_ele, dflt_dat_or_var_index, dflt_uni, &
                    dflt_eval_point, dflt_s_offset, dflt_orbit, datum)

type (tao_eval_stack1_struct) stack
type (tao_eval_stack1_struct), allocatable :: stack2(:)
type (tao_real_pointer_struct), allocatable :: re_array(:)
type (tao_data_array_struct), allocatable :: d_array(:)
type (tao_integer_array_struct), allocatable :: int_array(:)
type (tao_var_array_struct), allocatable :: v_array(:)
type (tao_data_struct), pointer :: d
type (tao_var_struct), pointer :: v
type (tao_lattice_struct), pointer :: tao_lat
type (ele_struct), pointer, optional :: dflt_ele_ref, dflt_ele_start, dflt_ele
type (coord_struct), optional :: dflt_orbit
type (tao_expression_info_struct), allocatable :: info(:)
type (tao_data_struct), optional :: datum

real(rp), optional :: dflt_s_offset
real(rp), allocatable :: value(:)

integer, optional :: dflt_uni, dflt_eval_point
integer ios, i, m, n, ix, ix2, ix_word, ix_uni

character(*) str, saved_prefix
character(*), optional :: dflt_source, dflt_component
character(*), optional :: dflt_dat_or_var_index

character(1) delim
character(16) s_str, source
character(60) name
character(200) str2, word2
character(*), parameter :: r_name = 'tao_param_value_routine'

logical use_good_user, err_flag, print_err, print_error, delim_found, valid_value, exterminate

! See if it is a constant like pi, etc.

err_flag = .false.

call match_word(str, physical_const_list%name, ix, .true., .false.)
if (ix > 0) then
  call re_allocate(stack%value, 1)
  stack%value(1) = physical_const_list(ix)%value
  call tao_re_allocate_expression_info (stack%info, 1)
  stack%info%good = .true.
  return
endif

! Named constants

if (allocated(s%com%symbolic_num)) then
  call match_word(str, s%com%symbolic_num%name, ix, .true., .false.)
  if (ix > 0) then
    call re_allocate(stack%value, 1)
    stack%value(1) = s%com%symbolic_num(ix)%value
    call tao_re_allocate_expression_info (stack%info, 1)
    stack%info%good = .true.
    return
  endif
endif

! An array "[...]", but not "[1,2]@ele::q[k1]"

if (str(1:1) == '[' .and. index(str, ']@') == 0) then
  n = len_trim(str)
  if (str(n:n) /= ']') then
    if (print_err) call out_io (s_error$, r_name, "Malformed array: " // str)
    err_flag = .true.
    return
  endif

  str2 = str(2:n-1)
  call re_allocate (stack%value, 0)

  do
    call word_read (str2, ',', word2, ix_word, delim, delim_found, str2, ignore_interior = .true.)
    call tao_evaluate_expression (word2, 1, use_good_user, value, err_flag, print_err, &
                         info, stack2, dflt_component, dflt_source, dflt_ele_ref, dflt_ele_start, &
                         dflt_ele, dflt_dat_or_var_index, dflt_uni, dflt_eval_point, dflt_s_offset, dflt_orbit)
    if (err_flag) return
    m = size(stack%value)
    call re_allocate(stack%value, m+1)
    call tao_re_allocate_expression_info (stack%info, m+1)
    stack%value(m+1) = value(1)
    stack%info(m+1)%good = .true.
    do i = 1, size(stack2)
      if (.not. allocated(stack2(i)%info)) cycle
      stack%info(m+1)%good = (stack%info(m+1)%good .and. stack2(i)%info(1)%good)
    enddo
    if (.not. delim_found) return
  enddo

endif

! Case where str represents a number.

if (is_real(str)) then
  call re_allocate(stack%value, 1)
  read (str, *, iostat = ios) stack%value(1)
  if (ios /= 0) then
    if (print_err) call out_io (s_error$, r_name, "This doesn't seem to be a number: " // str)
    err_flag = .true.
  endif
  call tao_re_allocate_expression_info (stack%info, 1)
  stack%info%good = (ios == 0)
  return
endif

! Case where str is a variable name.
! Remember the last string so 'orbit.x|meas-ref' translates to 'orbit.x|meas - orbit.x|ref.'

ix = index(str, '|')
stack%name = str

if (ix == 0) then
  select case (str)
  case ('model', 'design', 'base', 'meas', 'ref', 'old', 'fit')
    stack%name = trim(saved_prefix) // str
  end select
else
  saved_prefix = str(1:ix)
endif

name = stack%name

allocate (re_array(0))
allocate (d_array(0))
allocate (int_array(0))
allocate (v_array(0))

! Decide data source

source = dflt_source
ix = index(name, '::')
if (ix /= 0) then
  s_str = name(max(1,ix-4):ix-1)
  if (s_str(1:1) == '@') s_str = s_str(2:)
  if (s_str == 'lat' .or. s_str == 'ele' .or. s_str == 'data' .or. s_str == 'var') then
    source = s_str
  else if (name(max(1,ix-7):ix-1) == 'ele_mid') then
    source = 'ele'
  else if (name(max(1,ix-9):ix-1) == 'ele_begin') then
    source = 'ele'
  elseif (name(max(1, ix-4):ix-1) == 'beam') then
    source = 'beam'
  endif
endif

! source = "at_ele" is used for plotting. In this case, dflt_ele and dflt_orbit will be present.

if (source == 'at_ele') then
  call re_allocate(stack%value, 1)
  stack%value(1) = tao_param_value_at_s (str, dflt_ele, dflt_orbit, err_flag)
  call tao_re_allocate_expression_info (stack%info, 1)
  stack%info%good = (.not. err_flag)
  return
endif

! Look for a lat datum.

if (source == 'lat' .or. source == 'beam') then
  call tao_evaluate_lat_or_beam_data (err_flag, name, stack%value, print_err, dflt_source, dflt_ele_ref, &
                              dflt_ele_start, dflt_ele, dflt_component, dflt_uni, dflt_eval_point, dflt_s_offset)
  call tao_re_allocate_expression_info (stack%info, size(stack%value))
  stack%info%good = (.not. err_flag)
  stack%type = lat_num$
  return

! Look for a lattice element parameter 

elseif (source == 'ele') then
  call tao_evaluate_element_parameters (err_flag, name, stack%value, print_err, dflt_ele, &
                                                  dflt_source, dflt_component, dflt_uni, dflt_eval_point, stack%info)
  call tao_re_allocate_expression_info (stack%info, size(stack%value))
  stack%info%good = (.not. err_flag)
  stack%type = ele_num$
  return

! Look for variable or data values

else

  err_flag = .true.
  print_error = print_err
  if (source == '') print_error = .false. ! Don't generate unnecessary messages

  if (index(name, '|') == 0 .and. present(dflt_component)) name = trim(name) // '|' // trim(dflt_component)
  
  if (source == 'var' .or. source == '') then
    call tao_find_var (err_flag, name, v_array = v_array,  re_array = re_array, print_err = print_error, dflt_var_index = dflt_dat_or_var_index)
    stack%type = var_num$
  endif

  if (source == 'data' .or. (err_flag .and. source == '')) then
    call tao_find_data (err_flag, name, d_array = d_array, re_array = re_array, int_array = int_array, &
                        dflt_index = dflt_dat_or_var_index, print_err = print_error, ix_uni = dflt_uni)
    stack%type = data_num$
    ! Error if datum associated with the expression is evaluated before the datum evaluated here.
    if (present(datum) .and. .not. err_flag) then
      ! Only check if this is a user defined datum (not temp datum used for plotting).
      if (datum%ix_uni > 0 .and. datum%ix_data > 0 .and. substr(d_array(1)%d%data_type,1,11) == 'expression:') then
        if (datum%ix_uni < d_array(1)%d%ix_uni .or. (datum%ix_uni == d_array(1)%d%ix_uni .and. datum%ix_data < d_array(1)%d%ix_data)) then
          err_flag = .true.
          if (print_err) call out_io (s_error$, r_name, 'THE EXPRESSION ASSOCIATED WITH DATUM: ' // tao_datum_name(datum), &
                          'DEPENDS UPON AN EXPRESSION  DATUM (' // trim(tao_datum_name(d_array(1)%d)) // ') WHICH IS EVALUATED AFTER ' // tao_datum_name(datum))
          return
        endif
      endif
    endif
  endif

  if (err_flag) then
    if (source == '' .and. print_err) call out_io (s_error$, r_name, 'CANNOT EVALUATE AS DATUM OR VARIABLE VALUE: ' // name)
    return
  endif

endif

! Now transfer the information to the stack

if (size(re_array) /= 0) then
  n = size(re_array)
  if (allocated(stack%value_ptr)) then
    if (size(stack%value_ptr) /= n) deallocate (stack%value_ptr)
  endif

  if (.not. allocated(stack%value_ptr)) allocate (stack%value_ptr(n))
  call re_allocate (stack%value, n)
  call tao_re_allocate_expression_info (stack%info, n)

  stack%scale =  1

  if (index(name, 'ping_a.') /= 0 .and. index(name, 'ping_a.phase') == 0) then
    ix_uni = d_array(1)%d%d1%d2%ix_universe
    if (index(name, '|meas') /= 0) then
      stack%scale = s%u(ix_uni)%ping_scale%a_mode_meas
    elseif (index(name, '|ref') /= 0) then
      stack%scale = s%u(ix_uni)%ping_scale%a_mode_ref
    endif
  endif

  if (index(name, 'ping_b.') /= 0 .and. index(name, 'ping_b.phase') == 0) then
    ix_uni = d_array(1)%d%d1%d2%ix_universe
    if (index(name, '|meas') /= 0) then
      stack%scale = s%u(ix_uni)%ping_scale%b_mode_meas
    elseif (index(name, '|ref') /= 0) then
      stack%scale = s%u(ix_uni)%ping_scale%b_mode_ref
    endif
  endif

  do i = 1, n
    stack%value(i) = re_array(i)%r
    stack%value_ptr(i)%r => re_array(i)%r
    stack%info(i)%ele => null()

    ! good is only used with data and not variables
    if (associated(re_array(i)%good_value)) then
      stack%value_ptr(i)%good_user => re_array(i)%good_user
      stack%value_ptr(i)%good_value => re_array(i)%good_value
    else
      stack%info(i)%good = .true.
    endif

    select case (stack%type)
    case (var_num$)
      v => v_array(i)%v
      if (.not. v_array(i)%v%exists) cycle
      if (v%slave(1)%ix_ele >= 0) then
        stack%info(i)%ele     => s%u(v%slave(1)%ix_uni)%model%lat%branch(v%slave(1)%ix_branch)%ele(v%slave(1)%ix_ele)
      endif
      stack%info(i)%s      = v%s
      stack%info(i)%good   = v%exists
      stack%value_ptr(i)%good_user => v%good_user

    case (data_num$)
      d => d_array(i)%d
      if (.not. d%exists) cycle
      if (d%ix_ele >= 0) then
        stack%info(i)%ele  => s%u(d%ix_uni)%model%lat%branch(d%ix_branch)%ele(d%ix_ele)
      endif
      stack%info(i)%s      = d%s
      stack%info(i)%good   = d%exists
      stack%value_ptr(i)%good_user => d%good_user
    end select
  enddo

elseif (size(int_array) /= 0) then
  n = size(int_array)
  call re_allocate (stack%value, n)
  call tao_re_allocate_expression_info (stack%info, n)
  do i = 1, n
    stack%value(i) = int_array(i)%i
    stack%info(i)%good  = .true.
  enddo

else
  if (print_err) call out_io (s_error$, r_name, &
               'THIS IS NOT A DATUM OR A VARIABLE VALUE: ' // name, &
               '[PERHAPS MISSING "|<component>" SUFFIX.]')
  err_flag = .true.
  return
endif

end subroutine tao_param_value_routine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function tao_evaluate_datum_at_s (datum, tao_lat, ele, ele_ref, valid_value, err_str, bad_datum) result(value)
!
! Routine to evaluate a datum at a given s-position in the lattice
!
! Input:
!   datum         -- tao_data_struct: Datum to evaluate.
!   tao_lat       -- tao_lattice_struct: 
!   ele           -- ele_struct, pointer: Evaluation element.
!   ele_ref       -- ele_struct, pointer: Reference element.
!   valid_value   -- logical: True if evaluation was sucessful. False if not.
!
! Output:
!   err_str       -- character(*): Error string for printing an error message.
!   bad_datum     -- logical: True -> datum is malformed. False -> Could evaluate or evaluation problem was not due to the datum
!                               itself (EG: the lattice was unstable).
!   value         -- real(rp): Datum value.
!-

function tao_evaluate_datum_at_s (datum, tao_lat, ele, ele_ref, valid_value, err_str, bad_datum) result(value)

use twiss_and_track_mod, only: twiss_and_track_at_s
use transfer_map_mod, only: mat6_from_s_to_s

type (tao_data_struct) datum
type (tao_lattice_struct), target :: tao_lat
type (ele_struct), pointer :: ele, ele_ref
type (branch_struct), pointer :: branch
type (ele_struct) ele_at_s
type (coord_struct) :: orb_at_s, orb2
type (coord_struct), pointer :: orbit(:)

real(rp) s_offset, value
real(rp) s_eval_ref
integer ix_ref
logical valid_value, bad_datum, compute_floor, err
character(*) err_str

!

compute_floor = (datum%data_type(1:min(5,len(datum%data_type))) == 'floor')
branch => pointer_to_branch(ele)
orbit => tao_lat%tao_branch(ele%ix_branch)%orbit
bad_datum = .false.

valid_value = .false.
value = real_garbage$
s_eval_ref = branch%ele(0)%s
ix_ref = 0
if (associated(ele_ref)) ix_ref = ele_ref%ix_ele

if (.not. associated(ele)) then
  err_str = 'THERE MUST BE AN ASSOCIATED ELEMENT WHEN S_OFFSET IS NON-ZERO OR EVAL_POINT != END.'
  bad_datum = .true.
  return
endif

select case (datum%eval_point)
case (anchor_beginning$)
  datum%s = ele%s_start + datum%s_offset
  if (associated(ele_ref)) s_eval_ref = ele_ref%s_start
case (anchor_center$)
  datum%s = (ele%s_start + ele%s) / 2 + datum%s_offset
  if (associated(ele_ref)) s_eval_ref = (ele_ref%s_start + ele_ref%s) / 2
case (anchor_end$)
  datum%s = ele%s + datum%s_offset
  if (associated(ele_ref)) s_eval_ref = ele_ref%s
end select

!--------------------------------------------

if (substr(datum%data_type,1,2) == 'r.') then
  orb_at_s = orbit(ix_ref)
  call mat6_from_s_to_s (branch%lat, ele_at_s%mat6, ele_at_s%vec0, s_eval_ref, datum%s, &
                                                       orbit(ele%ix_ele), orb2, branch%ix_branch, .true.)
  value = tao_param_value_at_s (datum%data_type, ele_at_s, orb_at_s, err, bad_datum = bad_datum)
  if (err) then
    err_str = 'CANNOT EVALUATE DATUM AT OFFSET POSITION.'
    return
  endif

else
  call twiss_and_track_at_s (branch%lat, datum%s, ele_at_s, orbit, orb_at_s, branch%ix_branch, &
                                                                      err, compute_floor_coords = compute_floor)
  if (err) then
    err_str = 'CANNOT TRACK TO OFFSET POSITION.'
    return
  endif

  value = tao_param_value_at_s (datum%data_type, ele_at_s, orb_at_s, err, bad_datum = bad_datum)
  if (err) then
    err_str = 'CANNOT EVALUATE DATUM AT OFFSET POSITION.'
    return
  endif

  if (associated(ele_ref)) then
    call twiss_and_track_at_s (branch%lat, s_eval_ref, ele_at_s, orbit, orb_at_s, branch%ix_branch, &
                                                                      err, compute_floor_coords = compute_floor)
    if (err) then
      err_str = 'CANNOT TRACK TO REFERENCE POSITION.'
      return
    endif

    value = value - tao_param_value_at_s (datum%data_type, ele_at_s, orb_at_s, err, bad_datum = bad_datum)
    if (err) then
      err_str = 'CANNOT EVALUATE DATUM AT REFERENCE POSITION.'
      return
    endif
  endif
endif

valid_value = .true.

end function tao_evaluate_datum_at_s

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine tao_evaluate_stack (stack, n_size_in, use_good_user, value, info, err_flag, print_err, expression)
!
! Routine to evaluate an expression stack.
!
! Input:
!   stack(:)      -- tao_eval_stack1_struct: Expression stack
!   n_size_in     -- integer: Desired array size. If the expression evaluates to a
!                      a scalar, each value in the value array will get this value.
!                      If n_size = 0, the natural size is determined by the expression itself.
!   use_good_user -- logical: Use the good_user logical in evaluating good(:)
!   print_err     -- logical: If False then supress evaluation error messages.
!                      This does not affect syntax error messages. Default is True.
!   expression    -- character(*): Original expression. Used for error messages.
!
! Output:
!   value(:)      -- Real(rp), allocatable: Value of arithmetic expression.
!   info(:)       -- tao_expression_info_struct, allocatable: Is the value valid? 
!                      Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                      orbit.x[23]|good_user is False.
!   err_flag      -- Logical: True on error. False otherwise
!-

subroutine tao_evaluate_stack (stack, n_size_in, use_good_user, value, err_flag, print_err, expression, info_in)

use expression_mod

type (tao_eval_stack1_struct), target :: stack(:)
type (tao_eval_stack1_struct), pointer :: ss
type (tao_eval_stack1_struct), pointer :: s(:)
type (tao_eval_stack1_struct) stk2(20)
type (tao_expression_info_struct), allocatable, optional :: info_in(:)
type (tao_expression_info_struct), allocatable :: info(:)

real(rp), allocatable :: value(:)

integer n_size_in, species
integer i, i2, j, n, ns, ni, n_size

logical err_flag, use_good_user, print_err, info_allocated

character(*) expression
character(*), parameter :: r_name = 'tao_evaluate_stack'

! Calculate good

s => stack   ! For debugging purposes
err_flag = .true.
n_size = max(1, n_size_in)

do i = 1, size(stack)
  ss => stack(i)

  select case (ss%type)
  case (average$, sum$, rms$, min$, max$); n_size = 1
  end select

  if (allocated(ss%value)) then
    if (size(ss%value) > 1 .and. n_size > 1 .and. size(ss%value) /= n_size) then
      if (print_err) call out_io (s_error$, r_name, 'Array sizes mismatch in expression')
      err_flag = .true.
      return
    endif
    n_size = max(n_size, size(ss%value))
  endif

  if (allocated(ss%value_ptr)) then
    if (associated(ss%value_ptr(1)%good_value)) then    
      do j = 1, size(ss%value_ptr)
        if (use_good_user) then
          ss%info(j)%good = ss%value_ptr(j)%good_value .and. ss%value_ptr(j)%good_user
        else
          ss%info(j)%good = ss%value_ptr(j)%good_value
        endif
      enddo
    endif
  endif
enddo

call tao_re_allocate_expression_info(info, n_size)

! Go through the stack and perform the operations...

i2 = 0  ! stack pointer
do i = 1, size(stack)

  if (allocated(stack(i)%info)) then
    ns = size(stack(i)%info)
    if (.not. allocated(info)) then
      call tao_re_allocate_expression_info(info, ns)
      info = tao_expression_info_struct()
    endif
    ni = size(info)

    if (.not. this_size_check(ns, ni)) return
    if (ns == ni) then
      info%good = info%good .and. stack(i)%info%good
      do j = 1, size(info)
        if (stack(i)%info(j)%s /= real_garbage$) info(j)%s = stack(i)%info(j)%s
        info(j)%ele => stack(i)%info(j)%ele
      enddo 
    elseif (ns == 1) then
      info%good = info%good .and. stack(i)%info(1)%good
    elseif (ni == 1) then
      call tao_re_allocate_expression_info(info, ns)
      info%good = info(1)%good .and. stack(i)%info(1)%good
    endif
  endif

  !

  select case (stack(i)%type)
  case (arg_count$)
    cycle

  case (numeric$)
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (species_const$) 
    i2 = i2 + 1
    stk2(i2)%name = stack(i)%name
    call re_allocate(stk2(i2)%value, 1)

  case (lat_num$, ele_num$)
    !!! This needs to be fixed to include default stuff
    !!! call tao_param_value_routine (stack(i)%name, '', stack(i), err_flag, print_err)
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (var_num$, data_num$)
    do j = 1, size(stack(i)%value)
      stack(i)%value(j) = stack(i)%value_ptr(j)%r
    enddo
    stack(i)%value = stack(i)%value * stack(i)%scale
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (unary_minus$) 
    stk2(i2)%value = -stk2(i2)%value

  case (unary_plus$) 
    ! Nothing to do

  case (plus$)
    if (.not. this_size_check(size(stk2(i2)%value), size(stk2(i2-1)%value))) return
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value + stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) + stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value + stk2(i2)%value
    endif
    i2 = i2 - 1

  case (minus$)
    if (.not. this_size_check(size(stk2(i2)%value), size(stk2(i2-1)%value))) return
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value - stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) - stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value - stk2(i2)%value
    endif
    i2 = i2 - 1

  case (times$)
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value * stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) * stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value * stk2(i2)%value
    endif
    i2 = i2 - 1

  case (divide$)
    if (.not. this_size_check(size(stk2(i2)%value), size(stk2(i2-1)%value))) return
    n = size(stk2(i2)%value)
    do j = 1, n
      if (stk2(i2)%value(j) /= 0) cycle
      if (print_err) call out_io (s_error$, r_name, 'Divide by zero!')
      err_flag = .true.
      return
      ! Propably can get rid of this stuff...
      stk2(i2)%value(j) = 1
      if (n == 1) then
        info%good = .false.  ! All are false
      else
        info(j)%good = .false.
      endif
    enddo

    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value / stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) / stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value / stk2(i2)%value
    endif
    i2 = i2 - 1

  case (power$)
    if (.not. this_size_check(size(stk2(i2)%value), size(stk2(i2-1)%value))) return
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value ** stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) ** stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value ** stk2(i2)%value
    endif
    i2 = i2 - 1

  case (cot$)
    stk2(i2)%value = 1.0_rp / tan(stk2(i2)%value)

  case (csc$) 
    stk2(i2)%value = 1.0_rp / sin(stk2(i2)%value)

  case (sec$)
    stk2(i2)%value = 1.0_rp / cos(stk2(i2)%value)

  case (sin$) 
    stk2(i2)%value = sin(stk2(i2)%value)

  case (sinc$) 
    stk2(i2)%value = sinc(stk2(i2)%value)

  case (cos$) 
    stk2(i2)%value = cos(stk2(i2)%value)

  case (tan$) 
    stk2(i2)%value = tan(stk2(i2)%value)

  case (asin$) 
    stk2(i2)%value = asin(stk2(i2)%value)

  case (acos$) 
    stk2(i2)%value = acos(stk2(i2)%value)

  case (atan$) 
    stk2(i2)%value = atan(stk2(i2)%value)

  case (atan2$) 
    stk2(i2-1)%value = atan2(stk2(i2-1)%value, stk2(i2)%value)
    i2 = i2 - 1

  case (modulo$) 
    stk2(i2-1)%value = modulo(stk2(i2-1)%value, stk2(i2)%value)
    i2 = i2 - 1

  case (sinh$)
    stk2(i2)%value = sinh(stk2(i2)%value)

  case (cosh$) 
    stk2(i2)%value = cosh(stk2(i2)%value)

  case (tanh$) 
    stk2(i2)%value = tanh(stk2(i2)%value)

  case (coth$)
    stk2(i2)%value = 1.0_rp / tanh(stk2(i2)%value)

  case (asinh$) 
    stk2(i2)%value = asinh(stk2(i2)%value)

  case (acosh$) 
    stk2(i2)%value = acosh(stk2(i2)%value)

  case (atanh$) 
    stk2(i2)%value = atanh(stk2(i2)%value)

  case (acoth$)
    stk2(i2)%value = 1.0_rp / atanh(stk2(i2)%value)

  case (abs$) 
    stk2(i2)%value = abs(stk2(i2)%value)

  case (rms$)
    stk2(i2)%value(1) = rms_value(stk2(i2)%value, info%good)
    call re_allocate(stk2(i2)%value, 1)
    info(1)%good = any(info%good)
    call tao_re_allocate_expression_info(info, 1)

  case (average$)
    if (any(info%good)) then
      stk2(i2)%value(1) = sum(stk2(i2)%value, mask = info%good) / count(info%good)
    endif
    call re_allocate(stk2(i2)%value, 1)
    info(1)%good = any(info%good)
    call tao_re_allocate_expression_info(info, 1)

  case (sum$, min$, max$)
    select case (stack(i)%type)
    case (sum$); stk2(i2)%value(1) = sum(stk2(i2)%value, mask = info%good)
    case (min$); stk2(i2)%value(1) = minval(stk2(i2)%value, mask = info%good)
    case (max$); stk2(i2)%value(1) = maxval(stk2(i2)%value, mask = info%good)
    end select
    call re_allocate(stk2(i2)%value, 1)
    info(1)%good = .true.
    call tao_re_allocate_expression_info(info, 1)

  case (sqrt$) 
    stk2(i2)%value = sqrt(stk2(i2)%value)

  case (log$) 
    stk2(i2)%value = log(stk2(i2)%value)

  case (exp$) 
    stk2(i2)%value = exp(stk2(i2)%value)

  case (factorial$) 
    do n = 1, size(stk2(i2)%value)
      stk2(i2)%value(n) = factorial(nint(stk2(i2)%value(n)))
    enddo

  case (ran$) 
    i2 = i2 + 1
    call re_allocate(stk2(i2)%value, n_size)
    call ran_uniform(stk2(i2)%value)

    if (size(info) == 1) then
      call tao_re_allocate_expression_info(info, n_size)
      info%good = info(1)%good
    endif

  case (ran_gauss$) 
    if (nint(stack(i-1)%value(1)) == 0) then
      i2 = i2 + 1
      call re_allocate(stk2(i2)%value, n_size)
      call ran_gauss(stk2(i2)%value)
    else
      call re_allocate(value, n_size)
      call ran_gauss(value, sigma_cut = stk2(i2)%value(1))
      call re_allocate(stk2(i2)%value, n_size)
      stk2(i2)%value = value
    endif

    if (size(info) == 1) then
      call tao_re_allocate_expression_info(info, n_size)
      info%good = info(1)%good
    endif

  case (int$)
    stk2(i2)%value = int(stk2(i2)%value)

  case (sign$)
    stk2(i2)%value = sign_of(stk2(i2)%value)

  case (nint$)
    stk2(i2)%value = nint(stk2(i2)%value)

  case (floor$)
    stk2(i2)%value = floor(stk2(i2)%value)

  case (ceiling$)
    stk2(i2)%value = ceiling(stk2(i2)%value)

  case (mass_of$, charge_of$, anomalous_moment_of$)
    species = species_id(stk2(i2)%name)
    if (species == invalid$) then
      if (print_err) call out_io (s_error$, r_name, 'Not a valid species name: ' // stk2(i2)%name)
      err_flag = .true.
      return
    endif
    select case (stack(i)%type)
    case (mass_of$);              stk2(i2)%value = mass_of(species)
    case (charge_of$);            stk2(i2)%value = charge_of(species)
    case (anomalous_moment_of$);  stk2(i2)%value = anomalous_moment_of(species)
    end select

  case default
    call out_io (s_error$, r_name, 'INTERNAL ERROR')
    call err_exit
  end select
enddo

!

if (i2 /= 1) then
  call out_io (s_error$, r_name, 'INTERNAL ERROR')
  call err_exit
endif

if (size(stk2(1)%value) == 1 .and. n_size_in > 1) then
  call re_allocate(value, n_size_in)
  value = stk2(1)%value(1)
  if (.not. info(1)%good) value = 0
elseif (size(stk2(1)%value) > 1 .and. size(info) == 1) then
  call value_transfer (value, stk2(1)%value)
  if (.not. info(1)%good) value = 0
else
  call value_transfer (value, stk2(1)%value)
  where (.not. info%good) value = 0
endif

if (present(info_in)) then
  if (allocated(info_in)) deallocate(info_in)
  info_in = info
endif

n_size = size(value)
if (n_size_in /= 0) then
  if (n_size /= 1 .and. n_size_in /= n_size) then
    call out_io (s_error$, r_name, 'ARRAY SIZE MISMATCH FROM WHAT IS EXPECTED IN EXPRESSION: ' // expression)
    return
  endif
endif

err_flag = .false.

!-------------------------------------------------------------------------
contains

subroutine value_transfer (to_array, from_array)

real(rp), allocatable :: to_array(:)
real(rp) from_array(:)

!

call re_allocate (to_array, size(from_array))
to_array = from_array

end subroutine value_transfer

!-------------------------------------------------------------------------
! contains

function this_size_check (isize1, isize2) result (ok)
integer isize1, isize2
logical ok
ok = (isize1 == 1 .or. isize2 == 1 .or. isize1 == isize2)
if (.not. ok) then
  call out_io (s_error$, r_name, 'ARRAY SIZE MISMATCH IN EXPRESSION: ' // expression)
endif
end function this_size_check

end subroutine tao_evaluate_stack 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_to_int (str, i_int, err)
! 
! Converts a string to an integer
!
! If the string str is blank then i_int = 0
!-

subroutine tao_to_int (str, i_int, err)

character(*) str
integer ios, i_int
logical err
character(12) :: r_name = "tao_to_int"

!

call string_trim (str, str, ios)
if (ios .eq. 0) then
  i_int = 0
  return
endif
 
err = .false.
read (str, *, iostat = ios) i_int

if (ios /= 0) then
  call out_io (s_error$, r_name, 'EXPECTING INTEGER: ' // str)
  err = .true.
  return
endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_set_invalid (datum, message, why_invalid, exterminate, err_level)
!
! Routine to either print an error message to the terminal (if why_invalid
! is not present) or set the why_invalid string to the error message.
!
! Note: The exterminate argument should be set to False if the datum is invalid for
! reasons like the beam has been lost. In this case, the datum could possibly
! become valid in the future. Exterminate should be set to True if the datum could
! not possibly become valid. For example, the datum's reference element does not
! exist in the lattice.
!
! Input:
!   datum       -- tao_data_struct: Bad datum.
!   message     -- character(*): Error message.
!   exterminate -- logical, optional: Default is False. If True, set datum%exists
!                   to False so that Tao will ignore this datum from now on.
!   err_level   -- integer, optional: s_error$ (default), s_warn$, etc.
!
! Output:
!   why_invalid -- character(*), optional: Set to message if present.
!-

Subroutine tao_set_invalid (datum, message, why_invalid, exterminate, err_level)

type (tao_data_struct) datum
type (tao_universe_struct), pointer :: u

logical, optional :: exterminate
logical identified_err_found, found

integer, optional :: err_level
integer i

character(*) message
character(*), optional :: why_invalid
character(*), parameter :: r_name = 'tao_set_invalid'
character(40) :: err_str(2) = [character(40):: 'NO BEAM TRACKING HAS BEEN DONE', &
                                               'CANNOT EVALUATE DUE TO PARTICLE LOSS']


! The idea with err_str is to limit the number of error messages generated of a given type.

datum%why_invalid = message

if (logic_option(.false., exterminate)) then
  datum%exists = .false. 
endif

if (present(why_invalid)) then
  why_invalid = message
elseif (.not. datum%err_message_printed) then
  datum%err_message_printed = .true.
  identified_err_found = .false.
  do i = 1, size(err_str)
    if (index(message, trim(err_str(i))) == 0) cycle
    identified_err_found = .true.
    if (s%com%is_err_message_printed(i)) return
    s%com%is_err_message_printed(i) = .true.
    identified_err_found = .true.
  enddo

  s%com%n_err_messages_printed = s%com%n_err_messages_printed + 1
  if (s%com%n_err_messages_printed == s%global%datum_err_messages_max + 1) then
    call out_io (s_warn$, r_name, 'NUMBER OF ERROR MESSAGES EXCEEDS GLOBAL%DATUM_ERR_MESSAGES_MAX.', &
                                  'WILL NOT PRINT ANY MORE DATUM ERROR MESSAGES FOR THIS EVALUATION CYCLE.')
  endif
  if (s%com%n_err_messages_printed > s%global%datum_err_messages_max) return

  call out_io (integer_option(s_error$, err_level), r_name, message, &
                'FOR DATUM: ' // trim(tao_datum_name(datum)) // ' with data_type: ' // datum%data_type)

  if (index(upcase(message), 'UNSTABLE') /= 0) then
    u => tao_pointer_to_universe(datum%ix_uni)
    found = .false.
    do i = 1, size(u%data)
      if (substr(u%data(i)%data_type,1,8) == 'unstable') found = .true.
    enddo
    if (.not. found) call out_io(s_info$, r_name, &
              'NO unstable.orbit, unstable.ring, nor unstable.eigen FOR THE UNIVERSE WITH THE PROBLEM EXISTS.', &
              'YOU MIGHT WANT TO CONSIDER ADDING SUCH A DATUM.')
  endif

  if (identified_err_found) then
    call out_io (s_warn$, r_name, 'WILL NOT PRINT ANY MORE OF THIS KIND OF DATUM ERROR MESSAGE FOR THIS EVALUATION CYCLE.')
  endif
endif

end subroutine tao_set_invalid

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function tao_ele_geometry_with_misalignments (datum, ele, valid_value, why_invalid) result (value)
!
! Routine to evaluate a floor position with misalignments at a given element.
! This routine is private and not for general use.
!
! Input:
!   datum         -- tao_data_sturct: Datum info
!   ele           -- ele_struct: Lattice element to evaluate at.
!
! Output:
!   valid_value   -- logical: Was able to evalute the datum?
!   why_invalid   -- character(*): If not valid, why not.
!   value         -- real(rp): Datum value.
!-

function tao_ele_geometry_with_misalignments (datum, ele, valid_value, why_invalid) result (value)

type (tao_data_struct) datum
type (ele_struct) ele
type (floor_position_struct) position

real(rp) value
logical valid_value
character(*) why_invalid

!

valid_value = .false.

position = ele_geometry_with_misalignments(ele)

select case (datum%data_type)
case ('floor_actual.x');       value = position%r(1)
case ('floor_actual.y');       value = position%r(2)
case ('floor_actual.z');       value = position%r(3)
case ('floor_actual.theta');   value = position%theta
case ('floor_actual.phi');     value = position%phi
case ('floor_actual.psi');     value = position%psi
case default
  call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" IS NOT VALID', why_invalid, .true.)
  value = 0
  return
end select

valid_value = .true.

end function tao_ele_geometry_with_misalignments

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Function tao_eval_floor_orbit (datum, ele, orbit, bunch_params, valid_value, why_invalid) result (value)
!
! Routine to evaluate a floor_orbit datum at a given element.
! This routine is private and not for general use.
!
! Input:
!   datum         -- tao_data_sturct: Datum info
!   ele           -- ele_struct: Lattice element to evaluate at.
!   orbit         -- coord_struct: Particle orbit at element.
!   bunch_params  -- bunch_params_struct: Bunch parameters at element.
!
! Output:
!   valid_value   -- logical: Was able to evalute the datum?
!   why_invalid   -- character(*): If not valid, why not.
!   value         -- real(rp): Datum value.
!-

function tao_eval_floor_orbit (datum, ele, orbit, bunch_params, valid_value, why_invalid) result (value)

type (tao_data_struct) datum
type (ele_struct) ele
type (coord_struct) orbit
type (bunch_params_struct) bunch_params
type (floor_position_struct) position

real(rp) value, vec(6), p(3)
logical valid_value
character(*) why_invalid

!

valid_value = .false.

if (datum%data_source == 'lat') then
  vec = orbit%vec
else
  vec = bunch_params%centroid%vec
endif

position = orbit_to_local_curvilinear(orbit, ele, relative_to = downstream_end$)
position = coords_local_curvilinear_to_floor (position, ele, .false., relative_to = downstream_end$)

!

select case (datum%data_type)
case ('floor_orbit.x');   value = position%r(1)
case ('floor_orbit.y');   value = position%r(2)
case ('floor_orbit.z');   value = position%r(3)
case ('floor_orbit.theta');   value = position%theta
case ('floor_orbit.phi');     value = position%phi
case ('floor_orbit.psi');     value = position%psi
case default
  call tao_set_invalid (datum, 'DATA_TYPE = "' // trim(datum%data_type) // '" IS NOT VALID', why_invalid, .true.)
  value = 0
  return
end select

valid_value = .true.

end function tao_eval_floor_orbit

end module tao_data_and_eval_mod
