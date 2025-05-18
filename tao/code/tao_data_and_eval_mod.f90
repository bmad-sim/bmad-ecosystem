module tao_data_and_eval_mod

use tao_interface
use bmad_interface

implicit none

!! private tao_scratch_values_calc, tao_eval_floor_orbit, tao_ele_geometry_with_misalignments

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
            this_ele => dflt_ele%branch%ele(datum%ix_ele)
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

    err_flag = .not. tao_data_sanity_check(datum, .true., '', u)
    if (err_flag) return

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

    if (.not. valid) return
  enddo

  n_tot = n_tot + n_loc
enddo

if (n_tot == 0) then
  if (print_err) call out_io (s_error$, r_name, 'ELEMENT NOT FOUND: ' // ele_name)
  err = .true.
  return
endif

err = .false.

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

end subroutine tao_get_data

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

end subroutine tao_data_coupling_init

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
character(200) name
character(len(str)) :: str2, word2
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
  stack%value(1) = tao_param_value_at_s (str, dflt_ele, dflt_ele, dflt_orbit, err_flag)
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
  stack%type = ele_num$

  if (err_flag) then
    call tao_evaluate_lat_or_beam_data (err_flag, name, stack%value, print_err, "lat", dflt_ele_ref, &
                              dflt_ele_start, dflt_ele, dflt_component, dflt_uni, dflt_eval_point, dflt_s_offset)
    stack%type = lat_num$
  endif

  call tao_re_allocate_expression_info (stack%info, size(stack%value))
  stack%info%good = (.not. err_flag)
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
  value = tao_param_value_at_s (datum%data_type, ele_at_s, ele_at_s, orb_at_s, err, bad_datum = bad_datum)
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

  value = tao_param_value_at_s (datum%data_type, ele_at_s, ele_at_s, orb_at_s, err, bad_datum = bad_datum)
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

    value = value - tao_param_value_at_s (datum%data_type, ele_at_s, ele_at_s, orb_at_s, err, bad_datum = bad_datum)
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

  case (species_const$) ! Something like "electron". Just push on stack.
    i2 = i2 + 1
    stk2(i2)%name = stack(i)%name
    call re_allocate(stk2(i2)%value, 1)

  case (species$)
    stk2(i2)%value = species_id(stk2(i2)%name)

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
