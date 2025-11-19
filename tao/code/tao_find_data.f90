!+
! Subroutine tao_find_data (err, data_name, d2_array, d1_array, d_array, re_array, 
!                       log_array, str_array, int_array, ix_uni, dflt_index, print_err, component)
!
! Routine to set data pointers to the correct data structures. 
!
! The re_array will be used if the component is one of:
!   model, base, design, meas, ref, old, fit, weight
! The l_array will be used if the component is one of:
!   exists, good_meas, good_ref, good_user, good_opt, good_plot
!
! Also see: tao_find_datum_using_ele_name
! 
! Example:
!   data_name = '*@orbit.x'
! re_array & l_array will have size = 0 since there is no data component specified.
!
! Example:
!   data_name = 'orbit'
! In this case the default universe will be used. The d1_array will have two components
! pointing to orbit.x and orbit.y.
! re_array & l_array will have size = 0 since there is no data component specified.
!
! Example:
!   data_name = '2@orbit.x[3,7:9]|meas'
! The measured values for the 3rd, 7th, 8th and 9th elements of orbit.x in universe #2.
! r_arrray will be allocated and l_array will have size = 0.
!
! Example:
!   data_name = 'orbit.x'
!   dflt_index = '4'
! This is equivalent to:
!   data_name = 'orbit.x[4]'
! Notice that if dflt_index is not present, or is negative, 'orbit.x' will evaluate
! to an array of numbers.
!
! Input:
!   data_name    -- Character(*): The data name type. Eg: "3@orbit.x[2:5,10]|meas"
!   ix_uni       -- Integer, optional: Index of default universe to use.
!                     If ix_uni = 0 then "viewed" universe will be used.
!                     Also, if not present then the "viewed" universe will be used.
!   dflt_index   -- character, optional: If present and non-negative, and if no index is specified
!                     by the data_name argument, this index is used in the evaluation.
!   print_err    -- Logical, optional: Print error message if data is 
!                     not found? Default is True.
!
! Output:
!   err          -- Logical: Err condition
!   d2_array(:)  -- Tao_d2_data_array_struct, allocatable, optional: Array of pointers to all 
!                     the matching d2_data structure. Size(d2_array) = 0 if no structures found.
!   d1_array(:)  -- Tao_d1_data_array_struct, allocatable, optional: Array of pointers to all 
!                     the matching d1_data structures. Size(d1_array) = 0 if no structures found.
!   d_array(:)   -- Tao_data_array_struct, allocatable, optional: Array of pointers to all 
!                     the matching tao_data_structs.  Size(d_array) = 0 if no structures found.
!   re_array(:)  -- Tao_real_pointer_struct, allocatable, optional: Array of pointers to real 
!                     component values.  Size(re_array) = 0 if no structures found.
!   log_array(:) -- Tao_logical_array_struct, allocatable, optional: Array of pointers to
!                     logical component values.  Size(log_array) = 0 if no structures found.
!   str_array(:) -- Tao_string_array_struct, allocatable, optional: Array of pointers to 
!                     character component values.  Size(str_array) = 0 if no structures found.
!   int_array(:) -- Tao_integer_array_struct, allocatable, optional: Array of pointers to
!                     integer component values.  Size(int_array) = 0 if no structures found.
!   component    -- Character(*), optional: Name of the component. E.G: 'good_user'
!                     set to ' ' if no component present.
!-

subroutine tao_find_data (err, data_name, d2_array, d1_array, d_array, re_array, &
                           log_array, str_array, int_array, ix_uni, dflt_index, print_err, component)

use tao_interface, except_dummy => tao_find_data

implicit none

type (tao_d2_data_array_struct), allocatable, optional :: d2_array(:)
type (tao_d1_data_array_struct), allocatable, optional :: d1_array(:)
type (tao_data_array_struct), allocatable, optional    :: d_array(:)
type (tao_real_pointer_struct), allocatable, optional    :: re_array(:)
type (tao_integer_array_struct), allocatable, optional :: int_array(:)
type (tao_logical_array_struct), allocatable, optional :: log_array(:)
type (tao_string_array_struct), allocatable, optional  :: str_array(:)
type (tao_universe_struct), pointer :: u

character(*) :: data_name
character(*), optional :: component
character(*), optional :: dflt_index

character(20) :: r_name = 'tao_find_data'
character(80) dat_name, component_name
character(16), parameter :: real_components(20) = [character(16) :: &
             'model', 'base', 'design', 'meas', 'ref', 'old', &
             'model_value', 'base_value', 'design_value', 'meas_value', 'ref_value', 'old_value', &
             'weight', 'invalid', 'invalid_value', 's_offset', 'ref_s_offset', 'delta_merit', 'merit', 'error_rms']
character(16), parameter :: logic_components(9) = [ &
             'exists    ', 'good_meas ', 'good_ref  ', 'good_user ', 'good_opt  ', &
             'good_plot ', 'good_base ', 'useit_opt ', 'useit_plot']
character(16), parameter :: integer_components(6) = [character(16):: &
             'ix_ele', 'ix_ele_start', 'ix_ele_ref', 'ix_d1', 'ix_uni', 'eval_point']
character(16), parameter :: string_components(6) = [character(16) :: &
              'merit_type', 'ele_name', 'ele_start_name', 'ele_ref_name', 'data_type', 'data_source']

integer, optional :: ix_uni
integer :: data_num, ios, n_found
integer i, ix, iu

logical err, component_here, this_err, print_error, error, found_data, data_exists, explicit_uni
logical, optional :: print_err

! Init

print_error = logic_option(.true., print_err)

if (present(d2_array)) then
  if (allocated (d2_array)) then
    if (size(d2_array) /= 0) deallocate (d2_array)
  endif
  if (.not. allocated(d2_array)) allocate (d2_array(0))
endif

if (present(d1_array)) then
  if (allocated (d1_array)) then
    if (size(d1_array) /= 0) deallocate (d1_array)
  endif
  if (.not. allocated(d1_array)) allocate (d1_array(0))
endif

if (present(d_array)) then
  if (allocated (d_array)) then
    if (size(d_array) /= 0) deallocate (d_array)
  endif
  if (.not. allocated(d_array)) allocate (d_array(0))
endif

if (present(re_array)) then
  if (allocated (re_array)) then
    if (size(re_array) /= 0) deallocate (re_array)
  endif
  if (.not. allocated(re_array)) allocate (re_array(0))
endif

if (present(int_array)) then
  if (allocated (int_array)) then
    if (size(int_array) /= 0) deallocate (int_array)
  endif
  if (.not. allocated(int_array)) allocate (int_array(0))
endif

if (present(log_array)) then
  if (allocated (log_array)) then
    if (size(log_array) /= 0) deallocate (log_array)
  endif
  if (.not. allocated(log_array)) allocate (log_array(0))
endif

if (present(str_array)) then
  if (allocated (str_array)) then
    if (size(str_array) /= 0) deallocate (str_array)
  endif
  if (.not. allocated(str_array)) allocate (str_array(0))
endif

err = .true.
found_data = .false.

if (data_name == '') then
  if (print_error) call out_io (s_error$, r_name, 'DATA NAME IS BLANK')
  return
endif

! Select meas, ref, etc.

ix = index(data_name, '|')
if (ix == 0) then  ! not present
  component_here = .false.
  component_name = ''
  dat_name = data_name
else
  component_here = .true.
  component_name = data_name(ix+1:)
  dat_name = data_name(:ix-1)
endif
if (present(component)) component = component_name

call string_trim (dat_name, dat_name, ix)

if (component_here) then
  call string_trim (component_name, component_name, ix)
  if (.not. any(component_name == real_components) .and. &
      .not. any(component_name == logic_components) .and. &
      .not. any(component_name == integer_components) .and. &
      .not. any(component_name == string_components)) then
    if (print_error) call out_io (s_error$, r_name, "BAD COMPONENT NAME: " // data_name)
    return            
  endif
endif

! Select universe

call tao_pick_universe (dat_name, dat_name, scratch%picked, this_err, explicit_uni = explicit_uni)
if (this_err) return

! Trim 'data::' suffix if present

if (dat_name(1:5) == 'dat::') then
  call out_io (s_error$, r_name, 'NAME USES OLD "dat::" SYNTAX. PLEASE CHANGE TO "data::": ' // dat_name)
  dat_name = data_name(6:)
endif

if (dat_name(1:6) == 'data::') dat_name = dat_name(7:)

! Find the d2 data.

data_exists = .false.  ! Does *any* data exist in the universes searched?

if (present(ix_uni) .and. .not. explicit_uni) then
  u => tao_pointer_to_universe (ix_uni)
  if (.not. associated(u)) return
  call find_this_d2 (u, dat_name, this_err)
  if (u%n_d2_data_used > 0) data_exists = .true.
else
  do i = lbound(s%u, 1), ubound(s%u, 1)
    if (.not. scratch%picked(i)) cycle
    u => tao_pointer_to_universe (i)
    call find_this_d2 (u, dat_name, this_err)
    if (this_err) return
    if (u%n_d2_data_used > 0) data_exists = .true.
  enddo
endif

! error check

if (this_err .or. .not. found_data) then
  if (print_error) then
    if (.not. data_exists) then
      call out_io (s_error$, r_name, 'Cannot find: ' // dat_name, &
                                     'No data defined in universes searched. [Check your init file?!]')
    else
      call out_io (s_error$, r_name, 'Could not find data: ' // data_name)
    endif
  endif
  return
endif

err = .false.

!----------------------------------------------------------------------------
contains

subroutine find_this_d2 (uu, name, this_err)

type (tao_universe_struct), target :: uu
integer i, ix
character(*) name
character(80) d1_name, d2_name
logical this_err

! Everything before a period is the d2 name.
! if no period then must be something like name = "orbit" and everything is the d2 name.

this_err = .false.

ix = index(name, '.')
if (ix == 0) then
  ix = index(name, '[')
  if (ix /= 0) then
    d2_name = name(1:ix-1)
    d1_name = name(ix:)
  else
    d2_name = name
    d1_name = '*'
    if (present(dflt_index)) d1_name = dflt_index
  endif
else
  d2_name = name(1:ix-1)
  d1_name = name(ix+1:)
endif

! loop over matching d2 names

do i = 1, uu%n_d2_data_used
  if (.not. match_wild(uu%d2_data(i)%name, d2_name)) cycle
  call find_this_d1 (uu%d2_data(i), d1_name, .false., this_err)
  if (this_err) return
enddo

end subroutine find_this_d2

!----------------------------------------------------------------------------
! contains

subroutine find_this_d1 (d2, name, found_d1, this_err)

type (tao_d2_data_struct), target :: d2
type (tao_d2_data_array_struct), allocatable :: d2_temp(:)
type (tao_real_pointer_struct), allocatable :: ra(:)

integer i, j, ix, nd

character(*) name
character(80) d1_name, d_name

logical found_d1, this_err

! d2_array

if (present(d2_array)) then
  if (allocated(d2_array)) then
    nd = size(d2_array)
    call move_alloc(d2_array, d2_temp)
    allocate (d2_array(nd+1))
    d2_array(1:nd) = d2_temp
    deallocate(d2_temp)
    d2_array(nd+1)%d2 => d2
  else
    allocate (d2_array(1))
    d2_array(1)%d2 => d2
  endif
endif

! Special case: d2 level components

!! if (component_name == 'scale') then
!!   if (name /= '*') then
!!     call out_io (s_error$, r_name, 'Malformed datum: ' // data_name)
!!     this_err = .true.
!!     return
!!   endif
!! 
!!   if (.not. present(re_array)) return
!!   if (allocated(re_array)) then
!!     nd = size(re_array)
!!     call move_alloc(re_array, ra)
!!     allocate (re_array(nd+1))
!!     j = nd
!!     re_array(1:nd) = ra
!!     deallocate(ra)
!!   else
!!     allocate (re_array(1))
!!     j = 0
!!   endif
!!   re_array(nd+1)%r => d2%scale
!!   found_data = .true.
!!   return
!! endif

! Everything before a '[' is the d1 name.

ix = index(name, '[')

if (ix == 0) then
  d1_name = name
  d_name = '*'
  if (present(dflt_index)) d_name = dflt_index
else
  d1_name = name(1:ix-1)
  d_name = name(ix+1:)
  ix = index(d_name, ']')
  if (ix == 0) then
    if (print_error) call out_io (s_error$, r_name, "NO MATCHING ']': " // data_name)
    this_err = .true.
    return
  endif
  if (d_name(ix+1:) /= ' ') then
    if (print_error) call out_io (s_error$, r_name, "GARBAGE AFTER ']': " // data_name)
    this_err = .true.
    return
  endif
  d_name = d_name(:ix-1)
endif

if (size(d2%d1) == 1 .and. d1_name == '') d1_name = '*'

do i = 1, size(d2%d1)
  if (.not. match_wild(d2%d1(i)%name, d1_name)) cycle
  call find_this_data (d2%d1(i), d_name, this_err)
  if (this_err) return
enddo

end subroutine find_this_d1

!----------------------------------------------------------------------------
! contains

subroutine find_this_data (d1, name, this_err)

type (tao_d1_data_struct), target :: d1
type (tao_d1_data_array_struct), allocatable :: d1_temp(:)
type (tao_data_array_struct), allocatable :: da(:)
type (tao_real_pointer_struct), allocatable :: ra(:)
type (tao_integer_array_struct), allocatable :: ia(:)
type (tao_logical_array_struct), allocatable :: la(:)
type (tao_string_array_struct), allocatable  :: sa(:)

integer i, j, nd, nl, i1, i2, num

character(*) name
character(80) d1_name, d_name

logical this_err
logical, allocatable :: list(:)

! d1_array

if (present(d1_array)) then
  if (allocated(d1_array)) then
    nd = size(d1_array)
    call move_alloc(d1_array, d1_temp)
    allocate (d1_array(nd+1))
    d1_array(1:nd) = d1_temp
    deallocate(d1_temp)
    d1_array(nd+1)%d1 => d1
  else
    allocate (d1_array(1))
    d1_array(1)%d1 => d1
  endif
endif

! Special case: d1 level components

if (name == '*' .and. component == 'name') then
  if (.not. present(str_array)) return
  if (allocated(str_array)) then
    nd = size(str_array)
    call move_alloc(str_array, sa)
    allocate (str_array(nd+1))
    j = nd
    str_array(1:nd) = sa
    deallocate(sa)
  else
    allocate (str_array(1))
    j = 0
  endif
  str_array(nd+1)%s => d1%name
  found_data = .true.
  return
endif

!

if (allocated(list)) deallocate(list)
i1 = lbound(d1%d, 1)
i2 = ubound(d1%d, 1)
allocate (list(i1:i2))

if (name == '*') then
  list = .true.

else
  call location_decode (name, list, i1, num, can_abbreviate = .false., print_err = print_err)
  if (num <  0) then
    if (logic_option(.true., print_err)) call out_io (s_error$, r_name, "BAD DATA INDEX NUMBER(S): " // name)
    this_err = .true.
    return  
  endif
endif

err = .false.
nl = count(list)
if (nl > 0) found_data = .true.

! data array

if (present(d_array)) then

  if (allocated(d_array)) then
    nd = size(d_array)
    call move_alloc(d_array, da)
    allocate (d_array(nl+nd))
    j = nd
    d_array(1:nd) = da
    deallocate(da)
  else
    allocate (d_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      d_array(j)%d => d1%d(i)
    endif
  enddo

endif

! Integer component array

if (present(int_array) .and.  any(component_name == integer_components)) then

  if (allocated(int_array)) then
    nd = size(int_array)
    call move_alloc(int_array, ia)
    allocate (int_array(nl+nd))
    j = nd
    int_array(1:nd) = ia
    deallocate(ia)
  else
    allocate (int_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      select case (component_name)
      case ('ix_bunch')
        int_array(j)%i => d1%d(i)%ix_bunch
      case ('ix_branch')
        int_array(j)%i => d1%d(i)%ix_branch
      case ('ix_ele')
        int_array(j)%i => d1%d(i)%ix_ele
      case ('ix_ele_start')
        int_array(j)%i => d1%d(i)%ix_ele_start
      case ('ix_ele_ref')
        int_array(j)%i => d1%d(i)%ix_ele_ref
      case ('ix_d1')
        int_array(j)%i => d1%d(i)%ix_d1
      case ('ix_uni', 'ix_universe')
        int_array(j)%i => d1%d(i)%d1%d2%ix_universe
      case ('eval_point')
        int_array(j)%i => d1%d(i)%eval_point
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: INTEGER DATA")
        call err_exit
      end select
    endif
  enddo

endif

! real component array

if (present(re_array) .and.  any(component_name == real_components)) then

  if (allocated(re_array)) then
    nd = size(re_array)
    call move_alloc(re_array, ra)
    allocate (re_array(nl+nd))
    j = nd
    re_array(1:nd) = ra
    deallocate(ra)
  else
    allocate (re_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      re_array(j)%good_user  => d1%d(i)%good_user
      re_array(j)%good_value => forever_true$

      select case (component_name)
      case ('model', 'model_value')
        re_array(j)%r => d1%d(i)%model_value
        re_array(j)%good_value => d1%d(i)%good_model
      case ('base', 'base_value')
        re_array(j)%r => d1%d(i)%base_value
        re_array(j)%good_value => d1%d(i)%good_base
      case ('design', 'design_value')
        re_array(j)%r => d1%d(i)%design_value
        re_array(j)%good_value => d1%d(i)%good_model
      case ('meas', 'meas_value')
        re_array(j)%r => d1%d(i)%meas_value
        re_array(j)%good_value => d1%d(i)%good_meas
      case ('ref', 'ref_value')
        re_array(j)%r => d1%d(i)%ref_value
        re_array(j)%good_value => d1%d(i)%good_ref
      case ('old', 'old_value')
        re_array(j)%r => d1%d(i)%old_value
      case ('error_rms')
        re_array(j)%r => d1%d(i)%error_rms
      case ('invalid', 'invalid_value')
        re_array(j)%r => d1%d(i)%invalid_value
      case ('weight')
        re_array(j)%r => d1%d(i)%weight
      case ('merit')
        re_array(j)%r => d1%d(i)%merit
      case ('delta_merit')
        re_array(j)%r => d1%d(i)%delta_merit
      case ('s_offset')
        re_array(j)%r => d1%d(i)%s_offset
      case ('ref_s_offset')
        re_array(j)%r => d1%d(i)%ref_s_offset
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: REAL DATA")
        call err_exit
      end select
    endif
  enddo

endif

! logical component array

if (present(log_array) .and. any(component_name == logic_components)) then

  if (allocated(log_array) .and. component_here) then
    nd = size(log_array)
    allocate (la(nd))
    la = log_array
    deallocate(log_array)
    allocate (log_array(nl+nd))
    j = nd
    log_array(1:nd) = la
    deallocate(la)
  else
    allocate (log_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      select case (component_name)
      case ('exists')
        log_array(j)%l => d1%d(i)%exists
      case ('good_base')
        log_array(j)%l => d1%d(i)%good_base
      case ('good_meas')
        log_array(j)%l => d1%d(i)%good_meas
      case ('good_ref')
        log_array(j)%l => d1%d(i)%good_ref
      case ('good_user')
        log_array(j)%l => d1%d(i)%good_user
      case ('good_opt')
        log_array(j)%l => d1%d(i)%good_opt
      case ('good_plot')
        log_array(j)%l => d1%d(i)%good_plot
      case ('useit_opt')
        log_array(j)%l => d1%d(i)%useit_opt
      case ('useit_plot')
        log_array(j)%l => d1%d(i)%useit_plot
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: LOGIC DATA")
        call err_exit
      end select
    endif
  enddo

endif

! string component array

if (present(str_array) .and. any(component_name == string_components)) then

  if (allocated(str_array) .and. component_here) then
    nd = size(str_array)
    allocate (sa(nd))
    sa = str_array
    deallocate(str_array)
    allocate (str_array(nl+nd))
    j = nd
    str_array(1:nd) = sa
    deallocate(sa)
  else
    allocate (str_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      select case (component_name)
      case ('merit_type')
        str_array(j)%s => d1%d(i)%merit_type
      case ('ele_name')
        str_array(j)%s => d1%d(i)%ele_name
      case ('ele_start_name')
        str_array(j)%s => d1%d(i)%ele_start_name
      case ('ele_ref_name')
        str_array(j)%s => d1%d(i)%ele_ref_name
      case ('data_type')
        str_array(j)%s => d1%d(i)%data_type
      case ('data_source')
        str_array(j)%s => d1%d(i)%data_source
      case default
        call out_io (s_fatal$, r_name, "INTERNAL ERROR: STRING DATA")
        call err_exit
      end select
    endif
  enddo

endif

end subroutine find_this_data

end subroutine tao_find_data
