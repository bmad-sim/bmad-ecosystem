!+
! Subroutine tao_find_var (err, var_name, v1_array, v_array, re_array, log_array, str_array, print_err, component, dflt_var_index)
!
! Find a v1 variable type, and variable component then point to it.
!
! The re_array will be used if the component is one of:
!   model, base, design, meas, ref, old, step, weight, high_lim, low_lim 
! The log_array will be used if the component is one of:
!   exists, good_var, good_user, good_opt, good_plot
! 
! Note: Any of the output allocatable arrays will have size = 0 if not used.
! For example, if var_name is "model", then str_array, if present, have size = 0.
!
! Examples:
!   var_name = 'quad_k1[3]|design'
!   var_name = 'quad_k1[3]|slave[2]'
!
! Input:
!   var_name       -- Character(*): Name of the variable.
!   print_err      -- Logical, optional: Print error message if data is 
!                       not found? Default is True.
!   dflt_var_index -- chracter(*), optional: If present and "[...]" var selection substring is not present,
!                       then dflt_var_index will be used. [Do not include the brackets in this string.]
!
! Output:
!   err            -- Logical: err condition
!   v1_array(:)    -- Tao_v1_var_array_struct, allocatable, optional: Array of pointers to 
!                       all the v1_var structures.
!   v_array(:)     -- Tao_var_array_struct, allocatable, optional: Array of pointers to the 
!                       variable data point.
!   re_array(:)    -- Tao_real_pointer_struct, allocatable, optional: Array of pointers to 
!                       the real component values.
!   log_array(:)   -- Tao_logical_array_struct, allocatable, optional: Array of pointers to
!                       logical component values.
!   str_array(:)   -- Tao_string_array_struct, allocatable, optional: Array of pointers to 
!                       character component values.
!   component      -- Character(*), optional: Name of the component. E.G: 'good_user'
!                     set to ' ' if no component present.
!-

subroutine tao_find_var (err, var_name, v1_array, v_array, re_array, log_array, str_array, print_err, component, dflt_var_index)

use tao_interface, except_dummy => tao_find_var

implicit none

type (tao_v1_var_array_struct), allocatable, optional  :: v1_array(:)
type (tao_var_array_struct), allocatable, optional     :: v_array(:)
type (tao_real_pointer_struct), allocatable, optional  :: re_array(:)
type (tao_logical_array_struct), allocatable, optional :: log_array(:)
type (tao_string_array_struct), allocatable, optional  :: str_array(:)

integer i, ix, n_var, ios, n_v1

character(16), parameter :: real_components(19) = [character(16) :: &
             'model', 'base', 'design', 'meas', 'ref', 'old', &
             'model_value', 'base_value', 'design_value', 'meas_value', 'ref_value', 'old_value', &
             'step', 'weight', 'high_lim', 'low_lim', 'key_delta', 'merit', 'delta_merit']
character(16), parameter :: logic_components(8) = [character(16) :: &
             'exists', 'good_var', 'good_user', 'good_opt', 'good_plot', &
             'useit_opt', 'useit_plot', 'key_bound']
character(16), parameter :: string_components(3) = [character(16):: 'merit_type', 'ele_name', 'attrib_name']

character(*) :: var_name
character(*), optional :: component, dflt_var_index
character(20) :: r_name = 'tao_find_var'
character(80) v1_name, v_name, component_name

logical, optional :: print_err
logical err, component_here, this_err, print_error

! Init

print_error = logic_option(.true., print_err)

if (present(v1_array)) then
  if (allocated (v1_array)) then
    if (size(v1_array) /= 0) deallocate (v1_array)
  endif
  if (.not. allocated(v1_array)) allocate (v1_array(0))
endif

if (present(v_array)) then
  if (allocated (v_array)) then
    if (size(v_array) /= 0) deallocate (v_array)
  endif
  if (.not. allocated(v_array)) allocate (v_array(0))
endif

if (present(re_array)) then
  if (allocated (re_array)) then
    if (size(re_array) /= 0) deallocate (re_array)
  endif
  if (.not. allocated(re_array)) allocate (re_array(0))
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

! Error if no variables exist

if (s%n_var_used == 0) then
  if (print_error) call out_io (s_warn$, r_name, &
                        "NO VARIABLES HAVE BEEN DEFINED IN THE INPUT FILES!")
  return
endif

! Select meas, ref, etc.

ix = index(var_name, '|')
if (ix == 0) then  ! not present
  component_here = .false.
  component_name = ' '   ! garbage
  v1_name = var_name
else
  component_here = .true.
  component_name = var_name(ix+1:)
  v1_name = var_name(:ix-1)
endif
if (present(component)) component = component_name

call string_trim (v1_name, v1_name, ix)
call string_trim (component_name, component_name, ix)

if (component_here) then
  if (component_name(1:5) /= 'slave' .and. .not. any(component_name == real_components) .and. &
      .not. any(component_name == logic_components) .and. &
      .not. any(component_name == string_components)) then
    if (print_error) call out_io (s_error$, r_name, "BAD COMPONENT NAME: " // var_name)
    return            
  endif
endif

! Trim 'var::' suffix if present

if (v1_name(1:5) == 'var::') v1_name = v1_name(6:)

! split on '['

ix = index(v1_name, '[')
if (ix == 0) then
  v_name = '*'
  if (present(dflt_var_index)) v_name = dflt_var_index
else
  v_name  = v1_name(ix+1:)
  v1_name = v1_name(1:ix-1)
  ix = index(v_name, ']')
  if (ix == 0) then
    if (print_error) call out_io (s_error$, r_name, "NO MATCHING ']': " // var_name)
    return
  endif
  if (v_name(ix+1:) /= ' ') then
    if (print_error) call out_io (s_error$, r_name, "GARBAGE AFTER ']': " // var_name)
    return
  endif
  v_name = v_name(:ix-1)
endif

call string_trim(v1_name, v1_name, ix)
if (ix == 0) then
  if (print_error) call out_io (s_error$, r_name, 'VARIABLE NAME IS BLANK')
  return
endif

! Point to the correct v1 var type 

n_v1 = 0
do i = 1, s%n_v1_var_used
  if (.not. match_wild (s%v1_var(i)%name, v1_name)) cycle
  n_v1 = n_v1 + 1
enddo

if (n_v1 == 0) then
  if (print_error) call out_io (s_error$, r_name, "COULDN'T FIND V1 VARIABLES MATCHING: " // v1_name)
  return
endif

if (present(v1_array)) then
  deallocate (v1_array)
  allocate (v1_array(n_v1))
endif

n_v1 = 0
do i = 1, s%n_v1_var_used
  if (.not. match_wild (s%v1_var(i)%name, v1_name)) cycle
  n_v1 = n_v1 + 1
  if (present(v1_array)) v1_array(n_v1)%v1 => s%v1_var(i)
  call find_this_var (s%v1_var(i), v_name, this_err)
  if (this_err) return
enddo

! error check

if (err) then
  if (print_error) call out_io (s_error$, r_name, "COULDN'T FIND VARIABLE: " // var_name)
  return
endif

!----------------------------------------------------------------------------
contains

subroutine find_this_var (v1, name, this_err)

type (tao_v1_var_struct) :: v1
type (tao_var_array_struct), allocatable :: va(:)
type (tao_real_pointer_struct), allocatable :: ra(:)
type (tao_logical_array_struct), allocatable :: la(:)
type (tao_string_array_struct), allocatable  :: sa(:)

integer i, j, nd, nl, i1, i2, num, ix

character(*) name
character(40), allocatable :: names(:)
character(80) v1_name, v_name

logical this_err
logical, allocatable :: list(:)

!

if (allocated(list)) deallocate(list, names)
i1 = lbound(v1%v, 1)
i2 = ubound(v1%v, 1)
allocate (list(i1:i2), names(i1:i2))
this_err = .false.

if (name == '*') then
  list = .true.

else
  do i = i1, i2
    names(i) = v1%v(i)%ele_name
  enddo
  call location_decode (name, list, i1, num, names, can_abbreviate = .false., print_err = print_err)
  if (num <  0) then
    if (logic_option(.true., print_err)) call out_io (s_error$, r_name, "BAD VAR INDEX NUMBER(S): " // name)
    this_err = .true.
    return  
  endif
endif

err = .false.
nl = count(list)

! Record variable pointed to

if (present(v_array)) then

  if (allocated(v_array)) then
    nd = size(v_array)
    allocate (va(nd))
    va = v_array
    deallocate(v_array)
    allocate (v_array(nl+nd))
    j = nd
    v_array(1:nd) = va
    deallocate(va)
  else
    allocate (v_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (list(i)) then
      j = j + 1
      v_array(j)%v => v1%v(i)
    endif
  enddo

endif

! real component

if (present(re_array) .and.  any(component_name == real_components)) then

  if (allocated(re_array)) then
    nd = size(re_array)
    allocate (ra(nd))
    ra = re_array
    deallocate(re_array)
    allocate (re_array(nl+nd))
    j = nd
    re_array(1:nd) = ra
    deallocate(ra)
  else
    allocate (re_array(nl))
    j = 0
  endif

  do i = i1, i2
    if (.not. list(i)) cycle
    j = j + 1
    select case (component_name)
    case ('correction_value');        re_array(j)%r => v1%v(i)%correction_value
    case ('model', 'model_value');    re_array(j)%r => v1%v(i)%model_value
    case ('base', 'base_value');      re_array(j)%r => v1%v(i)%base_value
    case ('design', 'design_value');  re_array(j)%r => v1%v(i)%design_value
    case ('dmerit_dvar');             re_array(j)%r => v1%v(i)%dmerit_dvar
    case ('meas', 'meas_value');      re_array(j)%r => v1%v(i)%meas_value
    case ('ref', 'ref_value');        re_array(j)%r => v1%v(i)%ref_value
    case ('old', 'old_value');        re_array(j)%r => v1%v(i)%old_value
    case ('s');                       re_array(j)%r => v1%v(i)%s
    case ('step');                    re_array(j)%r => v1%v(i)%step
    case ('weight');                  re_array(j)%r => v1%v(i)%weight
    case ('high_lim');                re_array(j)%r => v1%v(i)%high_lim
    case ('low_lim');                 re_array(j)%r => v1%v(i)%low_lim
    case ('key_delta');               re_array(j)%r => v1%v(i)%key_delta
    case ('merit');                   re_array(j)%r => v1%v(i)%merit
    case ('delta_merit');             re_array(j)%r => v1%v(i)%delta_merit
    case default
      call out_io (s_fatal$, r_name, "INTERNAL ERROR: REAL VAR")
      this_err = .true.
      call err_exit
    end select
  enddo

endif

! logical component

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
    if (.not. list(i)) cycle
    j = j + 1
    select case (component_name)
    case ('exists');      log_array(j)%l => v1%v(i)%exists
    case ('good_var');    log_array(j)%l => v1%v(i)%good_var
    case ('good_user');   log_array(j)%l => v1%v(i)%good_user
    case ('good_opt');    log_array(j)%l => v1%v(i)%good_opt
    case ('good_plot');   log_array(j)%l => v1%v(i)%good_plot
    case ('useit_opt');   log_array(j)%l => v1%v(i)%useit_opt
    case ('useit_plot');  log_array(j)%l => v1%v(i)%useit_plot
    case ('key_bound');   log_array(j)%l => v1%v(i)%key_bound
    case default
      call out_io (s_fatal$, r_name, "INTERNAL ERROR: LOGIC VAR")
      this_err = .true.
      call err_exit
    end select
  enddo

endif

! string component

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
    if (.not. list(i)) cycle
    j = j + 1
    select case (component_name)
    case ('merit_type');    str_array(j)%s => v1%v(i)%merit_type
    case ('ele_name');      str_array(j)%s => v1%v(i)%merit_type
    case ('attrib_name');   str_array(j)%s => v1%v(i)%merit_type
    case default
      call out_io (s_fatal$, r_name, "INTERNAL ERROR: STRING VAR")
      this_err = .true.
      call err_exit
    end select
  enddo

endif

end subroutine find_this_var

end subroutine tao_find_var
