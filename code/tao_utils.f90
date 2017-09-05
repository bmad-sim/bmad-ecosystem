!+
! Module tao_utils
!
! helper subroutines available for communal usage.
!-

module tao_utils

use tao_pointer_to_universe_mod
use bmad

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_re_allocate_expression_info (info, n, exact, init_val)
!
! Routine to reallocate an array of tao_expression_info_struct structs.
! This is modeled after the reallocate functions in Numerical Recipes.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   info(:) -- tao_expression_info_struct, allocatable:
!   n       -- Integer: Size wanted.
!   exact   -- Logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   info(:) -- Real(rp), Allocatable: Allocated array with size(re) >= n.
!-

subroutine tao_re_allocate_expression_info (info, n, exact)

implicit none

type (tao_expression_info_struct), allocatable :: info(:), temp_info(:)

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact

!

if (allocated(info)) then
  n_old = size(info)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  call move_alloc (info, temp_info)
  allocate (info(n))
  n_save = min(n, n_old)
  info(1:n_save) = temp_info(1:n_save)
  deallocate (temp_info)  

else
  allocate (info(n))
endif

end subroutine tao_re_allocate_expression_info


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

subroutine tao_key_info_to_str (ix_key, ix_min_key, ix_max_key, key_str, header_str)

implicit none

type (tao_var_struct), pointer :: var

integer ix_key, ix_min_key, ix_max_key
integer i, n, m, p, ix_var
integer :: j_var1, j_att

real(rp) :: y_here, norm, v, x1, x2, y1, y2

character(*) key_str
character(*) header_str
character(60) fmt, fmt2
character(15) model_str, val0_str, delta_str
character(4) exp_str
character(24) :: r_name = 'tao_key_info_to_str'

! Compute widths of var1 and attrib fields.

j_var1 = 4
j_att = 5

do i = ix_min_key, ix_max_key
  if (i > ubound(s%key, 1)) cycle
  ix_var = s%key(i)
  if (ix_var < 1) cycle
  j_var1 = max(j_var1, len_trim(tao_var1_name(s%var(ix_var))))
  j_att  = max(j_att,  len_trim(tao_var_attrib_name(s%var(ix_var))))
enddo

! Write header 

write (fmt, '(a, i5, a, i2, a)') '(a, ', j_var1-2, 'x, a, ', j_att, 'x, a)'
write (header_str, fmt) 'Name', 'Attrib', 'Value     Value0      Delta Opt'

! Write key info

key_str = ''
if (ix_key > ubound(s%key, 1)) return
  
ix_var = s%key(ix_key)
if (ix_var < 1) return

var => s%var(ix_var)
v = max(abs(var%model_value), abs(var%key_val0), abs(var%key_delta))
if (v == 0) then
  n = 0
  m = 2
else
  m = 1.001 * log10(v)
  n = 3 * floor(m/3.0)
  p = 3 - (m - n)
endif

if (m >= -1 .and. m <= 1) then
  fmt2 = '(f11.4, a)'
  n = 0
elseif (m == 2) then
  fmt2 = '(f11.2, a)'
  n = 0
else
  write (fmt2, '(a, i1, a)') '(f7.', p, ', a)'
  write (exp_str, '(a, i3.2)') 'E', n
  if (exp_str(2:2) == ' ') exp_str(2:2) = '+'
endif

model_str = ''; val0_str = ''; delta_str = ''

write (model_str, fmt2) var%model_value / 10.0**n, exp_str
write (val0_str,  fmt2) var%key_val0 / 10.0**n, exp_str
write (delta_str, fmt2) var%key_delta / 10.0**n, exp_str

write (fmt, '(3(a, i2.2))') '(a', j_var1, ', 2x, a', j_att, ', 3a15, 3x, l1)'
write (key_str, fmt) tao_var1_name(var), tao_var_attrib_name(var), model_str, &
                        val0_str, delta_str, var%useit_opt

end subroutine tao_key_info_to_str

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_data_check (err)
!
! Routine to do some checking of data.
!-

subroutine tao_data_check (err)

implicit none

type (tao_data_struct), pointer :: datum
integer iu, id
logical err
character(16) :: r_name = 'tao_data_check'

!

err = .false.

do iu = lbound(s%u, 1), ubound(s%u, 1)
  do id = 1, size(s%u(iu)%data)
    datum => s%u(iu)%data(id)
    if (datum%merit_type(1:4) == 'int_' .and. &
            (s%global%opt_with_ref .or. s%global%opt_with_base)) then
      call out_io (s_error$, r_name, &
                        'BAD DATUM INTEGRATION FOR: ' // tao_datum_name(datum))
      err = .true.
    endif
  enddo
enddo

end subroutine tao_data_check

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_locate_elements (ele_list, ix_universe, eles, err, lat_type, ignore_blank, print_err, ix_dflt_branch) 
!
! Subroutine to find the lattice elements in the lattice
! corresponding to the ele_list argument. 
!
! Note: If ele_list can contain elements from different universes, use
! the routine:
!   tao_locate_all_elements 
!
! Input:
!   ele_list       -- Character(*): String with element names using element list format.
!   ix_universe    -- Integer: Universe to search. 0 => search s%com%default_universe.
!   lat_type       -- Integer, optional: model$ (default), design$, or base$.
!   ignore_blank   -- Logical, optional: If present and true then do nothing if
!                     ele_list is blank. otherwise treated as an error.
!   print_err      -- Logical, optional: If present and False then do not print error messages.
!   ix_dflt_branch -- Integer, optional: If present and positive then use this as the branch index 
!                       for elements specified using an integer index (EG: "43").
!                       If not present or -1 the default branch is branch 0.
!
! Output:
!   eles  -- ele_pointer_struct(:), allocatable: Array of elements in the model lat. 
!   err   -- Logical: Set true on error.
!-

subroutine tao_locate_elements (ele_list, ix_universe, eles, err, lat_type, ignore_blank, &
                                                                 print_err, above_ubound_is_err, ix_dflt_branch)

implicit none

type (tao_universe_struct), pointer :: u
type (tao_lattice_struct), pointer :: tao_lat
type (ele_pointer_struct), allocatable :: eles(:)

integer, optional :: lat_type, ix_dflt_branch
integer ios, ix, ix_universe, num, i, i_ix_ele, n_loc

character(*) ele_list
character(200) ele_name
character(20) :: r_name = 'tao_locate_elements'

logical err, printit
logical, optional :: ignore_blank, print_err, above_ubound_is_err

! 

err = .true.
printit = logic_option (.true., print_err)

call re_allocate_eles (eles, 0, exact = .true.)

call str_upcase (ele_name, ele_list)
call string_trim (ele_name, ele_name, ix)

if (ix == 0 .and. logic_option(.false., ignore_blank)) return

if (ix == 0) then
  if (printit) call out_io (s_error$, r_name, 'ELEMENT NAME IS BLANK')
  return
endif

u => tao_pointer_to_universe (ix_universe)
if (.not. associated(u)) return

tao_lat => tao_pointer_to_tao_lat (u, lat_type)

call lat_ele_locator (ele_name, tao_lat%lat, eles, n_loc, err, above_ubound_is_err, ix_dflt_branch)
if (err) return

if (n_loc == 0) then
  if (printit) call out_io (s_error$, r_name, 'ELEMENT(S) NOT FOUND: ' // ele_list)
  err = .true.
  return
endif

call re_allocate_eles (eles, n_loc, .true., .true.)

end subroutine tao_locate_elements

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_pointer_to_tao_lat (u, lat_type) result (tao_lat)
!
! Routine to set a pointer to a tao_lat.
!
! Also see:
!   tao_pointer_to_universe
!
! Input:
!   u         -- Tao_universe_struct: Universe to work with
!   lat_type  -- Integer, optional: model$ (default), design$, or base$.
!
! Output:
!   tao_lat   -- tao_lattice_struct, pointer: Tao_lat pointer. 
!                   Points to u%model, u%design, or u%base
!-

function tao_pointer_to_tao_lat (u, lat_type) result (tao_lat)

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct), pointer :: tao_lat
integer, optional :: lat_type
character(28) :: r_name = 'tao_pointer_to_tao_lat'

!

if (.not. present(lat_type)) then
  tao_lat => u%model
  return
endif

select case (lat_type)
case (design$)
  tao_lat => u%design
case (model$)
  tao_lat => u%model
case (base$)
  tao_lat => u%base
case default
  call err_exit
end select

end function tao_pointer_to_tao_lat

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_var_useit_plot_calc (graph, var)
!
! Subroutine to set the variables for plotting.
!
! Input:
!
! Output:
!   var     -- Tao_var_struct:
!     %useit_plot -- True if good for plotting.
!-

subroutine tao_var_useit_plot_calc (graph, var)

implicit none

type (tao_graph_struct) graph
type (tao_var_struct) var(:)

!

var%useit_plot = var%exists .and. var%good_plot .and. var%good_var .and. &
                 (var%good_user .or. .not. graph%draw_only_good_user_data_or_vars)

end subroutine tao_var_useit_plot_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_orbit_value (component, orbit, value, err)
!
! Routine to return the orbit component indicated by component
!
! Input:
!   component -- character(*): 'orbit.x', 'orbit.px', 'intensity.x', 'phase.y', 'energy', 'pc', etc.
!   orbit     -- coord_struct: Particle orbit.
!
! Output:
!   value -- real(rp): orbit component.
!   err   -- logical: Set True if component is not recognized. False otherwise.
!-

subroutine tao_orbit_value (component, orbit, value, err)

implicit none

type (coord_struct) orbit

real(rp) value
character(*) component
logical err

!

err = .true.

if (component == 'state') then
  value = orbit%state
  err = .false.
  return
endif

if (orbit%state /= alive$) return

select case (component)
case ('orbit_x', 'orbit.x')
  value = orbit%vec(1)
case ('orbit_px', 'orbit.px')
  value = orbit%vec(2)
case ('orbit_y', 'orbit.y')
  value = orbit%vec(3)
case ('orbit_py', 'orbit.py')
  value = orbit%vec(4)
case ('orbit_z', 'orbit.z')
  value = orbit%vec(5)
case ('orbit_pz', 'orbit.pz')
  value = orbit%vec(6)
case ('spin.x', 'spin_x')
  value = orbit%spin(1)
case ('spin.y', 'spin_y')
  value = orbit%spin(2)
case ('spin.z', 'spin_z')
  value = orbit%spin(3)
case ('intensity')
  value = orbit%field(1)**2 + orbit%field(2)**2
case ('intensity_x', 'intensity.x')
  value = orbit%field(1)**2
case ('intensity_y', 'intensity.y')
  value = orbit%field(2)**2
case ('phase_x', 'phase.x')
  value = orbit%phase(1)
case ('phase_y', 'phase.y')
  value = orbit%phase(2)
case ('t', 'time')
  value = orbit%t
case ('beta')
  value = orbit%beta
case ('energy')
  if (orbit%species == photon$) then
    value = orbit%p0c
  else
    call convert_pc_to(orbit%p0c * (1 + orbit%vec(6)), orbit%species, e_tot = value)
  endif
case ('pc')
  if (orbit%species == photon$) then
    value = orbit%p0c
  else
    value = orbit%p0c * (1 + orbit%vec(6))
  endif
case default
  return
end select

err = .false.

end subroutine tao_orbit_value

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_var_target_calc ()
! 
! Subroutine to calculate the variable target values (the values that they need
! to be set to to do a correction of the orbit, phase, etc.
!
! Input:
!
! Output:
!-

subroutine tao_var_target_calc ()

implicit none

type (tao_var_struct), pointer :: var

integer i, j

!

do j = 1, s%n_var_used
  var => s%var(j)
  var%correction_value = var%meas_value + (var%design_value - var%model_value)
enddo

end subroutine tao_var_target_calc

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_set_var_model_value (var, value, print_limit_warning)
!
! Subroutine to set the value for a model variable and do the necessary bookkeeping.
! If value is past the variable's limit, and s%global%var_limits_on = True, the
! variable will be set to the limit value.
!
! Input:
!   var   -- Tao_var_struct: Variable to set
!   value -- Real(rp): Value to set to
!   print_limit_warning
!         -- Logical, optional: Print a warning if the value is past the variable's limits.
!             Default is True.
!-

subroutine tao_set_var_model_value (var, value, print_limit_warning)

implicit none

type (tao_var_struct), target :: var
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (tao_var_slave_struct), pointer :: var_slave
type (ele_struct), pointer :: ele

real(rp) value
integer i
logical, optional :: print_limit_warning

!

if (.not. var%exists) return

! check if hit variable limit
if (s%global%var_limits_on .and. (.not. s%global%only_limit_opt_vars .or. var%useit_opt)) then
  if (value < var%low_lim) then
    if (logic_option (.true., print_limit_warning)) &
          call out_io (s_blank$, ' ', "Hit lower limit of variable: " // tao_var1_name(var))
    value = var%low_lim
  elseif (value > var%high_lim) then
    if (logic_option (.true., print_limit_warning)) &
          call out_io (s_blank$, ' ', "Hit upper limit of variable: " // tao_var1_name(var))
    value = var%high_lim
  endif
endif

var%model_value = value
do i = 1, size(var%slave)
  var_slave => var%slave(i)
  var_slave%model_value = value
  u => s%u(var_slave%ix_uni)

  if (s%com%common_lattice .and.  var_slave%ix_uni == ix_common_uni$) then
    s%u(:)%calc%lattice = .true.
  else
    u%calc%lattice = .true.
  endif

  lat => u%model%lat
  if (var%ele_name == 'BEAM_START') then
    u%model%tao_branch(0)%orb0%vec = lat%beam_start%vec
    u%model%tao_branch(0)%orb0%t   = lat%beam_start%t
    u%model%tao_branch(0)%orb0%p0c = lat%beam_start%p0c
  else
    ele => lat%branch(var_slave%ix_branch)%ele(var_slave%ix_ele)
    call set_flags_for_changed_attribute (ele, var_slave%model_value)
  endif
enddo

end subroutine tao_set_var_model_value

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! subroutine tao_count_strings (string, pattern, num)
! 
! Subroutine to count the number of a specific pattern in the string
!
! Input:    
!  string    -- character(*): the string to look at
!  pattern   -- character(*): the search pattern
!
! Output:
!  num       -- integer: number of occurances
!-

subroutine tao_count_strings (string, pattern, num)

implicit none

character(*) string, pattern

integer num, len_string, len_pattern, i

num = 0
len_pattern = len(pattern)
len_string  = len(string)

do i = 1, len(string)
  if (i+len_pattern-1 .gt. len_string) return
  if (string(i:i+len_pattern-1) .eq. pattern) num = num + 1
enddo

end subroutine tao_count_strings

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_lat_bookkeeper (u, tao_lat)
!
! This will make sure all bookkeeping is up to date.
!
! Input:
!  u            -- tao_universe_struct
!  lat_name     -- Integer: Which lattice
!
! Output:
!  tao_lat      -- lat_struct
!-

subroutine tao_lat_bookkeeper (u, tao_lat)

implicit none

type (tao_universe_struct), target :: u
type (tao_lattice_struct) :: tao_lat

integer i, j

character(20) :: r_name = "tao_lat_bookkeeper"

! Setup from common if it exists

if (associated(u%common)) then

  ! First put in the common values

  do i = 1, s%n_var_used
    do j = 1, size(s%var(i)%slave)
      s%var(i)%slave(j)%model_value = s%var(i)%common_slave%model_value
      s%var(i)%slave(j)%base_value  = s%var(i)%common_slave%base_value
    enddo
  enddo

  ! Then put in the values for this universe

  do i = 1, s%n_var_used
    do j = 1, size(s%var(i)%slave)
      if (s%var(i)%slave(j)%ix_uni /= u%ix_uni) cycle
      s%var(i)%slave(j)%model_value = s%var(i)%model_value
      s%var(i)%slave(j)%base_value = s%var(i)%base_value
    enddo
  enddo

endif

! Do bookkeeping

call lattice_bookkeeper (tao_lat%lat)

end subroutine tao_lat_bookkeeper

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function tao_read_this_index (name, ixc) result (ix)
!
! Returns the integer value in the array <name> at position <ixc>. This is used
! for finding a 6-dimensional index reference so any value at <ixc> greater than
! 6 returns an error.
!
! Input:
!  name     -- Character(*): character array holding the index
!  ixc      -- Integer: location within <name> to evaluate
!
! Output:
!  ix       -- Integer: Index at <name>(<ixc>:<ixc>)
!
! Example:
!      name = r:26
!      ixc  = 3
!
! Gives:
!      ix = 2
!
! Example:
!      name = mat_94
!      ixc  = 7
! Gives:
!      Error: "BAD INDEX CONSTRAINT: mat_94"
!-

function tao_read_this_index (name, ixc) result (ix)

character(*) name
integer ix, ixc
character(20) :: r_name = 'tao_read_this_index'

ix = index('123456', name(ixc:ixc))
if (ix == 0) then
  call out_io (s_abort$, r_name, 'BAD INDEX CONSTRAINT: ' // name)
  call err_exit
endif

end function tao_read_this_index


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_var1_name(var) result (var1_name)
!
! Function to return the variable name in the form:
!   var1_name[index]
! For example:
!   quad_k1[23]
!
! Input:
!   var -- Tao_var_struct: Variable
!
! Output:
!   var1_name -- Character(60): Appropriate name.
!-

function tao_var1_name(var) result (var1_name)

implicit none

type (tao_var_struct) var
character(60) var1_name

!

write (var1_name, '(2a, i0, a)') trim(var%v1%name), '[', var%ix_v1, ']'

end function tao_var1_name

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_var_attrib_name(var) result (var_attrib_name)
!
! Function to return a string encoding the attributes controlled by a 
! variable in the form:
!   {universe@}element[attribute]
! For example:
!   Q03W[k1]
!   [2,3]@V11E[HKICK]
!
! Input:
!   var -- Tao_var_struct: Variable
!
! Output:
!   var_attrib_name -- Character(60): Attribute list.
!-

function tao_var_attrib_name(var) result (var_attrib_name)

use location_encode_mod

implicit none

type (tao_var_struct) var

character(60) var_attrib_name
integer i

!

if (size(s%u) > 1 .and. size(var%slave) > 0) then
  call location_encode (var_attrib_name, var%slave%ix_uni, ',')
  if (index(var_attrib_name, ',') /= 0 .or. index(var_attrib_name, ':') /= 0) then  
    var_attrib_name = '[' // trim(var_attrib_name) // ']' 
  endif
  var_attrib_name = trim(var_attrib_name) // '@' // &
                    trim(var%ele_name) // '[' // trim(var%attrib_name) // ']'
else
  var_attrib_name = trim(var%ele_name) // '[' // trim(var%attrib_name) // ']'
endif

end function tao_var_attrib_name

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_constraint_type_name (datum) result (datum_name)
!
! Function to return the constraint type in the form:
!   data_type <merit_type>
! For example:
!   eta.x <target>
!
! Input:
!   datum      -- Tao_data_struct: Datum
!
! Output:
!   datum_name -- Character(60): Appropriate name.
!-

function tao_constraint_type_name(datum) result (datum_name)

implicit none

type (tao_data_struct) datum
character(200) datum_name
integer ix

! Expressions are too long so shorten the name

datum_name = trim(datum%data_type) // ' <' // trim(datum%merit_type) // '>'
if (datum_name(1:11) == 'expression:') call string_trim (datum_name(12:), datum_name, ix)

end function tao_constraint_type_name

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_datum_name (datum, show_universe) result (datum_name)
!
! Function to return the datum name in the form:
!   d2_name.d1_name[index]
! or (if show_universe is True and there is more than one universe):
!   universe@d2_name.d1_name[index]
! For example:
!   2@orbit.x[23]
!   
!
! Input:
!   datum         -- Tao_data_struct: Datum
!   show_universe -- Logical, optional: Show the datum's universe.
!                       Default is True.
!
! Output:
!   datum_name -- Character(60): Appropriate name.
!-

function tao_datum_name(datum, show_universe) result (datum_name)

implicit none

type (tao_data_struct) datum
character(60) datum_name
logical, optional :: show_universe

! If this datum is "isolated". That is, it does not have an associated d1_data 
! structure then just use it's data_type.
! This can happen if the datum is derived from a curve.

if (.not. associated(datum%d1)) then
  datum_name = datum%data_type
  return
endif

! Normal case

datum_name = tao_d2_d1_name (datum%d1, show_universe)
write (datum_name, '(2a, i0, a)') trim(datum_name), '[', datum%ix_d1, ']'

end function tao_datum_name

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_d2_d1_name (d1, show_universe) result (d2_d1_name)
!
! Function to return the datum name in the form:
!   d2_name.d1_name
! If there is only one d1_data array associated with the d2_data 
! array then the name is shortened to:
!   d2_name
! Additionally, if show_universe is True and there is more than one universe
! then "universe@" is prepended to the name.
! For example:
!   2@orbit.x
!   
!
! Input:
!   d1            -- Tao_d1_data_struct: Data array.
!   show_universe -- Logical, optional: Show the datum's universe.
!                       Default is True.
!
! Output:
!   d2_d1_name -- Character(60): Appropriate name.
!-

function tao_d2_d1_name(d1, show_universe) result (d2_d1_name)

implicit none

type (tao_d1_data_struct) d1
character(60) d2_d1_name, temp_str
logical, optional :: show_universe

! If there is only one d1 array associated with the d2_data array then
! drop the d1 name.

if (size(d1%d2%d1) == 1) then
  write (d2_d1_name, '(a)') trim(d1%d2%name)
else
  write (d2_d1_name, '(3a)') trim(d1%d2%name), '.', trim(d1%name)
endif

!

if (size(s%u) > 1 .and. logic_option(.true., show_universe)) then
  ! Stupid gfortran compiler requires a temp string for this
  temp_str = d2_d1_name
  write (d2_d1_name, '(i0, 2a)') d1%d2%ix_uni, '@', trim(temp_str)
endif

end function tao_d2_d1_name

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Function tao_curve_name (curve, use_region) result (curve_name)
!
! Function to return the curve name in the form:
!   plot_name.graph_name.curve_name
! For example:
!   orbit.x.c1
!
! Input:
!   curve      -- Tao_curve_struct: Curve
!   use_region -- Logical: If present and True then use the region 
!                  name instead of the plot name. Region name is
!                  'NULL_REGION' if there is no assocaited region.
!
! Output:
!   curve_name -- Character(60): Appropriate name.
!-

function tao_curve_name(curve, use_region) result (curve_name)

implicit none

type (tao_curve_struct) curve
character(60) curve_name
logical, optional :: use_region

!

curve_name = '.' // trim(curve%g%name) // '.' // trim(curve%name)

if (logic_option(.false., use_region)) then
  if (associated(curve%g%p%r)) then
    curve_name = trim(curve%g%p%r%name) // trim(curve_name)
  else
    curve_name = 'NULL_REGION' // trim(curve_name)
  endif
else
    curve_name = trim(curve%g%p%name) // trim(curve_name)
endif

end function tao_curve_name 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_parse_command_args (error, cmd_words)
!
! Subroutine to parse the command line arguments.
!
! Input:
!   cmd_words(:) -- Character(*), optional: If present then this is used
!                    in place of the command line.
! Output:
!   error -- Logical: Set True if there is an error. False otherwise.
!-

subroutine tao_parse_command_args (error, cmd_words)

implicit none

character(*), optional :: cmd_words(:)
character(80) arg0, base, switch
character(24) :: r_name = 'tao_parse_command_args'

integer n_arg, i_arg, ix
logical error

! Get command line input

error = .false.

call tao_hook_parse_command_args()
if (.not. s%com%parse_cmd_args) return

if (present(cmd_words)) then
  n_arg = size(cmd_words)
  if (cmd_words(1) == '') return
else
  n_arg = cesr_iargc()
  if (n_arg == 0) return
endif

! loop over all arguments

i_arg = 0

do 

  if (i_arg == n_arg) exit
  call get_next_arg (arg0)

  call match_word (arg0,            ['-?                       ', &
        '-init                    ', '-noinit                  ', '-beam_all                ', '-beam0                   ', &
        '-noplot                  ', '-lat                     ', '-log_startup             ', '-beam                    ', &
        '-var                     ', '-data                    ', '-building_wall           ', '-plot                    ', &
        '-startup                 ', 'help                     ', '-help                    ', '?                        ', &
        '-geometry                ', '-rf_on                   ', '-debug                   ', '-disable_smooth_line_calc', &
        '-color_prompt            ', '-no_stopping             ', '-hook_init_file          ', '-gui_mode                '], &
              ix, .true., matched_name=switch)

  select case (switch)

  case ('-init')
    call get_next_arg (s%com%init_tao_file)
    s%com%init_tao_file_arg_set = .true.
    if (s%com%init_tao_file == '') then
      call out_io (s_fatal$, r_name, 'NO TAO INIT FILE NAME ON COMMAND LINE.')
      call err_exit
    endif
    ix = SplitFileName(s%com%init_tao_file, s%com%init_tao_file_path, base)

  case ('-beam')
    call get_next_arg (s%com%beam_file)

  case ('-beam_all')
    call get_next_arg (s%com%beam_all_file)

  case ('-beam0')
    call get_next_arg (s%com%beam0_file)

  case ('-building_wall')
    call get_next_arg (s%com%building_wall_file)

  case ('-color_prompt')
    s%global%prompt_color = 'BLUE'

  case ('-data')
    call get_next_arg (s%com%data_file)

  case ('-disable_smooth_line_calc')
    s%global%disable_smooth_line_calc = .true.

  case ('-debug')
    s%global%debug_on = .true.

  case ('-geometry')
    call get_next_arg (s%com%plot_geometry)

  case ('help', '-help', '?', '-?')
    call tao_print_command_line_info
    stop

  case ('-hook_init_file')
    call get_next_arg (s%com%hook_init_file)

  case ('-lat')
    call get_next_arg (s%com%lat_file)

  case ('-log_startup')
    s%com%log_startup = .true.

  case ('-no_stopping')
    s%global%stop_on_error = .false.

  case ('-noinit')
    s%com%init_tao_file = ''

  case ('-noplot')
    s%com%noplot_arg_set = .true.

  case ('-rf_on')
    s%global%rf_on = .true.

  case ('-plot')
    call get_next_arg (s%com%plot_file)

  case ('-gui_mode')
    call out_io (s_error$, r_name, 'NOTE: -gui_mode NO LONGER DOES ANYTHING!')

  case ('-startup')
    call get_next_arg (s%com%startup_file)

  case ('-var')
    call get_next_arg (s%com%var_file)

  case default
    call out_io (s_error$, r_name, 'BAD COMMAND LINE ARGUMENT: ' // arg0)
    call tao_print_command_line_info
    error = .true.
    return
  end select

enddo

!-----------------------------
contains

subroutine get_next_arg(arg)

character(*) arg

!

if (i_arg == n_arg) then
  call out_io (s_error$, r_name, 'MISSING COMMAND LINE ARGUMENT FOR: ' // arg0)
  error = .true.
  return
endif

i_arg = i_arg + 1

if (present(cmd_words)) then
  arg = cmd_words(i_arg)
else
  call cesr_getarg(i_arg, arg)
endif

end subroutine get_next_arg

end subroutine tao_parse_command_args

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_print_command_line_info
!
! Routine to print a list of the command line options.
!-

subroutine tao_print_command_line_info

implicit none

character(40), parameter :: r_name = 'tao_print_command_line_info'

!

call out_io (s_blank$, r_name, [ &
        'Syntax:                                                                                  ', &
        '  <path-to-tao-exe-directory>/tao {OPTIONS}                                              ', &
        'Options are:                                                                             ', &
        '  -beam <beam_file>               # Beam init particle positions                         ', &
        '  -beam0 <beam0_file>             # Beam init params (beam size, etc.)                   ', &
        '  -beam_all <all_beam_file>       # Beam info from previous tracking                     ', &
        '  -building_wall <wall_file>      # Define the building tunnel wall                      ', &
        '  -color_prompt                   # Set color of prompt string to blue                   ', &
        '  -data <data_file>               # Define data for plotting and optimization            ', &
        '  -debug                          # Debug mode for Wizards                               ', &
        '  -disable_smooth_line_calc       # Disable the smooth line calc used in plotting        ', &
        '  -geometry <width>x<height>      # Plot window geometry                                 ', &
        '  -help                           # Display this list of command line options            ', &
        '  -hook_init_file <init_file>     # Init file for hook routines (Default = tao_hook.init)', &
        '  -init <tao_init_file>           # Tao init file                                        ', &
        '  -lat <bmad_lattice_file>        # Bmad lattice file                                    ', &
        '  -lat xsif::<xsif_lattice_file>  # XSIF lattice file                                    ', &
        '  -log_startup                    # Write startup debugging info                         ', &
        '  -no_stopping                    # For debugging: Prevents Tao from exiting on err      ', &
        '  -noinit                         # Do not use Tao init file                             ', &
        '  -noplot                         # Do not open a plotting window                        ', &
        '  -plot <plot_file>               # Define plot setup info                               ', &
        '  -rf_on                          # Keep RF on (Default is to turn off)                  ', &
        '  -startup <starup_command_file>  # Commands to run after parsing Tao init file          ', &
        '  -var <var_file>                 # Define variables for plotting and optimization       '])

end subroutine tao_print_command_line_info

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_ele_to_ele_track (ix_universe, ix_branch, ix_ele, ix_ele_track)
!
! Subroutine to compute ix_ele_track:
!   = ix_ele                 if ix_ele <= lat%branch(ix)branch)%n_ele_track
!   = ix_slave_at_exit_end   if ix_ele is a super_lord  
!   = -1                     otherwise
!
! Input:
!   ix_universe -- Integer: Universe index.
!   ix_branch   -- Integer: Branch index.
!   ix_ele      -- Integer: Element index
!
! Output:
!   ix_ele_track -- Integer: Corresponding element in the tracking 
!                         part of the lattice.
!-

subroutine tao_ele_to_ele_track (ix_universe, ix_branch, ix_ele, ix_ele_track)

implicit none

type (lat_struct), pointer :: lat
type (ele_struct), pointer :: slave

integer ix_universe, ix_branch, ix_ele, ix_ele_track
integer i_uni, ix_c

!

i_uni = tao_universe_number(ix_universe)
lat => s%u(i_uni)%model%lat

if (ix_ele < 0) then
  ix_ele_track = -1

elseif (ix_ele <= lat%branch(ix_branch)%n_ele_track) then
  ix_ele_track = ix_ele

elseif (lat%ele(ix_ele)%lord_status == super_lord$) then
  slave => pointer_to_slave (lat%ele(ix_ele), lat%ele(ix_ele)%n_slave)
  ix_ele_track = slave%ix_ele ! element at exit end.

else  ! overlays, multipass_lords, etc.
  ix_ele_track = -1
endif

end subroutine tao_ele_to_ele_track

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine tao_string_to_element_id (str, ix_class, ele_name, err, print_err)
!
! Routine to split a string in the form str = "xxx::yyy" into an element class
! and an element name. Example: 
!   str = "quad::q*".
! gives
!   ix_class = quadrupole$
!   ele_name = "Q*"
!
! If "::" is not found then ix_class is set to 0 (all classes).
! If str is of the form: "*::yyy" then ix_class is set to 0 (all classes).
! Class abbreviations and lower case names accepted 
! ele_name will be converted to upper case
!
! Input:
!   str       -- Character(*): Character string to parse.
!   print_err -- Logical, optional: If True then print an error message if 
!                 there is a problem. Default is True.
!
! Output:
!   ix_class  -- Integer: Element class. 0 => all classes. -1 => Not an element [EG: str = "var::..."].
!   ele_name  -- Character(*): Element name.
!   err       -- Set true if there is a problem translating the element class.
!-

subroutine tao_string_to_element_id (str, ix_class, ele_name, err, print_err)

implicit none

integer ix, ix_class

character(*) str, ele_name
character(40) :: r_name = 'tao_string_to_element_id'
character(20) class

logical, optional :: print_err
logical err

! 

if (index(str, 'dat::') /= 0) then
  call out_io (s_error$, r_name, 'NAME USES OLD "dat::" SYNTAX. PLEASE CHANGE TO "data::": ' // str)
  call err_exit
endif

!

err = .false.
ix_class = -1

if (str(1:6) == 'data::') return
if (str(1:5) == 'var::') return
if (str(1:5) == 'lat::') return
if (str(1:6) == 'wall::') return

ix = index(str, '::')

if (ix == 0) then
  ix_class = 0
  ele_name = str
  call str_upcase (ele_name, ele_name)
  return
endif

class = str(:ix-1)
ele_name = str(ix+2:)
call str_upcase (ele_name, ele_name)

if (class == '*') then
  ix_class = 0
  return
endif

ix_class = key_name_to_key_index (class, .true.)
if (ix_class < 1) then
  if (logic_option (.true., print_err)) call out_io (s_error$, r_name, 'BAD CLASS NAME: ' // class)
  err = .true.
endif

end subroutine tao_string_to_element_id

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine floor_to_screen_coords (graph, floor, screen)
!
! Routine to project a 3D floor coordinate onto a 2D projection plane
! along with projecting the orientation vector.
!
! Input:
!   graph  -- tao_graph_struct: Graph defining the projection plane.
!   floor  -- floor_position_struct: 3D coordinate.
!
! Output:
!   screen  -- floor_position_struct: Projected point
!     %r(3)   -- projected (x, y) = (%r(1), %r(2))
!     %theta  -- angle in (x, y) plane of projection of the orientation vector.
!-

subroutine floor_to_screen_coords (graph, floor, screen)

implicit none

type (tao_graph_struct) graph
type (floor_position_struct) floor, screen
real(rp) orient(3), theta, phi, x, y

! Get projection position

call floor_to_screen (graph, floor%r, screen%r(1), screen%r(2))

! screen%theta does not depend upon floor%psi

theta = floor%theta
phi = floor%phi
orient = [sin(theta) * cos(phi), sin(phi), cos(theta) * cos(phi)]  ! orientation vector

select case (graph%floor_plan_view(1:1))
case ('x')
  x = orient(1)
case ('y')
  x = orient(2)
case ('z')
  x = orient(3)
end select

select case (graph%floor_plan_view(2:2))
case ('x')
  y = orient(1)
case ('y')
  y = orient(2)
case ('z')
  y = orient(3)
end select

screen%theta = atan2(y, x) + twopi * graph%floor_plan_rotation

end subroutine floor_to_screen_coords

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine floor_to_screen (graph, floor, x_screen, y_screen)
!
! Routine to project a 3D floor coordinate onto a 2D projection plane.
!
! Input:
!   graph  -- tao_graph_struct: Graph defining the projection plane.
!   floor  -- floor_position_struct: 3D coordinate.
!
! Output:
!   x_screen -- real(rp): x-coordinate of projected point.
!   y_screen -- real(rp): y-coordinate of projected point.
!-

subroutine floor_to_screen (graph, r_floor, x_screen, y_screen)

implicit none

type (tao_graph_struct) graph

real(rp) r_floor(3), x_screen, y_screen
real(rp) x, y
real(rp), save :: t, old_t = 0
real(rp), save :: cc, ss

! 

select case (graph%floor_plan_view(1:1))
case ('x')
  x = r_floor(1)
case ('y')
  x = r_floor(2)
case ('z')
  x = r_floor(3)
end select

select case (graph%floor_plan_view(2:2))
case ('x')
  y = r_floor(1)
case ('y')
  y = r_floor(2)
case ('z')
  y = r_floor(3)
end select

t = graph%floor_plan_rotation
if (t == 0) then
  x_screen = x
  y_screen = y
else
  if (t /= old_t) then
    cc = cos(twopi * t)
    ss = sin(twopi * t)
    old_t = t
  endif
  x_screen =  x * cc - y * ss
  y_screen =  x * ss + y * cc 
endif

end subroutine floor_to_screen

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_split_component (comp_str, comp, err)
!
! Routine to split a component string.
!
! Input:
!   comp_str -- Character(*): Components. EG: 'meas - design'
!
! Output:
!   comp(:) -- Tao_data_var_component_struct, allocatable: Array of individual components.
!   err     -- Logical: Set True if there is an error, False otherwise.
!-

subroutine tao_split_component (comp_str, comp, err)

type (tao_data_var_component_struct), allocatable :: comp(:)

integer i, n, ix, ix1, ix2

character(*) comp_str
character(60) str
character(40) :: r_name = 'tao_split_component'

logical err

! Count number of components.

err = .true.
call string_trim (comp_str, str, ix)
if (.not. (str(1:1) == '+' .or. str(1:1) == '-')) str = '+' // trim(str)

n = 0
do i = 1, len_trim(str)
  if (str(i:i) == '+' .or. str(i:i) == '-') n = n + 1
enddo

! Allocate space and transfer info.

if (allocated(comp)) then
  if (size(comp) /= n) deallocate (comp)
endif
if (.not. allocated(comp)) allocate (comp(n))

n = 0
do n = 1, size(comp)
  if (str(1:1) == '+') then
    comp(n)%sign = 1
  elseif (str(1:1) == '-') then
    comp(n)%sign = -1
  else
    call out_io (s_error$, r_name, 'BAD COMPONENT LIST: ' // comp_str)
    return
  endif

  call string_trim(str(2:), str, ix)

  if (n == size(comp)) then
    comp(n)%name = str
    exit
  endif

  ix1 = index(str, "+")
  ix2 = index(str, "-")
  if (ix1 /= 0 .and. ix2 /= 0) then
    ix = min(ix1, ix2) - 1
  else
    ix = max(ix1, ix2) - 1
  endif

  comp(n)%name = str(1:ix)
  call string_trim(str(ix+1:), str, ix)

enddo

err = .false.

end subroutine tao_split_component

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_turn_on_special_calcs_if_needed_for_plotting ()
! 
! Routine to set u%dynch_rad_int_clac, etc to True if needed for a plot.
!-

subroutine tao_turn_on_special_calcs_if_needed_for_plotting ()

implicit none

type (tao_universe_struct), pointer :: u
type (tao_graph_struct), pointer :: graph
type (tao_curve_struct), pointer :: curve

integer i, j, k

character(*), parameter :: r_name = 'tao_turn_on_special_calcs_if_needed_for_plotting'

! Go through all the plots and find which universes need special_calcs.
! u%picked_uni = True => need calc.

s%u(:)%picked_uni  = s%u(:)%calc%rad_int_for_plotting

s%u(:)%calc%rad_int_for_plotting    = .false.
s%u(:)%calc%chrom_for_plotting      = .false.
s%u(:)%calc%beam_sigma_for_plotting = .false.

do i = 1, size(s%plot_page%region)
  if (.not. s%plot_page%region(i)%visible) cycle
  if (.not. allocated(s%plot_page%region(i)%plot%graph)) cycle

  do j = 1, size(s%plot_page%region(i)%plot%graph)
    graph => s%plot_page%region(i)%plot%graph(j)
    if (.not. allocated(graph%curve)) cycle

    do k = 1, size(graph%curve)
      curve => graph%curve(k)
      u => tao_pointer_to_universe(curve%ix_universe)

      if (.not. u%picked_uni .and. tao_rad_int_calc_needed(curve%data_type, curve%data_source)) then
        if (curve%ix_branch /= 0) then
          call out_io (s_fatal$, r_name, 'PLOTTING THIS: ' // curve%data_type, 'ON A BRANCH NOT YET IMPLEMENTED!')
          call err_exit
        endif
        u%calc%rad_int_for_plotting = .true.
        u%calc%lattice = .true.
        u%picked_uni = .true.
      endif

      if (tao_chrom_calc_needed(curve%data_type, curve%data_source)) then
        u%calc%chrom_for_plotting = .true.
        u%calc%lattice = .true.
      endif

      if (tao_beam_sigma_calc_needed(curve%data_type, curve%data_source)) then
        u%calc%beam_sigma_for_plotting = .true.
        u%calc%lattice = .true.
      endif

    enddo     ! graph%curve
  enddo     ! plot%graph
enddo     ! plotting%region

end subroutine tao_turn_on_special_calcs_if_needed_for_plotting

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_rad_int_calc_needed (data_type, data_source) result (do_rad_int)
! 
! Routine decide if a datum or plot curve needs the radiation integrals 
! to be evaluated.
!-

function tao_rad_int_calc_needed (data_type, data_source) result (do_rad_int)

implicit none

character(*) data_type, data_source
logical do_rad_int

!

do_rad_int = .false.

if (data_source /= 'lat') return

if (data_type  == 'sigma.pz') do_rad_int = .true. 
if (data_type(1:5)  == 'emit.') do_rad_int = .true. 
if (data_type(1:10) == 'norm_emit.') do_rad_int = .true. 
if (data_type(1:7)  == 'rad_int') do_rad_int = .true.
if (data_type(1:16)  == 'apparent_rad_int') do_rad_int = .true.

end function tao_rad_int_calc_needed

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_chrom_calc_needed (data_type, data_source) result (do_chrom)
! 
! Routine decide if a datum or plot curve needs the chromaticity calculation.
!-

function tao_chrom_calc_needed (data_type, data_source) result (do_chrom)

implicit none

character(*) data_type, data_source
logical do_chrom

!

do_chrom = .false.

if (data_source /= 'lat') return
if (data_type(1:6)  == 'chrom.') do_chrom = .true. 

end function tao_chrom_calc_needed

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_beam_sigma_calc_needed (data_type, data_source) result (do_beam_sigma)
! 
! Routine decide if a datum or plot curve needs a beam sigma calculation.
!-

function tao_beam_sigma_calc_needed (data_type, data_source) result (do_beam_sigma)

implicit none

character(*) data_type, data_source
logical do_beam_sigma

!

do_beam_sigma = .false.

if (data_source /= 'lat') return
if (data_type(1:6)  == 'sigma.') do_beam_sigma = .true. 

end function tao_beam_sigma_calc_needed

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_subin_uni_number (name_in, ix_uni, name_out) result (ok)
!
! Routine to mangle a file name based upon the universe number ix_uni.
! If there is a "#" in the name and the number of universes is more than one,
! the universe number will be substituted.
! If the number of universes is 1, the "#" will be stripped from the name.
!
! Input:
!   name_in -- Character(*): Input name with "#" character
!   ix_uni  -- Integer: Universe index.
!
! Output:
!   name_out -- Character(*): Output name.
!   ok       -- Logical: False if multiple universes and no "#" in name_in.
!                 True otherwise.
!-

function tao_subin_uni_number (name_in, ix_uni, name_out) result (ok)

character(*) name_in, name_out
character(28) :: r_name = 'tao_subin_uni_number'
integer ix, ix_uni
logical ok

!

ok = .true.
name_out = name_in

ix = index(name_out, '#')
if (size(s%u) > 1 .and. ix == 0) then
  call out_io (s_info$, r_name, 'FILE NAME DOES NOT HAVE A "#" CHARACTER!', &
    ' YOU NEED THIS TO GENERATE A UNIQUE FILE NAME FOR EACH UNIVERSE!')
  ok = .false.
endif

if (ix /= 0) write (name_out, '(a, i0, a)') &
                        name_out(1:ix-1), ix_uni, trim(name_out(ix+1:))

end function tao_subin_uni_number

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_lat_emit_calc (plane, emit_type, ele, modes) result (emit)
!
! Routine to calculate the emittance.
! The "projected" emittance is:
!   emit = sqrt(sigma_xx * sigma_pp - sigma_xp^2)
! Where the sigmas are calculated including coupling effects but assuming that sigma_pz = 0.
! The "apparent" emittance is:
!   emit = sqrt(sigma_xx - sigma_xp^2 / sigma_pp) / beta
!
! Input:
!   plane     -- Integer: x_plane$ or y_plane$.
!   emit_type -- Integer: Either projected_emit$ or apparent_emit$
!   ele       -- ele_struct: Element holding the Twiss and coupling parameters.
!   modes     -- normal_modes_struct: Structure holding the emittances
!
! Output:
!   emit -- Real(rp): emittance.
!-

function tao_lat_emit_calc (plane, emit_type, ele, modes) result (emit)

implicit none

type (ele_struct) ele
type (normal_modes_struct) modes

real(rp) s_mat(4,4), v_mat(4,4), emit
real(rp), save :: a_mat(4,4) = 0
integer plane, emit_type

!

a_mat(1,1) =  modes%a%emittance * ele%a%beta
a_mat(1,2) = -modes%a%emittance * ele%a%alpha
a_mat(2,2) =  modes%a%emittance * ele%a%gamma
a_mat(2,1) = a_mat(1,2)

a_mat(3,3) =  modes%b%emittance * ele%b%beta
a_mat(3,4) = -modes%b%emittance * ele%b%alpha
a_mat(4,4) =  modes%b%emittance * ele%b%gamma
a_mat(4,3) = a_mat(3,4)

call make_v_mats (ele, v_mat)
s_mat = matmul(matmul(v_mat, a_mat), transpose(v_mat))

if (plane == x_plane$) then
  if (emit_type == projected_emit$) then
    emit = sqrt(s_mat(1,1) * s_mat(2,2) - s_mat(1,2)**2)
  elseif (emit_type == apparent_emit$) then
    emit = s_mat(1,1) / ele%a%beta
  else
    call err_exit    
  endif

elseif (plane == y_plane$) then
  if (emit_type == projected_emit$) then
    emit = sqrt(s_mat(3,3) * s_mat(4,4) - s_mat(3,4)**2)
  elseif (emit_type == apparent_emit$) then
    emit = s_mat(3,3) / ele%b%beta
  else
    call err_exit    
  endif

else
  call err_exit
endif

end function tao_lat_emit_calc

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_beam_emit_calc (plane, emit_type, ele, bunch_params) result (emit)
!
! Routine to calculate the emittance from beam parameters.
! The "projected" emittance is:
!   emit = sqrt(sigma_xx * sigma_pp - sigma_xp^2)
! Where the sigmas are calculated including coupling effects but assuming that sigma_pz = 0.
! The "apparent" emittance is:
!   emit = sqrt(sigma_xx - sigma_xp^2 / sigma_pp) / beta
!
! Input:
!   plane        -- Integer: x_plane$ or y_plane$.
!   emit_type    -- Integer: Either projected_emit$ or apparent_emit$
!   ele          -- ele_struct: Element.
!   bunch_params -- bunch_params_struct: Bunch sigma matrix
!
! Output:
!   emit -- Real(rp): emittance.
!-

function tao_beam_emit_calc (plane, emit_type, ele, bunch_params) result (emit)

implicit none

type (ele_struct) ele
type (bunch_params_struct) bunch_params
real(rp) emit
integer plane, emit_type

!

if (plane == x_plane$) then
  if (emit_type == projected_emit$) then
    emit = bunch_params%x%emit
  elseif (emit_type == apparent_emit$) then
    if (bunch_params%sigma(6,6) == 0) then
      emit = bunch_params%sigma(1,1) / ele%a%beta
    else
      emit = (bunch_params%sigma(1,1) - bunch_params%sigma(1,6)**2 / bunch_params%sigma(6,6)) / ele%a%beta
    endif
  else
    call err_exit    
  endif

elseif (plane == y_plane$) then
  if (emit_type == projected_emit$) then
    emit = bunch_params%y%emit
  elseif (emit_type == apparent_emit$) then
    if (bunch_params%sigma(6,6) == 0) then
      emit = bunch_params%sigma(3,3) / ele%b%beta
    else
      emit = (bunch_params%sigma(3,3) - bunch_params%sigma(3,6)**2 / bunch_params%sigma(6,6)) / ele%b%beta
    endif
  else
    call err_exit    
  endif

else
  call err_exit
endif

end function tao_beam_emit_calc

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function tao_pointer_to_ele_shape (ix_uni, ele, ele_shape, dat_var_name, dat_var_value, ix_shape) result (e_shape)
!
! Routine to return the shape associated with a lattice element
!
! Input:
!   ix_uni        -- integer: Universe index.
!   ele           -- ele_struct: Lattice element.
!   ele_shape(:)  -- tao_ele_shape_struct: Array of shapes to search.
!
! Output:
!   e_shape       -- tao_ele_shape_struct, pointer: Associated shape. 
!                       Nullified if there is no associated shape.
!   dat_var_name  -- character(*), optional: Name of datum or variable associated with e_shape. 
!                       Will be set to "" if there is no associated datum or variable.
!   dat_var_value -- real(rp), optional: Value of datum or variable associated with e_shape.
!                       Will be set to zero if there is no associated datum or variable.
!   ix_shape      -- integer, optional: Index of associated shape in ele_shape(:).
!-

function tao_pointer_to_ele_shape (ix_uni, ele, ele_shape, dat_var_name, dat_var_value, ix_shape) result (e_shape)

implicit none

type (ele_struct) ele

type (tao_ele_shape_struct), target :: ele_shape(:)
type (tao_ele_shape_struct), pointer :: e_shape, es
type (tao_data_array_struct), allocatable, target :: d_array(:)
type (tao_logical_array_struct), allocatable :: logic_array(:)
type (tao_data_struct), pointer :: datum
type (tao_var_array_struct), allocatable, target :: v_array(:)
type (tao_var_struct), pointer :: var
type (tao_real_pointer_struct), allocatable :: re_array(:)

real(rp), optional :: dat_var_value

integer, optional :: ix_shape
integer ix_uni
integer j, j2, k, n_ele_track

character(*), optional :: dat_var_name
character(*), parameter :: r_name = 'tao_pointer_to_ele_shape'

logical err

!

nullify(e_shape)
if (present(dat_var_name)) dat_var_name = ''
if (present(dat_var_value)) dat_var_value = 0

if (ele%lord_status == group_lord$) return
if (ele%lord_status == overlay_lord$) return
if (ele%slave_status == super_slave$) return

do k = 1, size(ele_shape)
  es => ele_shape(k)
  if (present(ix_shape)) ix_shape = k
  if (.not. es%draw) cycle

  ! Data

  if (index(es%ele_id, 'dat::') /= 0) then
    call out_io (s_error$, r_name, 'SHAPE USES OLD "dat::" SYNTAX. PLEASE CHANGE TO "data::": ' // es%ele_id)
    call err_exit
  endif

  if (es%ele_id(1:6) == 'data::') then
    call tao_find_data (err, es%ele_id, d_array = d_array, log_array = logic_array, re_array = re_array)
    if (err) cycle
    do j = 1, size(d_array)
      datum => d_array(j)%d
      if (datum%d1%d2%ix_uni /= ix_uni) cycle
      if (size(logic_array) /= 0) then
        if (.not. logic_array(j)%l) cycle
      endif
      if (ele%ix_branch /= datum%ix_branch .or. ele%ix_ele /=  datum%ix_ele) cycle
      e_shape => es
      if (present(dat_var_name)) dat_var_name = tao_datum_name(datum)
      if (present(dat_var_value)) then
        if (size(re_array) > 0) then
          dat_var_value = re_array(j)%r
        else
          dat_var_value = datum%model_value
        endif
      endif
      return
    enddo
    cycle
  endif

  ! Variables

  if (es%ele_id(1:5) == 'var::') then
    call tao_find_var (err, es%ele_id, v_array = v_array, log_array = logic_array, re_array = re_array)
    if (err) cycle

    do j = 1, size(v_array)
      var => v_array(j)%v
      if (size(logic_array) /= 0) then
        if (.not. logic_array(j)%l) cycle
      endif
      do j2 = 1, size(var%slave)
        if (var%slave(j2)%ix_uni /= ix_uni) cycle
        if (var%ele_name == 'BEAM_START') then
          if (ele%ix_branch /= 0 .or. ele%ix_ele /= 0) cycle
        else
          if (ele%ix_branch /= var%slave(j2)%ix_branch .or. ele%ix_ele /=  var%slave(j2)%ix_ele) cycle
        endif
        e_shape => es
        if (present(dat_var_name)) dat_var_name = tao_var1_name(var)
        if (present(dat_var_value)) then
          if (size(re_array) > 0) then
            dat_var_value = re_array(j)%r
          else
            dat_var_value = var%slave(j2)%model_value
          endif
        endif
        return
      enddo
    enddo
    cycle
  endif

  ! All else

  if (es%ele_id == '') cycle
  if (es%ix_ele_key == -1) cycle
  if (es%ix_ele_key /= 0 .and. es%ix_ele_key /= ele%key) cycle
  if (.not. match_wild(ele%name, es%name_ele)) cycle

  e_shape => es
  return
enddo

end function tao_pointer_to_ele_shape 

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Put the variables marked by key_bound in the key table.

subroutine tao_setup_key_table ()

implicit none

integer i, j

! Key table setup

call re_allocate (s%key, count(s%var%key_bound))
s%key = -1

j = 0
do i = 1, s%n_var_used
  if (s%var(i)%key_bound .and. s%var(i)%exists) then
    j = j + 1
    s%key(j) = i
    s%var(i)%key_val0 = s%var(i)%model_value
    s%var(i)%ix_key_table = j
  else
    s%var(i)%key_bound = .false.
    s%var(i)%ix_key_table = -1
  endif
enddo

end subroutine tao_setup_key_table

end module tao_utils
