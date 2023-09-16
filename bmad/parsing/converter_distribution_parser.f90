!+
! Subroutine converter_distribution_parser(ele, delim, delim_found, err_flag)
!
! Routine to parse a distribution parameter of a converter element
! This routine is not for general use.
!
! Input:
!   ele         -- ele_struct: Converter element.
!
! Output:
!   ele         -- ele_struct: Converter element with %converter field set.
!   delim       -- character(1): Ending delimitor.
!   delim_found -- logical: Has a delimitor been found?
!   err_flag    -- logical: Set True if there is an error. False otherwise.
!-

subroutine converter_distribution_parser(ele, delim, delim_found, err_flag)

use bmad_parser_mod, dummy => converter_distribution_parser
use object_model_mod

implicit none

type subobj
  character(40) name
  integer n_min
  integer n_max
end type

type (ele_struct), target :: ele
type (converter_struct), pointer :: c
type (converter_distribution_struct), allocatable :: dtemp(:)
type (converter_distribution_struct), pointer :: dist
type (converter_sub_distribution_struct), pointer :: sub_d
type (converter_prob_pc_r_struct), pointer :: prob_pc_r
type (converter_direction_out_struct), pointer :: dir_out
type (converter_dir_coef_struct), pointer :: alpha, beta, c_x, dds
type (object_struct), target :: obj0
type (object_struct), pointer :: obj1, obj2, obj3, obj4

real(rp), allocatable :: r(:)
real(rp) val
integer i, j, k, n, ix_word, n_r, n_pc, ixe, n_sd, n_subd, ix_sd

character(*) delim
character(200) str
character(:), allocatable :: line, why_invalid

logical delim_found, valid, err_flag

!

valid = object_document_parse (obj0, 'DISTRIBUTION', line, why_invalid, get_more_text_func, parse_one = .true.)
bp_com%parse_line = trim(line) // ' ' // bp_com%parse_line
call get_next_word(str, ix_word, '}],', delim, delim_found, err_flag = err_flag)

if (.not. valid) then
  ! Why_invalid may be blank if it is set by get_more_text_func. 
  ! In this case, get_more_text_func has called parser_error so we do not have to call parser_error again.
  if (why_invalid /= '') call parser_error (why_invalid, 'FOR ELEMENT: ' // ele%name)
  err_flag = .true.
endif

if (err_flag) return

!

err_flag = .false.

c => ele%converter
if (.not. allocated(c%dist)) then
  allocate (c%dist(1))
  dist => c%dist(1)
else
  n = size(c%dist) + 1
  call move_alloc(c%dist, dtemp)
  allocate (c%dist(n))
  c%dist(1:n-1) = dtemp
  dist => c%dist(n)
endif

n_subd = subobject_number(obj0, 'SUB_DISTRIBUTION')
allocate (dist%sub_dist(n_subd))

if (.not. parser_check_subobjects(obj0, [subobj('THICKNESS',1,1), subobj('SUB_DISTRIBUTION',1,9999), &
                       subobj('MATERIAL',1,1), subobj('SPECIES_OUT',1,1)], ele)) return

if (.not. parser_read_subobj_real(obj0, 'THICKNESS', dist%thickness, ele)) return

ix_sd = 0
do n_sd = 1, n_subd
  obj1 => pointer_to_subobject(obj0, 'SUB_DISTRIBUTION', ix_sd)

  if (.not. parser_check_subobjects (obj1, [subobj('PC_IN',1,1), subobj('PROB_PC_R',1,1), &
                                 subobj('SPIN_Z_OUT',0,1), subobj('SPIN_IN',0,1), subobj('DIRECTION_OUT',1,1)], ele)) return
  sub_d => dist%sub_dist(n_sd)

  ! PC_IN
  obj2 => pointer_to_subobject(obj1, 'PC_IN')
  if (.not. parser_read_obj_real(obj2, sub_d%pc_in, ele)) return

  ! SPIN_IN
  obj2 => pointer_to_subobject(obj1, 'SPIN_IN')
  if (associated(obj2)) then
    if (.not. parser_read_obj_real_array(obj2, sub_d%spin_in, .true., ele)) return
  endif

  ! PROB_PC_R
  obj2 => pointer_to_subobject(obj1, 'PROB_PC_R')
  if (.not. parser_check_subobjects (obj2, [subobj('R_VALUES',1,1), subobj('ROW',1,9999)], ele)) return
  prob_pc_r => sub_d%prob_pc_r
  obj3 => pointer_to_subobject(obj2, 'R_VALUES')
  n_r = obj3%n_token
  n_pc = subobject_number(obj2, 'ROW')
  allocate(prob_pc_r%r(n_r))
  if (.not. parser_read_obj_real_array(obj3, sub_d%prob_pc_r%r, .true., ele)) return
  allocate (prob_pc_r%prob(n_pc, n_r), prob_pc_r%pc_out(n_pc))
  ixe = 0
  do k = 1, n_pc
    obj3 => pointer_to_subobject(obj2, 'ROW', ixe)
    if (.not. parser_check_subobjects (obj3, [subobj('PC_OUT',1,1), subobj('PROB',1,1)], ele)) return
    if (.not. parser_read_subobj_real(obj3, 'PC_OUT', prob_pc_r%pc_out(k), ele)) return
    if (.not. parser_read_subobj_real_array(obj3, 'PROB', prob_pc_r%prob(k,:), .true., ele)) return
  enddo
  if (.not. parser_check_obj_order (prob_pc_r%pc_out(:), '<', 'ROW->PC_OUT', obj2, ele)) return

  ! SPIN_Z_OUT
  obj2 => pointer_to_subobject(obj1, 'SPIN_Z_OUT')
  if (associated(obj2)) then
    if (.not. parser_check_subobjects (obj2, [subobj('R_VALUES',1,1), subobj('ROW',1,9999)], ele)) return
    obj3 => pointer_to_subobject(obj2, 'R_VALUES')
    n_r = obj3%n_token
    n_pc = subobject_number(obj2, 'ROW')
    call re_allocate(r, n_r)
    if (.not. parser_read_obj_real_array(obj3, r, .true., ele)) return
    if (n_pc /= size(prob_pc_r%pc_out) .or. n_r+1 /= size(prob_pc_r%r)) then
      call parser_error('ARRAY SIZE MISMATCH IN CONVERTER DISTRIBUTION BETWEEN PROB_PC_R AND SPIN_Z_OUT MATRICES.', &
                        'FOR ELEMENT: ' // ele%name)
    endif
    if (any (r /= prob_pc_r%r(2:))) then
      call parser_error('R VALUE MISMATCH IN CONVERTER DISTRIBUTION BETWEEN PROB_PC_R AND SPIN_Z_OUT MATRICES.', &
                        'FOR ELEMENT: ' // ele%name)
    endif
    allocate (prob_pc_r%spin_z(n_pc, n_r+1))
    ixe = 0
    do k = 1, n_pc
      obj3 => pointer_to_subobject(obj2, 'ROW', ixe)
      if (.not. parser_check_subobjects (obj3, [subobj('SPIN_Z',1,1), subobj('PC_OUT',1,1)], ele)) return
      if (.not. parser_read_subobj_real_array(obj3, 'SPIN_Z', prob_pc_r%spin_z(k,2:), .true., ele)) return
      ! Polarization must be symmetric so make this approximation.
      prob_pc_r%spin_z(k,1) = prob_pc_r%spin_z(k,2)
      if (.not. parser_read_subobj_real(obj3, 'PC_OUT', val, ele)) return
      if (val /= prob_pc_r%pc_out(k)) then
        call parser_error('PC_OUT VALUE MISMATCH IN CONVERTER DISTRIBUTION BETWEEN PROB_PC_R AND SPIN_Z_OUT MATRICES.', &
                          'FOR ELEMENT: ' // ele%name)
      endif
    enddo
  endif

  ! DIRECTION_OUT
  obj2 => pointer_to_subobject(obj1, 'DIRECTION_OUT')
  if (.not. parser_check_subobjects (obj2, [subobj('ALPHA_X',1,1), subobj('ALPHA_Y',1,1), &
                                 subobj('BETA',1,1), subobj('C_X',1,1), subobj('DXDS_MIN',1,1), &
                                 subobj('DXDS_MAX',1,1), subobj('DYDS_MAX',1,1)], ele)) return
  dir_out => sub_d%dir_out

  ! Alpha_x, etc

  if (.not. parse_this_dir_param('ALPHA_X', ele, obj2, dir_out%alpha_x, .false.))    return
  if (.not. parse_this_dir_param('ALPHA_Y', ele, obj2, dir_out%alpha_y, .false.))    return
  if (.not. parse_this_dir_param('BETA', ele, obj2, dir_out%beta, .false.))          return
  if (.not. parse_this_dir_param('C_X', ele, obj2, dir_out%c_x, .false.))            return
  if (.not. parse_this_dir_param('DXDS_MIN', ele, obj2, dir_out%dxds_min, .true.))   return
  if (.not. parse_this_dir_param('DXDS_MAX', ele, obj2, dir_out%dxds_max, .false.))  return
  if (.not. parse_this_dir_param('DYDS_MAX', ele, obj2, dir_out%dyds_max, .false.))  return

enddo

err_flag = .false.

!------------------------
contains

function get_more_text_func (line, end_of_document, why_invalid) result (valid)

integer ix_word
character(:), allocatable :: line, why_invalid
character(200) str
logical end_of_document, valid, err_flag

!

call get_next_word(str, ix_word, '}], ', delim, delim_found, err_flag = err_flag)
end_of_document = (len(str) /= 0 .and. .not. delim_found)
line = trim(line) // ' ' // trim(str) // delim
valid = .not. err_flag
if (err_flag) call str_set(why_invalid, 'ERROR WHILE PARSING CONVERTER DISTRIBUTION')

end function get_more_text_func

!------------------------
! contains

function parse_this_dir_param (who, ele, obj2, coef_head, has_c) result (is_ok)

type (ele_struct) ele
type (object_struct), pointer :: obj2, obj3, obj4
type (converter_dir_coef_struct) coef_head

integer k, n_r, ixe
character(*) who
logical has_c, is_ok

!

is_ok = .false.
obj3 => pointer_to_subobject(obj2, who)

if (has_c) then
  if (.not. parser_check_subobjects (obj3, [subobj('FIT_1D_R',0,9999), subobj('FIT_2D_PC',1,1), &
                                                 subobj('FIT_2D_R',1,1), subobj('C',1,1)], ele)) return
  if (.not. parser_read_subobj_real(obj3, 'C', coef_head%c0, ele)) return
else
  if (.not. parser_check_subobjects (obj3, [subobj('FIT_1D_R',1,9999), subobj('FIT_2D_PC',1,1), &
                                                 subobj('FIT_2D_R',1,1)], ele)) return
endif

n_r = subobject_number(obj3, 'FIT_1D_R')
allocate(coef_head%fit_1d_r(n_r))
ixe = 0

do k = 1, n_r
  obj4 => pointer_to_subobject(obj3, 'FIT_1D_R', ixe)
  if (.not. parser_check_subobjects (obj4, [subobj('PC_OUT',1,1), subobj('POLY',1,1)], ele)) return
  if (.not. parser_read_subobj_real(obj4, 'PC_OUT', coef_head%fit_1d_r(k)%pc_out, ele)) return
  if (.not. parser_read_subobj_real_array(obj4, 'POLY', coef_head%fit_1d_r(k)%poly, .false., ele)) return
enddo

if (.not. parser_check_obj_order (coef_head%fit_1d_r(:)%pc_out, '<', 'FIT_1D_R->PC_OUT', obj3, ele)) return

obj4 => pointer_to_subobject(obj3, 'FIT_2D_PC')
if (.not. parser_check_subobjects (obj4, [subobj('K',1,1), subobj('POLY',1,1)], ele)) return
if (.not. parser_read_subobj_real(obj4, 'K', coef_head%fit_2d_pc%k, ele)) return
if (.not. parser_read_subobj_real_array(obj4, 'POLY', coef_head%fit_2d_pc%poly, .false., ele)) return

obj4 => pointer_to_subobject(obj3, 'FIT_2D_R')
if (.not. parser_check_subobjects (obj4, [subobj('K',1,1), subobj('POLY',1,1)], ele)) return
if (.not. parser_read_subobj_real(obj4, 'K', coef_head%fit_2d_r%k, ele)) return
if (.not. parser_read_subobj_real_array(obj4, 'POLY', coef_head%fit_2d_r%poly, .false., ele)) return

is_ok = .true.

end function parse_this_dir_param

!------------------------
! contains

function parser_check_subobjects (obj, subobj_info, ele) result (valid)

type (object_struct) obj
type (subobj), target :: subobj_info(:)
type (subobj), pointer :: si
type (ele_struct) ele

integer i, n
logical valid

!

do i = 1, obj%n_child
  if (.not. any (obj%child(i)%name == subobj_info%name)) then
    valid = this_error_out ('UNKNOWN COMPONENT: ' // obj%child(i)%name, obj, ele)
    return
  endif
enddo

do i = 1, size(subobj_info)
  si => subobj_info(i)
  n = subobject_number(obj, si%name)

  if (n < si%n_min .and. n == 0) then
    valid = this_error_out ('COMPONENT MISSING: ' // si%name, obj, ele)
    return
  endif

  if (n < si%n_min) then
    valid = this_error_out ('NOT ENOUGH COMPONENTS NAMED: ' // si%name, obj, ele)
    return
  endif

  if (n > si%n_max .and. n == 2) then
    valid = this_error_out ('MULTIPLE COMPONENTS NAMED: ' // si%name, obj, ele)
    return
  endif

  if (n > si%n_max) then
    valid = this_error_out ('TOO MANY COMPONENTS NAMED: ' // si%name, obj, ele)
    return
  endif
enddo

valid = .true.

end function parser_check_subobjects

!---------------------------------
! contains

function this_error_out (string, obj, ele) result (valid)

type (object_struct) obj
type (ele_struct) ele
logical valid
character(*) string

!

valid = .false.
call parser_error (string, 'AT: ' // object_tree_name(obj), 'IN ELEMENT: ' // ele%name)

end function this_error_out

!---------------------------------
! contains

function parser_read_subobj_real (obj, name, var, ele) result (valid)

type (object_struct), target :: obj
type (object_struct), pointer :: sub
type (ele_struct) ele
real(rp) var
logical valid
character(*) name

!

sub => pointer_to_subobject(obj, name)
valid = parser_read_obj_real(sub, var, ele)

end function parser_read_subobj_real

!---------------------------------
! contains

function parser_read_obj_real (obj, var, ele) result (valid)

type (object_struct) obj
type (ele_struct) ele
real(rp) var
logical valid

!

if (obj%n_token == 0) then
  valid = this_error_out ('COMPONENT DOES NOT HAVE A VALUE', obj, ele)
  return  
endif

if (obj%n_token > 1) then
  valid = this_error_out ('COMPONENT HAS MULTIPLE VALUES. EXPECTING ONLY ONE.', obj, ele)
  return  
endif

valid = is_real(obj%token(1)%str, real_num = var)
if (.not. valid) then
  valid = this_error_out('COMPONENT VALUE: ' // obj%token(1)%str // &
                                    ' CANNOT BE EVALUATED AS A REAL NUMBER.', obj, ele)
endif

end function parser_read_obj_real

!---------------------------------
! contains

function parser_read_subobj_real_array (obj, name, var, exact_size, ele) result (valid)

type (object_struct), target :: obj
type (object_struct), pointer :: sub
type (ele_struct) ele
real(rp) var(:)
logical exact_size, valid
character(*) name

!

sub => pointer_to_subobject(obj, name)
valid = parser_read_obj_real_array (sub, var, exact_size, ele)

end function parser_read_subobj_real_array

!---------------------------------
! contains

function parser_read_obj_real_array (obj, var, exact_size, ele) result (valid)

type (object_struct) obj
type (ele_struct) ele
real(rp) var(:)
integer i
logical exact_size, valid

!

if (obj%n_token == 0) then
  valid = this_error_out ('COMPONENT DOES NOT HAVE A VALUE', obj, ele)
  return  
endif

if ((exact_size .and. obj%n_token /= size(var)) .or. obj%n_token > size(var)) then
  valid = this_error_out ('COMPONENT ARRAY SIZE IS ' // int_str(obj%n_token) // &
                              '. EXPECTED SIZE TO BE: ' // int_str(size(var)), obj, ele)
  return  
endif

var = 0
do i = 1, obj%n_token
  valid = is_real(obj%token(i)%str, real_num = var(i))
  if (.not. valid) then
    valid = this_error_out('COMPONENT VALUE: ' // obj%token(i)%str // &
                                    ' CANNOT BE EVALUATED AS A REAL NUMBER.', obj, ele)
    return
  endif
enddo

end function parser_read_obj_real_array

!---------------------------------
! contains

function parser_check_obj_order (array, order, who, obj, ele) result (valid)

type (object_struct) obj
type (ele_struct) ele

real(rp) array(:)
integer i
logical valid
character(*) who, order

!

valid = .true.

do i = 2, size(array)
  select case (order)
  case ('<')
    if (array(i-1) >= array(i)) valid = this_error_out (who // 'COMPONENT NOT IN STRICTLY ASSENDING ORDER.', obj, ele)
  case ('<=')
    if (array(i-1) > array(i)) valid = this_error_out (who // 'COMPONENT NOT IN ASSENDING ORDER.', obj, ele)
  case ('>')
    if (array(i-1) <= array(i)) valid = this_error_out (who // 'COMPONENT NOT IN STRICTLY DESENDING ORDER.', obj, ele)
  case ('>=')
    if (array(i-1) < array(i)) valid = this_error_out (who // 'COMPONENT NOT IN DESENDING ORDER.', obj, ele)
  case default
    call err_exit
  end select
enddo

end function

end subroutine converter_distribution_parser

