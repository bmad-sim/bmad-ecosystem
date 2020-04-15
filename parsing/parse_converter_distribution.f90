subroutine parse_converter_distribution(ele, delim, delim_found)

use bmad_parser_mod, dummy => parse_converter_distribution
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
type (converter_alpha_struct), pointer :: alpha_x, alpha_y
type (converter_beta_struct), pointer :: beta
type (converter_c_x_struct), pointer :: c_x
type (object_struct), target :: obj0
type (object_struct), pointer :: obj1, obj2, obj3, obj4

integer i, j, k, n, ix_word, id, n_r, n_pc, ixe, n_sd, n_subd, ix_sd

character(*) delim
character(200) str
character(:), allocatable :: line, why_invalid

logical delim_found, valid

!

valid = object_document_parse (obj0, 'DISTRIBUTION', line, why_invalid, get_more_text_func, parse_one = .true.)
bp_com%parse_line = trim(line) // ' ' // bp_com%parse_line
call get_next_word(str, ix_word, '}],', delim, delim_found)

if (.not. valid) then
  call parser_error (why_invalid, 'FOR ELEMENT: ' // ele%name)
  return
endif

!

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
                                                                  subobj('DXY_DS_MAX',1,1)], ele)) return

if (.not. parser_read_subobj_real(obj0, 'DXY_DS_MAX', dist%dxy_ds_max, ele)) return
if (.not. parser_read_subobj_real(obj0, 'THICKNESS', dist%thickness, ele)) return

ix_sd = 0
do n_sd = 1, n_subd
  obj1 => pointer_to_subobject(obj0, 'SUB_DISTRIBUTION', ix_sd)

  if (.not. parser_check_subobjects (obj1, [subobj('PC_IN',1,1), subobj('PROB_PC_R',1,1)], ele)) return
  id = id + 1
  sub_d => dist%sub_dist(id)
  do j = 1, obj1%n_child
    obj2 => obj1%child(j)
    select case (obj2%name)
    case ('PC_IN')
      if (.not. parser_read_obj_real(obj2, sub_d%pc_in, ele)) return

    case ('PROB_PC_R')
      if (.not. parser_check_subobjects (obj2, [subobj('R_VALUES',1,1), subobj('ROW',1,9999)], ele)) return
      prob_pc_r => sub_d%prob_pc_r
      obj3 = pointer_to_subobject(obj2, 'R_VALUES')
      n_r = obj3%n_token
      n_pc = subobject_number(obj2, 'ROW')
      allocate(prob_pc_r%r(n_r))
      if (.not. parser_read_obj_real_array(obj3, sub_d%prob_pc_r%r, ele)) return
      allocate (prob_pc_r%prob(n_pc, n_r), prob_pc_r%pc_out(n_pc))
      ixe = 0
      do k = 1, n_pc
        obj3 => pointer_to_subobject(obj2, 'ROW', ixe)
        if (.not. parser_check_subobjects (obj3, [subobj('PC_OUT',1,1), subobj('PROB',1,1)], ele)) return
        if (.not. parser_read_subobj_real(obj3, 'PC_OUT', prob_pc_r%pc_out(k), ele)) return
        if (.not. parser_read_subobj_real_array(obj3, 'PROB', prob_pc_r%prob(k,:), ele)) return
      enddo
      if (.not. parser_check_obj_order (prob_pc_r%pc_out(:), '<', 'ROW->PC_OUT', obj2, ele)) return

    case ('DIRECTION_OUT')
      if (.not. parser_check_subobjects (obj2, [subobj('ALPHA_X',1,1), subobj('ALPHA_Y',1,1), &
                                                subobj('BETA',1,1), subobj('C_X',1,1)], ele)) return
      dir_out => sub_d%dir_out

      ! Alpha_x
      obj3 => pointer_to_subobject(obj2, 'ALPHA_X')
      alpha_x => dir_out%alpha_x
      if (.not. parser_check_subobjects (obj3, [subobj('FIT_1D_R',1,9999), subobj('FIT_2D_PC',1,1), subobj('FIT_2D_R',1,1)], ele)) return
      n_r = subobject_number(obj3, 'FIT_1D_R')
      allocate(alpha_x%fit_1d_r(n_r))
      ixe = 0
      do k = 1, n_r
        obj4 => pointer_to_subobject(obj3, 'FIT_1D_R', ixe)
        if (.not. parser_check_subobjects (obj4, [subobj('PC_OUT',1,1), subobj('K',1,1), subobj('POLY',1,1)], ele)) return
        if (.not. parser_check_subobjects (obj4, [subobj('FIT_1D_R',1,9999), subobj('FIT_1D_PC',1,1)], ele)) return
        if (.not. parser_read_subobj_real(obj4, 'PC_OUT', alpha_x%fit_1d_r(k)%pc_out, ele)) return
        if (.not. parser_read_subobj_real(obj4, 'K', alpha_x%fit_1d_r(k)%k, ele)) return
        if (.not. parser_read_subobj_real_array(obj3, 'POLY', alpha_x%fit_1d_r(k)%poly, ele)) return
      enddo
      if (.not. parser_check_obj_order (alpha_x%fit_1d_r(:)%pc_out, '<', 'FIT_1D_R->PC_OUT', obj3, ele)) return
      obj4 => pointer_to_subobject(obj3, 'FIT_2D_PC')
      if (.not. parser_check_subobjects (obj4, [subobj('K',1,1), subobj('POLY',1,1)], ele)) return
      if (.not. parser_read_subobj_real(obj4, 'K', alpha_x%fit_2d_pc%k, ele)) return
      if (.not. parser_read_subobj_real_array(obj4, 'POLY', alpha_x%fit_2d_pc%poly, ele)) return
      obj4 => pointer_to_subobject(obj3, 'FIT_2D_R')
      if (.not. parser_check_subobjects (obj4, [subobj('K',1,1), subobj('POLY',1,1)], ele)) return
      if (.not. parser_read_subobj_real(obj4, 'K', alpha_x%fit_2d_r%k, ele)) return
      if (.not. parser_read_subobj_real_array(obj4, 'POLY', alpha_x%fit_2d_r%poly, ele)) return

      ! Alpha_y
      obj3 => pointer_to_subobject(obj2, 'ALPHA_Y')
      alpha_y => dir_out%alpha_y
      if (.not. parser_check_subobjects (obj3, [subobj('FIT_1D_R',1,9999), subobj('FIT_2D_PC',1,1), subobj('FIT_2D_R',1,1)], ele)) return
      n_r = subobject_number(obj3, 'FIT_1D_R')
      allocate(alpha_y%fit_1d_r(n_r))
      ixe = 0
      do k = 1, n_r
        obj4 => pointer_to_subobject(obj3, 'FIT_1D_R', ixe)
        if (.not. parser_check_subobjects (obj4, [subobj('PC_OUT',1,1), subobj('K',1,1), subobj('POLY',1,1)], ele)) return
        if (.not. parser_read_subobj_real(obj4, 'PC_OUT', alpha_y%fit_1d_r(k)%pc_out, ele)) return
        if (.not. parser_read_subobj_real(obj4, 'K', alpha_y%fit_1d_r(k)%k, ele)) return
        if (.not. parser_read_subobj_real_array(obj3, 'POLY', alpha_y%fit_1d_r(k)%poly, ele)) return
      enddo
      if (.not. parser_check_obj_order (alpha_y%fit_1d_r(:)%pc_out, '<', 'FIT_1D_R->PC_OUT', obj3, ele)) return
      obj4 => pointer_to_subobject(obj3, 'FIT_2D_PC')
      if (.not. parser_check_subobjects (obj4, [subobj('K',1,1), subobj('POLY',1,1)], ele)) return
      if (.not. parser_read_subobj_real(obj4, 'K', alpha_y%fit_2d_pc%k, ele)) return
      if (.not. parser_read_subobj_real_array(obj4, 'POLY', alpha_y%fit_2d_pc%poly, ele)) return
      obj4 => pointer_to_subobject(obj3, 'FIT_2D_R')
      if (.not. parser_check_subobjects (obj4, [subobj('K',1,1), subobj('POLY',1,1)], ele)) return
      if (.not. parser_read_subobj_real(obj4, 'K', alpha_y%fit_2d_r%k, ele)) return
      if (.not. parser_read_subobj_real_array(obj4, 'POLY', alpha_y%fit_2d_r%poly, ele)) return

      ! Beta
      obj3 => pointer_to_subobject(obj2, 'BETA')
      beta => dir_out%beta
      if (.not. parser_check_subobjects (obj3, [subobj('FIT_1D_R',1,9999), subobj('POLY_PC',1,1), subobj('POLY_R',1,1)], ele)) return
      n_r = subobject_number(obj3, 'FIT_1D_R')
      allocate(beta%fit_1d_r(n_r))
      ixe = 0
      do k = 1, n_r
        obj4 => pointer_to_subobject(obj3, 'FIT_1D_R', ixe)
        if (.not. parser_check_subobjects (obj4, [subobj('PC_OUT',1,1), subobj('POLY',1,1)], ele)) return
        if (.not. parser_read_subobj_real(obj4, 'PC_OUT', beta%fit_1d_r(k)%pc_out, ele)) return
        if (.not. parser_read_subobj_real_array(obj3, 'POLY', beta%fit_1d_r(k)%poly, ele)) return
      enddo
      if (.not. parser_check_obj_order (beta%fit_1d_r(:)%pc_out, '<', 'FIT_1D_R->PC_OUT', obj3, ele)) return
      if (.not. parser_read_subobj_real_array(obj3, 'POLY_PC', beta%poly_pc, ele)) return
      if (.not. parser_read_subobj_real_array(obj3, 'POLY_R', beta%poly_r, ele)) return

      ! C_x
      obj3 => pointer_to_subobject(obj2, 'C_X')
      c_x => dir_out%c_x
      if (.not. parser_check_subobjects (obj3, [subobj('POLY_R',1,1), subobj('POLY_PC',1,1)], ele)) return
      if (.not. parser_read_subobj_real_array(obj3, 'POLY_R', c_x%poly_r, ele)) return
      if (.not. parser_read_subobj_real_array(obj3, 'POLY_PC', c_x%poly_pc, ele)) return
    end select
  enddo
enddo

!print *
!print *, '----------------'
!print *, delim, ': ', trim(bp_com%parse_line)
!print *
!
!call object_print(obj)
!stop

!------------------------
contains

function get_more_text_func (line, end_of_document) result (valid)

integer ix_word
character(:), allocatable :: line
character(200) str
logical end_of_document, valid

!

call get_next_word(str, ix_word, '}], ', delim, delim_found)
end_of_document = (len(str) /= 0 .and. .not. delim_found)
valid = .true.
line = trim(line) // ' ' // trim(str) // delim


end function get_more_text_func

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
call parser_error (string, 'AT: ' // object_tree_name(obj), &
                           'IN ELEMENT: ' // ele%name)

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

function parser_read_subobj_real_array (obj, name, var, ele) result (valid)

type (object_struct), target :: obj
type (object_struct), pointer :: sub
type (ele_struct) ele
real(rp) var(:)
logical valid
character(*) name

!

sub => pointer_to_subobject(obj, name)
valid = parser_read_obj_real_array (sub, var, ele)

end function parser_read_subobj_real_array

!---------------------------------
! contains

function parser_read_obj_real_array (obj, var, ele) result (valid)

type (object_struct) obj
type (ele_struct) ele
real(rp) var(:)
integer i
logical valid

!

if (obj%n_token == 0) then
  valid = this_error_out ('COMPONENT DOES NOT HAVE A VALUE', obj, ele)
  return  
endif

if (obj%n_token /= size(var)) then
  valid = this_error_out ('COMPONENT ARRAY SIZE IS ' // int_str(obj%n_token) // &
                              '. EXPECTED SIZE TO BE: ' // int_str(size(var)), obj, ele)
  return  
endif

do i = 1, size(var)
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

end subroutine parse_converter_distribution

