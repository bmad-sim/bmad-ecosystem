!+
! Module for parsing a document with a tree format of key/value pairs.
!
! The document model is a tree of objects. 
! Each object is a key/value pair.
! An object will either have an array of one or more tokens or 
! an array of of one or more subobjects but not both.
! tokens are either quoted strings or non-quoted strings with no embeded blanks or tabs.
!
! Example:
!   a = {b = 7*3, d = "xyz", e = [1, "true", rats], f = {g = 3^34}}
!
! Notes:
!   Subobject arrays are delimited using curley braces "{...}".
!   A quoted token may use single or double quotation marks for delimitors.
!   To include a quotation mark in a quoted token, escape the mark with the backslash character.
!   A non-quoted token may contain the characters "{", "[", or "(" as long 
!      as these are not the first character in the token. 
!   In a non-quoted token, "{...}", "(...)", and "[...]" constructs must be properly matched.
!   Token arrays may be delimited using "[...]" or "(...)".
!   Something like "a = OK" is the same as "a = [OK]".
!   Mixing quoted and non-quoted tokens in an array is allowed.
!   Mixing tokens and objects in an array is forbidden. Example:
!     a = [3, b = 8]      ! Not allowed!
!   Having multiple keys of the same name is allowed. Example:
!     a = r, a = 7, a = r
!   When the document is parsed, order will be preserved.
!-

module object_model_mod

use sim_utils

implicit none

type object_model_struct
  character(:), allocatable :: name
  type (var_length_string_struct), allocatable :: token(:) ! Array of tokens
  type (object_model_struct), pointer :: obj(:) => null() ! Array of subobjects
  integer :: n_obj = 0
  integer :: n_token = 0
end type

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function object_model_parse(obj, line, why_invalid, get_more_text_func) result (valid)
!
! Object model parsing routine.
! The initial "{" is optional.
!-

recursive function object_model_parse(obj, line, why_invalid, get_more_text_func) result (valid)

interface 
  subroutine get_more_text_func(line, valid)
    character(:), allocatable :: line
    logical valid
  end subroutine
end interface

type (object_model_struct), target :: obj
type (object_model_struct), pointer :: optr

integer n, ia, ip, parse_state, token_state, n_len, nb1, nb2, n_item
integer, allocatable :: ixm(:)
integer, parameter :: before_obj_name$ = 1, in_obj_name$ = 2, after_obj_name$ = 3, after_equal_sign$ = 4
integer, parameter :: void$ = 11, in_quote$ = 12, before_token$ = 13, in_token$ = 14, after_token$ = 15

logical valid, abort

character(*) why_invalid
character(:), allocatable :: line
character(*), parameter :: name_chars = '0123456789_'
character(1), parameter :: tab = char(9)
character(1) ach, quote_mark, array_mark

!

valid = .true.
allocate(character(100):: line)
allocate(ixm(10))

ip = 0
if (line(1:1) == '{') ip = 1
n_len = len(line)

parse_state = before_obj_name$
token_state = void$

do while (.true.)
  ip = ip + 1
  if (ip > n_len) then
    call get_more_text_func(line, valid)
    if (.not. valid) then
      why_invalid = 'ABORT'
      return
    endif
    n_len = len(line)
  endif

  ach = line(ip:ip)

  ! Handle quote character

  if (ach == '"' .or. ach == "'") then
    select case (token_state)
    case (before_token$)
      quote_mark = ach
      token_state = in_quote$
    case (in_quote$)
      if (line(ip-1:ip-1) == '\') cycle      ! '
      if (ach /= quote_mark) cycle
      token_state = after_token$
    case default
      valid = .false.
      why_invalid = 'QUOTE MARK OUT OF PLACE'
      return
    end select
    cycle
  endif

  ! Handle non-quote character

  select case (parse_state)
  ! Before a variable name
  case (before_obj_name$)
    if (ach == ' ' .or. ach == tab) cycle
    if (.not. is_alphabetic(ach)) then
      valid = .false.
      why_invalid = 'DATUM NAME DOES NOT START WITH AN ALPHEBETIC CHARACTER.'
      return
    endif
    parse_state = in_obj_name$

  ! Parsing a variable name
  case (in_obj_name$)
    if (is_alphabetic(ach)) cycle
    if (index(name_chars, ach) /= 0) cycle
    parse_state = after_obj_name$

  ! After a variable name before the equal sign
  case (after_obj_name$)
    if (ach == ' ' .or. ach == tab) cycle
    if (ach /= '=') then
      valid = .false.
      why_invalid = 'EXPECTING EQUAL SIGN AFTER DATUM NAME'
      return
    endif
    parse_state = after_equal_sign$
    token_state = before_token$

  ! After an equal sign.
  case (after_equal_sign$)
    select case (ach)
    case (' ', tab)
      cycle

    case ('{')
      line = line(ip+1:)
      valid = object_model_parse(obj, line, why_invalid, get_more_text_func)
      return

    case default
      parse_state = in_token$
      token_state = before_token$
      if (ach == '(' .or. ach == '[') then
        array_mark = ach
      else
        array_mark = ' '
      endif
      nb1 = 0;  nb2 = 0
      n_item = 0
      ixm(1) = ip-1
    end select

  ! In an array of values.
  case (in_token$)
    select case (ach)
    case ('[');  nb1 = nb1 + 1
    case ('(');  nb2 = nb2 + 1
    case (']', ')')
      if (ach == array_mark .and. ach /= ' ' .and. nb1 == 0 .and. nb2 == 0) then
        ixm(n_item+1) = ip
        call re_allocate (obj%token, n_item)
        do ia = 1, n_item
          n = ixm(ia+1) - ixm(ia)
          if (.not. allocated(obj%token(ia)%str)) allocate (character(n):: obj%token(ia)%str)
          obj%token(ia)%str = line(ixm(ia)+1:ixm(ia+1)-1)
        enddo
        return

      elseif (ach == ']') then
        nb1 = nb1 - 1
        if (nb1 < 0) then
          valid = .false.
          why_invalid = 'MISPLACED "]" BRACKET.'
          return
        endif
      else
        nb2 = nb2 - 1
        if (nb2 < 0) then
          valid = .false.
          why_invalid = 'MISPLACED ")" BRACKET.'
          return
        endif
      endif
    case (',')
      if (nb1 /= 0 .or. nb2 /= 0) cycle
      if (array_mark == ' ') then
        call re_allocate(obj%token, 1)
        n = ip - ixm(1)
        if (.not. allocated(obj%token(1)%str)) allocate (character(n):: obj%token(1)%str)
        obj%token(1)%str = line(ixm(1)+1:ip-1)
        return
      else
        n_item = n_item + 1
        if (size(ixm) == n_item) call re_allocate(ixm, 2*n_item)
        ixm(n_item+1) = ia
      endif
    end select
  end select
enddo

end function object_model_parse

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine object_print(obj, n_indent, last)
!
! Routine to print the contents of an object.
!
! Input:
!   obj       -- object_model_struct: Object to be printed.
!   n_indent  -- integer, optional: Indentation when printing.
!   last      -- logical, optional: Is last datum in array? If True then no
!-

recursive subroutine object_print (obj, n_indent, last)

type (object_model_struct), target :: obj
type (object_model_struct), pointer :: optr
character(100) :: blank = ''
character(:), allocatable :: line
integer, optional :: n_indent
integer i, ni
logical, optional :: last

! Head object may not have a name.
! Note: In this case, the object must have subobjects and not a value

ni = integer_option(0, n_indent)

if (obj%n_obj > 0) then
  if (obj%name == '') then
    print '(2a)', blank(1:ni), '{'
  else
    print '(3a)', blank(1:ni), obj%name, ' = {'
  endif

  do i = 1, obj%n_obj
    call object_print(obj%obj(i), ni+2, last = (i == obj%n_obj))
  enddo

  if (logic_option(.false., last)) then
    print '()', blank(1:ni), '}'
  else
    print '()', blank(1:ni), '},'
  endif

elseif (obj%n_token == 1) then
  print '(4a)', blank(1:ni), obj%name, ' = ', obj%token(1)%str

else
  allocate (character(40):: line)
  line = blank(1:ni) // obj%name // ' = [' // obj%token(1)%str
  do i = 2, obj%n_token
    line = line // ', ' // obj%token(i)%str
  enddo
  if (logic_option(.false., last)) then
    print '(a)', line // '],'
  else
    print '(a)', line // ']'
  endif
endif

end subroutine object_print

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine subobject_reallocate(obj, n)
!
! Routine to reallocate the obj%obj(:) array.
!
! Input:
!   obj   -- object_model_struct: Object.
!   n     -- integer: Minimum size for obj%obj(:)
!-

subroutine subobject_reallocate(obj, n)

type (object_model_struct) obj, temp
integer n

!

if (.not. associated(obj%obj)) then
  allocate (obj%obj(n))
  obj%n_obj = 0
  return
endif

if (size(obj%obj) >= n) return

temp%obj => obj%obj
allocate(obj%obj(n))
obj%obj(1:size(temp%obj)) = temp%obj
deallocate(temp%obj)

end subroutine subobject_reallocate

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine object_deallocate(obj)
!
! Routine to deallocate an object
!
! Input:
!   obj   -- object_model_struct: Object.
!
! Output:
!   obj   -- object_model_struct: Deallocated object.
!-

recursive subroutine object_deallocate (obj)

type (object_model_struct) obj
integer i

!

if (obj%n_obj > 0) then
  do i = 1, obj%n_obj
    call object_decallocate (obj%obj(i))
  enddo
  deallocate(obj%obj)
endif 

end subroutine object_deallocate

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function subobject_number (obj, name) result (num)
!
! Routine to return how many subobjects match the name.
!
! Input:
!   obj     -- object_model_struct: Object.
!   name    -- character(*): Name of subobject
!
! Output:
!   num     -- integer: Number of subobjects whose name matches name argument.
!-

function subobject_number (obj, name) result (num)

type (object_model_struct) obj, temp
integer num, n
character(*) name

!

num = 0
do n = 1, obj%n_obj
  if (obj%obj(n)%name == name) num = num + 1
enddo

end function subobject_number

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function subobject_ptr (obj, name, indx) result (subobj)  
!
! Routine to return a pointer to a subobject with a given name.
! The indx argument is used when searching for multiple subobjects of the same name.
!
! Input:
!   obj     -- object_model_struct: Object.
!   name    -- character(*): Name of subobject
!   indx    -- integer, optional: If present, start subobject search from this index + 1.
!                 Default is 0.
!
! Output:
!   indx    -- integer, optional: Index of subobject found.
!   subobj  -- object_model_struct, pointer: Pointer to subobject. Null if not found.
!-

function subobject_ptr (obj, name, indx) result (subobj)  

type (object_model_struct), target :: obj
type (object_model_struct), pointer :: subobj
integer, optional :: indx
integer n, n0
character(*) name

!

subobj => null()

n0 = integer_option(0, indx) + 1
do n = n0, obj%n_obj
  if (obj%obj(n)%name /= name) cycle
  subobj => obj%obj(n)
  if (present(indx)) indx = n
  return
enddo

end function subobject_ptr

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function object_token_read (obj, ix_token, real_val, int_val, logic_val, str_val) result (valid)
!
! Routine to evaluate a given token in the token array of an object.
! One and only one of real_val, int_val, logic_val, and str_val should be present
! depending upon the expected type of the value.
!
! Input:
!   obj         -- object_model_struct: Object.
!   ix_token    -- integer: Index of token.
!
! Output:
!   real_val    -- real(rp), optional: Real value.
!   int_val     -- integer, optional: Integer value.
!   logic_val   -- logical, optional: Logical value. 
!   str_val     -- character(*), optional: String value.
!   valid       -- logical: Valid value found.
!-

function object_token_read (obj, ix_token, real_val, int_val, logic_val, str_val) result (valid)

type (object_model_struct), target :: obj
real(rp), optional :: real_val
integer, optional :: int_val
logical, optional :: logic_val
character(*), optional :: str_val

integer ix_token, ios
logical valid

!

valid = .false.
if (ix_token < 0 .or. ix_token > obj%n_token) return

if (present(real_val)) then
  read (obj%token(ix_token)%str, *, iostat = ios) real_val

elseif (present(int_val)) then
  read (obj%token(ix_token)%str, *, iostat = ios) int_val

elseif (present(logic_val)) then
  read (obj%token(ix_token)%str, *, iostat = ios) logic_val

elseif (present(str_val)) then
  str_val = obj%token(ix_token)%str
  return

else
  return
endif

valid = (ios == 0)

end function object_token_read

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function object_array_read (obj, real_arr, int_arr, logic_arr, exact_size) result (valid)
!
! Routine to evaluate the token array of an object.
! One and only one of real_arr, int_arr, or logic_arr should be present
! depending upon the expected type of the value.
!
! Input:
!   obj         -- object_model_struct: Object.
!   exact_size  -- logical, optional: Can the size of the array in obj be smaller than the array argument? 
!                    Default is False. It is always an error if the obj array size is larger.
!
! Output:
!   real_arr(:)   -- real(rp), optional: Real value.
!   int_arr(:)    -- integer, optional: Integer value.
!   logic_arr(:)  -- logical, optional: Logical value.
!   valid         -- logical: Valid value found.
!-

function object_array_read (obj, real_arr, int_arr, logic_arr, exact_size) result (valid)

type (object_model_struct), target :: obj
type (var_length_string_struct), allocatable :: array(:)

real(rp), optional :: real_arr(:)
integer, optional :: int_arr(:)
logical, optional :: logic_arr(:), exact_size
logical valid

integer i, arr_size

!

if (present(real_arr)) then
  arr_size = size(real_arr)
elseif (present(int_arr)) then
  arr_size = size(int_arr)
elseif (present(logic_arr)) then
  arr_size = size(logic_arr)
endif

valid = .false.
if (obj%n_token == 0) return
if (obj%n_token > arr_size) return
if (logic_option(.true., exact_size) .and. obj%n_token /= arr_size) return

do i = 1, arr_size
  if (present(real_arr)) then
    valid = object_token_read(obj, i, real_val = real_arr(i))
  elseif (present(int_arr)) then
    valid = object_token_read(obj, i, int_val = int_arr(i))
  elseif (present(logic_arr)) then
    valid = object_token_read(obj, i, logic_val = logic_arr(i))
  endif
  if (.not. valid) return
enddo

end function object_array_read

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function object_alloc_array_read (obj, real_arr, int_arr, logic_arr) result (valid)
!
! Routine to evaluate the array value of an object which.
! One and only one of real_arr, int_arr, or logic_arr should be present
! depending upon the expected type of the value.
!
! Input:
!   obj         -- object_model_struct: Object.
!   name        -- character(*): Name of subobject
!
! Output:
!   real_arr(:)   -- real(rp), allocatable, optional: Real value.
!   int_arr(:)    -- integer, allocatable, optional: Integer value.
!   logic_arr(:)  -- logical, allocatable, optional: Logical value. 
!   valid         -- logical: Valid value found.
!-

function object_alloc_array_read (obj, real_arr, int_arr, logic_arr) result (valid)

type (object_model_struct), target :: obj
type (object_model_struct), pointer :: subobj
type (var_length_string_struct), allocatable :: array(:)

real(rp), allocatable, optional :: real_arr(:)
integer, allocatable, optional :: int_arr(:)
logical, allocatable, optional :: logic_arr(:)
logical valid

!

valid = .false.
if (obj%n_token == 0) return

if (present(real_arr)) then
  call re_allocate (real_arr, obj%n_token)
elseif (present(int_arr)) then
  call re_allocate (int_arr, obj%n_token)
elseif (present(logic_arr)) then
  call re_allocate (logic_arr, obj%n_token)
endif

valid = object_array_read (obj, real_arr, int_arr, logic_arr) 

end function object_alloc_array_read 

end module





