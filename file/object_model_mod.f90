!+
! Module for parsing a document with a tree format of key/value pairs.
!
! A key/value pair is called an "object".
! The document model is a tree of objects. 
! The "value" of an object will either be an array of one or more tokens or 
! an array of of one or more subobjects (but not both).
!
! Example:
!   {a = {b = 7*atan(3, 4), d = "xyz", e = [1, "true", rats], f = {g = 3^34}}}
! In this example, the root object, which has a blank key name, has a single subobject
! whose (key) name is "a". Object "a" has three subobjects whose names are "b", "d", e", and "f".
! Object "b" has a single token "7*atan(3, 4)".
! Object "e" has three tokens. Etc.
!
! Notes:
!   Keys or Tokens are either quoted strings or non-quoted strings.
!   Leading and trailing whitespace in non-quoted keys and tokens is ignored.
!   To include a quotation mark in a name or token, escape the mark with the backslash character.
!   A key or token may use single or double quotation marks for delimitors.
!   In a non-quoted key or token, a tab character is converted to a single space character.
!   A non-quoted key or token may contain the "[" character as long as the first character is not a "[".
!   In a non-quoted key or token, parentheses or square brackets (that is "(...)", and "[...]") must
!      be present in matched left/right pairs.
!   Characters not allowed in non-quoted keys or tokens:
!     "{", "}", or "="
!   In a quoted key or token: Equal sign "=" is not allowed. [This rule is needed to catch missing quote marks.]
!   In a non-quoted key or token, commas are only allowed if the comma appears within
!      brackets. For example: 
!         z = atan2(3, 4)
!   For a non-quoted key, whitespace is only allowed within brackets.
!   For a non-quoted token, whitespace is allowed (but leading and trailing whitespace will be ignored). Example:
!         z = 2 * 3
!   Subobject arrays are delimited using curley braces "{...}".
!   Token arrays are delimited using square brackets "[...]".
!   A single token (not withing square backets) is equivalent to an array of a single token.
!      That is, the token
!         abc
!      is equivalent to
!         [abc]
!   Mixing quoted and non-quoted tokens in an array is allowed.
!   Mixing tokens and objects in an array is forbidden. Example:
!     a = [3, b = 8]      ! Not allowed. "3" is a token, "b = 8" is an object.
!   Having multiple keys of the same name is allowed. Example:
!     a = r, a = 7, a = r
!   The outer "{...}" braces in a document are optional.
!   When the document is parsed, key/value order will be preserved.
!   When the document is parsed, the quote marks of quoted keys or tokens will be preserved.
!     That is, quoted keys or tokens are distinct from nonquoted keys or tokens. 
!-

module object_model_mod

use sim_utils

implicit none

! One and only one of %n_obj and %n_token is nonzero.

type temp_struct
  character(40) :: str = ''
end type

type object_struct
  character(:), allocatable :: name
  type (var_length_string_struct), allocatable :: token(:) ! Array of tokens
  type (object_struct), pointer :: obj(:) => null() ! Array of subobjects
  integer :: n_obj = 0        ! Number of subobjects.
  integer :: n_token = 0      ! Number of tokens.
end type

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function object_document_parse(obj, line, why_invalid, get_more_text_func, parse_one, sub_call) result (valid)
!
! Object model parsing routine.
! Note: The initial "{" is optional.
! Note: At the end, the returned line argument will contain any unparsed text that was outside the tree.
!  For example, if the text was:
!     "{a = 5, b = 7} xxx"
!  Then the unparsed text will be " xxx".
! If the document does not contain the outer "{...}" braces:
!   if parse_one = False (the default), the entire document is parsed.
!   If parse_one = True, only one object is parsed. Example: "a = 5, b = 7" --> only "a = 5" is parsed.
!
! Input:
!   line          -- character(:), allocatable: Initial text to parse. If no text present then
!                       get_more_text_func will be called. 
!   why_invalid   -- character(*): String set to error message if there is an error.
!   get_more_text_func(line) result (valid)
!                 -- Function: Routine called to append more text to the string being parsed. 
!                     The interface is:
!                         function get_more_text_func(line, why_invalid) result (valid)
!                           character(*) why_invalid
!                           character(:), allocatable :: line
!                           logical valid
!   parse_one     -- logical, optional: If True, parse only one object and then return. Default = False.
!   sub_call      -- logical, optional: Used to signal if calling recursively. 
!                     For internal use only. Do not set this argument.
!
! output:
!   obj           -- object_struct: Root object for parsed tree.
!   line          -- character(:), allocatable: Unparsed text. 
!
!-

recursive function object_document_parse(obj, line, why_invalid, get_more_text_func, parse_one, sub_call) result (valid)

interface 
  function get_more_text_func(line, why_invalid) result (valid)
    character(:), allocatable :: line
    character(*) why_invalid
    logical valid
  end function
end interface

type (object_struct), target :: obj
type (object_struct), pointer :: optr

integer n, ns, ip, ip0, parse_state, string_state, nb1, nb2, n_item
integer, parameter :: at_start$ = 0, before_obj_name$ = 1, in_obj_name$ = 2, after_obj_name$ = 3
integer, parameter :: after_equal_sign$ = 4, parsing_tokens$ = 5, parsing_subobjects$ = 6
integer, parameter :: void$ = 11, before_string$ = 12, in_quote$ = 13, in_nonquote$ = 14, after_string$ = 15

logical, optional :: parse_one, sub_call
logical valid, abort, has_opening_bracket, in_array, deb

character(*) why_invalid
character(:), allocatable :: line
character(1), parameter :: tab = char(9)
character(1) ach, quote_mark

!

valid = .true.
deb = .true.
has_opening_bracket = .false.

if (.not. allocated(line)) then
  allocate(character(1):: line)
  line = ''
  valid = get_more_text_func(line, why_invalid)
  if (.not. valid) return
endif

ip = 0

if (present(sub_call)) then
  parse_state = before_obj_name$
else
  parse_state = at_start$
endif
string_state = before_string$

do while (.true.)
  ip = ip + 1
  if (ip > len(line)) then
    valid = get_more_text_func(line, why_invalid)
    if (.not. valid) return
  endif

  ach = line(ip:ip)
  if (deb) print '(i4, 2x, 3a)', ip, ach, ': ', line(ip+1:)

  ! Handle quote character

  if (ach == '"' .or. ach == "'") then
    select case (string_state)
    case (before_string$)
      quote_mark = ach
      string_state = in_quote$
      ip0 = ip
      select case (parse_state)
      case (before_obj_name$);   parse_state = in_obj_name$
      case (after_equal_sign$);  parse_state = parsing_tokens$
      end select
      cycle
        
    case (in_quote$)
      if (line(ip-1:ip-1) == '\') cycle      ! '
      if (ach /= quote_mark) cycle
      select case (parse_state)
      case (in_obj_name$)
        obj%name = line(ip0:ip)
        parse_state = after_obj_name$
      case (parsing_tokens$)
        obj%n_token = obj%n_token + 1
        n = obj%n_token
        call re_allocate (obj%token, n)
        if (.not. allocated(obj%token(n)%str)) allocate (character(ip-ip+1):: obj%token(n)%str)
        obj%token(n)%str = line(ip0:ip)
        line = line(ip+1:)
        ip = 0
      case default
        call err_exit
      end select
      string_state = after_string$
      quote_mark = ' '
      cycle

    case default
      valid = .false.
      why_invalid = 'QUOTE MARK OUT OF PLACE'
      return
    end select
    cycle
  endif

  if (quote_mark /= ' ') then
    if (ach == '=') then
      valid = .false.
      why_invalid = 'EQUAL SIGN "=" NOT ALLOWED IN QUOTED STRING.'
      return
    endif
    cycle
  endif

  ! Handle non quote mark character

  if (ach == tab) ach = ' '

  select case (parse_state)
  ! At start before possible '{' character. This character is considered optional.
  case (at_start$)
    if (ach == ' ') cycle
    if (ach == '{') then
      has_opening_bracket = .true.
    else
      ip = ip - 1  ! Reparse character
    endif
    line = line(ip+1:)
    call subobject_reallocate(obj, 5)
    valid = object_document_parse(obj%obj(1), line, why_invalid, get_more_text_func, sub_call = .true.)
    if (.not. valid) return
    obj%n_obj = 1
    if (logic_option(.false., parse_one) .and. .not. has_opening_bracket) return
    parse_state = parsing_subobjects$
    ip = 0

  ! Before a variable name
  case (before_obj_name$)
    if (ach == ' ') cycle
    select case (ach)
    case (',', '{', '}', '(', ')', '[', ']')
      valid = .false.
      why_invalid = 'MISPLACED ' // ach // ' CHARACTER'
      return
    end select
    nb1 = 0;  nb2 = 0
    parse_state = in_obj_name$
    ip0 = ip

  ! Parsing a variable name
  case (in_obj_name$)
    select case (ach)
    case ('[');  nb1 = nb1 + 1
    case ('(');  nb2 = nb2 + 1

    case (']')
      nb1 = nb1 - 1
      if (nb1 < 0) then
        valid = .false.
        why_invalid = 'MISPLACED "]" BRACKET.'
        return
      endif

    case (')')
      nb2 = nb2 - 1
      if (nb2 < 0) then
        valid = .false.
        why_invalid = 'MISPLACED ")" BRACKET.'
        return
      endif

    case (',')
      if (nb1 /= 0 .or. nb2 /= 0) then
        valid = .false.
        why_invalid = 'COMMA IN NAME'
        return
      endif

    case ('=')
      if (nb1 /= 0 .or. nb2 /= 0) then
        valid = .false.
        why_invalid = 'EQUAL SIGN IN NAME'
        return
      endif
      parse_state = after_obj_name$
      ip = ip - 1

    case (' ')
      if (nb1 /= 0 .or. nb2 /= 0) cycle
      parse_state = after_obj_name$
    end select

  ! After a variable name before the equal sign
  case (after_obj_name$)
    if (ach == ' ') cycle
    if (ach /= '=') then
      valid = .false.
      why_invalid = 'EXPECTING EQUAL SIGN AFTER DATUM NAME'
      return
    endif
    n = ip - ip0
    if (.not. allocated(obj%name)) allocate(character(n):: obj%name)
    obj%name = trim(line(ip0:ip-1))
    parse_state = after_equal_sign$
    string_state = before_string$

  ! After an equal sign.
  case (after_equal_sign$)
    select case (ach)
    case (' ')
      cycle

    case (',', ')', ']')
      valid = .false.
      why_invalid = 'MISPLACED "' // ach // '" CHARACTER.'
      return

    case ('{')
      line = line(ip+1:)
      call subobject_reallocate(obj, 5)
      valid = object_document_parse(obj%obj(1), line, why_invalid, get_more_text_func, sub_call = .true.)
      if (.not. valid) return
      obj%n_obj = 1
      if (logic_option(.false., parse_one) .and. .not. has_opening_bracket) return
      parse_state = parsing_subobjects$
      ip = 0

    case ('[')
      parse_state = parsing_tokens$
      string_state = before_string$
      in_array = .true.
      if (.not. allocated(obj%token)) allocate (obj%token(10))
      nb1 = 0;  nb2 = 0

    case default
      parse_state = parsing_tokens$
      string_state = in_nonquote$
      in_array = .false.
      ip0 = ip
      nb1 = 0;  nb2 = 0
      if (ach == '(') nb2 = 1
    end select

  ! In an array of values.
  case (parsing_tokens$)
    if (string_state == before_string$) then
      if (ach == ' ') cycle
      ip0 = ip
      string_state = in_nonquote$
    endif

    select case (ach)
    case ('[');  nb1 = nb1 + 1
    case ('(');  nb2 = nb2 + 1
    case (']')
      if (in_array .and. nb1 == 0 .and. nb2 == 0) then
        obj%n_token = obj%n_token + 1
        n = obj%n_token
        call re_allocate (obj%token, n)
        ns = ip - ip0
        if (.not. allocated(obj%token(n)%str)) allocate (character(ns):: obj%token(n)%str)
        obj%token(n)%str = trim(line(ip0:ip-1))
        line = line(ip+1:)
        return
      else
        nb1 = nb1 - 1
        if (nb1 < 0) then
          valid = .false.
          why_invalid = 'MISPLACED "]" BRACKET.'
          return
        endif
      endif

    case (')')
      nb2 = nb2 - 1
      if (nb2 < 0) then
        valid = .false.
        why_invalid = 'MISPLACED ")" PARENTHESES.'
        return
      endif

    case (',')
      if (nb1 /= 0 .or. nb2 /= 0) cycle
      obj%n_token = obj%n_token + 1
      n = obj%n_token
      if (in_array) then
        if (n > size(obj%token)) call re_allocate (obj%token, 2*n)
      else
        call re_allocate (obj%token, 1)
      endif
      ns = ip - ip0
      if (.not. allocated(obj%token(1)%str)) allocate (character(ns):: obj%token(n)%str)
      obj%token(n)%str = trim(line(ip0:ip-1))
      if (.not. in_array) then
        line = line(ip:)
        return
      endif
      string_state = before_string$

    case ('}')
      if (nb1 /= 0 .or. nb2 /= 0) then
        valid = .false.
        why_invalid = 'MISPLACED "}" CHARACTER.'
        return
      endif
      call re_allocate(obj%token, 1)
      n = ip - ip0
      if (.not. allocated(obj%token(1)%str)) allocate(character(n):: obj%token(1)%str)
      obj%token(1)%str = line(ip0:ip-1)
      obj%n_token = 1
      line = line(ip:)
      if (deb) print *, 'Out: ', line
      return
    end select

  ! In subobject
  case (parsing_subobjects$)
    select case (ach)
    case (' ') 
      cycle
    case (',')
      obj%n_obj = obj%n_obj + 1
      line = line(ip+1:)
      if (obj%n_obj > size(obj%obj)) call subobject_reallocate(obj,2*obj%n_obj)
      valid = object_document_parse(obj%obj(obj%n_obj), line, why_invalid, get_more_text_func, sub_call = .true.)
      if (.not. valid) return
      ip = 0
    case ('}')
      call subobject_reallocate(obj,obj%n_obj)
      line = line(ip+1:)
      return
    case default
      valid = .false.
      why_invalid = 'BAD CHARACTER: ' // ach
      return
    end select

  end select
enddo

end function object_document_parse

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine object_print(obj, n_indent, last)
!
! Routine to print the contents of an object.
!
! Input:
!   obj       -- object_struct: Object to be printed.
!   n_indent  -- integer, optional: Indentation when printing.
!   last      -- logical, optional: Is last datum in array? If True then no
!-

recursive subroutine object_print (obj, n_indent, last)

type (object_struct), target :: obj
type (object_struct), pointer :: optr
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
    print '(3a)', blank(1:ni), trim(obj%name), ' = {'
  endif

  do i = 1, obj%n_obj
    call object_print(obj%obj(i), ni+2, last = (i == obj%n_obj))
  enddo

  if (logic_option(.true., last)) then
    print '(2a)', blank(1:ni), '}'
  else
    print '(2a)', blank(1:ni), '},'
  endif

elseif (obj%n_token == 0) then
  print '(4a)', blank(1:ni), trim(obj%name), ' = ????'
  
elseif (obj%n_token == 1) then
  if (logic_option(.true., last)) then
    print '(4a)', blank(1:ni), trim(obj%name), ' = ', trim(obj%token(1)%str)
  else
    print '(5a)', blank(1:ni), trim(obj%name), ' = ', trim(obj%token(1)%str), ','
  endif

else
  allocate (character(40):: line)
  line = blank(1:ni) // trim(obj%name) // ' = [' // trim(obj%token(1)%str)
  do i = 2, obj%n_token
    line = line // ', ' // trim(obj%token(i)%str)
  enddo
  if (logic_option(.true., last)) then
    print '(a)', line // ']'
  else
    print '(a)', line // '],'
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
!   obj   -- object_struct: Object.
!   n     -- integer: Minimum size for obj%obj(:)
!-

subroutine subobject_reallocate(obj, n)

type (object_struct) obj, temp
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
!   obj   -- object_struct: Object.
!
! Output:
!   obj   -- object_struct: Deallocated object.
!-

recursive subroutine object_deallocate (obj)

type (object_struct) obj
integer i

!

if (obj%n_obj > 0) then
  do i = 1, obj%n_obj
    call object_deallocate (obj%obj(i))
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
!   obj     -- object_struct: Object.
!   name    -- character(*): Name of subobject
!
! Output:
!   num     -- integer: Number of subobjects whose name matches name argument.
!-

function subobject_number (obj, name) result (num)

type (object_struct) obj, temp
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
!   obj     -- object_struct: Object.
!   name    -- character(*): Name of subobject
!   indx    -- integer, optional: If present, start subobject search from this index + 1.
!                 Default is 0.
!
! Output:
!   indx    -- integer, optional: Index of subobject found.
!   subobj  -- object_struct, pointer: Pointer to subobject. Null if not found.
!-

function subobject_ptr (obj, name, indx) result (subobj)  

type (object_struct), target :: obj
type (object_struct), pointer :: subobj
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
!   obj         -- object_struct: Object.
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

type (object_struct), target :: obj
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
!   obj         -- object_struct: Object.
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

type (object_struct), target :: obj
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
!   obj         -- object_struct: Object.
!   name        -- character(*): Name of subobject
!
! Output:
!   real_arr(:)   -- real(rp), allocatable, optional: Real value.
!   int_arr(:)    -- integer, allocatable, optional: Integer value.
!   logic_arr(:)  -- logical, allocatable, optional: Logical value. 
!   valid         -- logical: Valid value found.
!-

function object_alloc_array_read (obj, real_arr, int_arr, logic_arr) result (valid)

type (object_struct), target :: obj
type (object_struct), pointer :: subobj
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





