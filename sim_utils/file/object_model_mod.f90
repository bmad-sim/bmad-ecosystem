!+
! Module for parsing a document with a tree format of key/value pairs.
!
! A key/value pair is called an "object".
! The document model is a tree of objects. 
! The "value" of an object will either be an array of one or more tokens or 
! an array of of one or more subobjects (children) (but not both).
!
! Example:
!   {a = {b = 7*atan(3, 4), d = "xyz", e = [1, "true", rats], f = {g = 3^34}}}
! In this example, the root object, which has a blank key name, has a single child subobject
! whose (key) name is "a". Object "a" has three children whose names are "b", "d", e", and "f".
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

! One and only one of %n_child and %n_token is nonzero.

type temp_struct
  character(40) :: str = ''
end type

type object_struct
  character(:), allocatable :: name
  type (object_struct), pointer :: parent => null()
  type (var_length_string_struct), allocatable :: token(:) ! Array of tokens
  type (object_struct), pointer :: child(:) => null() ! Array of subobjects
  integer :: n_child = 0        ! Number of subobjects.
  integer :: n_token = 0      ! Number of tokens.
end type

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function object_document_parse(obj, root_name, line, why_invalid, get_more_text_func, parse_one, sub_call) result (valid)
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
!   root_name     -- character(*): Name to use of the root object obj.
!   line          -- character(:), allocatable: Initial text to parse. If no text present then
!                       get_more_text_func will be called. 
!   why_invalid   -- character(*): String set to error message if there is an error.
!   get_more_text_func(line) result (valid)
!                 -- Function: Routine called to append more text to the string being parsed. 
!                     The interface is:
!                         function get_more_text_func(line, end_of_document, why_invalid) result (valid)
!                           character(:), allocatable :: line, why_invalid
!                           logical end_of_document
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

recursive function object_document_parse(obj, root_name, line, why_invalid, get_more_text_func, parse_one, sub_call) result (valid)

interface 
  function get_more_text_func(line, end_of_document, why_invalid) result (valid)
    character(:), allocatable :: line, why_invalid
    logical end_of_document
    logical valid
  end function
end interface

type (object_struct), target :: obj

integer n, ns, ip, ip0, parse_state, nb1, nb2, n_item
integer, parameter :: at_start$ = 0, before_obj_name$ = 1, in_obj_name$ = 2, in_obj_quoted_name$ = 3
integer, parameter :: after_obj_name$ = 4, after_equal_sign$ = 5
integer, parameter :: in_token_before_string$ = 10, in_token_in_string$ = 11
integer, parameter :: in_token_in_quoted_string$ = 12, in_token_after_string$ = 13
integer, parameter :: in_subobjects$ = 20

logical, optional :: parse_one, sub_call
logical valid, abort, has_opening_bracket, in_array, end_of_document, debug

character(*) root_name
character(:), allocatable :: line, why_invalid, why_func_invalid
character(1), parameter :: tab = char(9)
character(1) ach, quote_mark

!

valid = .true.
debug = .false.
if (.not. present(sub_call)) then
  if (.not. allocated(obj%name)) allocate (character(len_trim(root_name)):: obj%name)
  obj%name = trim(root_name)
endif
has_opening_bracket = .false.
quote_mark = ' '

if (.not. allocated(line)) then
  allocate(character(1):: line)
  line = ''
  valid = get_more_text_func(line, end_of_document, why_func_invalid)
  if (.not. valid) then
    valid = set_this_invalid (why_invalid, why_func_invalid, obj)
    return
  endif
endif

ip = 0

if (present(sub_call)) then
  parse_state = before_obj_name$
else
  parse_state = at_start$
endif

do while (.true.)
  ip = ip + 1
  if (ip > len(line)) then
    valid = get_more_text_func(line, end_of_document, why_func_invalid)

    if (.not. valid) then
      valid = set_this_invalid (why_invalid, why_func_invalid, obj)
      return
    endif

    if (end_of_document) then
      if (.not. has_opening_bracket .and. .not. in_array .and. &
                  (parse_state == in_token_in_string$ .or. parse_state == in_token_after_string$)) then
        obj%n_token = obj%n_token + 1
        n = obj%n_token
        call re_allocate (obj%token, n)
        ns = ip - ip0 - 1
        if (.not. allocated(obj%token(n)%str)) allocate (character(ns):: obj%token(n)%str)
        obj%token(n)%str = trim(line(ip0:ip-1))
        line = ''
        return
      else
        valid = set_this_invalid (why_invalid, 'MISSING END "}" DELIMITOR.', obj)
        return
      endif
    endif
  endif

  ach = line(ip:ip)
  if (debug) print '(i4, 2x, 3a)', ip, ach, ' :::', line(ip+1:)

  ! Handle quote character

  if (ach == '"' .or. ach == "'") then

    select case (quote_mark)
    case (' ')
      quote_mark = ach
      ip0 = ip
      select case (parse_state)
      case (before_obj_name$);   parse_state = in_obj_quoted_name$
      case (after_equal_sign$, in_token_before_string$);  parse_state = in_token_in_quoted_string$
      case default
        valid = set_this_invalid (why_invalid, 'QUOTE MARK OUT OF PLACE', obj)
        return
      end select
      cycle
        
    case ('"', "'")
      if (line(ip-1:ip-1) == '\') cycle      ! '
      if (ach /= quote_mark) cycle
      select case (parse_state)
      case (in_obj_quoted_name$)
        obj%name = line(ip0:ip)
        parse_state = after_obj_name$
      case (in_token_in_quoted_string$)
        parse_state = in_token_after_string$
      case default
        call err_exit
      end select
      quote_mark = ' '
      cycle
    end select
    cycle
  endif

  if (quote_mark /= ' ') then
    if (ach == '=') then
      valid = set_this_invalid (why_invalid, 'EQUAL SIGN "=" NOT ALLOWED IN QUOTED STRING.', obj)
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
    valid = object_document_parse(obj%child(1), '', line, why_invalid, get_more_text_func, sub_call = .true.)
    if (.not. valid) return
    obj%n_child = 1
    if (logic_option(.false., parse_one) .and. .not. has_opening_bracket) return
    parse_state = in_subobjects$
    ip = 0

  ! Before a variable name
  case (before_obj_name$)
    if (ach == ' ') cycle
    select case (ach)
    case (',', '{', '}', '(', ')', '[', ']')
      valid = set_this_invalid (why_invalid, 'MISPLACED ' // ach // ' CHARACTER', obj)
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
        valid = set_this_invalid (why_invalid, 'MISPLACED "]" BRACKET.', obj)
        return
      endif

    case (')')
      nb2 = nb2 - 1
      if (nb2 < 0) then
        valid = set_this_invalid (why_invalid, 'MISPLACED ")" BRACKET.', obj)
        return
      endif

    case (',')
      if (nb1 == 0 .and. nb2 == 0) then
        valid = set_this_invalid (why_invalid, 'COMMA IN NAME', obj)
        return
      endif

    case ('=')
      if (nb1 /= 0 .or. nb2 /= 0) then
        valid = set_this_invalid (why_invalid, 'EQUAL SIGN IN NAME', obj)
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
      valid = set_this_invalid (why_invalid, 'EXPECTING EQUAL SIGN AFTER DATUM NAME', obj)
      return
    endif
    n = ip - ip0
    if (.not. allocated(obj%name)) allocate(character(n):: obj%name)
    obj%name = trim(line(ip0:ip-1))
    parse_state = after_equal_sign$

  ! After an equal sign.
  case (after_equal_sign$)
    select case (ach)
    case (' ')
      cycle

    case (',', ')', ']')
      valid = set_this_invalid (why_invalid, 'MISPLACED "' // ach // '" CHARACTER.', obj)
      return

    case ('{')
      line = line(ip+1:)
      call subobject_reallocate(obj, 5)
      valid = object_document_parse(obj%child(1), '', line, why_invalid, get_more_text_func, sub_call = .true.)
      if (.not. valid) return
      obj%n_child = 1
      if (logic_option(.false., parse_one) .and. .not. has_opening_bracket) return
      parse_state = in_subobjects$
      ip = 0

    case ('[')
      parse_state = in_token_before_string$
      in_array = .true.
      if (.not. allocated(obj%token)) allocate (obj%token(10))
      nb1 = 0;  nb2 = 0

    case default
      parse_state = in_token_before_string$
      in_array = .false.
      ip0 = ip
      nb1 = 0;  nb2 = 0
      ip = ip - 1  ! Reeval character
    end select

  ! Before next token
  case (in_token_before_string$)
    if (ach == ' ') cycle
    ip0 = ip
    ip = ip - 1  ! Reeval
    parse_state = in_token_in_string$

  ! In token string.
  case (in_token_in_string$)
    select case (ach)
    case ('[');  nb1 = nb1 + 1
    case ('(');  nb2 = nb2 + 1
    case (']')
      if (in_array .and. nb1 == 0 .and. nb2 == 0) then
        parse_state = in_token_after_string$
        ip = ip - 1
      else
        nb1 = nb1 - 1
        if (nb1 < 0) then
          valid = set_this_invalid (why_invalid, 'MISPLACED "]" BRACKET.', obj)
          return
        endif
      endif

    case (')')
      nb2 = nb2 - 1
      if (nb2 < 0) then
        valid = set_this_invalid (why_invalid, 'MISPLACED ")" PARENTHESES.', obj)
        return
      endif

    case (',')
      if (nb1 /= 0 .or. nb2 /= 0) cycle
      ip = ip - 1
      parse_state = in_token_after_string$

    case ('}')
      if (nb1 /= 0 .or. nb2 /= 0) then
        valid = set_this_invalid (why_invalid, 'MISPLACED "}" CHARACTER.', obj)
        return
      endif
      ip = ip - 1
      parse_state = in_token_after_string$

    case (' ')
      if (nb1 /= 0 .or. nb2 /= 0) cycle
      parse_state = in_token_after_string$
    end select

  ! After token string.
  case (in_token_after_string$)
    select case (ach)
    case (',', ']', '}')
      obj%n_token = obj%n_token + 1
      n = obj%n_token
      if (n == 1) call re_allocate(obj%token, 2)
      if (n > size(obj%token)) call re_allocate (obj%token, max(n, 2*(n-1)))
      ns = ip - ip0
      if (.not. allocated(obj%token(1)%str)) allocate (character(ns):: obj%token(n)%str)
      obj%token(n)%str = trim(line(ip0:ip-1))
      if (debug) print *, 'Out: ', line

      if (in_array .and. ach == '}') then
        valid = set_this_invalid (why_invalid, 'expecting "," or "]", got: "}"', obj)
        return
      elseif (.not. in_array .and. ach == ']') then
        valid = set_this_invalid (why_invalid, 'expecting "," or "}", got: "]"', obj)
        return
      endif

      if (ach == ']') then
        line = line(ip+1:)
        return
      elseif (ach == '}' .or. .not. in_array) then
        line = line(ip:)
        return
      else
        parse_state = in_token_before_string$
      endif

    case default
      if (in_array) then
        valid = set_this_invalid (why_invalid, 'expecting "," or "]" but got: ' // quote(ach), obj)
      else
        valid = set_this_invalid (why_invalid, 'expecting "," or "}" but got: ' // quote(ach), obj)
      endif
      return
    end select

  ! In subobject
  case (in_subobjects$)
    select case (ach)
    case (' ') 
      cycle
    case (',')
      obj%n_child = obj%n_child + 1
      line = line(ip+1:)
      if (obj%n_child == 1) call subobject_reallocate(obj,2)
      if (obj%n_child > size(obj%child)) call subobject_reallocate(obj,2*obj%n_child)
      valid = object_document_parse(obj%child(obj%n_child), '', line, why_invalid, get_more_text_func, sub_call = .true.)
      if (.not. valid) return
      ip = 0
    case ('}')
      call subobject_reallocate(obj,obj%n_child)
      line = line(ip+1:)
      return
    case default
      valid = set_this_invalid (why_invalid, 'BAD CHARACTER: ' // ach, obj)
      return
    end select

  end select
enddo

!--------------------------------------------------
contains

function set_this_invalid (why_invalid, err_string, obj) result (valid)

type (object_struct) obj
character(*) err_string
character(:), allocatable :: why_invalid
logical valid

! If err_string is blank (can happen if err_string is set by get_more_text_func), set why_invalid to blank.

valid = .false.
if (err_string == '') then
  why_invalid = err_string
else
  call str_set(why_invalid, err_string // '.' // ' WHILE PARSING: ' // object_tree_name(obj))
endif

end function set_this_invalid

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
character(100) :: blank = ''
character(:), allocatable :: line
integer, optional :: n_indent
integer i, ni
logical, optional :: last

! Head object may not have a name.
! Note: In this case, the object must have subobjects and not a value

ni = integer_option(0, n_indent)

if (obj%n_child > 0) then
  if (obj%name == '') then
    print '(2a)', blank(1:ni), '{'
  else
    print '(3a)', blank(1:ni), trim(obj%name), ' = {'
  endif

  do i = 1, obj%n_child
    call object_print(obj%child(i), ni+2, last = (i == obj%n_child))
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
! Routine to reallocate the obj%child(:) array.
!
! Input:
!   obj   -- object_struct: Object.
!   n     -- integer: Minimum size for obj%child(:)
!-

subroutine subobject_reallocate(obj, n)

type (object_struct), target :: obj, temp
integer n, ix

!

if (.not. associated(obj%child)) then
  allocate (obj%child(n))
  obj%n_child = 0
else
  if (size(obj%child) >= n) return
  temp%child => obj%child
  allocate(obj%child(n))
  obj%child(1:size(temp%child)) = temp%child
  deallocate(temp%child)
endif

do ix = 1, size(obj%child)
  obj%child(ix)%parent => obj
enddo

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

if (obj%n_child > 0) then
  do i = 1, obj%n_child
    call object_deallocate (obj%child(i))
  enddo
  deallocate(obj%child)
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
do n = 1, obj%n_child
  if (obj%child(n)%name == name) num = num + 1
enddo

end function subobject_number

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function pointer_to_subobject (obj, name, indx) result (subobj)  
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

function pointer_to_subobject (obj, name, indx) result (subobj)  

type (object_struct), target :: obj
type (object_struct), pointer :: subobj
integer, optional :: indx
integer n, n0
character(*) name

!

subobj => null()

n0 = integer_option(0, indx) + 1
do n = n0, obj%n_child
  if (obj%child(n)%name /= name) cycle
  subobj => obj%child(n)
  if (present(indx)) indx = n
  return
enddo

end function pointer_to_subobject

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function object_tree_name(obj) result (tree_name)
!
! Routine to return the tree name of an object which is the path from the root to
! the object. For example, given the object tree
!     {a = {b = 1, b = {c = 7, c = 8}}}
! The tree name for "c = 7" object is
!     "a->b[2]->c[1]"
!
! Input:
!   obj       -- object_struct: Object
!
! Output:
!   tree_name -- character(:), allocatable: Tree name
!-

function object_tree_name(obj) result (tree_name)

type (object_struct), target :: obj
type (object_struct), pointer :: op1, op2

integer i, n
character(:), allocatable :: tree_name

!

op1 => obj
allocate(character(0):: tree_name)

do
  op2 => op1%parent
  if (.not. associated(op2)) then
    tree_name = op1%name // '->' // tree_name
    exit
  endif

  if (subobject_number(op2, op1%name) == 1) then
    tree_name = op1%name // '->' // tree_name
  else
    n = 0
    do i = 1, op2%n_child
      if (op2%child(i)%name /= op1%name) cycle
      n = n + 1
      if (associated(op1, op2%child(i))) exit
    enddo
    tree_name = op1%name // '[' // int_str(n) // ']->' // tree_name
  endif

  op1 => op2
  if (.not. associated(op1%parent) .and. op1%name == '') exit
enddo

n = len(tree_name)
tree_name = tree_name(1:n-2)  ! Trim off trailing "->"

end function object_tree_name

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function object_root (obj) result (root_obj_ptr)
!
! Routine to point to the root object of the tree obj is contained in.
!
! Input:
!   obj           -- object_struct: Object.
!
! Output:
!   root_obj_ptr  -- object_struct, pointer: Pointer to the root object
!-

function object_root (obj) result (root_obj_ptr)

type (object_struct), target :: obj
type (object_struct), pointer :: root_obj_ptr

root_obj_ptr => obj
do
  if (.not. associated(root_obj_ptr%parent)) return
  root_obj_ptr => obj%parent
enddo

end function object_root

end module
