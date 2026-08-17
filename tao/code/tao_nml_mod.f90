!+
! Module tao_nml_mod
!
! Reader for Tao initialization files.
!
! This replaces the use of Fortran namelist reads. The file syntax accepted is the existing
! namelist syntax, so all existing Tao init files continue to work unchanged. What changes is
! that the reading is done here rather than by the Fortran runtime, which buys:
!   * Error messages that name the file, the line, and the offending token.
!   * Values that may be expressions (see tao_nml_value_set).
!   * Array components that are sized from what is actually in the file rather than being
!     declared at some generous fixed maximum.
!   * No dependence on compiler specific namelist behavior.
!
! The parsing here is deliberately generic and the binding of a name to a variable is left to
! the caller, which does a select case on the item name. That is the same information the
! "namelist /group/ a, b, c" statement used to carry, and keeping it at the call site means the
! per group quirks (deprecated names, defaults that carry over between groups, etc.) stay
! visible instead of being hidden in a general mechanism.
!
! Three deliberate differences from a namelist read, all of which only ever accept more:
!   * A value that is not a plain number may be an expression. See tao_nml_set_ptr.
!   * A group with no terminating "/" is ended by end of file rather than discarded.
!   * A blank value is an error rather than being silently ignored. See tao_set_ptr_value.
!
! Syntax accepted, otherwise matching what gfortran's namelist read accepts:
!   &group_name  or  $group_name          Start of a group.
!   /  or  &end  or  $end                 End of a group.
!   name = value                          An item. Value runs to the next item or the end.
!   name(3) = value                       Subscripted item.
!   name%comp = value                     Structure component.
!   name(3)%comp = value                  Both.
!   a = 1, 2, 3                           Multiple values.
!   a = 3*1.5                             Repeat count.
!   'text'  or  "text"                    Quoted string. Double the quote mark to embed one.
!   ! comment                             Comment to end of line. A gfortran extension that
!                                         Tao init files rely on heavily.
!
! Typical use:
!   do
!     call tao_nml_group_read (iu, file_name, 'tao_d2_data', group, eof, err)
!     if (err .or. eof) exit
!     do i = 1, group%n_item
!       select case (group%item(i)%name)
!       case ('n_d1_data');  call tao_nml_value_set (group%item(i), n_d1_data, err)
!       ...
!-

module tao_nml_mod

use tao_attrib_ptr_mod
use tao_interface, only: tao_evaluate_expression

implicit none

integer, parameter :: tao_nml_name_len$ = 200

! Line number bookkeeping. A group read picks up where the previous one left off, so the line
! counter has to persist between calls to give a line number that matches the file.
!
! The counter is keyed on both the unit number and the file name. Keying on the unit number
! alone is not enough because a newunit= open can hand back a number that a previously closed
! file was using, which would make the new file inherit the old count. Call tao_nml_rewind
! rather than a bare rewind so that the count is reset when a unit is rewound.

integer, private :: nml_iu_last = -1
integer, private :: nml_line_count = 0
character(400), private :: nml_file_last = ''

! One "name = value" item from a group.

type tao_nml_item_struct
  character(tao_nml_name_len$) :: name = ''   ! Item name, down cased. Eg "datum(3)%ele_name".
  character(:), allocatable :: value          ! Raw value text as it appeared in the file.
  character(:), allocatable :: file           ! File the item came from.
  integer :: i_line = 0                       ! Line number of the item. For error messages.
end type

! An item name split into its parts, for dispatching at a call site.
! Eg for "datum(3)%ele_name": head = "datum", has_sub = True, isub = 3, rest = "ele_name".

type tao_nml_ref_struct
  character(:), allocatable :: head
  character(:), allocatable :: rest
  integer :: isub = 0        ! Subscript, or first subscript of a range.
  integer :: isub2 = 0       ! Last subscript of a range. Equal to isub if not a range.
  logical :: has_sub = .false.
end type

! One "&group ... /" group.

type tao_nml_group_struct
  character(60) :: name = ''                  ! Group name, down cased.
  type (tao_nml_item_struct), allocatable :: item(:)
  integer :: n_item = 0
  integer :: i_line = 0                       ! Line number the group starts on.
end type

private nml_find_group_start, nml_gather_group_text, nml_split_items, nml_is_name_char

! Set a variable from an item. The variable may be an intrinsic scalar or rank 1 array, or a
! component reached through a resolver in tao_attrib_resolve_mod. Errors are reported through
! out_io with the file and line, and err is set.

interface tao_nml_value_set
  module procedure tao_nml_set_real, tao_nml_set_int, tao_nml_set_logic, tao_nml_set_str
  module procedure tao_nml_set_real_arr, tao_nml_set_int_arr, tao_nml_set_logic_arr
  module procedure tao_nml_set_str_arr, tao_nml_set_ptr
end interface

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_nml_group_read (iu, file_name, want_name, group, eof, err, why)
!
! Read the next group from an open init file.
!
! If want_name is not blank then groups whose name does not match are skipped. This mirrors
! what a namelist read does: "read (iu, nml = tao_d2_data)" ignores intervening groups of a
! different name.
!
! Input:
!   iu        -- integer: Open unit number.
!   file_name -- character(*): Name of the file. Used only for error messages.
!   want_name -- character(*): Group name wanted, or '' to take the next group whatever it is.
!
! Output:
!   group     -- tao_nml_group_struct: The group read.
!   eof       -- logical: True if end of file was reached without finding a group.
!   err       -- logical: True on a syntax error.
!   why       -- character(*), optional: Explanation if err is True. Already includes the file
!                  name and line number.
!-

subroutine tao_nml_group_read (iu, file_name, want_name, group, eof, err, why)

type (tao_nml_group_struct), intent(out) :: group

character(*), intent(in) :: file_name, want_name
character(*), optional, intent(out) :: why
character(:), allocatable :: text, why_str

integer, intent(in) :: iu
integer, allocatable :: line_of(:)
integer i_line

logical, intent(out) :: eof, err

!

err = .false.
eof = .false.
why_str = ''

! A different unit or a different file means start counting again.

if (iu /= nml_iu_last .or. file_name /= nml_file_last) then
  nml_iu_last = iu
  nml_file_last = file_name
  nml_line_count = 0
endif

i_line = nml_line_count

do
  call nml_find_group_start (iu, file_name, group%name, i_line, eof, err, why_str)
  if (err .or. eof) exit
  group%i_line = i_line

  call nml_gather_group_text (iu, file_name, group%name, i_line, text, line_of, err, why_str)
  nml_line_count = i_line
  if (err) exit

  if (want_name == '' .or. group%name == downcase(want_name)) then
    call nml_split_items (text, line_of, file_name, group, err, why_str)
    exit
  endif
enddo

nml_line_count = i_line

if (err .and. present(why)) why = why_str
if (present(why) .and. .not. err) why = ''

end subroutine tao_nml_group_read

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_nml_rewind (iu)
!
! Rewind a unit and reset the line counter that tao_nml_group_read keeps for error messages.
! Use this instead of a bare rewind so that reported line numbers stay correct.
!-

subroutine tao_nml_rewind (iu)

integer, intent(in) :: iu

!

rewind (iu)
nml_iu_last = iu
nml_line_count = 0
nml_file_last = ''

end subroutine tao_nml_rewind

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine nml_find_group_start (iu, file_name, name, i_line, eof, err, why)
!
! Scan forward to the next line starting a group and return the group name.
!-

subroutine nml_find_group_start (iu, file_name, name, i_line, eof, err, why)

character(*), intent(in) :: file_name
character(*), intent(out) :: name
character(:), allocatable, intent(inout) :: why
character(1000) line
character(:), allocatable :: str

integer, intent(in) :: iu
integer, intent(inout) :: i_line
integer ios, ix

logical, intent(out) :: eof, err

!

err = .false.
eof = .false.
name = ''

do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) then
    eof = .true.
    return
  endif
  i_line = i_line + 1
  str = trim(adjustl(nml_strip_comment(line)))
  if (str == '') cycle
  if (str(1:1) /= '&' .and. str(1:1) /= '$') cycle

  str = adjustl(str(2:))
  ix = 1
  do while (ix <= len(str))
    if (.not. nml_is_name_char(str(ix:ix))) exit
    ix = ix + 1
  enddo

  if (ix == 1) then
    err = .true.
    why = nml_err_str(file_name, i_line, 'GROUP NAME MISSING AFTER "' // line(1:1) // '"')
    return
  endif

  name = downcase(str(1:ix-1))
  ! Anything after the group name on the same line is put back by the caller via the
  ! gather step, which re-reads from the current position. Tao init files always put the
  ! group name on a line of its own, but handle the general case anyway.
  backspace (iu)
  i_line = i_line - 1
  return
enddo

end subroutine nml_find_group_start

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine nml_gather_group_text (iu, file_name, name, i_line, text, line_of, err, why)
!
! Accumulate the text of a group, from just after the "&name" to the terminator, into a single
! string. line_of(i) is the line number that text(i:i) came from, for error messages.
!
! Quote state is carried across lines. A "/" outside of quotes terminates the group, which is
! the namelist rule and is why an unquoted file path in a namelist has always been a mistake.
!-

subroutine nml_gather_group_text (iu, file_name, name, i_line, text, line_of, err, why)

character(*), intent(in) :: file_name, name
character(:), allocatable, intent(out) :: text
character(:), allocatable, intent(inout) :: why
character(1000) line
character(:), allocatable :: str
character(1) quote_ch

integer, intent(in) :: iu
integer, intent(inout) :: i_line
integer, allocatable, intent(out) :: line_of(:)
integer ios, i, ix, n, n_alloc, depth

logical, intent(out) :: err
logical first_line, in_quote

!

err = .false.
text = ''
n = 0
n_alloc = 1024
allocate (line_of(n_alloc))
in_quote = .false.
quote_ch = ' '
depth = 0
first_line = .true.

do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) then
    ! End of file terminates the group rather than being an error.
    !
    ! This is a deliberate difference from a namelist read, not a match. A namelist read of a
    ! group with no terminating "/" aborts, or with iostat= reports end of file and throws away
    ! whatever it had read. Files like that do occur: the "write namelist -variable" command
    ! crashes partway through on some lattices, a pre-existing bug unrelated to the reader, and
    ! leaves a truncated file behind. Using what was read beats discarding it silently.
    !
    ! No well formed file is affected, since a well formed file terminates its groups.
    call nml_shrink (text, line_of, n)
    return
  endif
  i_line = i_line + 1

  str = nml_strip_comment(line)

  ! Skip past the "&name" on the first line.

  if (first_line) then
    first_line = .false.
    str = adjustl(str)
    str = adjustl(str(2:))          ! Drop the & or $
    ix = 1
    do while (ix <= len(str))
      if (.not. nml_is_name_char(str(ix:ix))) exit
      ix = ix + 1
    enddo
    str = str(ix:)
  endif

  ! Append, watching for the terminator.

  i = 1
  do while (i <= len(str))
    if (in_quote) then
      if (str(i:i) == quote_ch) then
        ! A doubled quote mark is an embedded quote, not the end of the string.
        if (i < len(str)) then
          if (str(i+1:i+1) == quote_ch) then
            call nml_append (text, line_of, n, n_alloc, str(i:i), i_line)
            call nml_append (text, line_of, n, n_alloc, str(i+1:i+1), i_line)
            i = i + 2
            cycle
          endif
        endif
        in_quote = .false.
      endif

    else
      select case (str(i:i))
      case ("'", '"')
        in_quote = .true.
        quote_ch = str(i:i)
      case ('(', '[')
        depth = depth + 1
      case (')', ']')
        depth = depth - 1
      case ('/')
        if (depth == 0) then
          call nml_shrink (text, line_of, n)
          return
        endif
      case ('&', '$')
        ! "&end" or "$end" terminator.
        if (depth == 0) then
          if (len(str) >= i+3) then
            if (downcase(str(i+1:i+3)) == 'end') then
              call nml_shrink (text, line_of, n)
              return
            endif
          endif
        endif
      end select
    endif

    call nml_append (text, line_of, n, n_alloc, str(i:i), i_line)
    i = i + 1
  enddo

  ! A line break acts as a value separator, same as a comma.

  call nml_append (text, line_of, n, n_alloc, ' ', i_line)
enddo

end subroutine nml_gather_group_text

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine nml_split_items (text, line_of, file_name, group, err, why)
!
! Split the gathered group text into name = value items.
!
! An item boundary is an "=" at top level (not inside quotes or parentheses) preceded by a
! name. So the value of an item runs from its "=" up to the name of the next item, which is
! how a namelist read decides where one value list ends and the next item begins.
!-

subroutine nml_split_items (text, line_of, file_name, group, err, why)

type (tao_nml_group_struct), intent(inout) :: group
type (tao_nml_item_struct), allocatable :: item(:)

character(*), intent(in) :: text, file_name
character(:), allocatable, intent(inout) :: why
character(:), allocatable :: nam, val

integer, intent(in) :: line_of(:)
integer ix_eq, ix_name, ix_next_eq, ix_next_name, i, n

logical, intent(out) :: err

!

err = .false.
group%n_item = 0
allocate (item(64))

ix_eq = nml_next_top_level_eq (text, 1)

do
  if (ix_eq == 0) exit

  ! The name ends just before the "=". Walk back over blanks then over name characters.

  i = ix_eq - 1
  do while (i >= 1)
    if (text(i:i) /= ' ') exit
    i = i - 1
  enddo
  ix_name = i
  do while (ix_name >= 1)
    if (.not. nml_is_name_char(text(ix_name:ix_name)) .and. text(ix_name:ix_name) /= '%' .and. &
        text(ix_name:ix_name) /= '(' .and. text(ix_name:ix_name) /= ')' .and. &
        text(ix_name:ix_name) /= ':' .and. text(ix_name:ix_name) /= ',') exit
    ix_name = ix_name - 1
  enddo
  ix_name = ix_name + 1

  if (ix_name > i) then
    err = .true.
    why = nml_err_str(file_name, line_of(ix_eq), 'NO NAME BEFORE "=" IN GROUP "' // &
                                                                    trim(group%name) // '"')
    return
  endif

  nam = downcase(text(ix_name:i))

  ! The value runs to just before the name of the next item.

  ix_next_eq = nml_next_top_level_eq (text, ix_eq + 1)

  if (ix_next_eq == 0) then
    val = text(ix_eq+1:)
  else
    i = ix_next_eq - 1
    do while (i >= 1)
      if (text(i:i) /= ' ') exit
      i = i - 1
    enddo
    ix_next_name = i
    do while (ix_next_name >= 1)
      if (.not. nml_is_name_char(text(ix_next_name:ix_next_name)) .and. &
          text(ix_next_name:ix_next_name) /= '%' .and. &
          text(ix_next_name:ix_next_name) /= '(' .and. &
          text(ix_next_name:ix_next_name) /= ')' .and. &
          text(ix_next_name:ix_next_name) /= ':' .and. &
          text(ix_next_name:ix_next_name) /= ',') exit
      ix_next_name = ix_next_name - 1
    enddo
    val = text(ix_eq+1:ix_next_name)
  endif

  ! Strip the trailing comma that separated this item from the next.

  val = trim(adjustl(val))
  n = len(val)
  do while (n > 0)
    if (val(n:n) /= ',' .and. val(n:n) /= ' ') exit
    n = n - 1
  enddo
  val = val(1:n)

  group%n_item = group%n_item + 1
  if (group%n_item > size(item)) call nml_reallocate_items (item, 2*size(item))
  item(group%n_item)%name = nam
  item(group%n_item)%value = val
  item(group%n_item)%file = file_name
  item(group%n_item)%i_line = line_of(min(ix_eq, size(line_of)))

  ix_eq = ix_next_eq
enddo

call nml_reallocate_items (item, max(group%n_item, 1))
call move_alloc (item, group%item)

end subroutine nml_split_items

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function nml_next_top_level_eq (text, ix_start) result (ix)
!
! Index of the next "=" that is not inside quotes or parentheses. Returns 0 if there is none.
!-

function nml_next_top_level_eq (text, ix_start) result (ix)

character(*), intent(in) :: text
character(1) quote_ch

integer, intent(in) :: ix_start
integer ix, i, depth

logical in_quote

!

in_quote = .false.
quote_ch = ' '
depth = 0

do i = ix_start, len(text)
  if (in_quote) then
    if (text(i:i) == quote_ch) in_quote = .false.
    cycle
  endif

  select case (text(i:i))
  case ("'", '"')
    in_quote = .true.
    quote_ch = text(i:i)
  case ('(', '[')
    depth = depth + 1
  case (')', ']')
    depth = depth - 1
  case ('=')
    if (depth == 0) then
      ix = i
      return
    endif
  end select
enddo

ix = 0

end function nml_next_top_level_eq

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function nml_strip_comment (line) result (str)
!
! Remove a trailing "!" comment, respecting quoted strings.
!
! Note: comments in namelist input are a gfortran extension rather than standard Fortran, but
! Tao init files use them everywhere so they must be supported.
!-

function nml_strip_comment (line) result (str)

character(*), intent(in) :: line
character(len(line)) str
character(1) quote_ch
character(1), parameter :: tab = achar(9)

integer i

logical in_quote

!

str = line
in_quote = .false.
quote_ch = ' '

do i = 1, len_trim(line)
  if (in_quote) then
    if (line(i:i) == quote_ch) in_quote = .false.
    cycle
  endif

  select case (line(i:i))
  case ("'", '"')
    in_quote = .true.
    quote_ch = line(i:i)
  case (tab)
    ! A tab outside of a quoted string is whitespace. Tao init files use tabs to line up
    ! value lists, and a namelist read treats a tab as a blank, so do the same here. Tabs
    ! inside a quoted string are left alone since they are part of the string.
    str(i:i) = ' '
  case ('!')
    ! Blank out from the comment on. Note: str, not line, since tabs converted above are in
    ! str only and reassigning from line here would throw those conversions away.
    str(i:) = ''
    return
  end select
enddo

end function nml_strip_comment

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function nml_is_name_char (ch) result (is_name)
!
! True if ch may appear in a name.
!-

function nml_is_name_char (ch) result (is_name)

character(1), intent(in) :: ch
logical is_name

!

is_name = ((ch >= 'a' .and. ch <= 'z') .or. (ch >= 'A' .and. ch <= 'Z') .or. &
           (ch >= '0' .and. ch <= '9') .or. ch == '_')

end function nml_is_name_char

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function nml_err_str (file_name, i_line, msg) result (str)
!
! Compose an error message that names the file and line.
!-

function nml_err_str (file_name, i_line, msg) result (str)

character(*), intent(in) :: file_name, msg
character(:), allocatable :: str

integer, intent(in) :: i_line

!

str = trim(msg) // '  [FILE: ' // trim(file_name) // ', LINE: ' // int_str(i_line) // ']'

end function nml_err_str

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Internal helpers for growing the gathered text and the item array.
!--------------------------------------------------------------------------

subroutine nml_append (text, line_of, n, n_alloc, ch, i_line)

character(:), allocatable :: text
character(1), intent(in) :: ch
integer, allocatable :: line_of(:)
integer, allocatable :: tmp(:)
integer n, n_alloc, i_line

!

n = n + 1
if (n > n_alloc) then
  n_alloc = 2 * n_alloc
  allocate (tmp(n_alloc))
  tmp(1:n-1) = line_of(1:n-1)
  call move_alloc (tmp, line_of)
endif
text = text // ch
line_of(n) = i_line

end subroutine nml_append

!--------------------------------------------------------------------------

subroutine nml_shrink (text, line_of, n)

character(:), allocatable :: text
integer, allocatable :: line_of(:)
integer, allocatable :: tmp(:)
integer n

!

if (n < 1) then
  text = ''
  allocate (tmp(1))
  tmp(1) = 1
  call move_alloc (tmp, line_of)
  return
endif

allocate (tmp(n))
tmp(1:n) = line_of(1:n)
call move_alloc (tmp, line_of)

end subroutine nml_shrink

!--------------------------------------------------------------------------

subroutine nml_reallocate_items (item, n)

type (tao_nml_item_struct), allocatable :: item(:), tmp(:)
integer n, n_copy

!

allocate (tmp(n))
n_copy = min(n, size(item))
tmp(1:n_copy) = item(1:n_copy)
call move_alloc (tmp, item)

end subroutine nml_reallocate_items

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_nml_item_split (item, head, has_sub, isub, rest, err)
!
! Split an item name into its leading component, subscript, and remainder, so that the caller
! can dispatch on the leading name. See tao_split_attrib_name.
!
! Errors are reported here with the file and line.
!-

subroutine tao_nml_item_split (item, head, has_sub, isub, rest, err)

type (tao_nml_item_struct), intent(in) :: item

character(:), allocatable, intent(out) :: head, rest
character(200) why
character(*), parameter :: r_name = 'tao_nml_item_split'

integer, intent(out) :: isub

logical, intent(out) :: has_sub, err

!

call tao_split_attrib_name (item%name, head, has_sub, isub, rest, err, why)
if (err) call tao_nml_err (item, why)

end subroutine tao_nml_item_split

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_nml_item_ref (item, ref, err)
!
! Same as tao_nml_item_split but collecting the parts into a tao_nml_ref_struct, which keeps
! the declarations at a call site down to one variable.
!-

subroutine tao_nml_item_ref (item, ref, err)

type (tao_nml_item_struct), intent(in) :: item
type (tao_nml_ref_struct), intent(out) :: ref

character(200) why

logical, intent(out) :: err

!

call tao_split_attrib_name (item%name, ref%head, ref%has_sub, ref%isub, ref%rest, err, why, ref%isub2)
if (err) call tao_nml_err (item, why)

end subroutine tao_nml_item_ref

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_nml_max_subscript (group, name) result (n_max)
!
! Largest subscript used for the named item anywhere in a group. 0 if the name never appears
! with a subscript.
!
! This is what lets an input array be sized from what is actually in the file. Under a namelist
! read there was no way to ask how much was coming, so the arrays these groups read into had to
! be declared at some generous fixed maximum, which is where var(-100:5000) and the 4000 element
! datum guess came from.
!
! An open ended range such as "var(3:)" is counted as reaching from its first subscript through
! as many values as the item supplies, which is how many elements a namelist read would set.
!
! Input:
!   group -- tao_nml_group_struct: The group.
!   name  -- character(*): Item name, without any subscript. Eg "datum".
!
! Output:
!   n_max -- integer: Largest subscript used.
!-

function tao_nml_max_subscript (group, name) result (n_max)

type (tao_nml_group_struct), intent(in) :: group
type (tao_nml_ref_struct) ref

character(*), intent(in) :: name
character(200) why
character(:), allocatable :: val(:)

integer n_max, i, n_val, i_hi

logical err

!

n_max = 0

do i = 1, group%n_item
  call tao_split_attrib_name (group%item(i)%name, ref%head, ref%has_sub, ref%isub, ref%rest, &
                                                                        err, why, ref%isub2)
  if (err) cycle
  if (ref%head /= downcase(name)) cycle
  if (.not. ref%has_sub) cycle

  i_hi = ref%isub2

  if (i_hi == int_garbage$) then
    ! Open ended: as far as the values reach.
    call tao_nml_split_values (group%item(i)%value, val, n_val, err, why)
    if (err) cycle
    if (ref%isub == int_garbage$) cycle          ! "(:)" says nothing about the size.
    i_hi = ref%isub + max(n_val, 1) - 1
  endif

  n_max = max(n_max, i_hi)
enddo

end function tao_nml_max_subscript

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_nml_int_item (group, name, dflt) result (n)
!
! Value of a plain integer item in a group, or dflt if the item is not present or does not
! decode. Used to size an array when the group also states a bound explicitly, as tao_var does
! with ix_max_var.
!-

function tao_nml_int_item (group, name, dflt) result (n)

type (tao_nml_group_struct), intent(in) :: group

character(*), intent(in) :: name

integer, intent(in) :: dflt
integer n, i, ios

!

n = dflt

do i = 1, group%n_item
  if (group%item(i)%name /= downcase(name)) cycle
  read (group%item(i)%value, *, iostat = ios) n
  if (ios /= 0) n = dflt
  return
enddo

end function tao_nml_int_item

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_nml_ref_bounds (item, ref, lo, hi, i1, i2, err)
!
! Resolve a possibly open ended subscript range against the bounds of the array referenced.
!
! An open ended bound such as the one in "var(1:)" comes back from tao_split_attrib_name as
! int_garbage$ because the splitter does not know the shape of what is being referenced. This
! fills it in and checks that the result is within bounds.
!
! Input:
!   item  -- tao_nml_item_struct: The item. Used for the error message.
!   ref   -- tao_nml_ref_struct: The split item name.
!   lo    -- integer: Lower bound of the array.
!   hi    -- integer: Upper bound of the array.
!
! Output:
!   i1    -- integer: First subscript.
!   i2    -- integer: Last subscript.
!   err   -- logical: Set True, with a message printed, if out of bounds.
!-

subroutine tao_nml_ref_bounds (item, ref, lo, hi, i1, i2, err)

type (tao_nml_item_struct), intent(in) :: item
type (tao_nml_ref_struct), intent(in) :: ref

integer, intent(in) :: lo, hi
integer, intent(out) :: i1, i2

logical, intent(out) :: err

!

err = .false.

i1 = ref%isub
i2 = ref%isub2
if (i1 == int_garbage$) i1 = lo
if (i2 == int_garbage$) i2 = hi

if (.not. ref%has_sub) then
  err = .true.
  call tao_nml_err (item, 'COMPONENT IS AN ARRAY SO A SUBSCRIPT IS NEEDED: ' // trim(ref%head))
  return
endif

if (i1 < lo .or. i2 > hi .or. i2 < i1) then
  err = .true.
  call tao_nml_err (item, 'SUBSCRIPT OUT OF THE RANGE [' // int_str(lo) // ', ' // &
                                int_str(hi) // '] FOR: ' // trim(ref%head))
endif

end subroutine tao_nml_ref_bounds

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_nml_split_values (str, val, n_val, err, why)
!
! Split a namelist value list into individual values.
!
! This is needed for a whole structure assignment such as
!     ele_shape(3) = "Quadrupole::*", "box", "blue", 0.2
! where the values are assigned to the components of the structure in declaration order. The
! generated tao_res_*_slot routines give the pointer for each position.
!
! Values are separated by commas and/or blanks. Two commas in a row give a null value, which a
! namelist read takes to mean "leave this component alone", and which is returned here as a
! blank string. A repeat count "n*value" expands to n copies, and a bare "n*" expands to n
! null values.
!
! Input:
!   str    -- character(*): The value list.
!
! Output:
!   val    -- character(*), allocatable: The individual values. Quote marks are kept so that
!               the setter can tell a quoted string from a bare word.
!   n_val  -- integer: Number of values.
!   err    -- logical: Set True on a malformed list.
!   why    -- character(*): Explanation if err is True.
!-

subroutine tao_nml_split_values (str, val, n_val, err, why)

character(*), intent(in) :: str
character(:), allocatable, intent(out) :: val(:)
character(*), intent(out) :: why
character(:), allocatable :: tok
character(1) quote_ch

integer, intent(out) :: n_val
integer i, j, n, i_rep, ix_star

logical, intent(out) :: err

!

err = .false.
why = ''
n_val = 0
allocate (character(max(len(str), 1)) :: val(max(32, len(str)/2 + 4)))

i = 1

do
  ! Skip blanks.
  do while (i <= len(str))
    if (str(i:i) /= ' ') exit
    i = i + 1
  enddo
  if (i > len(str)) exit

  ! A comma here means a null value.

  if (str(i:i) == ',') then
    call add_val ('')
    i = i + 1
    cycle
  endif

  ! Read one token.

  if (str(i:i) == '"' .or. str(i:i) == "'") then
    quote_ch = str(i:i)
    j = i + 1
    do
      if (j > len(str)) then
        err = .true.
        why = 'UNBALANCED QUOTE MARK IN VALUE: ' // trim(str)
        return
      endif
      if (str(j:j) == quote_ch) then
        if (j < len(str)) then
          if (str(j+1:j+1) == quote_ch) then
            j = j + 2
            cycle
          endif
        endif
        exit
      endif
      j = j + 1
    enddo
    tok = str(i:j)
    i = j + 1

  else
    j = i
    do while (j <= len(str))
      if (str(j:j) == ',' .or. str(j:j) == ' ') exit
      j = j + 1
    enddo
    tok = str(i:j-1)
    i = j
  endif

  ! A leading "n*" is a repeat count.

  i_rep = 1
  ix_star = index(tok, '*')
  if (ix_star > 1 .and. tok(1:1) /= '"' .and. tok(1:1) /= "'") then
    if (is_integer(tok(1:ix_star-1), n)) then
      if (n < 0) then
        err = .true.
        why = 'NEGATIVE REPEAT COUNT IN VALUE: ' // trim(tok)
        return
      endif
      i_rep = n
      tok = tok(ix_star+1:)
    endif
  endif

  do n = 1, i_rep
    call add_val (tok)
  enddo

  ! Consume one separating comma if present.

  do while (i <= len(str))
    if (str(i:i) /= ' ') exit
    i = i + 1
  enddo
  if (i <= len(str)) then
    if (str(i:i) == ',') i = i + 1
  endif
enddo

!--------------------------------------------------------------------------
contains

subroutine add_val (v)
character(*) v
character(len(val)), allocatable :: tmp(:)
n_val = n_val + 1
if (n_val > size(val)) then
  allocate (tmp(2*size(val)))
  tmp(1:n_val-1) = val(1:n_val-1)
  call move_alloc (tmp, val)
endif
val(n_val) = v
end subroutine add_val

end subroutine tao_nml_split_values

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_nml_err (item, msg)
!
! Report an error for an item, naming the file and line.
!-

subroutine tao_nml_err (item, msg)

type (tao_nml_item_struct), intent(in) :: item
character(*), intent(in) :: msg
character(*), parameter :: r_name = 'tao_nml_read'

!

call out_io (s_error$, r_name, trim(msg), &
             'IN FILE: ' // trim(item%file) // '    LINE: ' // int_str(item%i_line), &
             'ITEM: ' // trim(item%name) // ' = ' // trim(item%value))

end subroutine tao_nml_err

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_nml_unknown (item, group_name, err)
!
! Report an item name that the group does not have.
!
! Note: a namelist read of an unknown name is an error that says only "Cannot match namelist
! object name". Naming the item, the group, the file, and the line is the main day to day
! benefit of doing the reading here.
!-

subroutine tao_nml_unknown (item, group_name, err)

type (tao_nml_item_struct), intent(in) :: item
character(*), intent(in) :: group_name
logical, intent(out) :: err

!

err = .true.
call tao_nml_err (item, 'NO SUCH ITEM "' // trim(item%name) // '" IN GROUP "' // &
                                                                  trim(group_name) // '"')

end subroutine tao_nml_unknown

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! The tao_nml_set_* routines behind the tao_nml_value_set generic interface.
!
! Each builds a tao_ptr_struct for the variable and defers to tao_set_ptr_value so that the
! accepted value syntax is identical everywhere, including repeat counts and short value
! lists for arrays. See the notes in tao_set_ptr_value.
!--------------------------------------------------------------------------

subroutine tao_nml_set_ptr (item, ptr, err)

type (tao_nml_item_struct), intent(in) :: item
type (tao_ptr_struct), intent(in) :: ptr

real(rp), allocatable :: eval(:)

character(200) why, why2
character(:), allocatable :: val(:)

integer n_val

logical, intent(out) :: err
logical is_num_target, split_err

! A real or integer scalar may be given an expression rather than a plain number. This is the
! one place the reader accepts something a namelist read would not, and it only ever accepts
! more: a value that is a single plain number takes the ordinary path and behaves exactly as
! before, and the expression evaluator is only consulted for values that a namelist read would
! either have rejected or, worse, silently truncated. Eg "n_opti_cycles = 40 + 60" reads as
! just "40" under list directed input, so the expression has to be tried before the plain read
! rather than only after it fails.

is_num_target = (associated(ptr%r) .or. associated(ptr%i))

if (is_num_target) then
  call tao_nml_split_values (item%value, val, n_val, split_err, why2)
  if (.not. split_err .and. n_val /= 1) then
    if (eval_expression()) return
  endif
endif

! The ordinary path.

call tao_set_ptr_value (ptr, item%value, err, why)
if (.not. err) return

! A single token that is not a number may still be an expression. Eg "pi" or "1/3".

if (is_num_target) then
  if (eval_expression()) return
endif

call tao_nml_err (item, why)

!----------
contains

! Try the value as an expression. Returns True, with err cleared, if it evaluated.

function eval_expression () result (did_it)
logical did_it, eval_err
did_it = .false.
call tao_evaluate_expression (item%value, 1, .false., eval, eval_err, print_err = .false.)
if (eval_err) return
if (associated(ptr%r)) then
  ptr%r = eval(1)
else
  ptr%i = nint(eval(1))
endif
err = .false.
did_it = .true.
end function eval_expression

end subroutine tao_nml_set_ptr

!--------------------------------------------------------------------------

subroutine tao_nml_set_real (item, var, err)
type (tao_nml_item_struct), intent(in) :: item
real(rp), target :: var
type (tao_ptr_struct) ptr
logical, intent(out) :: err
ptr%r => var
call tao_nml_set_ptr (item, ptr, err)
end subroutine tao_nml_set_real

!--------------------------------------------------------------------------

subroutine tao_nml_set_int (item, var, err)
type (tao_nml_item_struct), intent(in) :: item
integer, target :: var
type (tao_ptr_struct) ptr
logical, intent(out) :: err
ptr%i => var
call tao_nml_set_ptr (item, ptr, err)
end subroutine tao_nml_set_int

!--------------------------------------------------------------------------

subroutine tao_nml_set_logic (item, var, err)
type (tao_nml_item_struct), intent(in) :: item
logical, target :: var
type (tao_ptr_struct) ptr
logical, intent(out) :: err
ptr%l => var
call tao_nml_set_ptr (item, ptr, err)
end subroutine tao_nml_set_logic

!--------------------------------------------------------------------------

subroutine tao_nml_set_str (item, var, err)
type (tao_nml_item_struct), intent(in) :: item
character(*), target :: var
type (tao_ptr_struct) ptr
logical, intent(out) :: err
ptr%str => var
call tao_nml_set_ptr (item, ptr, err)
end subroutine tao_nml_set_str

!--------------------------------------------------------------------------

subroutine tao_nml_set_real_arr (item, var, err)
type (tao_nml_item_struct), intent(in) :: item
real(rp), target :: var(:)
type (tao_ptr_struct) ptr
logical, intent(out) :: err
ptr%r1 => var
call tao_nml_set_ptr (item, ptr, err)
end subroutine tao_nml_set_real_arr

!--------------------------------------------------------------------------

subroutine tao_nml_set_int_arr (item, var, err)
type (tao_nml_item_struct), intent(in) :: item
integer, target :: var(:)
type (tao_ptr_struct) ptr
logical, intent(out) :: err
ptr%i1 => var
call tao_nml_set_ptr (item, ptr, err)
end subroutine tao_nml_set_int_arr

!--------------------------------------------------------------------------

subroutine tao_nml_set_logic_arr (item, var, err)
type (tao_nml_item_struct), intent(in) :: item
logical, target :: var(:)
type (tao_ptr_struct) ptr
logical, intent(out) :: err
ptr%l1 => var
call tao_nml_set_ptr (item, ptr, err)
end subroutine tao_nml_set_logic_arr

!--------------------------------------------------------------------------

subroutine tao_nml_set_str_arr (item, var, err)
type (tao_nml_item_struct), intent(in) :: item
character(*), target :: var(:)
type (tao_ptr_struct) ptr
logical, intent(out) :: err
ptr%str1 => var
call tao_nml_set_ptr (item, ptr, err)
end subroutine tao_nml_set_str_arr

end module tao_nml_mod
