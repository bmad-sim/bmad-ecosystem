!+
! Module tao_attrib_ptr_mod
!
! Runtime support for the generated structure component resolvers in tao_attrib_resolve_mod.
!
! The resolvers map a component name string such as
!     "floor_plan%ele_shape(3)%ele_id"
! onto a pointer to the corresponding component of a structure instance.  This replaces the use
! of Fortran namelist reads for the purpose of setting structure components from a string.
!
! Note: sim_utils has all_pointer_struct which serves a similar purpose but has no character
! pointer component and so cannot be used here.  tao_ptr_struct is the Tao local equivalent.
!
! The closure of structures reachable from the Tao input structures contains only rank 0 and
! rank 1 components of type real, integer, logical, and character.  tao_ptr_struct therefore
! has no rank 2 or complex components.  If that ever changes, generate_attrib_tables.py will
! flag the offending component rather than silently ignore it.
!-

module tao_attrib_ptr_mod

use sim_utils

implicit none

! Pointer to a structure component.
! At most one component of a given tao_ptr_struct instance is associated.

type tao_ptr_struct
  real(rp), pointer :: r => null()
  integer, pointer :: i => null()
  logical, pointer :: l => null()
  character(:), pointer :: str => null()
  real(rp), pointer :: r1(:) => null()
  integer, pointer :: i1(:) => null()
  logical, pointer :: l1(:) => null()
  character(:), pointer :: str1(:) => null()
end type

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_split_attrib_name (name, head, has_sub, isub, rest, err, why)
!
! Split off the leading component of a structure component reference.
!
! Example: name = "ele_shape(3)%ele_id" gives
!   head    = "ele_shape"
!   has_sub = True
!   isub    = 3
!   rest    = "ele_id"
!
! A subscript may be a range, as in "var(1:4)%ele_name", which a namelist read takes to mean
! that the value list is spread over those elements. A range is only accepted if the caller
! asks for it by supplying isub2. The generated resolvers do not, since a resolver returns a
! pointer to one component, so a range reaching a resolver is reported as an error rather than
! being silently taken as its first element.
!
! Input:
!   name    -- character(*): Component reference. Case insensitive.
!
! Output:
!   head    -- character(:), allocatable: Leading component name, down cased.
!   has_sub -- logical: True if a subscript was present.
!   isub    -- integer: Subscript value, or the first value of a range. 0 if has_sub is False.
!   rest    -- character(:), allocatable: Remainder of the reference. Blank if none.
!   err     -- logical: Set True on a malformed reference.
!   why     -- character(*): Explanation if err is True.
!   isub2   -- integer, optional: Last value of a subscript range. Equal to isub if the
!                subscript is a single value. If not present then a range is an error.
!-

subroutine tao_split_attrib_name (name, head, has_sub, isub, rest, err, why, isub2)

character(*), intent(in) :: name
character(:), allocatable, intent(out) :: head, rest
character(*), intent(out) :: why
character(:), allocatable :: nam, sub_str, s1, s2

integer, intent(out) :: isub
integer, optional, intent(out) :: isub2
integer ix_pct, ix_lp, ix_rp, ix_colon, i2

logical, intent(out) :: has_sub, err

!

err = .false.
why = ''
has_sub = .false.
isub = 0
if (present(isub2)) isub2 = 0
head = ''
rest = ''

nam = trim(adjustl(downcase(name)))

if (nam == '') then
  err = .true.
  why = 'BLANK COMPONENT NAME'
  return
endif

! Split at the first "%" that is not inside a subscript.

ix_pct = index_outside_paren(nam, '%')

if (ix_pct == 0) then
  head = nam
else
  head = nam(1:ix_pct-1)
  rest = nam(ix_pct+1:)
  if (rest == '') then
    err = .true.
    why = 'TRAILING "%" IN COMPONENT NAME: ' // trim(name)
    return
  endif
endif

! Pull off a subscript if present.

ix_lp = index(head, '(')

if (ix_lp /= 0) then
  ix_rp = index(head, ')', back = .true.)
  if (ix_rp /= len(head)) then
    err = .true.
    why = 'MALFORMED SUBSCRIPT IN: ' // trim(name)
    return
  endif
  sub_str = trim(adjustl(head(ix_lp+1:ix_rp-1)))
  head = trim(head(1:ix_lp-1))
  if (head == '') then
    err = .true.
    why = 'SUBSCRIPT WITH NO COMPONENT NAME IN: ' // trim(name)
    return
  endif
  ! A subscript range "lo:hi".

  ix_colon = index(sub_str, ':')

  if (ix_colon /= 0) then
    if (.not. present(isub2)) then
      err = .true.
      why = 'A SUBSCRIPT RANGE IS NOT ALLOWED HERE: ' // trim(name)
      return
    endif
    s1 = trim(adjustl(sub_str(1:ix_colon-1)))
    s2 = trim(adjustl(sub_str(ix_colon+1:)))

    ! An open ended bound, as in "var(1:)" or "var(:4)" or "var(:)", is returned as
    ! int_garbage$ for the caller to replace with the array bound. tao_split_attrib_name does
    ! not know the shape of the thing being referenced.

    if (s1 == '') then
      isub = int_garbage$
    elseif (.not. is_integer(s1, isub)) then
      err = .true.
      why = 'SUBSCRIPT RANGE BOUND IS NOT AN INTEGER: ' // trim(name)
      return
    endif

    if (s2 == '') then
      i2 = int_garbage$
    elseif (.not. is_integer(s2, i2)) then
      err = .true.
      why = 'SUBSCRIPT RANGE BOUND IS NOT AN INTEGER: ' // trim(name)
      return
    endif

    if (isub /= int_garbage$ .and. i2 /= int_garbage$ .and. i2 < isub) then
      err = .true.
      why = 'SUBSCRIPT RANGE IS EMPTY: ' // trim(name)
      return
    endif

    isub2 = i2
    has_sub = .true.

  else
    if (.not. is_integer(sub_str, isub)) then
      err = .true.
      why = 'SUBSCRIPT IS NOT AN INTEGER: ' // trim(name)
      return
    endif
    if (present(isub2)) isub2 = isub
    has_sub = .true.
  endif
endif

if (head == '') then
  err = .true.
  why = 'BLANK COMPONENT NAME IN: ' // trim(name)
  return
endif

end subroutine tao_split_attrib_name

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function index_outside_paren (str, ch) result (ix)
!
! Index of the first occurrence of ch in str that is not nested inside parentheses.
! Returns 0 if there is none.
!-

function index_outside_paren (str, ch) result (ix)

character(*), intent(in) :: str, ch
integer ix, i, depth

!

depth = 0

do i = 1, len(str)
  if (str(i:i) == '(') then
    depth = depth + 1
  elseif (str(i:i) == ')') then
    depth = depth - 1
  elseif (depth == 0 .and. str(i:i) == ch) then
    ix = i
    return
  endif
enddo

ix = 0

end function index_outside_paren

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_ptr_is_associated (ptr) result (is_assoc)
!
! True if any component of ptr is associated.
!-

function tao_ptr_is_associated (ptr) result (is_assoc)

type (tao_ptr_struct), intent(in) :: ptr
logical is_assoc

!

is_assoc = (associated(ptr%r) .or. associated(ptr%i) .or. associated(ptr%l) .or. &
            associated(ptr%str) .or. associated(ptr%r1) .or. associated(ptr%i1) .or. &
            associated(ptr%l1) .or. associated(ptr%str1))

end function tao_ptr_is_associated

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Function tao_ptr_type_name (ptr) result (name)
!
! Name of the type that ptr points to. Used for error messages.
!-

function tao_ptr_type_name (ptr) result (name)

type (tao_ptr_struct), intent(in) :: ptr
character(20) name

!

if (associated(ptr%r)) then;         name = 'REAL'
elseif (associated(ptr%i)) then;     name = 'INTEGER'
elseif (associated(ptr%l)) then;     name = 'LOGICAL'
elseif (associated(ptr%str)) then;   name = 'CHARACTER'
elseif (associated(ptr%r1)) then;    name = 'REAL ARRAY'
elseif (associated(ptr%i1)) then;    name = 'INTEGER ARRAY'
elseif (associated(ptr%l1)) then;    name = 'LOGICAL ARRAY'
elseif (associated(ptr%str1)) then;  name = 'CHARACTER ARRAY'
else;                                name = 'NOT ASSOCIATED'
endif

end function tao_ptr_type_name

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine tao_set_ptr_value (ptr, value_str, err, why)
!
! Set the component pointed to by ptr from a string.
!
! Note on back compatibility: the conversion is done with an internal list directed read, which
! is the same mechanism a namelist read uses. This is deliberate. It means the set of value
! strings accepted here is the same as the set the namelist scratch file read accepted
! (T, .true., TRUE, 1.5d0, 1e-3, quoted and unquoted strings, comma separated array values,
! repeat counts such as 3*0.0, etc.) without having to enumerate them.
!
! A plain list directed read is not by itself equivalent to a namelist read in two respects,
! both of which are handled explicitly below:
!
!   1) Excess values. A namelist read of a scalar rejects "x = 2*5" because the repeat count
!      supplies two values for a scalar. A plain list directed read quietly takes the first
!      and ignores the rest. So the item count is checked: a scalar target requires exactly
!      one value and an array target requires no more than size(array) values.
!
!   2) Short value lists for an array. A namelist read of "d_orb = 1e-6, 2e-6" sets the first
!      two elements and leaves the remaining four unchanged. A plain list directed read into
!      the whole array instead hits end of record and fails. So an array is filled only as far
!      as the values provided, counted by reading successively longer leading sections.
!
! The one intentional difference from a namelist read: a blank value string is an error here.
! A namelist read silently ignores a blank value field, so "set bmad_com sr_wakes_on = " used
! to do nothing at all and report success. No previously working input is rejected.
!
! Input:
!   ptr       -- tao_ptr_struct: Pointer to the component to set.
!   value_str -- character(*): Value to set the component to.
!
! Output:
!   err       -- logical: Set True on error.
!   why       -- character(*): Explanation if err is True.
!-

subroutine tao_set_ptr_value (ptr, value_str, err, why)

type (tao_ptr_struct), intent(in) :: ptr

character(*), intent(in) :: value_str
character(*), intent(out) :: why

integer ios, n, lb

logical, intent(out) :: err

! Scalars. The second read of one extra item must fail, otherwise too many values were given.

err = .true.
ios = 0

if (value_str == '') then
  why = 'VALUE IS BLANK'
  return
endif

if (associated(ptr%r)) then
  block
    real(rp) v1, v2
    read (value_str, *, iostat = ios) v1
    if (ios /= 0) goto 9000
    read (value_str, *, iostat = ios) v1, v2
    if (ios == 0) goto 9100
    ptr%r = v1
  end block

elseif (associated(ptr%i)) then
  block
    integer v1, v2
    read (value_str, *, iostat = ios) v1
    if (ios /= 0) goto 9000
    read (value_str, *, iostat = ios) v1, v2
    if (ios == 0) goto 9100
    ptr%i = v1
  end block

elseif (associated(ptr%l)) then
  block
    logical v1, v2
    read (value_str, *, iostat = ios) v1
    if (ios /= 0) goto 9000
    read (value_str, *, iostat = ios) v1, v2
    if (ios == 0) goto 9100
    ptr%l = v1
  end block

elseif (associated(ptr%str)) then
  block
    character(len(ptr%str)) v1, v2
    read (value_str, *, iostat = ios) v1
    if (ios /= 0) goto 9000
    read (value_str, *, iostat = ios) v1, v2
    if (ios == 0) goto 9100
    ptr%str = v1
  end block

! Arrays. Read successively longer leading sections to find how many values were supplied,
! then assign only that many so that trailing elements keep their present value.

elseif (associated(ptr%r1)) then
  block
    real(rp) v(size(ptr%r1)+1)
    ! Seed with the present values so that a null value in the list, as in "a = 1.0,,3.0",
    ! leaves that element unchanged, which is what a namelist read does.
    v(1:size(ptr%r1)) = ptr%r1
    n = n_items_read_r(value_str, v, size(ptr%r1))
    if (n < 0) goto 9100
    if (n == 0) goto 9000
    lb = lbound(ptr%r1, 1)
    ptr%r1(lb:lb+n-1) = v(1:n)
  end block

elseif (associated(ptr%i1)) then
  block
    integer v(size(ptr%i1)+1)
    ! Seed with the present values so that a null value in the list, as in "a = 1.0,,3.0",
    ! leaves that element unchanged, which is what a namelist read does.
    v(1:size(ptr%i1)) = ptr%i1
    n = n_items_read_i(value_str, v, size(ptr%i1))
    if (n < 0) goto 9100
    if (n == 0) goto 9000
    lb = lbound(ptr%i1, 1)
    ptr%i1(lb:lb+n-1) = v(1:n)
  end block

elseif (associated(ptr%l1)) then
  block
    logical v(size(ptr%l1)+1)
    ! Seed with the present values so that a null value in the list, as in "a = 1.0,,3.0",
    ! leaves that element unchanged, which is what a namelist read does.
    v(1:size(ptr%l1)) = ptr%l1
    n = n_items_read_l(value_str, v, size(ptr%l1))
    if (n < 0) goto 9100
    if (n == 0) goto 9000
    lb = lbound(ptr%l1, 1)
    ptr%l1(lb:lb+n-1) = v(1:n)
  end block

elseif (associated(ptr%str1)) then
  block
    character(len(ptr%str1)) v(size(ptr%str1)+1)
    ! See the note above about seeding with the present values.
    v(1:size(ptr%str1)) = ptr%str1
    n = n_items_read_str(value_str, v, size(ptr%str1))
    if (n < 0) goto 9100
    if (n == 0) goto 9000
    lb = lbound(ptr%str1, 1)
    ptr%str1(lb:lb+n-1) = v(1:n)
  end block

else
  why = 'INTERNAL ERROR: POINTER NOT ASSOCIATED'
  return
endif

err = .false.
why = ''
return

!

9000 continue
why = 'CANNOT DECODE THE VALUE ' // trim(value_str) // ' AS TYPE ' // &
                                                        trim(tao_ptr_type_name(ptr))
return

9100 continue
why = 'TOO MANY VALUES GIVEN FOR ' // trim(tao_ptr_type_name(ptr)) // &
                                                        ' COMPONENT: ' // trim(value_str)
return

end subroutine tao_set_ptr_value

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! The n_items_read_* functions return the number of values that a list directed read of str
! supplies, up to a maximum of n_max. A return of 0 means the first value could not be decoded
! and a return of -1 means more than n_max values were supplied.
!
! The count is found by reading successively longer leading sections of the target array. This
! is O(n_max^2) reads but n_max is the size of a structure component array, which is small, and
! it gets repeat counts such as "3*0.0" right for free since the list directed reader expands
! them.
!--------------------------------------------------------------------------

function n_items_read_r (str, v, n_max) result (n)
character(*), intent(in) :: str
real(rp), intent(out) :: v(:)
integer, intent(in) :: n_max
integer n, i, ios
n = 0
do i = 1, n_max + 1
  read (str, *, iostat = ios) v(1:i)
  if (ios /= 0) return
  n = i
enddo
n = -1
end function n_items_read_r

!--------------------------------------------------------------------------

function n_items_read_i (str, v, n_max) result (n)
character(*), intent(in) :: str
integer, intent(out) :: v(:)
integer, intent(in) :: n_max
integer n, i, ios
n = 0
do i = 1, n_max + 1
  read (str, *, iostat = ios) v(1:i)
  if (ios /= 0) return
  n = i
enddo
n = -1
end function n_items_read_i

!--------------------------------------------------------------------------

function n_items_read_l (str, v, n_max) result (n)
character(*), intent(in) :: str
logical, intent(out) :: v(:)
integer, intent(in) :: n_max
integer n, i, ios
n = 0
do i = 1, n_max + 1
  read (str, *, iostat = ios) v(1:i)
  if (ios /= 0) return
  n = i
enddo
n = -1
end function n_items_read_l

!--------------------------------------------------------------------------

function n_items_read_str (str, v, n_max) result (n)
character(*), intent(in) :: str
character(*), intent(out) :: v(:)
integer, intent(in) :: n_max
integer n, i, ios
n = 0
do i = 1, n_max + 1
  read (str, *, iostat = ios) v(1:i)
  if (ios /= 0) return
  n = i
enddo
n = -1
end function n_items_read_str

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! The following tao_res_*_bad functions are the guard clauses used by the generated
! resolvers in tao_attrib_resolve_mod. Each returns True, and sets err and why, if the
! reference is not valid for the kind of component involved.
!--------------------------------------------------------------------------
!+
! Function tao_res_scalar_bad (head, has_sub, rest, err, why) result (bad)
!
! Guard for a scalar component of intrinsic type.
!-

function tao_res_scalar_bad (head, has_sub, rest, err, why) result (bad)

character(*), intent(in) :: head, rest
character(*), intent(out) :: why
logical, intent(in) :: has_sub
logical, intent(out) :: err
logical bad

!

bad = .true.
err = .true.

if (has_sub) then
  why = 'COMPONENT IS NOT AN ARRAY SO A SUBSCRIPT IS NOT ALLOWED: ' // head
  return
endif

if (rest /= '') then
  why = 'COMPONENT IS NOT A STRUCTURE SO IT HAS NO SUBCOMPONENT ' // trim(rest) // ': ' // head
  return
endif

bad = .false.
err = .false.
why = ''

end function tao_res_scalar_bad

!--------------------------------------------------------------------------
!+
! Function tao_res_array_bad (head, rest, has_sub, isub, lo, hi, err, why) result (bad)
!
! Guard for a rank 1 array component of intrinsic type.
! A missing subscript is allowed and means the whole array.
!-

function tao_res_array_bad (head, rest, has_sub, isub, lo, hi, err, why) result (bad)

character(*), intent(in) :: head, rest
character(*), intent(out) :: why
integer, intent(in) :: isub, lo, hi
logical, intent(in) :: has_sub
logical, intent(out) :: err
logical bad

!

bad = .true.
err = .true.

if (rest /= '') then
  why = 'COMPONENT IS NOT A STRUCTURE SO IT HAS NO SUBCOMPONENT ' // trim(rest) // ': ' // head
  return
endif

if (has_sub) then
  if (isub < lo .or. isub > hi) then
    why = 'SUBSCRIPT ' // int_str(isub) // ' IS OUT OF RANGE [' // int_str(lo) // ', ' // &
                                                    int_str(hi) // '] FOR COMPONENT: ' // head
    return
  endif
endif

bad = .false.
err = .false.
why = ''

end function tao_res_array_bad

!--------------------------------------------------------------------------
!+
! Function tao_res_struct_bad (head, has_sub, rest, err, why) result (bad)
!
! Guard for a scalar component of derived type. A subcomponent must be named.
!-

function tao_res_struct_bad (head, has_sub, rest, err, why) result (bad)

character(*), intent(in) :: head, rest
character(*), intent(out) :: why
logical, intent(in) :: has_sub
logical, intent(out) :: err
logical bad

!

bad = .true.
err = .true.

if (has_sub) then
  why = 'COMPONENT IS NOT AN ARRAY SO A SUBSCRIPT IS NOT ALLOWED: ' // head
  return
endif

if (rest == '') then
  why = 'COMPONENT IS A STRUCTURE SO A SUBCOMPONENT MUST BE GIVEN. EG: ' // &
                                                            trim(head) // '%<SUBCOMPONENT>'
  return
endif

bad = .false.
err = .false.
why = ''

end function tao_res_struct_bad

!--------------------------------------------------------------------------
!+
! Function tao_res_struct_array_bad (head, rest, has_sub, isub, lo, hi, err, why) result (bad)
!
! Guard for a rank 1 array component of derived type.
! Both a subscript and a subcomponent are required.
!-

function tao_res_struct_array_bad (head, rest, has_sub, isub, lo, hi, err, why) result (bad)

character(*), intent(in) :: head, rest
character(*), intent(out) :: why
integer, intent(in) :: isub, lo, hi
logical, intent(in) :: has_sub
logical, intent(out) :: err
logical bad

!

bad = .true.
err = .true.

if (rest == '') then
  why = 'COMPONENT IS A STRUCTURE SO A SUBCOMPONENT MUST BE GIVEN. EG: ' // &
                                                      trim(head) // '(1)%<SUBCOMPONENT>'
  return
endif

if (.not. has_sub) then
  why = 'COMPONENT IS AN ARRAY OF STRUCTURES SO A SUBSCRIPT MUST BE GIVEN: ' // head
  return
endif

if (isub < lo .or. isub > hi) then
  why = 'SUBSCRIPT ' // int_str(isub) // ' IS OUT OF RANGE [' // int_str(lo) // ', ' // &
                                                  int_str(hi) // '] FOR COMPONENT: ' // head
  return
endif

bad = .false.
err = .false.
why = ''

end function tao_res_struct_array_bad

!--------------------------------------------------------------------------
!+
! Function tao_res_struct_bad_for_slot (head, has_sub, err, why) result (bad)
!
! Guard for a scalar structure component reached during a positional assignment. Unlike
! tao_res_struct_bad, a blank remainder is fine here: it means the value list is to be spread
! over the components of that structure.
!-

function tao_res_struct_bad_for_slot (head, has_sub, err, why) result (bad)

character(*), intent(in) :: head
character(*), intent(out) :: why
logical, intent(in) :: has_sub
logical, intent(out) :: err
logical bad

!

bad = .false.
err = .false.
why = ''

if (has_sub) then
  bad = .true.
  err = .true.
  why = 'COMPONENT IS NOT AN ARRAY SO A SUBSCRIPT IS NOT ALLOWED: ' // head
endif

end function tao_res_struct_bad_for_slot

!--------------------------------------------------------------------------
!+
! Function tao_res_struct_array_bad_for_slot (head, has_sub, isub, lo, hi, err, why) result (bad)
!
! Guard for an array of structures reached during a positional assignment.
!-

function tao_res_struct_array_bad_for_slot (head, has_sub, isub, lo, hi, err, why) result (bad)

character(*), intent(in) :: head
character(*), intent(out) :: why
integer, intent(in) :: isub, lo, hi
logical, intent(in) :: has_sub
logical, intent(out) :: err
logical bad

!

bad = .true.
err = .true.

if (.not. has_sub) then
  why = 'COMPONENT IS AN ARRAY OF STRUCTURES SO A SUBSCRIPT MUST BE GIVEN: ' // head
  return
endif

if (isub < lo .or. isub > hi) then
  why = 'SUBSCRIPT ' // int_str(isub) // ' IS OUT OF RANGE [' // int_str(lo) // ', ' // &
                                                  int_str(hi) // '] FOR COMPONENT: ' // head
  return
endif

bad = .false.
err = .false.
why = ''

end function tao_res_struct_array_bad_for_slot

!--------------------------------------------------------------------------
!+
! Function tao_res_slot_leaf_bad (head, rest, err, why) result (bad)
!
! Guard for a component of intrinsic type reached during a positional assignment. There is
! nothing below it to name.
!-

function tao_res_slot_leaf_bad (head, rest, err, why) result (bad)

character(*), intent(in) :: head, rest
character(*), intent(out) :: why
logical, intent(out) :: err
logical bad

!

bad = .false.
err = .false.
why = ''

if (rest /= '') then
  bad = .true.
  err = .true.
  why = 'COMPONENT IS NOT A STRUCTURE SO IT HAS NO SUBCOMPONENT ' // trim(rest) // ': ' // head
endif

end function tao_res_slot_leaf_bad

!--------------------------------------------------------------------------
!+
! Function tao_res_alloc_bad (head, is_alloc, err, why) result (bad)
!
! Guard for an allocatable component. It must be allocated to be referenced.
!-

function tao_res_alloc_bad (head, is_alloc, err, why) result (bad)

character(*), intent(in) :: head
character(*), intent(out) :: why
logical, intent(in) :: is_alloc
logical, intent(out) :: err
logical bad

!

if (is_alloc) then
  bad = .false.
  err = .false.
  why = ''
  return
endif

bad = .true.
err = .true.
why = 'COMPONENT IS NOT ALLOCATED SO IT CANNOT BE SET: ' // head

end function tao_res_alloc_bad

end module tao_attrib_ptr_mod
