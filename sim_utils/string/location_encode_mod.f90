module location_encode_mod

!+
! Subroutine location_encode
!
! Subroutine to encode a list of locations.
! This routine is an overloaded name for:
!  location_encode1 (string, loc, exists, ix_min, names, separator, err_flag)
!  location_encode2 (string, loc_array, separator, err_flag)
!
! Input (for location_encode1):
!   loc(ix_min:)    -- logcial: Array of locations.
!                        loc(I) = True means include in the location string
!   exists(ix_min:) -- logical: Array of locations that are relavent.
!   ix_min          -- Integer: Array minimum
!   names(ix_min:)  -- Character(*), optional: Array of names to be used in place
!                       of numbers. 
!   separator       -- Character(*), optional: Separator string between numbers.
!                       Default is a space: " ".
!
! Input (for location_encode2):
!   loc_array(:) -- Integer: Array of locations. 
!                     It is assumed that loc_array(i) is in increasing order.
!   separator       -- Character(*), optional: Separator string between numbers.
!                       Default is a space: " ".
!
! Output:
!   string    -- Character(*): encoded location string. ":" is used 
!                  for a reange of locations
!   err_flag  -- Logical, optional: Set True if encoding error. EG: string too short.
!
! Example: given
!   exists(6) = .false.
!   exists(i) = .true. otherwise
!   loc(-3) = .true.
!   loc(5)  = .true.
!   loc(7)  = .true.
!   loc(57) = .true.
!   loc(j) = .false. for j /= above
!
! With:
!   ix_min = -10
!
! Then:
!   location_encode (string, loc, ix_min)
! Produces:
!   string = '-3 5:7 57'
!-

interface location_encode
  module procedure location_encode1
  module procedure location_encode2
end interface

contains

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine location_encode2 (string, loc_array, separator, err_flag)
!
! Routine to encode a list of locations. 
! This routine is overloaded by the routine: location_encode.
! See location_encode for more details.
!-

subroutine location_encode2 (string, loc_array, separator, err_flag)

implicit none

integer i, range, ix_str
integer loc_array(:)

character(*) string
character(*), optional :: separator

logical err
logical, optional :: err_flag

!

if (present(err_flag)) err_flag = .false.
err = .false.

string = ''
if (size(loc_array) == 0) return

ix_str = 1
call write_to_string ('', loc_array(1))
if (err) return
ix_str = len_trim(string) + 1

range = 0

do i = 2, size(loc_array)

  if (loc_array(i) == loc_array(i-1)) cycle ! Skip duplicates

  if (loc_array(i) == loc_array(i-1) + 1) then
    range = 1
    cycle
  endif

  if (range == 1) then
    call write_to_string (':', loc_array(i-1))
    if (err) return
    ix_str = len_trim(string) + 1
    range = -1
  endif

  if (present(separator)) then
    string(ix_str:) = separator
    ix_str = ix_str + len(separator)
  else
    ix_str = ix_str + 1
  endif

  call write_to_string ('', loc_array(i))
  if (err) return
  ix_str = len_trim(string) + 1
  range = 0

enddo

if (range == 1) call write_to_string (':', loc_array(size(loc_array)))

!------------------------------
contains

subroutine write_to_string (prefix, loc)

character(*) prefix
character(10) str
integer loc

!

write (str, '(a, i0)') prefix, loc
if (len(string) < ix_str + len_trim(str) - 1) then
  err = .true.
  if (present(err_flag)) err_flag = .true.
endif

string(ix_str:) = str 

end subroutine write_to_string

end subroutine location_encode2

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!+
! Subroutine location_encode1 (string, loc, exists, ix_min, names, separator, err_flag)
!
! Routine to encode a list of locations. 
! This routine is overloaded by the routine: location_encode.
! See location_encode for more details.
!-

subroutine location_encode1 (string, loc, exists, ix_min, names, separator, err_flag)

implicit none

integer ix_min
integer n, ix_str, str_max, range

logical loc(ix_min:), exists(ix_min:)
logical separator_needed, err
logical, optional :: err_flag

character(*) string
character(*), optional :: names(ix_min:), separator
character(20) last_name

! init
  
str_max = len(string)
ix_str = 1                 ! Next character in string to write to
string = ''
range = -1                 ! In the middle of a "nn:mm" construct?
                           ! -1 -> no, 0 -> First loc, 1 = in range
separator_needed = .false. ! Need a separator?
err = .false.
if (present(err_flag)) err_flag = .false.

! loop over all locations

do n = lbound(loc, 1), ubound(loc, 1)

  ! ignore locations that do not exist.

  if (.not. exists(n)) cycle

  ! if location is not to be encoded then finish off any range.

  if (.not. loc(n)) then
    if (range == 1) then
      call add_to_string (':' // last_name)
      if (err) return
      ix_str = len_trim(string) + 1 
      separator_needed = .true.
    endif
    range = -1
    cycle
  endif

  ! Construct the name

  if (present(names)) then
    last_name = names(n)
  else
    write (last_name, '(i0)') n
  endif

  ! Range detected? Then just skip to the next loc.

  if (range == 0 .or. range == 1) then
    range = 1
    cycle
  endif

  ! If no range then just add the loc to the list

  if (separator_needed) then
    if (present(separator)) then
      call add_to_string (separator)
      if (err) return
      ix_str = ix_str + len(separator)
    else
      ix_str = ix_str + 1
    endif
    call add_to_string (last_name)
    if (err) return
  endif

  call add_to_string (last_name)
  if (err) return
  ix_str = len_trim(string) + 1
  separator_needed = .true.
  range = 0  

enddo

! cleanup if we are in a range

if (range == 1) call add_to_string (':' // last_name)

!-----------------------------------
contains

subroutine add_to_string (str)

character(*) str

!

if (len(string) < ix_str + len_trim(str) - 1) then
  err = .true.
  if (present(err_flag)) err_flag = .true.
  return
endif

string(ix_str:) = str

end subroutine add_to_string

end subroutine location_encode1
      
end module
