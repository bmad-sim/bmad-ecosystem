!+
! Subroutine pointer_to_locations (string, array, num, ix_min, ix_max, names, exact_case, print_err)
!
! Routine to return an integer array that corresponds to the array encoded in a string.
!
! Locations are encoded in the string argument using the following rules:
!   A "," or a space delimits location numbers 
!   "i1:i2" specifies a range of locations between i1 and i2 inclusive.
!     If "i1" is not present then the beginning of the range is taken to be the beginning of the array.
!     If "i2" is not present then the end of the range is taken to be the end of the array.
!   If string = "*", array will be set to [ix_min:ix_max]
! Example:
!     string = '3:5 98, 50' -> array = [3, 4, 5, 98, 50]
!
! Also see:
!   location_encode
!   location_decode
!
! Input:
!   string         -- Character(*): Array of locations.
!   ix_min         -- Integer: Locations lower bound.
!   ix_max         -- Integer: Locations upper bound.
!   names(ix_min:) -- Character(*), optional: Array of location names that can be used
!                       instead used of numbers. Names cannot contain blanks, or ":" 
!   exact_case     -- Logical, optional: Name matching is case sensitive? 
!                       Default is False.
!   print_err      -- logical, optional: If present and False, do not print error messages.
!
! Output:
!   array(:)       -- integer, allocatable: Array of locations.
!   num            -- Integer: Number of locations.
!                        Set to -1 if there is a decode error.
!-

subroutine pointer_to_locations(string, array, num, ix_min, ix_max, names, exact_case, print_err)

use sim_utils, dummy => pointer_to_locations

implicit none

integer, allocatable :: array(:)
integer :: num, ix_min, ix_max

integer indx, start_loc, end_loc, ix_word, ix_next
integer i, ios, n_num, n_add

character(*) string
character(len(string)) :: str 
character(2) :: delimitors = ':,'
character(1) delim
character(*), optional :: names(:)
character(*), parameter :: r_name = 'pointer_to_locations'

logical :: found_colon
logical found
logical, optional :: exact_case, print_err

! initialize

found_colon = .false.
ix_next = 1
num = -1
n_num = 0

call string_trim (string, str, i)

!

if (string(1:1) == '*') then
  call re_allocate(array, ix_max-ix_min+1, .false.)
  do i = 1, size(array)
    array(i) = i - ix_min + 1
  enddo

  if (str == '') then
    num = ix_max - ix_min + 1
  else
    if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'EXTRA SUFF AFTER "*"')
  endif
  return
endif

!

do
  if (ix_next == 0) then
    if (delim == ',') then
      if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'MISPLACED COMMA')
      return          
    endif
    str = '';  ix_word = 0;  delim = ''
  else
    call string_trim2 (str(ix_next:), delimitors, str, ix_word, delim, ix_next)
  endif
  if (ix_word == 0 .and. delim == '' .and. .not. found_colon) exit

  ! Look for a name match

  found = .false.
  if (present(names) .and. ix_word /= 0) then
    call match_word (str(:ix_word), names, indx, exact_case)
    if (indx < 0) then
      if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'NAME MATCHES TO MULTIPLE LOCATIONS: ' // str(:ix_word))
      return
    elseif (indx > 0) then
      found = .true.
      ! Correction since match_word assumes that names(:) has lower bound of 1.
      indx = indx + (ix_min - 1)  
    endif
  endif

  ! If there is no name match then assume it is a number

  if (.not. found) then
    if (ix_word == 0) then  ! Use default
      if (delim == ':') then
        if (found_colon) then
          indx = ix_max
        else
          indx = ix_min
        endif

      else                   ! delim /= ':'
        if (found_colon) then
          indx = ix_max
        else
          if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'MISPLACED COMMA')
          return          
        endif
      endif

    else
      read (str(:ix_word), *, iostat = ios) indx
      if (ios /= 0) then
        if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'BAD LOCATION: ' // str(:ix_word))
        return
      endif
    endif
  endif

  ! A negative number means count backwards from the end of the array

  if (indx < 0 .and. ix_min >= 0) indx = ix_max + 1 + indx

  ! Check for an error

  if (indx < ix_min .or. indx > ix_max) then
    if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'LOCATION OUT OF BOUNDS: ' // str(:ix_word))
    return
  endif

  !

  if (delim == ':') then
    if (found_colon) then
      if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'ERROR: EXTRA COLON FOUND ' // str)
      return
    else
      found_colon = .true.
      start_loc = indx
    endif

  else
    if (found_colon) then
      n_add = indx - start_loc + 1
      if (n_num == 0) allocate (array(n_add + 10))
      if (n_num+n_add > size(array)) call re_allocate(array, n_num+n_add+10, .false.)
      do i = 1, n_add
        array(i+n_num) = start_loc + i - 1
      enddo
      n_num = n_num + n_add
      found_colon = .false.
    else
      if (n_num == 0) allocate (array(10))
      if (n_num+1 > size(array)) call re_allocate(array, n_num+10, .false.)
      n_num = n_num + 1
      array(n_num) = indx
    endif
  endif

enddo

! count number of elements in arrray
                  
num = n_num

end subroutine
                          
