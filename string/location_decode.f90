!+
! Subroutine location_decode (string, array, ix_min, num, names, exact_case, can_abbreviate, print_err)
!
! Subroutine to set a list of locations in a logical array to True.
!
! Locations are encoded in the string argument using the following rules:
!    A "," or a space delimits location numbers 
!   "i1:i2" specifies a range of locations between i1 and i2 inclusive.
!     If "i1" is not present then the beginning of the range is taken to be the beginning of the array.
!     If "i2" is not present then the end of the range is taken to be the end of the array.
!   "i1:i2:iskip" specifies a range of locations between i1 and i2 skipping every "iskip" elements.
!   If ix_min is non-negative, a negative number indicates a location that is determined by counting 
!     backwards from the end of the array with "-1" indicating the end of the array.
!   If string is "ALL", "all", or "*" then entire array is set True.
! Examples:
!     string = '3:37 98, 101:115:2'
!
! Also see:
!   location_encode
!   pointer_to_locations
!
! Input:
!   string         -- Character(*): Array of locations.
!   ix_min         -- Integer: Array minimum.
!   names(ix_min:) -- Character(*), optional: Array of location names that can be used
!                       instead used of numbers. Names cannot contain blanks, or ":" 
!   exact_case     -- Logical, optional: Name matching is case sensitive? 
!                       Default is False.
!   can_abbreviate -- Logical, optional: If present and False then abbreviations
!                       are not permitted. Default is True.
!   print_err      -- logical, optional: If present and False, do not print error messages.
!
! Output:
!   array(ix_min:) -- Logical: Array of locations.
!                        ARRAY(I) Set True if location found in STRING.
!   num            -- Integer: Number of locations.
!                        Set to -1 if there is a decode error.
!
! Example:
!     string = '3:37 98, 101:115:2'
!     ix_min = 0
! Then:
!     call location_decode (string, array, ix_min, num)
! Gives:
!     array(I) = True      for i = 3 to 37, and for 98, and
!                           for the odd indexes from 101 to 115
!              = False     otherwise
!     num = 44
!-

subroutine location_decode(string, array, ix_min, num, names, exact_case, can_abbreviate, print_err)

use sim_utils, dummy => location_decode

implicit none

integer index, start_loc, end_loc, step, ix_word, ix_next, num
integer ix_min, ios, where, ix_max
integer, parameter :: no_range$ = 0, found_colon1$ = 1, found_colon2$ = 2

character(*) string
character(len(string)) :: str 
character(2) :: delimitors = ':,'
character(1) delim
character(*), optional :: names(ix_min:)
character(*), parameter :: r_name = 'location_decode'

logical array(ix_min:), found
logical, optional :: exact_case, can_abbreviate, print_err

! initialize array

str = string

array(:) = .false.
where = no_range$
ix_next = 1
ix_max = ubound(array, 1)
num = -1

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
  if (ix_word == 0 .and. delim == '' .and. where == no_range$) exit

  if (str(:ix_word) == 'ALL' .or. str(:ix_word) == 'all' .or. str(:ix_word) == '*') then
    array(:) = .true.
    cycle
  endif

  ! Look for a name match

  found = .false.
  if (present(names) .and. ix_word /= 0) then
    call match_word (str(:ix_word), names, index, exact_case, can_abbreviate)
    if (index < 0) then
      if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'NAME MATCHES TO MULTIPLE LOCATIONS: ' // str(:ix_word))
      return
    elseif (index > 0) then
      found = .true.
      ! Correction since match_word assumes that names(:) has lower bound of 1.
      index = index + (ix_min - 1)  
    endif
  endif

  ! If there is no name match then assume it is a number

  if (.not. found) then
    if (ix_word == 0) then  ! Use default
      if (delim == ':') then
        select case (where)
        case (found_colon2$)
          if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'MISPLACED COLON')
          return
        case (found_colon1$)
          index = ix_max
        case (no_range$)
          index = ix_min
        end select

      else                   ! delim /= ':'
        select case (where)
        case (found_colon1$)
          index = ix_max
        case (no_range$, found_colon2$)
          if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'MISPLACED COMMA')
          return          
        end select
      endif

    else
      read (str(:ix_word), *, iostat = ios) index
      if (ios /= 0) then
        if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'BAD LOCATION: ' // str(:ix_word))
        return
      endif
    endif
  endif

  ! A negative number means count backwards from the end of the array

  if (index < 0 .and. ix_min >= 0) index = ix_max + 1 + index

  ! Check for an error

  if (index < ix_min .or. index > ix_max) then
    if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'LOCATION OUT OF BOUNDS: ' // str(:ix_word))
    return
  endif

  !

  if (where == found_colon2$) step = index
    
  if (delim == ':') then
    select case (where)
    case (found_colon2$)
      if (logic_option(.true., print_err)) call out_io (s_error$, r_name, 'ERROR: BAD STEP(S) ' // str(1:20))
      return
    case (found_colon1$)
      where = found_colon2$
      end_loc = index
    case (no_range$)
      where = found_colon1$
      start_loc = index
    end select

  else
    select case (where)
    case (found_colon2$)
      array(start_loc:end_loc:step) = .true.
      where = no_range$
    case (found_colon1$)
      array(start_loc:index) = .true.
      where = no_range$
    case (no_range$)
      array(index:index) = .true.
    end select
  endif

enddo

! count number of elements in arrray
                  
num = count(array)

end subroutine
                          
