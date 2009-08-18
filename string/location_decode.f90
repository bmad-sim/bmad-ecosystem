!+
! Subroutine location_decode (string, array, ix_min, num, names, exact_case)
!
! Subroutine to set a list of locations in a logical array to True.
!
! Input:
!   string         -- Character*(*): Array of locations.
!                        "," or a space delimits location numbers 
!                        A ":" is used for a range of locations. 
!                        A second ":" specifies a step increment.
!                        "*" or "ALL" Sets everything to True.
!   ix_min         -- Integer: Array minimum.
!   names(ix_min:) -- Character(*), optional: Array of location names that can be used
!                       instead used of numbers. Names cannot contain blanks, or ":" 
!   exact_case     -- Logical, optional: Name matching is case sensitive? 
!                       Default is False.
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

#include "CESR_platform.inc"

subroutine location_decode(string, array, ix_min, num, names, exact_case)

use word_mod

implicit none

integer index, start_loc, end_loc, step, ix_word, ix_next, num
integer ix_min, ios

character(*) string
character(len(string)) :: str 
character(2) :: delimitors = ':,'
character(1) delim
character(*), optional :: names(ix_min:)

logical array(ix_min:), range_found, step_found
logical, optional :: exact_case

! initialize array

str = string

array(:) = .false.
range_found = .false.
step_found = .false.
ix_next = 1

do

  call string_trim2 (str(ix_next:), delimitors, str, ix_word, delim, ix_next)

  if (ix_word == 0) exit

  if (str(:ix_word) == 'ALL' .or. str(:ix_word) == 'all' .or. str(:ix_word) == '*') then
    array(:) = .true.
    if (ix_next == 0) exit
    cycle
  endif

  ! Look for a name match

  index = 0
  if (present(names)) then
    call match_word (str(:ix_word), names, index, exact_case)
    if (index < 0) then
      print *, 'ERROR: NAME MATCHES TO MULTIPLE LOCATIONS: ' // str(:ix_word)
      num = -1
      return
    endif
  endif

  ! If there is no name match then assume it is a number

  if (index == 0) then
    read (str(:ix_word), *, iostat = ios) index
    if (ios /= 0) then
      print *, 'ERROR: BAD LOCATION: ', str(:ix_word)
      num = -1
      return
    endif
  endif

  ! Check for an error

  if (index > ubound(array,1)) then
    print *, 'ERROR: LOCATION OUT OF BOUNDS: ', str(:ix_word)
    num = -1
    return
  endif

  if (step_found) step = index
    
  if (delim == ':') then
    if (range_found) then
      if (step_found) then
        print *, 'ERROR: BAD STEP(S) ', str(1:20)
        num = -1
        return
      else
        step_found = .true.
        end_loc = index
      endif
    else
      range_found = .true.
      start_loc = index
    endif
  else
    if (range_found) then
      if (step_found) then
        array(start_loc:end_loc:step) = .true.
        range_found = .false.
        step_found = .false.
      else
        array(start_loc:index) = .true.
        range_found = .false.
      endif
    else
      array(index:index) = .true.
    endif
  endif

  if (ix_next == 0) exit

enddo

!--------

if (range_found) then
  print *, 'ERROR IN LOCATION_DECODE: OPEN RANGE'
  num = -1
  return
endif

! count number of elements in arrray
                  
num = count(array)

end subroutine
                          
