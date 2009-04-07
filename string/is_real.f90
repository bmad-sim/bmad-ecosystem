!+
! Function is_real (string, ignore) result (good)
!
! Function to test if a string represents a real number.
! If the ignore argument is present and True then only the first "word" 
! will be considered and the rest of the line will be ignored. 
! For example:
!   print *, is_real('12.3 45.7', .true)   ! Result: True
!   print *, is_real('12.3 45.7')          ! Result: False
!
! Input:
!   string -- Character(*): Character string to check
!   ignore -- Logical, optional: Ignore everything after the first word?
!               Default is False.
!
! Output:
!   good -- Logical: Set True if string represents a real number. 
!                    Set False otherwise.
!-

function is_real (string, ignore) result (good)

use utilities_mod

implicit none

character(*) string
logical good, digit_found, point_found, exponent_found
logical, optional :: ignore
integer i

! first skip beginning white space

good = .false.

i = 1
do
  if (string(i:i) /= ' ') exit
  i = i + 1
  if (i > len(string)) return
enddo

! look for beginning "+" or "-" sign

if (string(i:i) == '+' .or. string(i:i) == '-') then
  i = i + 1
  if (i > len(string)) return
endif

! look for a digit, '.', or 'e'

digit_found = .false.
point_found = .false.
exponent_found = .false.

do
  if (index ('1234567890', string(i:i)) /= 0) then
    digit_found = .true.
  elseif (string(i:i) == '.') then
    if (point_found) return  ! cannot have two of '.'
    point_found = .true.
  elseif (string(i:i) == 'e' .or. string(i:i) == 'E') then
    exponent_found = .true.
  elseif (string(i:i) == ' ') then
    exit
  else
    return
  endif

  i = i + 1
  if (i > len(string)) then
    good = digit_found .and. .not. exponent_found
    return
  endif

  if (exponent_found) exit

enddo

if (.not. digit_found) return

! Parse the rest of the exponent if needed

if (exponent_found) then

  digit_found = .false.

  if (string(i:i) == '+' .or. string(i:i) == '-') then
    i = i + 1
    if (i > len(string)) return
  endif

  do
    if (index ('1234567890', string(i:i)) /= 0) then
      digit_found = .true.
    elseif (string(i:i) == ' ') then
      exit
    else
      return
    endif

    i = i + 1
    if (i > len(string)) then
      good = digit_found
      return
    endif
  enddo

  if (.not. digit_found) return

endif

! look for something more

good = .true.

if (.not. logic_option(.false., ignore)) then ! if not ignore
  if (string(i:) /= ' ') then
    good = .false.
    return
  endif
endif

end function
