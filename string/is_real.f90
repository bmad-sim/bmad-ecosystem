!+
! Function is_real (string, ignore, real_num) result (good)
!
! Function to test if a string represents a real number. 
! Both 'E' and 'D' are accepted as signifying an exponent.
!
! If the ignore argument is present and True then only the first "word" 
! will be considered and the rest of the line will be ignored. 
! For example:
!   print *, is_real('12.3 45.7', .true)   ! Result: True
!   print *, is_real('12.3 45.7')          ! Result: False
!
! Input:
!   string -- character(*): Character string to check
!   ignore -- logical, optional: Ignore everything after the first word?
!               Default is False.
!
! Output:
!   real_num  -- real(rp), optional: Real number represented by string.
!   good      -- logical: Set True if string represents a real number. 
!                         Set False otherwise.
!-

function is_real (string, ignore, real_num) result (good)

use utilities_mod

implicit none

character(*) string
real(rp), optional :: real_num
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

! look for a digit, '.' or 'E' or 'D'

digit_found = .false.
point_found = .false.
exponent_found = .false.

do
  if (index ('1234567890', string(i:i)) /= 0) then
    digit_found = .true.
  elseif (string(i:i) == '.') then
    if (point_found) goto 8000  ! cannot have two of '.'
    point_found = .true.
  elseif (string(i:i) == 'e' .or. string(i:i) == 'E' .or. string(i:i) == 'd' .or. string(i:i) == 'D') then
    exponent_found = .true.
  elseif (string(i:i) == ' ') then
    exit
  else
    goto 8000
  endif

  i = i + 1
  if (i > len(string)) then
    good = digit_found .and. .not. exponent_found
    goto 8000
  endif

  if (exponent_found) exit

enddo

if (.not. digit_found) goto 8000

! Parse the rest of the exponent if needed

if (exponent_found) then

  digit_found = .false.

  if (string(i:i) == '+' .or. string(i:i) == '-') then
    i = i + 1
    if (i > len(string)) goto 8000
  endif

  do
    if (index ('1234567890', string(i:i)) /= 0) then
      digit_found = .true.
    elseif (string(i:i) == ' ') then
      exit
    else
      goto 8000
    endif

    i = i + 1
    if (i > len(string)) then
      good = digit_found
      goto 8000
    endif
  enddo

  if (.not. digit_found) goto 8000

endif

! look for something more

good = .true.

if (.not. logic_option(.false., ignore)) then ! if not ignore
  if (string(i:) /= ' ') then
    good = .false.
    goto 8000
  endif
endif

!-----------
8000 continue

if (.not. present(real_num) .or. .not. good) return
i = min(i, len(string))
read (string(1:i), *) real_num

end function
