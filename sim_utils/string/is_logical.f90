!+
! Function is_logical (string, ignore) result (good)
!
! Function to test if a string represents a logical.
! Accepted possibilities are (individual characters can be either case):
!   .TRUE.  .FALSE. 
!    TRUE    FALSE
!    T       F
! If the ignore argument is present and True then only the first "word" 
! will be considered and the rest of the line will be ignored. 
! For example:
!   print *, is_logical('F F', .true.)  ! Result: True
!   print *, is_logical('F F')          ! Result: False
!
! Input:
!   string -- Character(*): Character string to check
!   ignore -- Logical, optional: Ignore everything after the first word?
!               Default is False.
!
! Output:
!   good -- Logical: Set True if string represents a logical. 
!                    Set False otherwise.
!-

function is_logical (string, ignore) result (good)

use sim_utils, except => is_logical

implicit none

character(*) string
character(8) tf
logical good
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

! check first word

tf = string(i:)
call str_upcase (tf, tf)

if (tf == '.TRUE. ') then
  i = i + 6
elseif (tf == 'TRUE ') then
  i = i + 4
elseif (tf == 'T ') then
  i = i + 1
elseif (tf == '.FALSE. ') then
  i = i + 7
elseif (tf == 'FALSE ') then
  i = i + 5
elseif (tf == 'F ') then
  i = i + 1
else
  return
endif

good = .true.
if (i > len(string)) return

! check for garbage after the first word

if (.not. logic_option(.false., ignore)) then ! if not ignore
  if (string(i:) /= ' ') then
    good = .false.
    return
  endif
endif

end function

