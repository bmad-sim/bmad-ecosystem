!+
! Function is_integer (string, int) result (valid)
!
! Function to tell if the first word in a string is a valid integer.
! A blank string is considered to NOT be a valid integer.
!
! Input:
!   string  -- Character(*): Character string.
!
! Output:
!   int     -- integer, optional: Integer value. Default is -999 if not valid.
!   valid   -- Logical: Is a valid integer 
!-

function is_integer (string, int) result (valid)

implicit none

character(*) string
integer i, i1
integer, optional :: int
logical valid

! ignore beginning spaces

valid = .false.
if (present(int)) int = -999

if (string == '') return

do i1 = 1, len(string)
  if (string(i1:i1) /= ' ') exit
enddo

! first char may be a + or -

if (index('+-', string(i1:i1)) /= 0) i1 = i1 + 1
if (i1 == len(string)+1) return

! need at least 1 digit

if (index('1234567890', string(i1:i1)) == 0) return
i1 = i1 + 1
valid = .true.

! check for a non-digit

do i = i1, len(string)
  if (string(i:i) == ' ') exit
  if (index('1234567890', string(i:i)) == 0) then
    valid = .false.
    return
  endif
enddo

if (present(int)) read(string, *) int

end function
