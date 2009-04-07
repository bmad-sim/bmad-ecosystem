!+
! Function is_integer (string)
!
! Function to tell if the first word in a string is a valid integer.
!
! Modules needed:
!   use dcslib
!
! Input:
!   string  -- Character(*): Character string.
!
! Output:
!   is_integer -- Logical: Is a valid integer 
!-

#include "CESR_platform.inc"

function is_integer (string) result (valid)

  character(*) string
  integer i, i1
  logical valid

! ignore beginning spaces

  valid = .false.

  do i1 = 1, len(string)
    if (string(i1:i1) /= ' ') exit
  enddo

  if (i1 == len(string)+1) return

! first char may be a + or -

  if (index('+-', string(i1:i1)) /= 0) i1 = i1 + 1
  if (i1 == len(string)+1) return

! need at least 1 digit

  if (index('1234567890', string(i1:i1)) == 0) return
  i1 = i1 + 1
  valid = .true.

! check for a non-digit

  do i = i1, len(string)
    if (string(i:i) == ' ') return
    if (index('1234567890', string(i:i)) == 0) then
      valid = .false.
      return
    endif
  enddo

end function
